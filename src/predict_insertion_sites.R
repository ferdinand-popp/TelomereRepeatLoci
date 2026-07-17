# Author: Lina Sieverling

# Usage:
#   R --no-save --slave --args --candidate_region_file <file> --clipped_reads_file <file> --discordant_read_file <file> \
#     --outfile <file> --function_file <file> --bamfile_tumor <bam> [--bamfile_control <bam> --clipped_reads_control_file <file>]
# Description: trys to predict a telomere insertion site for each candidate region from clipped reads of tumor sample
#              - takes the position where most clipped sequences start/end (if this is not unique it returns NA)
#              - add the result to the extended candidate region table
#
#              For each predicted site we additionally count, for tumor and control, the total reads at the site
#              and the subset of unique clipped reads whose clipped interval overlaps the 1-bp site window.
#              The clipped-read counts are deduplicated by read_name to avoid double counting soft-/hard-clipped
#              evidence for the same physical read.
#

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("--candidate_region_file"), type="character", help="candidate region table", metavar="file"),
  make_option(c("--clipped_reads_file"), type="character", help="tumor clipped reads table", metavar="file"),
  make_option(c("--discordant_read_file"), type="character", help="tumor discordant reads table", metavar="file"),
  make_option(c("--outfile"), type="character", help="output table", metavar="file"),
  make_option(c("--function_file"), type="character", help="helper function file", metavar="file"),
  make_option(c("--bamfile_tumor"), type="character", help="tumor BAM file", metavar="bam"),
  make_option(c("--bamfile_control"), type="character", default=NA, help="control BAM file", metavar="bam"),
  make_option(c("--clipped_reads_control_file"), type="character", default=NA, help="control clipped reads table", metavar="file"),
  make_option(c("--min_site_support"), type="integer", default=3, help="minimum reads_supporting_insertion_pos to pass [default %default]", metavar="int"),
  make_option(c("--max_control_tel_ratio"), type="double", default=0.10, help="maximum control_telomeric_clip_ratio_all to pass [default %default]", metavar="float"),
  make_option(c("--max_control_tel_reads"), type="integer", default=4, help="maximum control_telomeric_clipped_reads_at_site to pass [default %default]", metavar="int")
)

opt_parser = OptionParser(option_list=option_list)
args = parse_args(opt_parser)

candidate_region_file = args$candidate_region_file
clipped_reads_file = args$clipped_reads_file
discordant_read_file = args$discordant_read_file
outfile = args$outfile
function_file = args$function_file
bamfile_tumor = args$bamfile_tumor
bamfile_control = args$bamfile_control
min_site_support = args$min_site_support
max_control_tel_ratio = args$max_control_tel_ratio
max_control_tel_reads = args$max_control_tel_reads

# review-only thresholds (flag for manual review but do not affect pass/fail) -- edit
# here directly rather than via Snakemake config, since these are secondary signals
# rather than the primary pass/fail decision
REVIEW_CONTROL_TEL_READS = 2
REVIEW_CONTROL_TEL_RATIO = 0.02
REVIEW_CONTROL_CLIP_RATIO = 0.30
REVIEW_TUMOR_CLIP_RATIO = 0.30
REVIEW_TUMOR_TELCLIP_FRACTION = 0.50
REVIEW_TUMOR_NONTEL_RATIO = 0.25
REVIEW_LOW_SITE_SUPPORT_MAX = 4
REVIEW_LOW_SITE_SUPPORT_NONTEL_RATIO = 0.20
clipped_reads_control_file = args$clipped_reads_control_file

count_control = !is.na(bamfile_control) && !is.na(clipped_reads_control_file)

source(function_file)

#don't use exponential notation of numbers (e.g. 600000 instead of 6e+05)
options(scipen = 999)

##########################################################################################################################

candidate_regions = read.table(candidate_region_file, header=TRUE, sep = "\t", stringsAsFactors=FALSE)
row.names(candidate_regions) = candidate_regions$window

clipped_reads_all = read.table(clipped_reads_file, header=TRUE, sep = "\t", stringsAsFactors=FALSE, comment.char='')
discordant_read_table = read.table(discordant_read_file, header=TRUE, sep = "\t", stringsAsFactors=FALSE, comment.char='')

if(count_control){
  clipped_reads_control_all = read.table(clipped_reads_control_file, header=TRUE, sep = "\t", stringsAsFactors=FALSE, comment.char='')
}

print(paste("candidate_regions rows:", dim(candidate_regions)[1]))
print(paste("clipped_reads_all rows:", dim(clipped_reads_all)[1]))
print(paste("discordant_read_table rows:", dim(discordant_read_table)[1]))
if(count_control){
  print(paste("clipped_reads_control_all rows:", dim(clipped_reads_control_all)[1]))
}

#--------------------------------------------------------------------------------------------------
# helper functions
#--------------------------------------------------------------------------------------------------

count_reads_at_site = function(bamfile, chrom, site_pos){
  if(is.na(site_pos)){
    return(NA)
  }
  view_cmd = paste0("samtools view -F 1024 ", bamfile, " ", chrom, ":", site_pos, "-", site_pos)
  sam_lines = system(view_cmd, intern=TRUE)
  return(length(sam_lines))
}

clip_interval_overlaps_site = function(start, end, site_pos){
  if(is.na(start) || is.na(end) || is.na(site_pos)){
    return(FALSE)
  }
  return(start <= site_pos && end >= site_pos)
}

ratio_or_na = function(num, den){
  if(is.na(num) || is.na(den) || den == 0){
    return(NA)
  }
  return(num / den)
}

#------------------------------------------------------------------------------
# helper functions for clipped-read site counting
#------------------------------------------------------------------------------

cigar_tokens = function(cigar){
  if(is.na(cigar) || cigar == ""){
    return(list(lengths=integer(0), ops=character(0)))
  }

  m = regmatches(cigar, gregexpr("[0-9]+[MIDNSHP=X]", cigar, perl=TRUE))[[1]]
  if(length(m) == 0){
    return(list(lengths=integer(0), ops=character(0)))
  }

  lens = as.integer(sub("([0-9]+).*", "\\1", m))
  ops = sub("[0-9]+", "", m)
  return(list(lengths=lens, ops=ops))
}

# Leading (left) and trailing (right) soft/hard-clip length from a CIGAR string.
# These are the clip lengths at the very start/end of the CIGAR, i.e. the parts of
# the read that fall *outside* the aligned [start, end] reference span.
clip_lengths = function(cigar){
  if(is.na(cigar) || cigar == ""){
    return(list(left=0L, right=0L))
  }

  tok = cigar_tokens(cigar)
  if(length(tok$ops) == 0){
    return(list(left=0L, right=0L))
  }

  left = 0L
  right = 0L

  if(tok$ops[1] %in% c("S", "H")){
    left = tok$lengths[1]
  }

  n = length(tok$ops)
  if(tok$ops[n] %in% c("S", "H")){
    right = tok$lengths[n]
  }

  return(list(left=left, right=right))
}

# A read is "clipped at site_pos" if the genomic footprint of one of its soft/hard
# clips overlaps the 1-bp site window. The clip footprint is not part of the aligned
# [start, end] span (soft/hard clips don't consume reference), so we place it
# immediately adjacent to that span:
#   left clip  -> [start - left_len,  start - 1]
#   right clip -> [end + 1,           end + right_len]
read_clipped_at_site = function(start, end, cigar, site_pos){
  if(is.na(start) || is.na(end) || is.na(site_pos) || is.na(cigar)){
    return(FALSE)
  }

  cl = clip_lengths(cigar)

  if(cl$left > 0){
    left_clip_start = start - cl$left
    left_clip_end = start - 1
    if(clip_interval_overlaps_site(left_clip_start, left_clip_end, site_pos)){
      return(TRUE)
    }
  }

  if(cl$right > 0){
    right_clip_start = end + 1
    right_clip_end = end + cl$right
    if(clip_interval_overlaps_site(right_clip_start, right_clip_end, site_pos)){
      return(TRUE)
    }
  }

  return(FALSE)
}

count_unique_clipped_reads_at_site = function(clipped_reads, site_pos){
  if(is.na(site_pos) || dim(clipped_reads)[1] == 0){
    return(data.frame(clipped_reads_at_site=NA, telomeric_clipped_reads_at_site=NA))
  }

  needed_cols = c("read_name", "start", "end", "cigar", "part_telomere")
  if(any(!(needed_cols %in% colnames(clipped_reads)))){
    return(data.frame(clipped_reads_at_site=NA, telomeric_clipped_reads_at_site=NA))
  }

  clipped_reads$start = as.numeric(clipped_reads$start)
  clipped_reads$end = as.numeric(clipped_reads$end)

  keep = sapply(seq_len(dim(clipped_reads)[1]), function(i){
    read_clipped_at_site(
      clipped_reads$start[i],
      clipped_reads$end[i],
      clipped_reads$cigar[i],
      site_pos
    )
  })

  clipped_at_site = clipped_reads[keep, ]

  if(dim(clipped_at_site)[1] == 0){
    return(data.frame(clipped_reads_at_site=0, telomeric_clipped_reads_at_site=0))
  }

  clipped_at_site = clipped_at_site[!is.na(clipped_at_site$read_name) & clipped_at_site$read_name != "", ]
  if(dim(clipped_at_site)[1] == 0){
    return(data.frame(clipped_reads_at_site=0, telomeric_clipped_reads_at_site=0))
  }

  clipped_at_site_unique = clipped_at_site[!duplicated(clipped_at_site$read_name), ]
  clipped_count = dim(clipped_at_site_unique)[1]
  telomeric_count = sum(!is.na(clipped_at_site_unique$part_telomere) & clipped_at_site_unique$part_telomere)

  return(data.frame(
    clipped_reads_at_site=clipped_count,
    telomeric_clipped_reads_at_site=telomeric_count
  ))
}

##########################################################################################################################

#--------------------------------------------------------------------------------------------------
# predict insertion site per candidate region
#--------------------------------------------------------------------------------------------------

for(window in unique(clipped_reads_all$window)){

  #-------------------------------------------------------------------------------------------------------------------
  # get most likely insertion site (where do the most clipped sequences start or end?)
  #-------------------------------------------------------------------------------------------------------------------

  clipped_reads = clipped_reads_all[clipped_reads_all$window==window, ]

  # defaults
  candidate_regions[window, "insertion_site"] = NA
  candidate_regions[window, "pos_telomeres_from_insertion"] = NA
  candidate_regions[window, "reads_supporting_insertion_pos"] = NA
  candidate_regions[window, "sum_TTAGGG_count"] = NA
  candidate_regions[window, "sum_CCCTAA_count"] = NA
  candidate_regions[window, "repeat_forward"] = NA
  candidate_regions[window, "ambiguous_insertion_site"] = FALSE

  if(sum(clipped_reads$part_telomere, na.rm=TRUE) == 0){
    next
  }

  clipped_reads_filtered = clipped_reads[clipped_reads$part_telomere, ]


  #-----------------------------------------------------------------------------------
  # discordant reads on plus strand => clipped reads should end at same position
  #
  # 5' ------ chromosome ------ telomere 3'
  #
  #-----------------------------------------------------------------------------------

  if (compareNA(candidate_regions[window, "strand"], "+")){
    clipped_expected_pos_fusion = "downstream"
    clipped_start_end = "end"
    site_offset = 1
  }


  #-----------------------------------------------------------------------------------
  # discordant reads on minus strand => clipped reads should start at same position
  #
  # 5' telomere ------ chromosome ------ 3'
  #
  #-----------------------------------------------------------------------------------

  if (compareNA(candidate_regions[window, "strand"], "-")){
    clipped_expected_pos_fusion = "upstream"
    clipped_start_end = "start"
    site_offset = -1
  }


  #-------------------------------------------------------------------------
  # only keep clipped reads that match orientation of discordant reads
  #-------------------------------------------------------------------------

  #attention: does not consider reads that are clipped on both ends (these are unlikely to be indicative of telomere insertion)
  clipped_reads_filtered_matching_discordant = clipped_reads_filtered[compareNA(clipped_reads_filtered$expected_pos_fusion, clipped_expected_pos_fusion), ]

  #-------------------------------------------------------------------------
  # only keep clipped reads that match position of discordant reads
  #-------------------------------------------------------------------------

  # currently taking median to prevent missing insertions where there is another discordant read nearby
  chrom = candidate_regions[window, "chrom"]
  window_start = candidate_regions[window, "chromStart"]
  window_end = candidate_regions[window, "chromEnd"]

  discordant_reads = discordant_read_table[discordant_read_table$mate_chr==chrom &
                                           discordant_read_table$mate_position>= window_start &
                                           discordant_read_table$mate_position<= window_end &
                                           discordant_read_table$mate_strand==candidate_regions[window, "strand"], ]

  discordant_read_pos_median = median(discordant_reads$mate_position) + 50   #plus 50 to get middle of read => also not ideal

  if (compareNA(candidate_regions[window, "strand"], "+")){
    clipped_reads_filtered_matching_discordant = clipped_reads_filtered_matching_discordant[clipped_reads_filtered_matching_discordant$end>discordant_read_pos_median, ]

  }else if (compareNA(candidate_regions[window, "strand"], "-")){
    clipped_reads_filtered_matching_discordant = clipped_reads_filtered_matching_discordant[clipped_reads_filtered_matching_discordant$start<discordant_read_pos_median, ]
  }

  #------------------------------------------------------
  # where do most reads start/end?
  #------------------------------------------------------

  table_pos_insertion = as.data.frame(table(clipped_reads_filtered_matching_discordant[, clipped_start_end]),
                                        stringsAsFactors=FALSE)

  if (dim(table_pos_insertion)[1] != 0){
    colnames(table_pos_insertion) = c("pos", "Freq")

    # go through insertion positions again and count number of unique cigars (this prevents that we only count reads that map at exactly the same position)
    for(pos in table_pos_insertion$pos){
      unique_cigars = length(unique(clipped_reads_filtered_matching_discordant[clipped_reads_filtered_matching_discordant[, clipped_start_end]==pos, "cigar"]))
      table_pos_insertion[table_pos_insertion$pos==pos, "unique_cigars"] = unique_cigars
    }

  } else{
    table_pos_insertion = data.frame(pos=NA, Freq=NA, unique_cigars=NA)
  }


  insertion_pos = table_pos_insertion[table_pos_insertion$unique_cigars==max(table_pos_insertion$unique_cigars), "pos"]

  # ties happen when two candidate positions have equal unique_cigars support (e.g. two
  # nearby real insertion sites collapsed into the same window). Rather than dropping the
  # window entirely, break the tie using the discordant-read median position already
  # computed above as an independent signal for where the breakpoint should be, and flag
  # the call for manual review instead of silently discarding it.
  if(length(insertion_pos) > 1){
    insertion_pos = insertion_pos[which.min(abs(as.numeric(insertion_pos) - discordant_read_pos_median))]
    candidate_regions[window, "ambiguous_insertion_site"] = TRUE
  }

  if(length(insertion_pos)==1){
    candidate_regions[window, "insertion_site"] = insertion_pos
    candidate_regions[window, "pos_telomeres_from_insertion"] = clipped_expected_pos_fusion
    candidate_regions[window, "reads_supporting_insertion_pos"] = max(table_pos_insertion$unique_cigars)

    #-----------------------------------------------------------------------------------
    # get total TTAGGG and CCCTAA counts in telomere fusion reads at insertion site
    # and determine most likely repeat on forward strand
    #-----------------------------------------------------------------------------------

    clipped_reads_filtered_at_insertion = clipped_reads_filtered_matching_discordant[clipped_reads_filtered_matching_discordant[, clipped_start_end]==insertion_pos,]

    sum_TTAGGG_count = sum(clipped_reads_filtered_at_insertion$TTAGGG_count)
    sum_CCCTAA_count = sum(clipped_reads_filtered_at_insertion$CCCTAA_count)

    candidate_regions[window, "sum_TTAGGG_count"] = sum_TTAGGG_count
    candidate_regions[window, "sum_CCCTAA_count"] = sum_CCCTAA_count

    if(sum_TTAGGG_count > sum_CCCTAA_count){
      candidate_regions[window, "repeat_forward"] = "TTAGGG"
    }else if (sum_CCCTAA_count > sum_TTAGGG_count){
      candidate_regions[window, "repeat_forward"] = "CCCTAA"
    }else{
      candidate_regions[window, "repeat_forward"] = NA
    }
  }

  #--------------------------------------------------------------------------------------------------
  # site-level read counts for tumor and control
  #--------------------------------------------------------------------------------------------------

  if(is.na(candidate_regions[window, "insertion_site"])){
    candidate_regions[window, "tumor_all_reads_at_site"] = NA
    candidate_regions[window, "tumor_clipped_reads_at_site"] = NA
    candidate_regions[window, "tumor_telomeric_clipped_reads_at_site"] = NA
    candidate_regions[window, "tumor_telomeric_clip_ratio_all"] = NA
    candidate_regions[window, "tumor_clipped_ratio_all"] = NA
    candidate_regions[window, "tumor_telomeric_clip_ratio_clipped"] = NA

    candidate_regions[window, "control_all_reads_at_site"] = NA
    candidate_regions[window, "control_clipped_reads_at_site"] = NA
    candidate_regions[window, "control_telomeric_clipped_reads_at_site"] = NA
    candidate_regions[window, "control_telomeric_clip_ratio_all"] = NA
    candidate_regions[window, "control_clipped_ratio_all"] = NA
    candidate_regions[window, "control_telomeric_clip_ratio_clipped"] = NA
    next
  }

  # Count at insertion_site +/- 1 depending on strand/direction.
  # The 1-bp site window is evaluated from the read's clip genomic footprint (see
  # read_clipped_at_site() above) rather than by walking the CIGAR's reference-consuming ops.
  site_pos = as.numeric(candidate_regions[window, "insertion_site"]) + site_offset
  candidate_regions[window, "site_pos_used"] = site_pos

  candidate_regions[window, "tumor_all_reads_at_site"] = count_reads_at_site(bamfile_tumor, chrom, site_pos)
  tumor_clipped_counts = count_unique_clipped_reads_at_site(clipped_reads, site_pos)
  candidate_regions[window, "tumor_clipped_reads_at_site"] = tumor_clipped_counts$clipped_reads_at_site
  candidate_regions[window, "tumor_telomeric_clipped_reads_at_site"] = tumor_clipped_counts$telomeric_clipped_reads_at_site
  candidate_regions[window, "tumor_telomeric_clip_ratio_all"] = ratio_or_na(
    candidate_regions[window, "tumor_telomeric_clipped_reads_at_site"],
    candidate_regions[window, "tumor_all_reads_at_site"]
  )
  candidate_regions[window, "tumor_clipped_ratio_all"] = ratio_or_na(
    candidate_regions[window, "tumor_clipped_reads_at_site"],
    candidate_regions[window, "tumor_all_reads_at_site"]
  )
  candidate_regions[window, "tumor_telomeric_clip_ratio_clipped"] = ratio_or_na(
    candidate_regions[window, "tumor_telomeric_clipped_reads_at_site"],
    candidate_regions[window, "tumor_clipped_reads_at_site"]
  )
  candidate_regions[window, "tumor_nontelomeric_ratio_all"] = ratio_or_na(
    candidate_regions[window, "tumor_clipped_reads_at_site"] - candidate_regions[window, "tumor_telomeric_clipped_reads_at_site"],
    candidate_regions[window, "tumor_all_reads_at_site"] - candidate_regions[window, "tumor_telomeric_clipped_reads_at_site"]
  )
  candidate_regions[window, "tumor_nontelomeric_clip_ratio_all_reads"] = ratio_or_na(
    candidate_regions[window, "tumor_clipped_reads_at_site"] - candidate_regions[window, "tumor_telomeric_clipped_reads_at_site"],
    candidate_regions[window, "tumor_all_reads_at_site"]
  )

  if(count_control){
    candidate_regions[window, "control_all_reads_at_site"] = count_reads_at_site(bamfile_control, chrom, site_pos)
    control_clipped_counts = count_unique_clipped_reads_at_site(clipped_reads_control_all, site_pos)
    candidate_regions[window, "control_clipped_reads_at_site"] = control_clipped_counts$clipped_reads_at_site
    candidate_regions[window, "control_telomeric_clipped_reads_at_site"] = control_clipped_counts$telomeric_clipped_reads_at_site
    candidate_regions[window, "control_telomeric_clip_ratio_all"] = ratio_or_na(
      candidate_regions[window, "control_telomeric_clipped_reads_at_site"],
      candidate_regions[window, "control_all_reads_at_site"]
    )
    candidate_regions[window, "control_clipped_ratio_all"] = ratio_or_na(
      candidate_regions[window, "control_clipped_reads_at_site"],
      candidate_regions[window, "control_all_reads_at_site"]
    )
    candidate_regions[window, "control_telomeric_clip_ratio_clipped"] = ratio_or_na(
      candidate_regions[window, "control_telomeric_clipped_reads_at_site"],
      candidate_regions[window, "control_clipped_reads_at_site"]
    )
  }else{
    candidate_regions[window, "control_all_reads_at_site"] = NA
    candidate_regions[window, "control_clipped_reads_at_site"] = NA
    candidate_regions[window, "control_telomeric_clipped_reads_at_site"] = NA
    candidate_regions[window, "control_telomeric_clip_ratio_all"] = NA
    candidate_regions[window, "control_clipped_ratio_all"] = NA
    candidate_regions[window, "control_telomeric_clip_ratio_clipped"] = NA
  }
}

#--------------------------------------------------------------------------------------------------
# final pass / review classification
#--------------------------------------------------------------------------------------------------

candidate_regions[, "passed"] = rep(FALSE, nrow(candidate_regions))
candidate_regions[, "flagged_for_review"] = rep(FALSE, nrow(candidate_regions))
candidate_regions[, "filter_reason"] = rep(NA, nrow(candidate_regions))
candidate_regions[, "flagged_reason"] = rep(NA, nrow(candidate_regions))

for(window in candidate_regions$window){
  insertion_site = candidate_regions[window, "insertion_site"]

  if(is.na(insertion_site)){
    candidate_regions[window, "passed"] = FALSE
    candidate_regions[window, "flagged_for_review"] = FALSE
    candidate_regions[window, "filter_reason"] = "no_insertion_site"
    next
  }

  reasons = c()
  review_reasons = c()
  review_hits = 0

  site_support = candidate_regions[window, "reads_supporting_insertion_pos"]
  control_tel_ratio = candidate_regions[window, "control_telomeric_clip_ratio_all"]
  control_tel_reads = candidate_regions[window, "control_telomeric_clipped_reads_at_site"]
  control_clip_ratio = candidate_regions[window, "control_clipped_ratio_all"]
  tumor_clip_ratio = candidate_regions[window, "tumor_clipped_ratio_all"]
  tumor_telclip_within_clipped = candidate_regions[window, "tumor_telomeric_clip_ratio_clipped"]
  tumor_nontel_ratio_all = candidate_regions[window, "tumor_nontelomeric_ratio_all"]
  tumor_nontelclip_ratio_all_reads = candidate_regions[window, "tumor_nontelomeric_clip_ratio_all_reads"]
  tumor_site_support_ratio_all = ratio_or_na(site_support, candidate_regions[window, "tumor_all_reads_at_site"])

  pass_ok = TRUE

  if(compareNA(candidate_regions[window, "blacklisted"], "yes")){
    pass_ok = FALSE
    reasons = c(reasons, "blacklisted")
  }

  if(is.na(site_support) || site_support < min_site_support){
    pass_ok = FALSE
    reasons = c(reasons, "insufficient_site_support")
  }

  if(is.na(control_tel_ratio) || control_tel_ratio >= max_control_tel_ratio){
    pass_ok = FALSE
    reasons = c(reasons, "high_control_telomeric_ratio")
  }

  if(is.na(control_tel_reads) || control_tel_reads > max_control_tel_reads){
    pass_ok = FALSE
    reasons = c(reasons, "high_control_telomeric_count")
  }

  if(!is.na(control_tel_reads) && control_tel_reads >= REVIEW_CONTROL_TEL_READS){
    review_hits = review_hits + 1
    review_reasons = c(review_reasons, "control_telomeric_count")
  }
  if(!is.na(control_tel_ratio) && control_tel_ratio > REVIEW_CONTROL_TEL_RATIO){
    review_hits = review_hits + 1
    review_reasons = c(review_reasons, "control_telomeric_ratio")
  }
  if(!is.na(control_clip_ratio) && control_clip_ratio > REVIEW_CONTROL_CLIP_RATIO){
    review_hits = review_hits + 1
    review_reasons = c(review_reasons, "control_clip_ratio")
  }
  if(!is.na(tumor_clip_ratio) && tumor_clip_ratio > REVIEW_TUMOR_CLIP_RATIO && !is.na(tumor_telclip_within_clipped) && tumor_telclip_within_clipped < REVIEW_TUMOR_TELCLIP_FRACTION){
    review_hits = review_hits + 1
    review_reasons = c(review_reasons, "low_tumor_telomeric_clip_fraction")
  }
  if(!is.na(tumor_nontel_ratio_all) && tumor_nontel_ratio_all > REVIEW_TUMOR_NONTEL_RATIO){
    review_hits = review_hits + 1
    review_reasons = c(review_reasons, "high_tumor_nontelomeric_ratio")
  }
  if(!is.na(site_support) && site_support >= min_site_support && site_support <= REVIEW_LOW_SITE_SUPPORT_MAX && !is.na(tumor_nontelclip_ratio_all_reads) && tumor_nontelclip_ratio_all_reads > REVIEW_LOW_SITE_SUPPORT_NONTEL_RATIO){
    review_hits = review_hits + 1
    review_reasons = c(review_reasons, "low_site_support_high_nontelomeric_clips")
  }
  if(!is.na(tumor_site_support_ratio_all) && tumor_site_support_ratio_all < 0.05){
    review_hits = review_hits + 1
    review_reasons = c(review_reasons, "low_tumor_site_support_ratio")
  }

  if(compareNA(candidate_regions[window, "ambiguous_insertion_site"], TRUE)){
    review_hits = review_hits + 1
    review_reasons = c(review_reasons, "ambiguous_insertion_site_tiebreak")
  }

  candidate_regions[window, "passed"] = pass_ok
  candidate_regions[window, "flagged_for_review"] = review_hits >= 1

  if(length(reasons) == 0){
    candidate_regions[window, "filter_reason"] = NA
  }else{
    candidate_regions[window, "filter_reason"] = paste(unique(reasons), collapse = ";")
  }

  if(length(review_reasons) == 0){
    candidate_regions[window, "flagged_reason"] = NA
  }else{
    candidate_regions[window, "flagged_reason"] = paste(unique(review_reasons), collapse = ";")
  }
}

if (dim(candidate_regions)[1]==0){
  candidate_regions = data.frame(
    PID=NA, window=NA, chrom=NA, chromStart=NA, chromEnd=NA, strand=NA,
    tumor_discordant_read_count=NA, control_discordant_read_count=NA, blacklisted=NA,
    insertion_site=NA, pos_telomeres_from_insertion=NA, reads_supporting_insertion_pos=NA,
    sum_TTAGGG_count=NA, sum_CCCTAA_count=NA, repeat_forward=NA, ambiguous_insertion_site=NA,
    tumor_all_reads_at_site=NA, tumor_clipped_reads_at_site=NA, tumor_telomeric_clipped_reads_at_site=NA,
    tumor_telomeric_clip_ratio_all=NA, tumor_clipped_ratio_all=NA, tumor_telomeric_clip_ratio_clipped=NA,
    tumor_nontelomeric_ratio_all=NA, tumor_nontelomeric_clip_ratio_all_reads=NA,
    control_all_reads_at_site=NA, control_clipped_reads_at_site=NA, control_telomeric_clipped_reads_at_site=NA,
    control_telomeric_clip_ratio_all=NA, control_clipped_ratio_all=NA, control_telomeric_clip_ratio_clipped=NA,
    passed=NA, flagged_for_review=NA, filter_reason=NA, flagged_reason=NA
  )[numeric(0), ]
}

write.table(candidate_regions, file=outfile, quote=FALSE, row.names = FALSE, sep="\t")