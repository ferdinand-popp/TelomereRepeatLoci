# Author: Lina Sieverling

# Usage: R --no-save --slave --args <candidate_region_file> <clipped_reads_file> <discordant_read_file> <outfile> <function_file> <bamfile_tumor> [<bamfile_control> <clipped_reads_control_file>] < ...>
# Description: trys to predict a telomere insertion site for each candidate region from clipped reads of tumor sample
#              - takes the position where most clipped sequences start/end (if this is not unique it returns NA)
#              - add the result to the extended candidate region table
#
#              For each predicted site we additionally count, for tumor and control, the total reads at the site
#              and the subset of unique clipped reads whose clipped interval overlaps the 1-bp site window.
#              The clipped-read counts are deduplicated by read_name to avoid double counting soft-/hard-clipped
#              evidence for the same physical read.

# get commandline arguments
commandArgs = commandArgs()
candidate_region_file = commandArgs[5]
clipped_reads_file = commandArgs[6]
discordant_read_file = commandArgs[7]
outfile = commandArgs[8]
function_file = commandArgs[9]
bamfile_tumor = commandArgs[10]

# optional control BAM and control clipped-read file (2-sample mode only)
if(length(commandArgs) >= 13){
  bamfile_control = commandArgs[11]
  clipped_reads_control_file = commandArgs[12]
  count_control = TRUE
}else{
  count_control = FALSE
}

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

#--------------------------------------------------------------------------------------------------
# helper functions
#--------------------------------------------------------------------------------------------------

count_reads_at_site = function(bamfile, chrom, site_pos){
  if(is.na(site_pos)){
    return(NA)
  }
  view_cmd = paste0("samtools view ", bamfile, " ", chrom, ":", site_pos, "-", site_pos)
  sam_lines = system(view_cmd, intern=TRUE)
  return(length(sam_lines))
}

clip_interval_overlaps_site = function(start, end, site_pos){
  if(is.na(start) || is.na(end) || is.na(site_pos)){
    return(FALSE)
  }
  return(start <= site_pos && end >= site_pos)
}

count_unique_clipped_reads_at_site = function(clipped_reads, site_pos){
  if(is.na(site_pos) || dim(clipped_reads)[1] == 0){
    return(data.frame(clipped_reads_at_site=NA, telomeric_clipped_reads_at_site=NA))
  }

  if(!("read_name" %in% colnames(clipped_reads))){
    return(data.frame(clipped_reads_at_site=NA, telomeric_clipped_reads_at_site=NA))
  }

  # In find_fusion_reads.R:
  # - soft-clipped reads get genomic start/end from the alignment coordinates
  # - hard-clipped supplementary alignments are also represented with start/end
  # For both, a read overlaps the 1-bp site if start <= site_pos <= end.
  clipped_reads$start = as.numeric(clipped_reads$start)
  clipped_reads$end = as.numeric(clipped_reads$end)

  overlapping = clipped_reads[sapply(seq_len(dim(clipped_reads)[1]), function(i) clip_interval_overlaps_site(clipped_reads$start[i], clipped_reads$end[i], site_pos)), ]

  if(dim(overlapping)[1] == 0){
    return(data.frame(clipped_reads_at_site=0, telomeric_clipped_reads_at_site=0))
  }

  overlapping = overlapping[!is.na(overlapping$read_name) & overlapping$read_name != "", ]
  if(dim(overlapping)[1] == 0){
    return(data.frame(clipped_reads_at_site=0, telomeric_clipped_reads_at_site=0))
  }

  # Deduplicate by read_name so a single read only contributes once per site.
  overlapping_unique = overlapping[!duplicated(overlapping$read_name), ]
  clipped_count = dim(overlapping_unique)[1]
  telomeric_count = sum(!is.na(overlapping_unique$part_telomere) & overlapping_unique$part_telomere)

  return(data.frame(clipped_reads_at_site=clipped_count, telomeric_clipped_reads_at_site=telomeric_count))
}

ratio_or_na = function(num, den){
  if(is.na(num) || is.na(den) || den == 0){
    return(NA)
  }
  return(num / den)
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
  # The 1-bp site window is evaluated by overlap against clipped-read intervals.
  site_pos = as.numeric(candidate_regions[window, "insertion_site"]) + site_offset

  candidate_regions[window, "tumor_all_reads_at_site"] = count_reads_at_site(bamfile_tumor, chrom, site_pos)
  tumor_clipped_counts = count_unique_clipped_reads_at_site(clipped_reads, site_pos)
  candidate_regions[window, "tumor_clipped_reads_at_site"] = tumor_clipped_counts$clipped_reads_at_site
  candidate_regions[window, "tumor_telomeric_clipped_reads_at_site"] = tumor_clipped_counts$telomeric_clipped_reads_at_site
  candidate_regions[window, "tumor_telomeric_clip_ratio_all"] = ratio_or_na(candidate_regions[window, "tumor_telomeric_clipped_reads_at_site"], candidate_regions[window, "tumor_all_reads_at_site"])
  candidate_regions[window, "tumor_clipped_ratio_all"] = ratio_or_na(candidate_regions[window, "tumor_clipped_reads_at_site"], candidate_regions[window, "tumor_all_reads_at_site"])
  candidate_regions[window, "tumor_telomeric_clip_ratio_clipped"] = ratio_or_na(candidate_regions[window, "tumor_telomeric_clipped_reads_at_site"], candidate_regions[window, "tumor_clipped_reads_at_site"])

  if(count_control){
    candidate_regions[window, "control_all_reads_at_site"] = count_reads_at_site(bamfile_control, chrom, site_pos)
    control_clipped_counts = count_unique_clipped_reads_at_site(clipped_reads_control_all, site_pos)
    candidate_regions[window, "control_clipped_reads_at_site"] = control_clipped_counts$clipped_reads_at_site
    candidate_regions[window, "control_telomeric_clipped_reads_at_site"] = control_clipped_counts$telomeric_clipped_reads_at_site
    candidate_regions[window, "control_telomeric_clip_ratio_all"] = ratio_or_na(candidate_regions[window, "control_telomeric_clipped_reads_at_site"], candidate_regions[window, "control_all_reads_at_site"])
    candidate_regions[window, "control_clipped_ratio_all"] = ratio_or_na(candidate_regions[window, "control_clipped_reads_at_site"], candidate_regions[window, "control_all_reads_at_site"])
    candidate_regions[window, "control_telomeric_clip_ratio_clipped"] = ratio_or_na(candidate_regions[window, "control_telomeric_clipped_reads_at_site"], candidate_regions[window, "control_clipped_reads_at_site"])
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

candidate_regions[, "passed"] = FALSE
candidate_regions[, "flagged_for_review"] = FALSE
candidate_regions[, "filter_reason"] = NA

for(window in candidate_regions$window){
  insertion_site = candidate_regions[window, "insertion_site"]

  if(is.na(insertion_site)){
    candidate_regions[window, "passed"] = FALSE
    candidate_regions[window, "flagged_for_review"] = FALSE
    candidate_regions[window, "filter_reason"] = "no_insertion_site"
    next
  }

  reasons = c()
  review_hits = 0

  site_support = candidate_regions[window, "reads_supporting_insertion_pos"]
  control_tel_ratio = candidate_regions[window, "control_telomeric_clip_ratio_all"]
  control_tel_reads = candidate_regions[window, "control_telomeric_clipped_reads_at_site"]
  control_clip_ratio = candidate_regions[window, "control_clipped_ratio_all"]
  tumor_clip_ratio = candidate_regions[window, "tumor_clipped_ratio_all"]
  tumor_telclip_within_clipped = candidate_regions[window, "tumor_telomeric_clip_ratio_clipped"]

  pass_ok = TRUE

  if(is.na(site_support) || site_support < 3){
    pass_ok = FALSE
    reasons = c(reasons, "insufficient_site_support")
  }

  if(is.na(control_tel_ratio) || control_tel_ratio >= 0.10){
    pass_ok = FALSE
    reasons = c(reasons, "high_control_telomeric_ratio")
  }

  if(is.na(control_tel_reads) || control_tel_reads > 4){
    pass_ok = FALSE
    reasons = c(reasons, "high_control_telomeric_count")
  }

  if(!is.na(control_tel_reads) && control_tel_reads >= 2){
    review_hits = review_hits + 1
  }
  if(!is.na(control_tel_ratio) && control_tel_ratio > 0.02){
    review_hits = review_hits + 1
  }
  if(!is.na(control_clip_ratio) && control_clip_ratio > 0.30){
    review_hits = review_hits + 1
  }
  if(!is.na(tumor_clip_ratio) && tumor_clip_ratio > 0.30 && !is.na(tumor_telclip_within_clipped) && tumor_telclip_within_clipped < 0.50){
    review_hits = review_hits + 1
  }

  candidate_regions[window, "passed"] = pass_ok
  candidate_regions[window, "flagged_for_review"] = review_hits >= 2

  if(length(reasons) == 0){
    candidate_regions[window, "filter_reason"] = NA
  }else{
    candidate_regions[window, "filter_reason"] = paste(unique(reasons), collapse = ";")
  }
}

if (dim(candidate_regions)[1]==0){
  candidate_regions = data.frame(
    PID=NA, window=NA, chrom=NA, chromStart=NA, chromEnd=NA, strand=NA,
    tumor_discordant_read_count=NA, control_discordant_read_count=NA, blacklisted=NA,
    insertion_site=NA, pos_telomeres_from_insertion=NA, reads_supporting_insertion_pos=NA,
    sum_TTAGGG_count=NA, sum_CCCTAA_count=NA, repeat_forward=NA,
    tumor_all_reads_at_site=NA, tumor_clipped_reads_at_site=NA, tumor_telomeric_clipped_reads_at_site=NA,
    tumor_telomeric_clip_ratio_all=NA, tumor_clipped_ratio_all=NA, tumor_telomeric_clip_ratio_clipped=NA,
    control_all_reads_at_site=NA, control_clipped_reads_at_site=NA, control_telomeric_clipped_reads_at_site=NA,
    control_telomeric_clip_ratio_all=NA, control_clipped_ratio_all=NA, control_telomeric_clip_ratio_clipped=NA,
    passed=NA, flagged_for_review=NA, filter_reason=NA
  )[numeric(0), ]
}

write.table(candidate_regions, file=outfile, quote=FALSE, row.names = FALSE, sep="\t")