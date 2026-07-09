# Author: Lina Sieverling

# Usage: R --no-save --slave --args <candidate_regions_file> <outfile1> <outfile2> <outfile3> <PID> < make_bed_for_visualization.R
# Description: makes bed files for the determined telomere insertion sites
# outfile3 is the subset of outfile2 (100bp zoomed-in window) restricted to sites that are
# both passed and flagged_for_review, so downstream plotting can additionally copy those PDFs
# into a review_passed folder.

# get commandline arguments
commandArgs = commandArgs()
candidate_region_file = commandArgs[5]
outfile1 = commandArgs[6]
outfile2 = commandArgs[7]
outfile3 = commandArgs[8]
pid = commandArgs[9]

#don't use exponential notation of numbers (e.g. 600000 instead of 6e+05)
options(scipen = 999)

candidate_regions = read.table(candidate_region_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Only visualize candidate regions that passed the new post-prediction filter.
# Rows with insertion_site = NA are already excluded by passed == TRUE.
if("passed" %in% colnames(candidate_regions)){
  candidate_regions = candidate_regions[candidate_regions$passed %in% c(TRUE, "TRUE", "True", 1), ]
}else{
  candidate_regions = candidate_regions[!is.na(candidate_regions$insertion_site), ]
  candidate_regions = candidate_regions[candidate_regions$reads_supporting_insertion_pos > 2, ]
}

##############################################################################
### make bed file with 500 bp surrounding telomere insertion on either side
##############################################################################
if (dim(candidate_regions)[1] != 0){
  candidate_regions$windowStart1 = candidate_regions$insertion_site - 500
  candidate_regions$windowStart1[is.na(candidate_regions$windowStart1)] = candidate_regions$chromStart[is.na(candidate_regions$windowStart1)]
  candidate_regions$windowStart1[candidate_regions$windowStart1 < 0] = 0
  
  candidate_regions$windowEnd1 = candidate_regions$insertion_site + 500
  candidate_regions$windowEnd1[is.na(candidate_regions$windowEnd1)] = candidate_regions$chromEnd[is.na(candidate_regions$windowEnd1)]
  
  bed_file = data.frame(chrom = candidate_regions$chrom,
                        chromStart = candidate_regions$windowStart1,
                        chromEnd = candidate_regions$windowEnd1,
                        pos = candidate_regions$insertion_site,
                        pid = pid,
                        stringsAsFactors = FALSE)
  
  bed_file[is.na(bed_file$pos), "pos"] = bed_file[is.na(bed_file$pos), "chromStart"] + ((bed_file[is.na(bed_file$pos), "chromEnd"]-bed_file[is.na(bed_file$pos), "chromStart"]) /2)
  
}else{
  bed_file = data.frame(chrom=NA, chromStart=NA, chromEnd=NA, pos=NA, pid=NA)[numeric(0), ]
}

colnames(bed_file)[colnames(bed_file)=="chrom"] = "#chrom"
# Sort by #chrom (as integer) and chromStart (numeric) before writing
bed_file = bed_file[order(as.integer(bed_file$`#chrom`), bed_file$chromStart), ]
write.table(bed_file, file=outfile1, quote=FALSE, row.names = FALSE, sep="\t") 





##############################################################################
### make bed file with 100 bp surrounding telomere insertion on either side
##############################################################################
if (dim(candidate_regions)[1] != 0){
  candidate_regions$windowStart2 = candidate_regions$insertion_site - 100
  candidate_regions$windowStart2[is.na(candidate_regions$windowStart2)] = candidate_regions$chromStart[is.na(candidate_regions$windowStart2)]
  
  candidate_regions$windowEnd2 = candidate_regions$insertion_site + 100
  candidate_regions$windowEnd2[is.na(candidate_regions$windowEnd2)] = candidate_regions$chromEnd[is.na(candidate_regions$windowEnd2)]
  
  bed_file = data.frame(chrom = candidate_regions$chrom,
                        chromStart = candidate_regions$windowStart2,
                        chromEnd = candidate_regions$windowEnd2,
                        pos = candidate_regions$insertion_site,
                        pid = pid,
                        stringsAsFactors = FALSE)
  
  bed_file[is.na(bed_file$pos), "pos"] = bed_file[is.na(bed_file$pos), "chromStart"] + ((bed_file[is.na(bed_file$pos), "chromEnd"]-bed_file[is.na(bed_file$pos), "chromStart"]) /2)

  if("flagged_for_review" %in% colnames(candidate_regions)){
    review_bed_file = bed_file[candidate_regions$flagged_for_review %in% c(TRUE, "TRUE", "True", 1), ]
  }else{
    review_bed_file = bed_file[numeric(0), ]
  }
}else{
  bed_file = data.frame(chrom=NA, chromStart=NA, chromEnd=NA, pos=NA, pid=NA)[numeric(0), ]
  review_bed_file = bed_file
}


colnames(bed_file)[colnames(bed_file)=="chrom"] = "#chrom"
# Sort by #chrom (as integer) and chromStart (numeric) before writing
bed_file = bed_file[order(as.integer(bed_file$`#chrom`), bed_file$chromStart), ]
write.table(bed_file, file=outfile2, quote=FALSE, row.names = FALSE, sep="\t")


##############################################################################
### make bed file (subset of the 100bp window bed) for sites that are both
### passed and flagged_for_review, used to additionally copy their PDFs into
### a review_passed folder
##############################################################################
colnames(review_bed_file)[colnames(review_bed_file)=="chrom"] = "#chrom"
review_bed_file = review_bed_file[order(as.integer(review_bed_file$`#chrom`), review_bed_file$chromStart), ]
write.table(review_bed_file, file=outfile3, quote=FALSE, row.names = FALSE, sep="\t")


#############################################################################################################
# > sessionInfo()
# R version 3.2.2 (2015-08-14)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: openSUSE 13.1 (Bottle) (x86_64)
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base  