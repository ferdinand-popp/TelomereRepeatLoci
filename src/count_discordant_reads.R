# Author: Lina Sieverling
# (debug-instrumented version -- extra status/cat() prints added throughout)

# usage: R --no-save --slave --args -t <discordantReadFileTumor> -c <discordantReadFileControl> -b <blacklist_file> -o <outFile> -f <function_file> < count_discordant_reads.R
# description: - uses the output tables of add_mate_mapq.py
#              - gets the number of discordant reads in 1 kb windows
#              - merges adjacent windows with discordant reads in the tumor sample
#              - output: table for each pid containing raw discordant read counts for tumor and control


library("optparse")
library("data.table")

option_list <- list(
  make_option(c("-t", "--discordantReadFileTumor"), type = "character", default = NULL,
              help = "Discordant read file for tumor", metavar = "character"),
  make_option(c("-c", "--discordantReadFileControl"), type = "character", default = NULL,
              help = "Discordant read file for control", metavar = "character"),
  make_option(c("-b", "--blacklist_file"), type = "character", default = NULL,
              help = "Blacklist file", metavar = "character"),
  make_option(c("-o", "--outFile"), type = "character", default = NULL,
              help = "Output file", metavar = "character"),
  make_option(c("-f", "--function_file"), type = "character",
              default = "/abi/data/sieverling/global_scripts/functions.R",
              help = "File with R functions [default= %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

cat("opts:", deparse(opt), "\n")

# Extract arguments using the *exact* names defined in make_option
discordantReadFileTumor  <- opt$discordantReadFileTumor
discordantReadFileControl <- opt$discordantReadFileControl
blacklist_file           <- opt$blacklist_file   # <-- fixed here
outFile                  <- opt$outFile
function_file            <- opt$function_file

# Source helper functions
source(function_file)

# Disable scientific notation
options(scipen = 999)


#######################################################################################################################################################
cat("========== count_discordant_reads.R : DEBUG START ==========\n")
cat("discordantReadFileTumor:  ", discordantReadFileTumor, "\n")
cat("discordantReadFileControl:", discordantReadFileControl, "\n")
cat("blacklist_file:           ", blacklist_file, "\n")
cat("outFile:                  ", outFile, "\n")
cat("function_file:            ", function_file, "\n")
cat("==============================================================\n\n")

#------------------------------------------------------------------------------------
# get 1 kb windows with discordant read count
#------------------------------------------------------------------------------------

windowCountList = list()

discordantReadFileList = list(tumor=discordantReadFileTumor,
                            control=discordantReadFileControl)

for (sample in c("tumor", "control")){

  cat("---- processing sample:", sample, "----\n")

  # get discordant reads
  discordantReadFile = discordantReadFileList[[sample]]
  cat("  file path:", discordantReadFile, "\n")
  cat("  file exists:", file.exists(discordantReadFile), "\n")

  # check if file exists, if not: skip sample
  if(!file.exists(discordantReadFile)){
    cat("  -> SKIPPING sample", sample, ": file does not exist\n\n")
    windowCount = data.frame(window=NA)
    windowCount[, paste0(sample, '_discordant_read_count')] = NA
    windowCountList[[sample]] = windowCount[numeric(0), ]
    next
    }

  cat("  file size (bytes):", file.info(discordantReadFile)$size, "\n")

  discordantReads_all = read.table(discordantReadFile, header=TRUE, sep="\t", comment.char='')

  cat("  ncol(discordantReads_all):", ncol(discordantReads_all), "\n")
  cat("  colnames(discordantReads_all):", paste(colnames(discordantReads_all), collapse=", "), "\n")
  cat("  nrow(discordantReads_all):", nrow(discordantReads_all), "\n")
  cat("  class(mate_mapq):", class(discordantReads_all$mate_mapq), "\n")
  cat("  num NA in mate_mapq:", sum(is.na(discordantReads_all$mate_mapq)), "\n")
  cat("  head(discordantReads_all):\n")
  print(head(discordantReads_all))

  # only keep those with mate mapping quality > 30
  discordantReads = subset(discordantReads_all, mate_mapq>30)
  cat("  nrow after mate_mapq>30 filter:", nrow(discordantReads), "\n")

  if (nrow(discordantReads) == 0) {
    cat("  -> WARNING: 0 rows remain for sample", sample, "after mate_mapq>30 filter!\n")
  }

  # get 1kb window
  discordantReads$mate_position_1kb = floor(discordantReads$mate_position/1000) * 1000
  #discordantReads$window = paste0(discordantReads$mate_chr, '_', as.character(discordantReads$mate_position_1kb))
  discordantReads$window = paste0(discordantReads$mate_chr, '_', as.character(discordantReads$mate_position_1kb), '_', discordantReads$mate_strand)


  # get read count per window
  windowCount = as.data.frame(table(discordantReads$window), stringsAsFactors=FALSE)
  colnames(windowCount) = c('window', paste0(sample, '_discordant_read_count'))
  cat("  nrow(windowCount) for sample", sample, ":", nrow(windowCount), "\n")
  cat("  head(windowCount):\n")
  print(head(windowCount))
  cat("\n")

  #save in list
  windowCountList[[sample]] = windowCount
}

cat("---- after per-sample loop ----\n")
cat("nrow(windowCountList$tumor):  ", nrow(windowCountList$tumor), "\n")
cat("nrow(windowCountList$control):", nrow(windowCountList$control), "\n\n")

# merge window counts into 1 table and set missing values to 0
windowCountMerged = merge(windowCountList$tumor, windowCountList$control, by="window", all=TRUE)
windowCountMerged[is.na(windowCountMerged)] = 0

cat("nrow(windowCountMerged) after merge:", nrow(windowCountMerged), "\n")
cat("head(windowCountMerged) after merge:\n")
print(head(windowCountMerged))
cat("\n")

# -----------------------------------------------------------------------
# Handle the edge case where neither tumor nor control has any discordant
# read windows (e.g. mate_mapq > 30 filter removed everything, or one/both
# input files were empty). Rather than crashing when trying to assign PID
# to a 0-row data frame, write a valid empty output table and exit cleanly
# so this PID doesn't take down the whole Snakemake run.
# -----------------------------------------------------------------------
if (nrow(windowCountMerged) == 0) {
  cat("No discordant read windows found for this PID (tumor and control ",
      "both empty after the mate_mapq > 30 filter) -- writing an empty output table.\n")
  pid_value = gsub("_discordant_reads_1_kb_windows.tsv", "", basename(outFile))
  emptyTable = data.frame(
    PID = character(0),
    window = character(0),
    chrom = character(0),
    chromStart = numeric(0),
    chromEnd = numeric(0),
    strand = character(0),
    tumor_discordant_read_count = numeric(0),
    control_discordant_read_count = numeric(0),
    blacklisted = character(0),
    stringsAsFactors = FALSE
  )
  quit(status = 1)

  # write.table(emptyTable, file=outFile, quote=FALSE, row.names = FALSE, sep="\t")
  # cat("Wrote empty output table to:", outFile, "\n")
  # cat("========== count_discordant_reads.R : DEBUG END (empty exit) ==========\n")
}

# add chromosome, start and end coordinates
windowCountMerged$chrom = gsub("_.*", "", windowCountMerged$window)
# windowCountMerged$chromStart = as.numeric(gsub(".*_", "", windowCountMerged$window))
windowCountMerged$chromStart = as.numeric(gsub("_", "", regmatches(windowCountMerged$window,regexpr("_.*_",windowCountMerged$window))))
windowCountMerged$chromEnd = windowCountMerged$chromStart + 1000
windowCountMerged$strand = gsub(".*_", "", windowCountMerged$window)

#sort colnames and set rownames
windowCountMerged$PID = gsub("_discordant_reads_1_kb_windows.tsv", "", basename(outFile))
windowCountMerged = windowCountMerged[ , c("PID", "window", "chrom", "chromStart", "chromEnd", "strand", "tumor_discordant_read_count", "control_discordant_read_count")]
row.names(windowCountMerged) = windowCountMerged$window

cat("nrow(windowCountMerged) after adding chrom/coords/PID:", nrow(windowCountMerged), "\n\n")


#------------------------------------------------------------------------------------
# add blacklist
#------------------------------------------------------------------------------------

cat("---- blacklist step ----\n")
cat("blacklist_file:", blacklist_file, "\n")
cat("blacklist_file exists:", file.exists(blacklist_file), "\n")

if(file.exists(blacklist_file)){
  blacklist = read.table(blacklist_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
  cat("nrow(blacklist):", nrow(blacklist), "\n")
  windowCountMerged[windowCountMerged$window %in% blacklist$window, "blacklisted"] = "yes"
  windowCountMerged[! windowCountMerged$window %in% blacklist$window, "blacklisted"] = "no"
  cat("num windows flagged blacklisted=yes:", sum(windowCountMerged$blacklisted == "yes"), "\n")
}else{
  print("No blacklist file provided, or specified file does not exist. Continuing without blacklist")
  windowCountMerged$blacklisted = NA
}
cat("\n")


#------------------------------------------------------------------------------------
# merge adjacent 1 kb window with discordant reads in the tumor sample
#------------------------------------------------------------------------------------

cat("---- merging adjacent tumor windows ----\n")
cat("nrow(windowCountMerged) before adjacent-window merging:", nrow(windowCountMerged), "\n")
n_tumor_nonzero_windows = sum(windowCountMerged$tumor_discordant_read_count != 0)
cat("num windows with tumor_discordant_read_count != 0:", n_tumor_nonzero_windows, "\n\n")

for (window1 in windowCountMerged[windowCountMerged$tumor_discordant_read_count!=0, "window"]){

  #skip window if row has already been removed
  if (! window1 %in% windowCountMerged$window){next}

  while(1){

    chromEnd = windowCountMerged[window1, "chromEnd"]
    chrom = windowCountMerged[window1, "chrom"]
    strand = windowCountMerged[window1, "strand"]

    # get window where chromosome is the same and whose starting point is the same as the end point of the first window
    window2 = windowCountMerged[windowCountMerged$chrom==chrom &
                                  windowCountMerged$chromStart == chromEnd &
                                  windowCountMerged$strand == strand, "window"]

    if(identical(window2, character(0))){break}

    #skip if read counts in window 1 or 2 are zero in tumor sample
    if(windowCountMerged[window1, "tumor_discordant_read_count"]==0 | windowCountMerged[window2, "tumor_discordant_read_count"]==0){break}

    # add counts
    windowCountMerged[window1, "chromEnd"] = windowCountMerged[window2, "chromEnd"]
    windowCountMerged[window1, "tumor_discordant_read_count"] = windowCountMerged[window1, "tumor_discordant_read_count"] + windowCountMerged[window2, "tumor_discordant_read_count"]
    windowCountMerged[window1, "control_discordant_read_count"] = windowCountMerged[window1, "control_discordant_read_count"] + windowCountMerged[window2, "control_discordant_read_count"]

    # if any of the windows are blacklisted, merged window is also blacklisted

    if (is.na(windowCountMerged[window1, "blacklisted"]) && is.na(windowCountMerged[window2, "blacklisted"])){
      windowCountMerged[window1, "blacklisted"] = NA
    }else if(compareNA(windowCountMerged[window1, "blacklisted"],"yes") || compareNA(windowCountMerged[window2, "blacklisted"],"yes")){
      windowCountMerged[window1, "blacklisted"] = "yes"
    }else{
      windowCountMerged[window1, "blacklisted"] = "no"
    }


    # remove row with window 2
    windowCountMerged = windowCountMerged[windowCountMerged$window!=window2, ]

  }
}

cat("nrow(windowCountMerged) after adjacent-window merging:", nrow(windowCountMerged), "\n\n")

#------------------------------------------------------------------------------------
# save results
#------------------------------------------------------------------------------------
cat("---- writing output ----\n")
cat("final nrow(windowCountMerged):", nrow(windowCountMerged), "\n")
cat("head(windowCountMerged):\n")
print(head(windowCountMerged))
cat("writing to:", outFile, "\n")

write.table(windowCountMerged, file=outFile, quote=FALSE, row.names = FALSE, sep="\t")

cat("========== count_discordant_reads.R : DEBUG END ==========\n")

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
#
# other attached packages:
#   [1] optparse_1.3.2   data.table_1.9.6
#
# loaded via a namespace (and not attached):
#   [1] getopt_1.20.0 chron_2.3-47