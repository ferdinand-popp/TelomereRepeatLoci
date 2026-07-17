# Author: Lina Sieverling
# usage: R --no-save --slave --args -t <discordantReadFileTumor> -c <discordantReadFileControl> -b <blacklist_file> -o <outFile> -f <function_file> < count_discordant_reads.R
# description: - uses the output tables of add_mate_mapq.py
#              - gets the number of discordant reads in 1 kb windows
#              - merges adjacent windows with discordant reads in the tumor sample
#              - annotates merged windows with ENCODE blacklist overlap fraction
#              - output: table for each pid containing raw discordant read counts for tumor and control

library("optparse")
library("data.table")

option_list <- list(
  make_option(c("-t", "--discordantReadFileTumor"), type = "character", default = NULL,
              help = "Discordant read file for tumor", metavar = "character"),
  make_option(c("-c", "--discordantReadFileControl"), type = "character", default = NULL,
              help = "Discordant read file for control", metavar = "character"),
  make_option(c("-b", "--blacklist_file"), type = "character", default = NULL,
              help = "Blacklist file (BED/BED.GZ or legacy window table)", metavar = "character"),
  make_option(c("-o", "--outFile"), type = "character", default = NULL,
              help = "Output file", metavar = "character"),
  make_option(c("-f", "--function_file"), type = "character",
              default = "/abi/data/sieverling/global_scripts/functions.R",
              help = "File with R functions [default= %default]", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

discordantReadFileTumor   <- opt$discordantReadFileTumor
discordantReadFileControl <- opt$discordantReadFileControl
blacklist_file            <- opt$blacklist_file
outFile                   <- opt$outFile
function_file             <- opt$function_file

source(function_file)
options(scipen = 999)

# ------------------------------------------------------------------------------------
# helpers for ENCODE blacklist (hg19-blacklist.v2.bed.gz)
# ------------------------------------------------------------------------------------

normalize_chrom <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x
}

read_blacklist_as_bed <- function(path) {
  if (is.null(path) || is.na(path) || path == "" || path == "no_file" || !file.exists(path)) {
    return(NULL)
  }

  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con))

  bl <- tryCatch(
    read.table(con, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "", quote = ""),
    error = function(e) NULL
  )

  if (is.null(bl) || nrow(bl) == 0 || ncol(bl) < 3) return(NULL)

  # BED-like check: col2/col3 numeric
  s2 <- suppressWarnings(as.numeric(bl[[2]]))
  s3 <- suppressWarnings(as.numeric(bl[[3]]))
  if (all(is.na(s2)) || all(is.na(s3))) return(NULL)

  out <- data.frame(
    chrom = normalize_chrom(bl[[1]]),
    start = s2,
    end   = s3,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$start) & !is.na(out$end) & out$end > out$start, ]
  if (nrow(out) == 0) return(NULL)
  out
}

read_blacklist_as_window_table <- function(path) {
  # legacy format with column "window"
  if (is.null(path) || is.na(path) || path == "" || path == "no_file" || !file.exists(path)) {
    return(NULL)
  }
  x <- tryCatch(
    read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = ""),
    error = function(e) NULL
  )
  if (is.null(x) || !("window" %in% colnames(x))) return(NULL)
  x
}

compute_blacklist_fraction <- function(chrom, start, end, bl_by_chr) {
  width <- end - start
  if (is.na(width) || width <= 0) return(0)

  chr_key <- as.character(chrom)
  if (!(chr_key %in% names(bl_by_chr))) return(0)

  bl <- bl_by_chr[[chr_key]]
  if (is.null(bl) || nrow(bl) == 0) return(0)

  ov_start <- pmax(start, bl$start)
  ov_end <- pmin(end, bl$end)
  overlap <- pmax(0, ov_end - ov_start)
  covered <- sum(overlap)   # safe for Boyle v2 (non-overlapping intervals)

  covered / width
}

# ------------------------------------------------------------------------------------
# get 1 kb windows with discordant read count
# ------------------------------------------------------------------------------------

windowCountList <- list()
discordantReadFileList <- list(
  tumor   = discordantReadFileTumor,
  control = discordantReadFileControl
)

for (sample in c("tumor", "control")) {

  discordantReadFile <- discordantReadFileList[[sample]]

  # check if file exists, if not: skip sample
  if (!file.exists(discordantReadFile)) {
    windowCount <- data.frame(window=NA, stringsAsFactors=FALSE)
    windowCount[, paste0(sample, '_discordant_read_count')] <- NA
    windowCountList[[sample]] <- windowCount[numeric(0), ]
    next
  }

  discordantReads_all <- read.table(discordantReadFile, header=TRUE, sep="\t", comment.char='', stringsAsFactors=FALSE)

  # only keep those with mate mapping quality > 30
  discordantReads <- subset(discordantReads_all, mate_mapq > 30)

  if (nrow(discordantReads) == 0) {
    windowCount <- data.frame(window=character(0), stringsAsFactors=FALSE)
    windowCount[, paste0(sample, '_discordant_read_count')] <- numeric(0)
    windowCountList[[sample]] <- windowCount
    next
  }

  # get 1kb window
  discordantReads$mate_position_1kb <- floor(discordantReads$mate_position / 1000) * 1000
  discordantReads$window <- paste0(discordantReads$mate_chr, "_", as.character(discordantReads$mate_position_1kb), "_", discordantReads$mate_strand)

  # get read count per window
  windowCount <- as.data.frame(table(discordantReads$window), stringsAsFactors=FALSE)
  colnames(windowCount) <- c("window", paste0(sample, "_discordant_read_count"))
  windowCountList[[sample]] <- windowCount
}

# merge tumor/control window counts and set missing to 0
windowCountMerged <- merge(windowCountList$tumor, windowCountList$control, by="window", all=TRUE)
windowCountMerged[is.na(windowCountMerged)] <- 0

# edge case: no windows at all -> write valid empty output and exit success
if (nrow(windowCountMerged) == 0) {
  emptyTable <- data.frame(
    PID = character(0),
    window = character(0),
    chrom = character(0),
    chromStart = numeric(0),
    chromEnd = numeric(0),
    strand = character(0),
    tumor_discordant_read_count = numeric(0),
    control_discordant_read_count = numeric(0),
    blacklist_overlap_frac = numeric(0),
    blacklist_excluded = logical(0),
    blacklisted = character(0),
    stringsAsFactors = FALSE
  )
  write.table(emptyTable, file=outFile, quote=FALSE, row.names = FALSE, sep="\t")
  quit(save="no", status=0)
}

# add chromosome/start/end/strand
windowCountMerged$chrom <- gsub("_.*", "", windowCountMerged$window)
windowCountMerged$chromStart <- as.numeric(gsub("_", "", regmatches(windowCountMerged$window, regexpr("_.*_", windowCountMerged$window))))
windowCountMerged$chromEnd <- windowCountMerged$chromStart + 1000
windowCountMerged$strand <- gsub(".*_", "", windowCountMerged$window)

# add PID and order columns
windowCountMerged$PID <- gsub("_discordant_reads_1_kb_windows.tsv", "", basename(outFile))
windowCountMerged <- windowCountMerged[, c(
  "PID", "window", "chrom", "chromStart", "chromEnd", "strand",
  "tumor_discordant_read_count", "control_discordant_read_count"
)]
row.names(windowCountMerged) <- windowCountMerged$window

# ------------------------------------------------------------------------------------
# merge adjacent 1 kb windows with discordant reads in tumor sample
# (kept as in original logic, now before blacklist annotation)
# ------------------------------------------------------------------------------------

for (window1 in windowCountMerged[windowCountMerged$tumor_discordant_read_count != 0, "window"]) {

  # skip if row removed during prior merge
  if (!window1 %in% windowCountMerged$window) next

  while (1) {

    chromEnd <- windowCountMerged[window1, "chromEnd"]
    chrom <- windowCountMerged[window1, "chrom"]
    strand <- windowCountMerged[window1, "strand"]

    # adjacent downstream window on same chrom and strand
    window2 <- windowCountMerged[
      windowCountMerged$chrom == chrom &
        windowCountMerged$chromStart == chromEnd &
        windowCountMerged$strand == strand, "window"
    ]

    if (identical(window2, character(0))) break

    # require both windows to have non-zero tumor support
    if (windowCountMerged[window1, "tumor_discordant_read_count"] == 0 ||
        windowCountMerged[window2, "tumor_discordant_read_count"] == 0) break

    # merge counts + extend interval
    windowCountMerged[window1, "chromEnd"] <- windowCountMerged[window2, "chromEnd"]
    windowCountMerged[window1, "tumor_discordant_read_count"] <-
      windowCountMerged[window1, "tumor_discordant_read_count"] + windowCountMerged[window2, "tumor_discordant_read_count"]
    windowCountMerged[window1, "control_discordant_read_count"] <-
      windowCountMerged[window1, "control_discordant_read_count"] + windowCountMerged[window2, "control_discordant_read_count"]

    # remove second row
    windowCountMerged <- windowCountMerged[windowCountMerged$window != window2, ]
  }
}

# ------------------------------------------------------------------------------------
# annotate blacklist overlap
# ------------------------------------------------------------------------------------

# defaults
# (use rep() rather than a bare scalar: assigning a scalar to a new column on a
# 0-row data.frame errors in base R - "replacement has 1 row, data has 0")
windowCountMerged$blacklist_overlap_frac <- rep(0, nrow(windowCountMerged))
windowCountMerged$blacklist_excluded <- rep(FALSE, nrow(windowCountMerged))
windowCountMerged$blacklisted <- rep("no", nrow(windowCountMerged))  # legacy-compatible column

# try ENCODE BED(.gz) first
blacklist_bed <- read_blacklist_as_bed(blacklist_file)

if (!is.null(blacklist_bed) && nrow(blacklist_bed) > 0) {

  cand_chr_norm <- normalize_chrom(windowCountMerged$chrom)
  blacklist_bed <- blacklist_bed[blacklist_bed$chrom %in% unique(cand_chr_norm), , drop=FALSE]

  if (nrow(blacklist_bed) > 0) {
    bl_by_chr <- split(blacklist_bed, blacklist_bed$chrom)

    frac <- mapply(
      compute_blacklist_fraction,
      chrom = cand_chr_norm,
      start = windowCountMerged$chromStart,
      end   = windowCountMerged$chromEnd,
      MoreArgs = list(bl_by_chr = bl_by_chr)
    )

    windowCountMerged$blacklist_overlap_frac <- as.numeric(frac)
    windowCountMerged$blacklist_excluded <- windowCountMerged$blacklist_overlap_frac >= 0.5
    windowCountMerged$blacklisted <- ifelse(windowCountMerged$blacklist_excluded, "yes", "no")
  }

} else {
  # fallback: legacy blacklist table containing "window"
  legacy_blacklist <- read_blacklist_as_window_table(blacklist_file)
  if (!is.null(legacy_blacklist)) {
    windowCountMerged$blacklisted <- ifelse(windowCountMerged$window %in% legacy_blacklist$window, "yes", "no")
    windowCountMerged$blacklist_excluded <- windowCountMerged$blacklisted == "yes"
    windowCountMerged$blacklist_overlap_frac <- ifelse(windowCountMerged$blacklist_excluded, 1, 0)
  } else {
    # no usable blacklist provided
    windowCountMerged$blacklisted <- rep(NA, nrow(windowCountMerged))
    windowCountMerged$blacklist_excluded <- rep(FALSE, nrow(windowCountMerged))
    windowCountMerged$blacklist_overlap_frac <- rep(0, nrow(windowCountMerged))
  }
}

# ------------------------------------------------------------------------------------
# save results
# ------------------------------------------------------------------------------------
write.table(windowCountMerged, file=outFile, quote=FALSE, row.names = FALSE, sep="\t")