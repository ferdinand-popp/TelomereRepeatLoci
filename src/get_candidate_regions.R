# Author: Lina Sieverling
# usage: R --no-save --slave --args <window_file> <candidate_region_file> <tumor_discordant_read_lower_limit> <control_discordant_read_upper_limit> <consider_blacklist> <function_file> < get_candidate_regions.R
# description: filters windows to get telomere insertion candidate regions:
#              - window must contain at least <tumor_discordant_read_lower_limit> discordant reads in tumor sample
#                and at most <control_discordant_read_upper_limit> in control sample
#              - applies ENCODE blacklist overlap filter (>=50% window coverage)
#              output: table with filtered results (= telomere insertion candidate regions)

# get commandline arguments
commandArgs = commandArgs()
window_file = commandArgs[5]
candidate_region_file = commandArgs[6]
tumor_discordant_read_lower_limit = as.numeric(commandArgs[7])
control_discordant_read_upper_limit = as.numeric(commandArgs[8])
consider_blacklist = commandArgs[9]
function_file = commandArgs[10]

source(function_file)

#don't use exponential notation of numbers (e.g. 600000 instead of 6e+05)
options(scipen = 999)

#########################################################################################################################

windowTable = read.table(window_file, header=TRUE, sep="\t", comment.char='', stringsAsFactors=FALSE)

n_windows_total = nrow(windowTable)
cat(sprintf("[get_candidate_regions] read %d windows from %s\n", n_windows_total, window_file))
cat(sprintf("[get_candidate_regions] thresholds: tumor_discordant_read_count >= %s, control_discordant_read_count <= %s\n",
            tumor_discordant_read_lower_limit, control_discordant_read_upper_limit))

#------------------------------------------------------------------------------------
# filter by tumor and control discordant read thresholds
#------------------------------------------------------------------------------------
pass_tumor = windowTable$tumor_discordant_read_count >= tumor_discordant_read_lower_limit
pass_control = windowTable$control_discordant_read_count <= control_discordant_read_upper_limit

cat(sprintf("[get_candidate_regions] %d/%d windows pass tumor threshold\n", sum(pass_tumor), n_windows_total))
cat(sprintf("[get_candidate_regions] %d/%d windows pass control threshold\n", sum(pass_control), n_windows_total))

candidate_regions = windowTable[pass_tumor & pass_control, ]

cat(sprintf("[get_candidate_regions] %d/%d windows pass both thresholds\n", nrow(candidate_regions), n_windows_total))

if (nrow(candidate_regions) == 0) {
  cat(sprintf("[get_candidate_regions] max tumor_discordant_read_count in input: %s\n",
              ifelse(n_windows_total > 0, max(windowTable$tumor_discordant_read_count, na.rm=TRUE), NA)))
}

# Ensure required output columns exist even if blacklist filtering is disabled
# (use rep() rather than a bare scalar: assigning a scalar to a new column on a
# 0-row data.frame errors in base R - "replacement has 1 row, data has 0")
candidate_regions$blacklist_overlap_frac = rep(0, nrow(candidate_regions))
candidate_regions$blacklist_excluded = rep(FALSE, nrow(candidate_regions))
if (!("exclusion_reason" %in% colnames(candidate_regions))) {
  candidate_regions$exclusion_reason = rep("", nrow(candidate_regions))
}

#------------------------------------------------------------------------------------
# filter by blacklist (>=50% window coverage, OR-combined style via blacklist_excluded)
#------------------------------------------------------------------------------------
if (consider_blacklist == "True" && ("blacklist_overlap_frac" %in% colnames(candidate_regions))) {
  candidate_regions$blacklist_overlap_frac[is.na(candidate_regions$blacklist_overlap_frac)] = 0
  candidate_regions$blacklist_excluded = candidate_regions$blacklist_overlap_frac >= 0.5

  # append exclusion reason where excluded
  idx_excl = which(candidate_regions$blacklist_excluded)
  if (length(idx_excl) > 0) {
    reason_add = paste0("blacklist>=50%(cov=", sprintf("%.2f", candidate_regions$blacklist_overlap_frac[idx_excl]), ")")
    existing = candidate_regions$exclusion_reason[idx_excl]
    existing[is.na(existing)] = ""
    candidate_regions$exclusion_reason[idx_excl] = ifelse(
      existing == "",
      reason_add,
      paste(existing, reason_add, sep=";")
    )
  }

  cat(sprintf("[get_candidate_regions] %d/%d candidates excluded by blacklist (>=50%% coverage)\n",
              length(idx_excl), nrow(candidate_regions)))

  # final candidate output excludes blacklisted candidates
  candidate_regions = candidate_regions[!candidate_regions$blacklist_excluded, ]
} else if (consider_blacklist == "True" && ("blacklisted" %in% colnames(candidate_regions))) {
  # Backward-compatible fallback: old yes/no blacklist column
  candidate_regions$blacklist_excluded = compareNA(candidate_regions$blacklisted, "yes")
  candidate_regions$blacklist_overlap_frac = ifelse(candidate_regions$blacklist_excluded, 1, 0)

  idx_excl = which(candidate_regions$blacklist_excluded)
  if (length(idx_excl) > 0) {
    reason_add = "blacklist>=50%(cov=1.00)"
    existing = candidate_regions$exclusion_reason[idx_excl]
    existing[is.na(existing)] = ""
    candidate_regions$exclusion_reason[idx_excl] = ifelse(
      existing == "",
      reason_add,
      paste(existing, reason_add, sep=";")
    )
  }

  cat(sprintf("[get_candidate_regions] %d/%d candidates excluded by blacklist (legacy yes/no column)\n",
              length(idx_excl), nrow(candidate_regions)))

  candidate_regions = candidate_regions[!candidate_regions$blacklist_excluded, ]
}

cat(sprintf("[get_candidate_regions] %d final candidate regions (from %d input windows)\n",
            nrow(candidate_regions), n_windows_total))

#------------------------------------------------------------------------------------
# save raw and filtered results
#------------------------------------------------------------------------------------
write.table(candidate_regions, candidate_region_file, quote=FALSE, row.names = FALSE, sep="\t")