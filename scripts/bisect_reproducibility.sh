#!/usr/bin/env bash
# Reruns the telomere-repeat-loci pipeline at a series of checkpoint commits between the
# copilot/create-readme-and-update-setup divergence point and current main, to find which
# change broke reproducibility against an external cohort's expected results.
#
# Each checkpoint is checked out into its own git worktree (your working tree is untouched),
# run against the same input BAMs, and the resulting candidate-region table is diffed against
# both the previous checkpoint's table and (if provided) EXPECTED_TABLE.
#
# Usage:
#   TUMOR_BAM=/path/tumor.bam \
#   TEL_TUMOR_BAM=/path/tumor_filtered_intratelomeric.bam \
#   CONTROL_BAM=/path/control.bam \
#   TEL_CONTROL_BAM=/path/control_filtered_intratelomeric.bam \
#   REFERENCE_FASTA=/path/reference.fa \
#   EXPECTED_TABLE=/path/external_cohort_expected_candidate_regions.tsv \
#   ./scripts/bisect_reproducibility.sh
#
# All variables except TUMOR_BAM and TEL_TUMOR_BAM are optional.

set -uo pipefail

TUMOR_BAM="${TUMOR_BAM:?set TUMOR_BAM to the tumor BAM/CRAM path}"
TEL_TUMOR_BAM="${TEL_TUMOR_BAM:?set TEL_TUMOR_BAM to the filtered intratelomeric tumor BAM path}"
CONTROL_BAM="${CONTROL_BAM:-}"
TEL_CONTROL_BAM="${TEL_CONTROL_BAM:-}"
REFERENCE_FASTA="${REFERENCE_FASTA:-}"
EXTRA_ARGS="${EXTRA_ARGS:---skip-visualization}"
EXPECTED_TABLE="${EXPECTED_TABLE:-}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKTREE_ROOT="${WORKTREE_ROOT:-$REPO_ROOT/../telomere-bisect-worktrees}"
RESULTS_ROOT="${RESULTS_ROOT:-$REPO_ROOT/../telomere-bisect-results}"

mkdir -p "$WORKTREE_ROOT" "$RESULTS_ROOT"

# Ordered checkpoints: hash:label. Each one is a commit where pipeline *behavior* changed
# (not docs/formatting/tests) between the R/Snakemake divergence point (2440157) and HEAD.
CHECKPOINTS=(
  "622023a:01-python-cli-entrypoint"
  "918cf36:02-candidate-fusion-rework"
  "60cf627:03-human-chrom-filter"
  "7f2857c:04-distinct-read-site-pick"
  "8ff2da0:05-confidence-diagnostics-added"
  "c95d04c:06-confidence-filter-optin"
  "f00842d:07-confidence-filter-always-on"
  "03829f4:08-shrunk-mate-lookup-window"
  "088113d:09-window-fusion-capped"
  "a715514:10-clustered-discordant-regions-HEAD"
)

log() { printf '\n=== %s ===\n' "$*"; }
warn() { printf 'WARN: %s\n' "$*" >&2; }

find_entrypoint() {
  # Echoes the pipeline entrypoint script path (relative to worktree root) for this checkpoint.
  local root="$1"
  if [[ -f "$root/src/telomererepeatloci/main.py" ]]; then
    echo "src/telomererepeatloci/main.py"
  elif [[ -f "$root/src/main.py" ]]; then
    echo "src/main.py"
  elif [[ -f "$root/run_telomere_repeat_loci.py" ]]; then
    echo "run_telomere_repeat_loci.py"
  else
    return 1
  fi
}

pick_table() {
  # Prefer the most-processed candidate-region table that exists at this checkpoint.
  local out_dir="$1" pattern f
  for pattern in \
    "*_candidate_regions_extended_with_confidence_filtered.tsv" \
    "*_candidate_regions_extended_with_confidence.tsv" \
    "*_candidate_regions_extended_with_consensus.tsv"; do
    f=$(find "$out_dir" -name "$pattern" 2>/dev/null | head -n1)
    if [[ -n "$f" ]]; then
      echo "$f"
      return 0
    fi
  done
  return 1
}

prev_table=""
prev_label=""

for entry in "${CHECKPOINTS[@]}"; do
  hash="${entry%%:*}"
  label="${entry##*:}"
  wt_dir="$WORKTREE_ROOT/$label"
  out_dir="$RESULTS_ROOT/$label/pipeline_out"

  log "Checkpoint $label ($hash)"

  if [[ ! -d "$wt_dir" ]]; then
    if ! git -C "$REPO_ROOT" worktree add -f --detach "$wt_dir" "$hash"; then
      warn "worktree add failed for $hash, skipping"
      continue
    fi
  fi

  entrypoint=$(find_entrypoint "$wt_dir") || { warn "no known entrypoint found at $label, skipping"; continue; }

  (
    cd "$wt_dir" || exit 1
    if ! uv sync --quiet; then
      warn "uv sync failed at $label, trying editable pip install fallback"
      uv pip install -e . --quiet || exit 1
    fi
  ) || { warn "environment setup failed at $label, skipping run"; continue; }

  help_output=$(cd "$wt_dir" && uv run python "$entrypoint" --help 2>&1)
  if ! grep -q -- "--tumor-bam" <<<"$help_output"; then
    warn "$label ($entrypoint) does not expose --tumor-bam — its CLI shape differs from current main (likely uses --telomerehunter-dir/--results-per-pid-dir/--bam-suffix instead)."
    warn "Skipping automated run for $label; inspect manually with: (cd $wt_dir && uv run python $entrypoint --help)"
    continue
  fi

  mkdir -p "$out_dir"
  cmd=(uv run python "$entrypoint" --tumor-bam "$TUMOR_BAM" --tel-tumor-bam "$TEL_TUMOR_BAM" --output-dir "$out_dir")
  [[ -n "$CONTROL_BAM" ]] && cmd+=(--control-bam "$CONTROL_BAM")
  [[ -n "$TEL_CONTROL_BAM" ]] && cmd+=(--tel-control-bam "$TEL_CONTROL_BAM")
  [[ -n "$REFERENCE_FASTA" ]] && cmd+=(--reference-fasta "$REFERENCE_FASTA")
  # shellcheck disable=SC2206
  extra=($EXTRA_ARGS)
  cmd+=("${extra[@]}")

  echo "+ (cd $wt_dir && ${cmd[*]})"
  if ! (cd "$wt_dir" && "${cmd[@]}"); then
    warn "pipeline run failed at $label, skipping diff"
    continue
  fi

  table=$(pick_table "$out_dir") || { warn "no candidate-region table produced at $label"; continue; }
  cp "$table" "$RESULTS_ROOT/$label/candidate_regions.tsv"
  rows=$(( $(wc -l < "$table") - 1 ))
  echo "-> $label: $rows candidate region(s) [$(basename "$table")]"

  if [[ -n "$prev_table" ]]; then
    changed=$(diff <(sort "$prev_table") <(sort "$table") | grep -c '^[<>]' || true)
    echo "-> diff vs $prev_label: $changed differing line(s)"
  fi

  if [[ -n "$EXPECTED_TABLE" ]]; then
    changed_ref=$(diff <(sort "$EXPECTED_TABLE") <(sort "$table") | grep -c '^[<>]' || true)
    echo "-> diff vs EXPECTED_TABLE: $changed_ref differing line(s)"
  fi

  prev_table="$RESULTS_ROOT/$label/candidate_regions.tsv"
  prev_label="$label"
done

log "Done."
echo "Per-checkpoint outputs: $RESULTS_ROOT/<label>/candidate_regions.tsv"
echo "Worktrees left in place at $WORKTREE_ROOT for manual follow-up. Clean up with:"
for entry in "${CHECKPOINTS[@]}"; do
  label="${entry##*:}"
  echo "  git -C \"$REPO_ROOT\" worktree remove \"$WORKTREE_ROOT/$label\" --force"
done
echo "  git -C \"$REPO_ROOT\" worktree prune"
