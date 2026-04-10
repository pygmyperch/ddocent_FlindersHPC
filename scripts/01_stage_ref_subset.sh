#!/bin/bash
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=dd_ref_subset
#SBATCH --output=/scratch/user/brau0037/logs/%x.%j.out
#SBATCH --error=/scratch/user/brau0037/logs/%x.%j.err

set -euo pipefail

#################################
# USER SETTINGS
#################################

module load parallel || true

export MAMBA_ROOT_PREFIX=/scratch/user/brau0037/local/miniforge3
source /scratch/user/brau0037/local/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/user/brau0037/local/miniforge3/envs/kelp_ddocent

echo "conda: $(which conda)"
echo "python: $(which python)"
echo "dDocent: $(which dDocent)"
echo "ReferenceOpt.sh: $(which ReferenceOpt.sh)"
echo "RefMapOpt.sh: $(which RefMapOpt.sh)"

WD="/scratch/user/brau0037/kelp/ddocent"
INDIR="${WD}/raw"
SUBSET_LIST="${WD}/subset_samples.txt"
OUTFINAL="${WD}/ref_subset_stage"

# Dry run mode: DRYRUN=1 sbatch 01_stage_ref_subset.slurm
DRYRUN="${DRYRUN:-0}"

#################################
# JOB WORKDIR (fast IO)
#################################

WORK="${BGFS:-${TEMP:-/tmp/${USER}/dd_ref_subset_${SLURM_JOB_ID:-$$}}}"
INWORK="$WORK/input_symlinks"
REFWORK="$WORK/RefOpt"
LOGWORK="$WORK/log"
TMPWORK="$WORK/tmp"

mkdir -p "$WORK" "$INWORK" "$REFWORK" "$LOGWORK" "$TMPWORK" \
         "$OUTFINAL" "$OUTFINAL/logs"

echo "Working dir:   $WORK"
echo "Input dir:     $INDIR"
echo "Subset list:   $SUBSET_LIST"
echo "Output dir:    $OUTFINAL"
echo "DRYRUN:        $DRYRUN"

#################################
# VALIDATION
#################################

if [[ ! -d "$INDIR" ]]; then
  echo "ERROR: input directory not found: $INDIR" >&2
  exit 1
fi

if [[ ! -f "$SUBSET_LIST" ]]; then
  echo "ERROR: subset sample list not found: $SUBSET_LIST" >&2
  exit 1
fi

NSAMPLES=$(grep -v '^[[:space:]]*$' "$SUBSET_LIST" | wc -l)

echo "Samples listed in subset file: $NSAMPLES"

if [[ "$NSAMPLES" -eq 0 ]]; then
  echo "ERROR: subset sample list is empty: $SUBSET_LIST" >&2
  exit 1
fi

#################################
# STAGE SUBSET INPUTS INTO FAST IO AREA
#################################

echo "Staging subset symlinks into working directory"

MISSING="$LOGWORK/missing_files.tsv"
STAGED="$LOGWORK/staged_files.tsv"
: > "$MISSING"
: > "$STAGED"

while IFS= read -r sample; do
  [[ -z "$sample" ]] && continue

  f="${INDIR}/${sample}.F.fq.gz"
  r="${INDIR}/${sample}.R.fq.gz"

  missing_flag=0

  if [[ ! -e "$f" ]]; then
    echo -e "${sample}\tF\t$f" >> "$MISSING"
    missing_flag=1
  fi

  if [[ ! -e "$r" ]]; then
    echo -e "${sample}\tR\t$r" >> "$MISSING"
    missing_flag=1
  fi

  if [[ "$missing_flag" -eq 0 ]]; then
    ln -s "$(realpath "$f")" "$INWORK/$(basename "$f")"
    ln -s "$(realpath "$r")" "$INWORK/$(basename "$r")"
    echo -e "${sample}\t$(basename "$f")\t$(basename "$r")" >> "$STAGED"
  fi
done < "$SUBSET_LIST"

NMISSING=$(wc -l < "$MISSING")
NSTAGED_F=$(find "$INWORK" -maxdepth 1 -type l -name '*.F.fq.gz' | wc -l)
NSTAGED_R=$(find "$INWORK" -maxdepth 1 -type l -name '*.R.fq.gz' | wc -l)

echo "Missing entries: $NMISSING"
echo "Staged forward files: $NSTAGED_F"
echo "Staged reverse files: $NSTAGED_R"

if [[ "$NMISSING" -gt 0 ]]; then
  echo "ERROR: some subset files were missing, see $MISSING" >&2
  cp -f "$MISSING" "$OUTFINAL/" || true
  exit 1
fi

if [[ "$NSTAGED_F" -ne "$NSAMPLES" || "$NSTAGED_R" -ne "$NSAMPLES" ]]; then
  echo "ERROR: staged counts do not match expected sample count" >&2
  exit 1
fi

#################################
# CREATE DDOCENT REFERENCE WORKDIR
#################################

echo "Creating RefOpt working directory"

find "$INWORK" -maxdepth 1 -type l -name '*.fq.gz' -print0 | \
while IFS= read -r -d '' f; do
  ln -s "$(realpath "$f")" "$REFWORK/$(basename "$f")"
done

cd "$REFWORK"

#################################
# MANIFESTS AND SUMMARIES
#################################

echo "Writing manifests and summaries"

{
  echo -e "metric\tvalue"
  echo -e "samples_in_subset\t$NSAMPLES"
  echo -e "missing_entries\t$NMISSING"
  echo -e "staged_forward\t$NSTAGED_F"
  echo -e "staged_reverse\t$NSTAGED_R"
} > "$LOGWORK/subset_counts.tsv"

find "$INWORK" -maxdepth 1 -type l -printf '%f\t%l\n' | sort > "$LOGWORK/input_symlinks_manifest.tsv"
find "$REFWORK" -maxdepth 1 -type l -printf '%f\t%l\n' | sort > "$LOGWORK/refopt_symlinks_manifest.tsv"

# Follow symlinks and report target file sizes
find -L "$REFWORK" -maxdepth 1 -type f -name '*.fq.gz' -printf '%f\t%s\n' | sort > "$LOGWORK/refopt_file_sizes.tsv"

# Paired summary table
{
  echo -e "sample\tF_size_bytes\tR_size_bytes\ttotal_size_bytes"
  for f in "$REFWORK"/*.F.fq.gz; do
    sample=$(basename "$f" .F.fq.gz)
    fsize=$(stat -Lc '%s' "$REFWORK/${sample}.F.fq.gz")
    rsize=$(stat -Lc '%s' "$REFWORK/${sample}.R.fq.gz")
    total=$((fsize + rsize))
    printf "%s\t%s\t%s\t%s\n" "$sample" "$fsize" "$rsize" "$total"
  done | sort -k4,4n
} > "$LOGWORK/refopt_pair_sizes.tsv"

du -h "$WORK" > "$LOGWORK/work_disk_usage.txt"

#################################
# COPY OUTPUTS BACK TO SCRATCH
#################################

echo "Copying subset-stage outputs back to $OUTFINAL"

if [[ "$DRYRUN" -eq 1 ]]; then
  echo "DRYRUN copy skipped"
else
  rsync -av "$LOGWORK"/ "$OUTFINAL"/
fi

#################################
# SUMMARY
#################################

echo "Done"
echo "Final outputs: $OUTFINAL"
echo "Counts file:   $OUTFINAL/subset_counts.tsv"
echo "Pair sizes:    $OUTFINAL/refopt_pair_sizes.tsv"

# Optional cleanup once confident:
# rm -rf "$WORK"
