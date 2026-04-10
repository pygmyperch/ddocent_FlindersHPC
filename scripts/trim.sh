#!/bin/bash
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=7-00:00:00
#SBATCH --job-name=trim_er
#SBATCH --output=/scratch/user/brau0037/logs/%x.%j.out
#SBATCH --error=/scratch/user/brau0037/logs/%x.%j.err

set -euo pipefail

#################################
# USER SETTINGS
#################################

module load Miniconda3
module load parallel || true

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate kelp_ddocent

PIPE="/scratch/user/brau0037/kelp/SNPcallingPipe/SNPcallPipe.pl"
INDIR="/scratch/user/brau0037/kelp/demux_all_symln"
OUTFINAL="/scratch/user/brau0037/kelp/Er_trim"

# Trim settings SNC=62
FM="n"
SNC=62
EXF=".F.fq.gz,.R.fq.gz"
RBAR="n"

# Optional barcode file only needed if RBAR="y"
# BARCODES="/scratch/user/brau0037/kelp/path/to/barcodes.txt"

# Dry run mode: DRYRUN=1 sbatch trim_er.slurm
DRYRUN="${DRYRUN:-0}"

#################################
# JOB WORKDIR (fast IO)
#################################

WORK="${BGFS:-${TEMP:-/tmp/${USER}/trim_${SLURM_JOB_ID:-$$}}}"
INWORK="$WORK/input_symlinks"
OUTWORK="$WORK/trim_out"
LOGWORK="$WORK/log"
TMPWORK="$WORK/tmp"

mkdir -p "$WORK" "$INWORK" "$OUTWORK" "$LOGWORK" "$TMPWORK" "$OUTFINAL" "$OUTFINAL/logs"

echo "Working dir: $WORK"
echo "Input dir:   $INDIR"
echo "Output dir:  $OUTFINAL"
echo "DRYRUN:      $DRYRUN"

#################################
# VALIDATION
#################################

if [[ ! -d "$INDIR" ]]; then
  echo "ERROR: input directory not found: $INDIR" >&2
  exit 1
fi

if [[ ! -f "$PIPE" ]]; then
  echo "ERROR: SNPcallPipe.pl not found: $PIPE" >&2
  exit 1
fi

# Count .fq.gz inputs, including regular files and symlinks
NINPUT=$(find "$INDIR" -maxdepth 1 \( -type f -o -type l \) -name '*.fq.gz' | wc -l)

echo "Detected .fq.gz inputs in $INDIR: $NINPUT"

if [[ "$NINPUT" -eq 0 ]]; then
  echo "ERROR: no .fq.gz files or symlinks found in input directory: $INDIR" >&2
  exit 1
fi

#################################
# STAGE INPUT INTO FAST IO AREA
#################################

echo "Staging input symlinks into working directory"

# Recreate symlinks in BGFS/TEMP so the pipeline works from the fast filesystem.
find "$INDIR" -maxdepth 1 \( -type f -o -type l \) -name '*.fq.gz' -print0 | \
while IFS= read -r -d '' f; do
  ln -s "$(realpath "$f")" "$INWORK/$(basename "$f")"
done

NSTAGED=$(find "$INWORK" -maxdepth 1 -type l -name '*.fq.gz' | wc -l)

echo "Staged files: $NSTAGED"

if [[ "$NSTAGED" -eq 0 ]]; then
  echo "ERROR: staging produced no .fq.gz symlinks in $INWORK" >&2
  exit 1
fi

#################################
# RUN TRIM STEP
#################################

CMD=(
  perl "$PIPE"
  -stp trim
  -i "$INWORK"
  -o "$OUTWORK"
  -fm "$FM"
  -snc "$SNC"
  -exf "$EXF"
  -rbar "$RBAR"
)

# If you later switch to barcode removal:
# if [[ "$RBAR" == "y" ]]; then
#   CMD+=(-bf "$BARCODES")
# fi

echo "Trim command:"
printf ' %q' "${CMD[@]}"
echo

if [[ "$DRYRUN" -eq 1 ]]; then
  echo "DRYRUN enabled, not running trim"
else
  (
    "${CMD[@]}"
  ) 2>&1 | tee "$LOGWORK/trim.log"
fi

#################################
# COPY OUTPUTS BACK TO SCRATCH
#################################

echo "Copying trimmed outputs back to $OUTFINAL"

if [[ "$DRYRUN" -eq 1 ]]; then
  echo "DRYRUN copy skipped"
else
  # rsync preserves structure and is safer than mv across filesystems
  rsync -av --info=progress2 "$OUTWORK"/ "$OUTFINAL"/
  cp -f "$LOGWORK/trim.log" "$OUTFINAL/logs/" 2>/dev/null || true
fi

#################################
# SUMMARY
#################################

echo "Done"
echo "Final outputs: $OUTFINAL"
echo "Main log:      $OUTFINAL/logs/trim.log"

# Optional cleanup once confident:
# rm -rf "$WORK"
