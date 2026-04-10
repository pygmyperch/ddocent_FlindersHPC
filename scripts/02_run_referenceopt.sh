#!/bin/bash
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --job-name=dd_refopt
#SBATCH --output=/scratch/user/brau0037/logs/%x.%j.out
#SBATCH --error=/scratch/user/brau0037/logs/%x.%j.err

set -euo pipefail

#################################
# ENVIRONMENT
#################################

module load parallel || true

export MAMBA_ROOT_PREFIX=/scratch/user/brau0037/local/miniforge3
source /scratch/user/brau0037/local/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/user/brau0037/local/miniforge3/envs/kelp_ddocent

echo "python: $(which python)"
echo "dDocent: $(which dDocent)"
echo "ReferenceOpt.sh: $(which ReferenceOpt.sh)"
echo "RefMapOpt.sh: $(which RefMapOpt.sh)"

#################################
# USER SETTINGS
#################################

WD="/scratch/user/brau0037/kelp/ddocent"
INDIR="${WD}/raw"
SUBSET_LIST="${WD}/subset_samples.txt"
OUTFINAL="${WD}/refopt"

DRYRUN="${DRYRUN:-0}"

# ReferenceOpt settings
ASSEMBLY_TYPE="PE"
NPROC="${SLURM_CPUS_PER_TASK}"

MIN_K1=3
MAX_K1=9
MIN_K2=3
MAX_K2=9

MIN_SIM=0.80
MAX_SIM=0.98
SIM_INC=0.02

#################################
# JOB WORKDIR (fast IO)
#################################

WORK="${BGFS:-${TEMP:-/tmp/${USER}/dd_refopt_${SLURM_JOB_ID:-$$}}}"
INWORK="$WORK/input_symlinks"
RUNWORK="$WORK/RefOpt"
LOGWORK="$WORK/log"
TMPWORK="$WORK/tmp"

mkdir -p "$WORK" "$INWORK" "$RUNWORK" "$LOGWORK" "$TMPWORK" \
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

if ! command -v ReferenceOpt.sh >/dev/null 2>&1; then
  echo "ERROR: ReferenceOpt.sh not found in PATH" >&2
  exit 1
fi

NSAMPLES=$(grep -v '^[[:space:]]*$' "$SUBSET_LIST" | wc -l)
echo "Samples listed in subset file: $NSAMPLES"

if [[ "$NSAMPLES" -eq 0 ]]; then
  echo "ERROR: subset sample list is empty" >&2
  exit 1
fi

#################################
# STAGE INPUT INTO FAST IO AREA
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
  echo "ERROR: some subset files were missing" >&2
  cp -f "$MISSING" "$OUTFINAL/" || true
  exit 1
fi

if [[ "$NSTAGED_F" -ne "$NSAMPLES" || "$NSTAGED_R" -ne "$NSAMPLES" ]]; then
  echo "ERROR: staged counts do not match expected sample count" >&2
  exit 1
fi

#################################
# CREATE REFERENCEOPT WORKDIR
#################################

echo "Creating RefOpt working directory"

find "$INWORK" -maxdepth 1 -type l -name '*.fq.gz' -print0 | \
while IFS= read -r -d '' f; do
  ln -s "$(realpath "$f")" "$RUNWORK/$(basename "$f")"
done

cd "$RUNWORK"

#################################
# INVENTORY BEFORE RUN
#################################

{
  echo -e "metric\tvalue"
  echo -e "samples_in_subset\t$NSAMPLES"
  echo -e "missing_entries\t$NMISSING"
  echo -e "staged_forward\t$NSTAGED_F"
  echo -e "staged_reverse\t$NSTAGED_R"
  echo -e "hostname\t$(hostname)"
  echo -e "date\t$(date '+%F %T')"
  echo -e "assembly_type\t$ASSEMBLY_TYPE"
  echo -e "nproc\t$NPROC"
  echo -e "min_k1\t$MIN_K1"
  echo -e "max_k1\t$MAX_K1"
  echo -e "min_k2\t$MIN_K2"
  echo -e "max_k2\t$MAX_K2"
  echo -e "min_sim\t$MIN_SIM"
  echo -e "max_sim\t$MAX_SIM"
  echo -e "sim_inc\t$SIM_INC"
} > "$LOGWORK/refopt_counts.tsv"

find "$RUNWORK" -maxdepth 1 -type l -printf '%f\t%l\n' | sort > "$LOGWORK/refopt_symlinks_manifest.tsv"
find -L "$RUNWORK" -maxdepth 1 -type f -name '*.fq.gz' -printf '%f\t%s\n' | sort > "$LOGWORK/refopt_file_sizes.tsv"
du -h "$WORK" > "$LOGWORK/work_disk_usage_start.txt"

#################################
# RUN REFERENCEOPT
#################################

CMD=(
  ReferenceOpt.sh
  "$MIN_K1" "$MAX_K1"
  "$MIN_K2" "$MAX_K2"
  "$ASSEMBLY_TYPE"
  "$NPROC"
  "$MIN_SIM" "$MAX_SIM" "$SIM_INC"
)

echo "ReferenceOpt command:"
printf ' %q' "${CMD[@]}"
echo

if [[ "$DRYRUN" -eq 1 ]]; then
  echo "DRYRUN enabled, not running ReferenceOpt.sh"
else
  (
    /usr/bin/time -v "${CMD[@]}"
  ) 2>&1 | tee "$LOGWORK/ReferenceOpt.log"
fi

#################################
# POST-RUN SUMMARIES
#################################

du -h "$WORK" > "$LOGWORK/work_disk_usage_end.txt"

{
  for f in kopt.data uniqseq.data dDocent_main.LOG dDocent.runs; do
    [[ -e "$f" ]] && printf "%s\tpresent\n" "$f"
  done
} > "$LOGWORK/key_outputs.tsv"

find "$RUNWORK" -maxdepth 1 \( -type f -o -type l \) | sed 's#^\./##' | sort > "$LOGWORK/runwork_listing.txt"

#################################
# COPY OUTPUTS BACK TO SCRATCH
#################################

echo "Copying ReferenceOpt outputs back to $OUTFINAL"

if [[ "$DRYRUN" -eq 1 ]]; then
  echo "DRYRUN copy skipped"
else
  [[ -f kopt.data ]] && cp -f kopt.data "$OUTFINAL/"
  [[ -f uniqseq.data ]] && cp -f uniqseq.data "$OUTFINAL/" || true
  rsync -av "$LOGWORK"/ "$OUTFINAL/logs"/
  [[ -f dDocent_main.LOG ]] && cp -f dDocent_main.LOG "$OUTFINAL/logs/" || true
  [[ -f dDocent.runs ]] && cp -f dDocent.runs "$OUTFINAL/logs/" || true
fi

#################################
# SUMMARY
#################################

echo "Done"
echo "Final outputs: $OUTFINAL"
echo "Main log:      $OUTFINAL/logs/ReferenceOpt.log"

if [[ -f "$OUTFINAL/kopt.data" ]]; then
  echo "Key result:    $OUTFINAL/kopt.data"
else
  echo "Key result:    kopt.data was not created"
fi

# Optional cleanup once confident:
# rm -rf "$WORK"
