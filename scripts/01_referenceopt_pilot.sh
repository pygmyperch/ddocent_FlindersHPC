#!/bin/bash
#SBATCH --partition=high-capacity
#SBATCH --qos=hc-concurrent-jobs
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=dd_refopt
#SBATCH --output=/scratch/user/brau0037/logs/%x.%j.out
#SBATCH --error=/scratch/user/brau0037/logs/%x.%j.err

# =============================================================================
# 01_referenceopt_pilot.sh
#
# Run ReferenceOpt.sh on the 6-sample pilot to find optimal K1, K2, similarity
# parameters for de novo RADseq assembly.
#
# Input:  raw reads at $RAWDIR/$sample.{F,R}.fq.gz (6 samples in $SUBSET_LIST)
# Output: kopt.data, uniqseq.data, logs -> $OUTFINAL/
#
# Usage:
#   sbatch 01_referenceopt_pilot.sh
#   DRYRUN=1 sbatch 01_referenceopt_pilot.sh   (validate without running)
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# ENVIRONMENT
# -----------------------------------------------------------------------------

source /scratch/user/brau0037/local/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/user/brau0037/local/miniforge3/envs/kelp_ddocent

echo "=== Tool versions ==="
echo "ReferenceOpt.sh: $(which ReferenceOpt.sh)"
echo "dDocent:         $(which dDocent)"
echo "fastp:           $(fastp --version 2>&1 | head -1)"
echo "cd-hit-est:      $(cd-hit-est -h 2>&1 | head -1 || true)"
echo "====="

# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------

RAWDIR="/scratch/user/brau0037/kelp/ddocent/raw"
SUBSET_LIST="/scratch/user/brau0037/kelp/ddocent/subset_samples.txt"
OUTFINAL="/scratch/user/brau0037/kelp/ddocent/refopt"

DRYRUN="${DRYRUN:-0}"

# ReferenceOpt parameter grid (K2 max capped at 6 = pilot sample count)
ASSEMBLY_TYPE="PE"
NPROC=16
MIN_K1=2
MAX_K1=6
MIN_K2=2
MAX_K2=6
MIN_SIM=0.80
MAX_SIM=0.98
SIM_INC=0.02

# -----------------------------------------------------------------------------
# DIRECTORIES
# -----------------------------------------------------------------------------

# BGFS is set by the cluster to /cluster/jobs/$USER/$SLURM_JOB_ID
if [[ -z "${BGFS:-}" ]]; then
    echo "ERROR: BGFS environment variable is not set. This script must run under SLURM on a BGFS-enabled node." >&2
    exit 1
fi

WORK="${BGFS%/}/dd_refopt_${SLURM_JOB_ID}"
RUNWORK="${WORK}/RefOpt"
LOGDIR="${WORK}/logs"

mkdir -p "$RUNWORK" "$LOGDIR" "$OUTFINAL" "$OUTFINAL/logs"

# Process log — opened now, appended throughout
PROCLOG="$LOGDIR/process.log"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$PROCLOG"; }

log "=== Job start ==="
log "Job ID:       ${SLURM_JOB_ID:-unknown}"
log "Hostname:     $(hostname)"
log "BGFS workdir: $WORK"
log "Raw reads:    $RAWDIR"
log "Subset list:  $SUBSET_LIST"
log "Output dir:   $OUTFINAL"
log "DRYRUN:       $DRYRUN"
log "Parameters:   K1=${MIN_K1}-${MAX_K1} K2=${MIN_K2}-${MAX_K2} sim=${MIN_SIM}-${MAX_SIM}+${SIM_INC} type=${ASSEMBLY_TYPE} nproc=${NPROC}"

# -----------------------------------------------------------------------------
# VALIDATION
# -----------------------------------------------------------------------------

log "--- Validation ---"

for tool in ReferenceOpt.sh mawk samtools rainbow gnuplot seqtk cd-hit-est parallel fastp pearRM; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        log "ERROR: required tool not in PATH: $tool"
        exit 1
    fi
done
log "All required tools found in PATH"

if [[ ! -d "$RAWDIR" ]]; then
    log "ERROR: raw read directory not found: $RAWDIR"
    exit 1
fi

if [[ ! -f "$SUBSET_LIST" ]]; then
    log "ERROR: subset sample list not found: $SUBSET_LIST"
    exit 1
fi

NSAMPLES=$(grep -cv '^[[:space:]]*$' "$SUBSET_LIST" || true)
log "Samples in subset list: $NSAMPLES"

if [[ "$NSAMPLES" -eq 0 ]]; then
    log "ERROR: subset sample list is empty"
    exit 1
fi

if [[ "$MAX_K2" -gt "$NSAMPLES" ]]; then
    log "WARNING: MAX_K2=$MAX_K2 exceeds sample count $NSAMPLES; capping MAX_K2 at $NSAMPLES"
    MAX_K2="$NSAMPLES"
fi

# -----------------------------------------------------------------------------
# STAGE INPUT READS (LOG ARTIFACT 3 — input manifest)
# -----------------------------------------------------------------------------

log "--- Staging reads into BGFS ---"

INPUT_MANIFEST="$LOGDIR/input_manifest.tsv"
echo -e "sample\tsymlink_F\tsymlink_R\ttarget_F\ttarget_R\tsize_F_bytes\tsize_R_bytes" > "$INPUT_MANIFEST"

MISSING=0
while IFS= read -r sample; do
    [[ -z "$sample" || "$sample" =~ ^# ]] && continue

    src_f="${RAWDIR}/${sample}.F.fq.gz"
    src_r="${RAWDIR}/${sample}.R.fq.gz"
    dst_f="${RUNWORK}/${sample}.F.fq.gz"
    dst_r="${RUNWORK}/${sample}.R.fq.gz"

    ok=1
    for f in "$src_f" "$src_r"; do
        if [[ ! -e "$f" ]]; then
            log "ERROR: missing input: $f"
            MISSING=$((MISSING + 1))
            ok=0
        fi
    done

    if [[ "$ok" -eq 1 ]]; then
        tgt_f=$(realpath "$src_f")
        tgt_r=$(realpath "$src_r")
        sz_f=$(stat -c '%s' "$tgt_f")
        sz_r=$(stat -c '%s' "$tgt_r")
        ln -sf "$tgt_f" "$dst_f"
        ln -sf "$tgt_r" "$dst_r"
        echo -e "${sample}\t${dst_f}\t${dst_r}\t${tgt_f}\t${tgt_r}\t${sz_f}\t${sz_r}" >> "$INPUT_MANIFEST"
    fi
done < "$SUBSET_LIST"

NSTAGED=$(grep -c '/' "$INPUT_MANIFEST" || true)  # count data rows
log "Staged $NSTAGED sample pairs (missing: $MISSING)"

if [[ "$MISSING" -gt 0 ]]; then
    log "ERROR: $MISSING input files missing; aborting"
    cp -f "$INPUT_MANIFEST" "$OUTFINAL/logs/" || true
    exit 1
fi

# Initial disk usage snapshot
du -sh "${WORK}"/*/  > "$LOGDIR/bgfs_du_start.txt" 2>/dev/null || du -sh "$WORK" > "$LOGDIR/bgfs_du_start.txt"
log "BGFS disk usage at start written to $LOGDIR/bgfs_du_start.txt"

# -----------------------------------------------------------------------------
# RUN REFERENCEOPT (LOG ARTIFACT 5 — process log, time -v)
# -----------------------------------------------------------------------------

log "--- Running ReferenceOpt ---"

CMD=(
    ReferenceOpt.sh
    "$MIN_K1" "$MAX_K1"
    "$MIN_K2" "$MAX_K2"
    "$ASSEMBLY_TYPE"
    "$NPROC"
    "$MIN_SIM" "$MAX_SIM" "$SIM_INC"
)

log "Command: ${CMD[*]}"

# ReferenceOpt.sh must run from the directory containing the reads
cd "$RUNWORK"

REFOPT_LOG="$LOGDIR/ReferenceOpt.log"
EXIT_CODE=0

if [[ "$DRYRUN" -eq 1 ]]; then
    log "DRYRUN enabled — skipping ReferenceOpt.sh execution"
else
    log "ReferenceOpt start: $(date '+%F %T')"
    /usr/bin/time -v "${CMD[@]}" 2>&1 | tee "$REFOPT_LOG" || EXIT_CODE=$?
    log "ReferenceOpt end: $(date '+%F %T')  exit_code=$EXIT_CODE"

    if [[ "$EXIT_CODE" -ne 0 ]]; then
        log "ERROR: ReferenceOpt.sh exited with code $EXIT_CODE"
        # Still proceed to collect diagnostics and copy what we have
    fi
fi

# -----------------------------------------------------------------------------
# POST-RUN DIAGNOSTICS
# -----------------------------------------------------------------------------

log "--- Post-run diagnostics ---"

# LOG ARTIFACT 1 — BGFS inventory
find "$WORK" -maxdepth 3 -ls > "$LOGDIR/bgfs_listing.txt" 2>/dev/null || true
log "BGFS inventory written to $LOGDIR/bgfs_listing.txt"

# LOG ARTIFACT 2 — disk usage snapshot (end)
du -sh "${WORK}"/*/  > "$LOGDIR/bgfs_du_end.txt" 2>/dev/null || du -sh "$WORK" > "$LOGDIR/bgfs_du_end.txt"
log "BGFS disk usage at end written to $LOGDIR/bgfs_du_end.txt"

# Summarise kopt.data if present
if [[ -f "$RUNWORK/kopt.data" ]]; then
    NROWS=$(wc -l < "$RUNWORK/kopt.data")
    log "kopt.data rows: $NROWS"
    echo "--- kopt.data head ---" >> "$PROCLOG"
    head -20 "$RUNWORK/kopt.data" >> "$PROCLOG" || true
else
    log "WARNING: kopt.data not found in $RUNWORK"
fi

# -----------------------------------------------------------------------------
# COPY OUTPUTS TO SCRATCH (LOG ARTIFACT 4 — output manifest with md5sums)
# -----------------------------------------------------------------------------

log "--- Copying outputs to scratch ---"

OUTPUT_MANIFEST="$LOGDIR/output_manifest.tsv"
echo -e "filename\tsize_bytes\tmd5sum\tdestination" > "$OUTPUT_MANIFEST"

copy_out() {
    local src="$1" dst_dir="$2" do_md5="${3:-0}"
    if [[ -f "$src" ]]; then
        local bname sz
        bname=$(basename "$src")
        sz=$(stat -c '%s' "$src")
        local md5="-"
        if [[ "$do_md5" -eq 1 ]]; then
            md5=$(md5sum "$src" | cut -d' ' -f1)
        fi
        if [[ "$DRYRUN" -eq 0 ]]; then
            cp -f "$src" "${dst_dir}/${bname}"
        fi
        echo -e "${bname}\t${sz}\t${md5}\t${dst_dir}/${bname}" >> "$OUTPUT_MANIFEST"
        log "  copied: $bname ($sz bytes)"
    fi
}

copy_out "$RUNWORK/kopt.data"        "$OUTFINAL"      1
copy_out "$RUNWORK/uniqseq.data"     "$OUTFINAL"      1
copy_out "$RUNWORK/dDocent_main.LOG" "$OUTFINAL/logs" 0
copy_out "$RUNWORK/dDocent.runs"     "$OUTFINAL/logs" 0

# Summary stats from intermediate sequence files (not copied — too large)
# Per-sample unique sequence counts
SAMPLE_SEQS_SUMMARY="$LOGDIR/sample_uniqseqs_counts.tsv"
echo -e "sample\tuniq_seq_count" > "$SAMPLE_SEQS_SUMMARY"
for f in "$RUNWORK"/*.uniq.seqs; do
    [[ -f "$f" ]] || continue
    sample=$(basename "$f" .uniq.seqs)
    count=$(wc -l < "$f")
    echo -e "${sample}\t${count}" >> "$SAMPLE_SEQS_SUMMARY"
done
log "Per-sample unique sequence counts written to $LOGDIR/sample_uniqseqs_counts.tsv"

# Sequences surviving each K1/similarity filter
FILTER_COUNTS="$LOGDIR/filter_seqs_counts.tsv"
echo -e "filename\tseq_count" > "$FILTER_COUNTS"
for f in "$RUNWORK"/uniq.k.*.seqs; do
    [[ -f "$f" ]] || continue
    count=$(wc -l < "$f")
    echo -e "$(basename "$f")\t${count}" >> "$FILTER_COUNTS"
done
log "Per-filter sequence counts written to $LOGDIR/filter_seqs_counts.tsv"

# Sync all log artifacts to scratch (always — logs are the output of a dry run too)
log "Syncing logs to $OUTFINAL/logs/"
rsync -a "$LOGDIR/" "$OUTFINAL/logs/"

# -----------------------------------------------------------------------------
# SACCT SUMMARY (LOG ARTIFACT 6)
# -----------------------------------------------------------------------------

log "--- sacct summary ---"
{
    echo "=== sacct for job ${SLURM_JOB_ID:-unknown} ==="
    sacct -j "${SLURM_JOB_ID:-0}" \
        --format=JobID,AllocCPUS,Elapsed,TotalCPU,MaxRSS,MaxVMSize,State,ExitCode \
        --noheader 2>/dev/null || echo "(sacct unavailable)"
} | tee -a "$PROCLOG"

# Final rsync to flush process.log with sacct content
rsync -a "$LOGDIR/" "$OUTFINAL/logs/"

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------

log "=== Job complete ==="
log "Key outputs:"
log "  kopt.data:    $OUTFINAL/kopt.data"
log "  uniqseq.data: $OUTFINAL/uniqseq.data"
log "  logs:         $OUTFINAL/logs/"
log "  process log:  $OUTFINAL/logs/process.log"

if [[ "$EXIT_CODE" -ne 0 ]]; then
    log "NOTICE: ReferenceOpt.sh exit code was $EXIT_CODE — review logs before proceeding"
fi

# Uncomment once confident in output integrity:
# rm -rf "$WORK"
