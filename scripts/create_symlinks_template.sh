#!/usr/bin/env bash
# =============================================================================
# create_symlinks_template.sh
#
# Create dDocent-format symlinks from demultiplexed read files.
#
# dDocent requires reads named as:
#   PopID_IndID.F.fq.gz  (forward)
#   PopID_IndID.R.fq.gz  (reverse)
#
# This script reads a tab-delimited sample map and creates symlinks with the
# correct dDocent names in a target directory.
#
# Usage:
#   bash create_symlinks_template.sh <source_dir> <target_dir> <sample_map>
#
# Arguments:
#   source_dir   Directory containing the original read files
#   target_dir   Directory where symlinks will be created (created if absent)
#   sample_map   Tab-delimited file with two columns:
#                  col 1: source base name (without extension)
#                  col 2: dDocent base name (without extension)
#
# The script appends forward/reverse suffixes to both source and target names.
# Set FWD_IN / REV_IN to match your source file naming, e.g.:
#   .F.fq.gz / .R.fq.gz   (dDocent default)
#   .1.fq.gz / .2.fq.gz   (Illumina R1/R2 style with number suffix)
#   _R1.fastq.gz / _R2.fastq.gz
#
# Sample map format (no header):
#   BBCS01    BBCS_01
#   CH01      CH_01
#   GAB01     GAB_01
#
# Example:
#   bash create_symlinks_template.sh \
#     /scratch/user/$USER/project/demultiplexed \
#     /scratch/user/$USER/project/raw \
#     sample_map.txt
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Configure source file suffixes to match your input naming convention
# -----------------------------------------------------------------------------
FWD_IN=".F.fq.gz"     # forward suffix on source files  — edit as needed
REV_IN=".R.fq.gz"     # reverse suffix on source files  — edit as needed

# dDocent output suffixes — do not change
FWD_OUT=".F.fq.gz"
REV_OUT=".R.fq.gz"

# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------
SOURCE_DIR="${1:?Usage: $0 <source_dir> <target_dir> <sample_map>}"
TARGET_DIR="${2:?Usage: $0 <source_dir> <target_dir> <sample_map>}"
SAMPLE_MAP="${3:?Usage: $0 <source_dir> <target_dir> <sample_map>}"

if [[ ! -d "$SOURCE_DIR" ]]; then
    echo "ERROR: source directory not found: $SOURCE_DIR" >&2
    exit 1
fi

if [[ ! -f "$SAMPLE_MAP" ]]; then
    echo "ERROR: sample map not found: $SAMPLE_MAP" >&2
    exit 1
fi

mkdir -p "$TARGET_DIR"

# -----------------------------------------------------------------------------
# Create symlinks
# -----------------------------------------------------------------------------
CREATED=0
MISSING=0
SKIPPED=0

while IFS=$'\t' read -r src_name ddocent_name; do
    # Skip blank lines and comments
    [[ -z "$src_name" || "$src_name" =~ ^# ]] && continue

    src_f="${SOURCE_DIR}/${src_name}${FWD_IN}"
    src_r="${SOURCE_DIR}/${src_name}${REV_IN}"
    dst_f="${TARGET_DIR}/${ddocent_name}${FWD_OUT}"
    dst_r="${TARGET_DIR}/${ddocent_name}${REV_OUT}"

    ok=1

    for src in "$src_f" "$src_r"; do
        if [[ ! -e "$src" ]]; then
            echo "MISSING: $src" >&2
            MISSING=$((MISSING + 1))
            ok=0
        fi
    done

    [[ "$ok" -eq 0 ]] && continue

    for pair in "$src_f:$dst_f" "$src_r:$dst_r"; do
        src="${pair%%:*}"
        dst="${pair##*:}"
        if [[ -e "$dst" ]]; then
            echo "SKIP (exists): $dst"
            SKIPPED=$((SKIPPED + 1))
        else
            ln -s "$(realpath "$src")" "$dst"
            CREATED=$((CREATED + 1))
        fi
    done

done < "$SAMPLE_MAP"

echo ""
echo "Done."
echo "  Symlinks created: $CREATED"
echo "  Skipped (exist):  $SKIPPED"
echo "  Missing source:   $MISSING"

if [[ "$MISSING" -gt 0 ]]; then
    echo "WARNING: $MISSING source files were not found — check paths and FWD_IN/REV_IN settings" >&2
    exit 1
fi
