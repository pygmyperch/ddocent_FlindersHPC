# dDocent RADseq Pipeline on deepthought — User Manual

## Overview

This manual describes how to run the dDocent de novo RADseq pipeline on the
Flinders University HPC (deepthought). It covers reference assembly parameter
optimisation, de novo reference building, read mapping, and SNP calling.

The pipeline is structured as a series of SLURM batch scripts. Each script handles
one stage, runs on the high-speed BeeGFS filesystem (BGFS), and copies its
essential outputs back to scratch before exiting.

### Two key constraints of deepthought

**1. BGFS is ephemeral.**
The high-speed BeeGFS filesystem (`/cluster/jobs/$USER/$SLURM_JOB_ID/`) is
allocated when a job starts and deleted when it ends. All computation must happen
there; only selected outputs are copied to scratch before the job exits.

**2. Scratch space is shared and limited.**
`/scratch/user/$USER/` is an NFS filesystem shared across all users. It is
frequently near capacity. Intermediate files should be deleted just-in-time to
free space only when you are ready to use it.

---

## Cluster reference

### Filesystems

| Filesystem | Path | Notes |
|------------|------|-------|
| Home | `/home/$USER/` | Small quota; not for data |
| Scratch | `/scratch/user/$USER/` | Persistent but shared and currently woefully inadequate |
| BGFS | `/cluster/jobs/$USER/$SLURM_JOB_ID/` | Fast, ephemeral — compute here |

### SLURM partitions

| Partition | Max walltime | Notes |
|-----------|-------------|-------|
| `melfu` | infinite | Reserved for MELFU lab, see Yuma for access |
| `high-capacity` | 7 days | QOS `hc-concurrent-jobs` required |
| `general` | 14 days | Use for longer mapping and SNP calling jobs |


### Conda / software environment

dDocent and its dependencies should be installed in a conda environment.
I recommend installing your own Miniforge on and using mamba to install dDocent,
rather than the outdated system Miniconda3 module.
The example below assumes you did this on scratch/local/

A typical activation block for SLURM scripts:

```bash
source /scratch/user/$USER/local/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/user/$USER/local/miniforge3/envs/<env_name>
```

---

## Before you start

### 1. Prepare your reads

dDocent requires reads to be named in the format:

```
PopID_IndID.F.fq.gz   (forward)
PopID_IndID.R.fq.gz   (reverse)
```

Where `PopID` is a population identifier and `IndID` is an individual identifier.
If your demultiplexed reads use a different naming convention, you can just create symlinks with
the correct names in a new dir without having to edit the original files.
An example bash script to create symlinks is in <path/to/create_symlinks.sh>

Place (or symlink) all reads into a single directory on scratch, e.g.:

```
/scratch/user/$USER/<project>/raw/
```

### 2. Prepare sample lists

Create plain-text files listing sample base names (without `.F.fq.gz`), one per line:

```
POP1_001
POP1_002
POP2_001
```

You will need at minimum:
- A small **pilot list** (4–6 samples) for parameter optimisation
- A **reference subset list** (representative geographic/population coverage,
  ~100–200 samples) for reference assembly

### 3. Set up your scratch directory

Recommended layout:

```
/scratch/user/$USER/<project>/
├── raw/                    # read symlinks (dDocent naming)
├── pilot_samples.txt       # pilot sample list
├── subset_samples.txt      # reference-building sample list
├── refopt/                 # Script 01 outputs
├── assembly/               # Script 02 outputs
├── mapping/                # Script 03 outputs
└── snps/                   # Script 04 outputs
```

### 4. Manage scratch space proactively

Check available space before each major stage:

```bash
df -h /scratch
```

Scratch space is shared. Do not clear large files until you are ready to use the
space — another user may fill it in the meantime. The main opportunities to free
space in this pipeline are:

- After reference assembly is confirmed: delete trimmed reads from any prior runs
- Before mapping: delete pre-demultiplexing raw sequencing files if no longer needed
- After the VCF is verified: delete BAMs

---

## BGFS execution model

Every script in this pipeline follows the same pattern:

1. Activate conda environment
2. Validate required inputs exist on scratch
3. Create a working directory on BGFS: `$BGFS/<stage>_$SLURM_JOB_ID/`
4. Symlink required inputs from scratch into the BGFS working directory
5. Run the analysis stage inside BGFS
6. Write logs, manifests, and summary files
7. Copy retained outputs back to scratch
8. (Optional) Delete BGFS working directory — uncomment cleanup line once confident

Inputs are symlinked rather than copied so that large read files are not duplicated.
The analysis tools follow the symlinks and read directly from scratch or wherever
the real files live.

---

## Standard log artifacts

Every script produces the following log files before it exits, regardless of
whether `DRYRUN` mode is enabled:

| File | Contents |
|------|----------|
| `process.log` | Timestamped start/end, commands run, exit codes, `/usr/bin/time -v` resource summary, `sacct` output |
| `bgfs_listing.txt` | Full BGFS inventory at job end (`find -maxdepth 3 -ls`) |
| `bgfs_du_start.txt` | BGFS directory sizes at job start |
| `bgfs_du_end.txt` | BGFS directory sizes at job end |
| `input_manifest.tsv` | Every staged input: sample name, symlink path, resolved target, file size |
| `output_manifest.tsv` | Every file copied to scratch: filename, size, md5sum |

All log files are copied to a `logs/` subdirectory of the stage output directory
on scratch (e.g. `refopt/logs/`). These logs are the permanent audit trail for
the job since the BGFS working directory is deleted when the job ends.

---

## Testing scripts before a full run

All scripts support a dry run mode that validates inputs and environment without
running the analysis:

```bash
sbatch --export=ALL,DRYRUN=1 scripts/01_referenceopt_pilot.sh
```

In dry run mode the script will:
- activate the conda environment and check all tools are in PATH
- validate input directories and sample lists
- stage symlinks into BGFS
- write all log artifacts and sync them to scratch
- skip the actual analysis command

Review the logs from the dry run before submitting the real job.

---

## Script 01 — `01_referenceopt_pilot.sh`

**Purpose:** Run `ReferenceOpt.sh` on a small pilot sample set (4–6 samples) to
survey the K1, K2, and similarity parameter space and identify a suitable range
for full reference assembly.

**Partition:** `high-capacity` | **CPUs:** 16 | **RAM:** 32 GB | **Walltime:** 2 days

### Why raw reads?

dDocent's de novo assembly requires **raw, untrimmed reads**. The `ReferenceOpt.sh`
script generates per-sample unique sequence counts directly from the raw reads.
Do not use trimmed reads for this step.

### Parameters to set in the script

| Variable | Example | Description |
|----------|---------|-------------|
| `RAWDIR` | `/scratch/user/$USER/<project>/raw` | Directory containing raw reads |
| `SUBSET_LIST` | `/scratch/user/$USER/<project>/pilot_samples.txt` | Pilot sample list |
| `OUTFINAL` | `/scratch/user/$USER/<project>/refopt` | Output directory on scratch |
| `MIN_K1` | 2 | Minimum per-individual unique sequence count cutoff |
| `MAX_K1` | 6 | Maximum per-individual unique sequence count cutoff |
| `MIN_K2` | 2 | Minimum number of individuals sharing a sequence |
| `MAX_K2` | 6 | Maximum (should not exceed sample count) |
| `MIN_SIM` | 0.80 | Minimum CD-HIT similarity threshold |
| `MAX_SIM` | 0.98 | Maximum CD-HIT similarity threshold |
| `SIM_INC` | 0.02 | Similarity step size |
| `NPROC` | 16 | Number of threads (match `--cpus-per-task`) |
| `ASSEMBLY_TYPE` | PE | Assembly type: PE, SE, OL, HYB, RPE, or ROL |

**Note on K2:** `MAX_K2` is automatically capped at the sample count. Setting it
higher than the number of pilot samples is harmless but produces `NA` results for
those combinations.

### Outputs copied to scratch

| File | Description |
|------|-------------|
| `refopt/kopt.data` | K1, K2, similarity, contig count for each parameter combination |
| `refopt/uniqseq.data` | Unique sequence counts |
| `refopt/logs/ReferenceOpt.log` | Full ReferenceOpt stdout including resource usage |
| `refopt/logs/sample_uniqseqs_counts.tsv` | Unique sequence count per sample |
| `refopt/logs/filter_seqs_counts.tsv` | Sequence counts surviving each K1/similarity filter |

### Interpreting kopt.data

`kopt.data` has four space-separated columns:

```
K1  K2  similarity  n_contigs
3   4   0.86        12453
3   4   0.88        11987
...
```

- **K1** — minimum number of times a sequence must appear within an individual
- **K2** — minimum number of individuals a sequence must appear in
- **similarity** — CD-HIT clustering threshold
- **n_contigs** — number of assembled reference contigs

Plot `n_contigs` across the parameter grid. Look for the region where contig count
stabilises (a plateau). Avoid the extremes: very low K values include noise; very
high values discard real loci. Rows with `NA` indicate combinations that produced
no sequences — these are expected when K2 exceeds the number of pilot samples or
when similarity is very high.

`sample_uniqseqs_counts.tsv` is useful for flagging samples with unusually low
unique sequence counts before scaling to a larger sample set.

`filter_seqs_counts.tsv` shows how many sequences survive each K1/similarity
filter, helping you understand why certain parameter combinations produce very few
or very many contigs.

### Decision point

Review `kopt.data` and choose K1, K2, and similarity (c) before running Script 02.
Record your chosen values and the reasoning in your project notes.

---

## Script 02 — `02_ddocent_assemble.sh` *(planned)*

**Purpose:** Run dDocent (trim + assemble only) on a representative subset of
samples to build the de novo reference.

**Partition:** `high-capacity` | **CPUs:** 32 | **RAM:** 128 GB | **Walltime:** 5 days

### Inputs

- Raw reads for the reference-building subset (e.g. 100–200 samples)
- K1, K2, and c values chosen from Script 01

### How dDocent is run non-interactively

dDocent does not have a command-line flag interface. It is driven by piping a
config file to stdin:

```bash
dDocent << 'EOF'
32          # number of processors
yes         # perform trimming
yes         # perform assembly
PE          # assembly type
<K1>
<K2>
<c>
60          # BWA match score (not used when mapping=no)
11          # BWA mismatch penalty
5,5         # BWA gap costs
no          # map reads
no          # call SNPs
user@institution.edu.au
EOF
```

dDocent automatically detects existing output files (e.g. `reference.fasta`) and
skips completed stages, which makes staged execution possible.

### Key outputs to scratch

| File | Description |
|------|-------------|
| `assembly/reference.fasta` | De novo reference — input for all downstream scripts |
| `assembly/uniqseq.data` | Unique sequence counts |
| `assembly/logs/` | dDocent logs and manifests |

---

## Script 03 — `03_ddocent_map.sh` *(planned)*

**Purpose:** Trim and map all samples to the assembled reference.

**Partition:** `general` | **CPUs:** 32 | **RAM:** 256 GB | **Walltime:** 14 days

### Before submitting

Ensure sufficient scratch space is available for BAM output (~200–300 MB per
sample). Free space just-in-time by deleting files you no longer need (e.g.
pre-demultiplexing archives, trimmed reads from earlier runs). Do not clear space
too early as the shared filesystem fills quickly.

### Inputs

- All raw read pairs in dDocent naming format
- `reference.fasta` from Script 02

### Key outputs to scratch

| File | Description |
|------|-------------|
| `mapping/*.bam` | Aligned reads, one BAM per sample |
| `mapping/*.bam.bai` | BAM indexes |
| `mapping/flagstats/` | Per-sample alignment summary statistics |

---

## Script 04 — `04_ddocent_snpcall.sh` *(planned)*

**Purpose:** Call variants from the mapped BAMs to produce a raw SNP dataset.

**Partition:** `general` | **CPUs:** 32 | **RAM:** 128 GB | **Walltime:** 3–5 days

### Inputs

- BAMs and `.bai` files from Script 03
- `reference.fasta` from Script 02

### Key outputs to scratch

| File | Description |
|------|-------------|
| `snps/TotalRawSNPs.vcf.gz` | Raw variant calls across all samples |
| `snps/TotalRawSNPs.vcf.gz.tbi` | Tabix index |

After verifying the VCF, the BAMs can be deleted to free scratch space.

---

## Scratch space planning

Estimate your space requirements before each stage and check `df -h /scratch`
regularly. Rough estimates for a ~2000-sample project:

| Stage | BGFS peak | Scratch output |
|-------|-----------|----------------|
| Script 01 (pilot) | ~1–5 GB | < 1 MB |
| Script 02 (123 samples) | ~500–700 GB | ~50 MB |
| Script 03 (2000 samples) | ~1.5 TB | ~400–600 GB (BAMs) |
| Script 04 (SNP calling) | ~50–100 GB | ~5–10 GB (VCF) |

BGFS has a large shared pool (check with `df -h /cluster/jobs`) and no per-job
quota. Scratch is the bottleneck.

---

## Known issues with dDocent v2.9.5

The following bugs affect `ReferenceOpt.sh` and have been patched in the version
maintained in this repository (`kelp_ddocent/bin/ReferenceOpt.sh`).

### 1. fastp version check rejects fastp ≥ 1.0

The upstream version check uses `$FASTP1 -lt 2`, which incorrectly rejects
fastp 1.x (because 1 < 2 is true). The patch changes the logic so that only
`FASTP1 == 0` triggers a version check, correctly accepting any 1.x release.

### 2. Empty `totaluniqseq` causes crash

When K2 exceeds the sample count, or when similarity is very high, the filtered
sequence file can be empty. The upstream script passes this empty file to
downstream tools which then crash. The patch detects an empty file, writes `NA`
to `kopt.data` for that combination, and continues to the next combination.

### 3. Empty `LENGTH` / `MaxLen` causes arithmetic error

If `totaluniqseq` is empty or the FASTA produced from it contains no sequences,
the `LENGTH` and `MaxLen` variables become empty strings. The upstream script
then attempts `$(( $LENGTH * 3 / 4 ))` which fails with a shell error. The patch
validates both variables before use and skips the combination with `NA` if either
is invalid.

### 4. `NA` rows in `plot.kopt.data` break gnuplot

The upstream script extracts column 4 of `kopt.data` with `cut` and passes all
values including `NA` strings to gnuplot, which then errors. The patch uses
`awk '$4 ~ /^[0-9]+$/'` to pass only numeric rows.

### Deploying the patched ReferenceOpt.sh

```bash
rsync -av --perms kelp_ddocent/bin/ReferenceOpt.sh \
  $USER@deepthought.flinders.edu.au:/scratch/user/$USER/local/miniforge3/envs/<env_name>/bin/ReferenceOpt.sh

# On HPC:
chmod +x /scratch/user/$USER/local/miniforge3/envs/<env_name>/bin/ReferenceOpt.sh
```
