# dDocent on HPC with ephemeral BGFS

## Development workflow

Scripts are authored locally in this project directory and deployed to the HPC.
Claude Code reads/edits scripts locally; the user deploys and submits jobs.

HPC address: `brau0037@deepthought.flinders.edu.au`

### Deploy scripts to HPC

```bash
./deploy.sh
# or with dry-run preview:
./deploy.sh --dry-run
```

### Collect logs after a job finishes

```bash
# SLURM stdout/stderr
rsync -avz brau0037@deepthought.flinders.edu.au:/scratch/user/brau0037/logs/dd_refopt.*.{out,err} logs/

# Stage-specific outputs (kopt.data, manifests, etc.)
rsync -avz brau0037@deepthought.flinders.edu.au:/scratch/user/brau0037/kelp/ddocent/refopt/ logs/refopt/

# Staging-step outputs
rsync -avz brau0037@deepthought.flinders.edu.au:/scratch/user/brau0037/kelp/ddocent/ref_subset_stage/ logs/ref_subset_stage/
```

### Collaboration loop

1. Claude Code writes/edits script in `scripts/`
2. User reviews, then runs `./deploy.sh`
3. User submits: `sbatch scripts/NN_scriptname.sh`
4. User rsyncs logs back to `logs/`
5. Claude Code reads logs, diagnoses, iterates

---

## Project overview

This project is a staged implementation of the dDocent de novo RADseq reference-building and SNP-calling workflow on the Flinders HPC.

The goals are to:

- run dDocent on HPC in a way that works with the cluster's transient high-speed filesystem
- start from **raw demultiplexed reads** formatted for dDocent
- build a de novo reference using a representative subset of samples
- retain only the key outputs needed for downstream work because scratch space is limited
- add strong logging so transient job-local files can still be audited after the job finishes

The working directory for this project is:

```bash
/scratch/user/brau0037/kelp/ddocent
```

Important subpaths:

```bash
/scratch/user/brau0037/kelp/ddocent/raw
/scratch/user/brau0037/kelp/ddocent/subset_samples.txt
/scratch/user/brau0037/kelp/ddocent/ref_subset_stage
/scratch/user/brau0037/kelp/ddocent/refopt
```

The `raw/` directory contains correctly named dDocent input files like:

```text
AB_01.F.fq.gz
AB_01.R.fq.gz
```

These are symlinks to the real source files, but from dDocent's perspective they behave as the raw demultiplexed inputs.

---

## HPC / BGFS execution model

A key constraint of this system is that **BGFS is ephemeral**. It only exists during the job and disappears when the job ends.

That means:

- the BGFS working directory must be created **inside each SLURM job**
- inputs are staged into BGFS inside the job, usually as symlinks back to scratch
- the actual analysis runs inside BGFS because it is the fast job-local filesystem
- only selected outputs are copied back to scratch before the job exits
- everything left behind in BGFS should be treated as disposable

This is different from a persistent project filesystem. Each stage must therefore be **self-contained**.

A typical pattern is:

1. activate environment
2. validate required inputs on scratch
3. create BGFS working directory
4. symlink required inputs from scratch into BGFS
5. run one analysis stage
6. write logs, manifests, and summary files
7. copy retained outputs back to scratch

This model is already reflected in the trimmed-read SLURM template used previously.

---

## Environment

The jobs should use the user's own Miniforge/mamba installation rather than the system `Miniconda3` module.

Environment root:

```bash
/scratch/user/brau0037/local/miniforge3
```

dDocent environment:

```bash
/scratch/user/brau0037/local/miniforge3/envs/kelp_ddocent
```

Recommended activation block for SLURM scripts:

```bash
export MAMBA_ROOT_PREFIX=/scratch/user/brau0037/local/miniforge3
source /scratch/user/brau0037/local/miniforge3/etc/profile.d/conda.sh
conda activate /scratch/user/brau0037/local/miniforge3/envs/kelp_ddocent
```

Useful validation lines:

```bash
which python
which dDocent
which ReferenceOpt.sh
which RefMapOpt.sh
which fastp
fastp -v
```

---

## Retention policy

Because scratch space is limited, not all intermediates should be retained.

### Keep on scratch

- raw demultiplexed input symlinks in `raw/`
- curated subset list
- trimmed full-dataset reads
- final `reference.fasta`
- `uniqseq.data`
- canonical BAMs and `.bai`
- `TotalRawSNPs.vcf.gz` and index
- `dDocent_main.LOG`
- `dDocent.runs`
- stage manifests, disk-usage files, and summary logs

### Do not retain unless needed for debugging

- large temporary assembly intermediates
- rainbow/CD-HIT helper files
- temporary reference indexes
- raw scatter-gather VCF chunks
- other files that can be regenerated from retained checkpoints

---

## Current status

### Subset strategy

A curated list of **123** samples was created for RefOpt.

This subset was chosen to:

- preserve broad geographic coverage
- avoid obvious poor/very small samples
- provide a better starting point than a naive one-per-prefix selection

A smaller **6-sample pilot** is being used to test the mechanics of `ReferenceOpt.sh` and script behavior before scaling to the full 123-sample subset.

### Resource observations from pilot

The `ReferenceOpt.sh` pilot job:

- ran for about 10.5 hours
- used approximately 132.5 CPU-hours total
- therefore averaged roughly 12.7 CPUs during the run
- grew the BGFS working directory to about 11 GB at the end

This suggests that the original request of 64 CPUs / 128 GB / 7 days is over-allocated for pilot jobs.

Recommended initial resources:

#### Pilot jobs

```bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
```

#### First full RefOpt run

```bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
```

Memory accounting from `sacct` did not appear reliable, so adding `/usr/bin/time -v` around the main command is recommended for better memory reporting.

---

## Script 01: `01_stage_ref_subset.sh`

### Purpose

This job is a validation/staging job only.

It does **not** run `ReferenceOpt.sh`.

It exists to confirm that:

- the conda environment activates correctly inside SLURM
- the curated subset list resolves correctly to raw input files
- both `.F.fq.gz` and `.R.fq.gz` files exist for every subset sample
- BGFS job-local staging works
- manifests and summaries copy back to scratch correctly

### Inputs

- `/scratch/user/brau0037/kelp/ddocent/raw`
- `/scratch/user/brau0037/kelp/ddocent/subset_samples.txt`

### Main actions

- create BGFS work directory
- symlink subset raw reads into BGFS
- count and validate staged files
- generate manifests and pair-size summaries
- copy logs and summaries back to:

```bash
/scratch/user/brau0037/kelp/ddocent/ref_subset_stage
```

### Expected outputs

- `subset_counts.tsv`
- `input_symlinks_manifest.tsv`
- `refopt_symlinks_manifest.tsv`
- `refopt_file_sizes.tsv`
- `refopt_pair_sizes.tsv`
- `staged_files.tsv`
- `missing_files.tsv` if there is a problem

### Outcome so far

This step completed successfully and confirmed that all 123 selected samples staged correctly.

---

## Script 02: `02_run_referenceopt.sh`

### Purpose

This is the first real dDocent optimization job.

It:

- stages a subset of raw reads into BGFS
- runs `ReferenceOpt.sh`
- saves optimization outputs and logs back to scratch

### Inputs

- `/scratch/user/brau0037/kelp/ddocent/raw`
- `/scratch/user/brau0037/kelp/ddocent/subset_samples.txt`

### Outputs

Written to:

```bash
/scratch/user/brau0037/kelp/ddocent/refopt
```

Expected key outputs:

- `kopt.data`
- `uniqseq.data` if produced
- `ReferenceOpt.log`
- manifests and disk usage summaries
- any dDocent logs created by the script

### Important implementation details

- the script creates its BGFS workdir inside the job
- inputs are staged as symlinks
- `ReferenceOpt.sh` requires explicit arguments:

```bash
ReferenceOpt.sh minK1 maxK1 minK2 maxK2 Assembly_Type Number_of_Processors [minSim maxSim increment]
```

### Current parameter approach

The current pilot is being narrowed toward:

```bash
MIN_K1=3
MAX_K1=6
MIN_K2=3
MAX_K2=6
MIN_SIM=0.80
MAX_SIM=0.98
SIM_INC=0.02
```

This broader-but-coarser similarity sweep aligns better with dDocent's own defaults and tutorials than the earlier `0.90-0.99 by 0.01` test.

---

## dDocent / script issues already encountered

### 1. `fastp` version check bug

The installed `fastp` is version `1.3.1`, which is newer than the minimum requirement.

However, dDocent helper scripts were using brittle logic that incorrectly rejected `fastp 1.x`.

This affected:

- `ReferenceOpt.sh`
- `RefMapOpt.sh`

These version checks were patched so that all `1.x` versions are correctly accepted.

### 2. `ReferenceOpt.sh` LENGTH bug

Pilot runs revealed failures at higher `K2` values due to this type of line:

```bash
LENGTH=$(( $LENGTH * 3 / 4 ))
```

In some parameter combinations, `LENGTH` became empty or non-numeric, causing shell arithmetic errors.

The script was patched to:

- check that `totaluniqseq` is non-empty
- validate `LENGTH`
- validate `MaxLen`
- return `NA` for failed parameter combinations instead of crashing

A second change was made so that only numeric values from `kopt.data` are sent to `plot.kopt.data`, avoiding issues when some combinations are skipped.

---

## Meaning of the K parameters

In dDocent, the K parameters are **k-mer sizes** used during de novo reference construction.

Very roughly:

- `K1` affects the earlier clustering/assembly behavior
- `K2` affects a later refinement/filtering stage

Smaller K values are more permissive and sensitive.
Larger K values are more stringent and can fail when coverage or support is low.

That fits the pilot behavior observed so far: higher `K2` values are more likely to produce empty intermediates or edge-case failures.

---

## Plan for remaining steps

### Step 03: inspect `kopt.data`

After a clean pilot run:

- inspect `kopt.data`
- visualize contig counts across tested `K1`, `K2`, and `c`
- identify the region where contig counts stabilize or show an inflection
- choose a narrower final search region if needed

### Step 04: run RefOpt on full curated 123-sample subset

Once the pilot behaves cleanly:

- switch from 6-sample pilot list to the full 123-sample curated subset
- rerun `ReferenceOpt.sh` with the chosen coarse grid or refined grid
- retain `kopt.data`, logs, and manifests

### Step 05: run `RefMapOpt.sh`

After selecting the preferred similarity region:

- stage the same subset into BGFS
- run `RefMapOpt.sh`
- inspect output to refine the final choice of `c`

### Step 06: build the final reference on the subset

Run dDocent on the curated subset using the chosen `K1`, `K2`, and `c`, while skipping later stages.

Retain:

- `reference.fasta`
- `uniqseq.data`
- `dDocent_main.LOG`
- `dDocent.runs`
- compact assembly summaries

### Step 07: full dataset trimming

Run trimming once on the full raw dataset.

Retain:

- trimmed reads
- trim logs
- trim summaries

### Step 08: full mapping to BAMs

Use full trimmed reads plus `reference.fasta`.

Retain:

- BAMs
- `.bai`
- mapping summaries
- `flagstat`
- `idxstats`

### Step 09: BAM filtering

Produce and retain filtered BAMs as the canonical BAM set.

### Step 10: raw variant calling

Use filtered BAMs and the final reference to produce:

- `TotalRawSNPs.vcf.gz`
- index
- variant summary tables

No need to retain `Final.recode.vcf`, because downstream filtering will be done separately.

---

## Logging recommendations

For all future jobs, keep:

- stdout/stderr via SLURM log files
- stage manifests
- disk-usage summaries
- command-line capture for the main command
- `/usr/bin/time -v` for process-level resource use
- `sacct` output after job completion for CPU accounting

Useful post-job command:

```bash
sacct -j <jobid> --format=JobID,JobName%20,AllocCPUS,Elapsed,TotalCPU,MaxRSS,MaxVMSize,State
```

---

## Immediate next actions

1. let the current pilot finish
2. inspect the new `kopt.data`
3. decide whether the pilot grid is now clean enough
4. if yes, move to full 123-sample RefOpt
5. if needed, further refine `ReferenceOpt.sh` or narrow parameter ranges before scaling up

