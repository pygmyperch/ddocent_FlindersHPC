#!/bin/bash
# Deploy scripts to HPC scratch.
# Usage: ./deploy.sh
# Add --dry-run to preview without transferring.

REMOTE="brau0037@deepthought.flinders.edu.au:/scratch/user/brau0037/kelp/ddocent/scripts/"

rsync -avz --progress "$@" scripts/ "$REMOTE"
