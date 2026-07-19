#!/usr/bin/env bash
set -euo pipefail

# (A) If your cluster uses modules for conda, load it here (optional)
# module purge
# module load Miniconda3  # or Anaconda3 / Miniforge (use module spider to find exact name)

# (B) Initialize conda for NON-interactive shells (this is the key line)
# Try common locations—use the one that exists on your system.

if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/mambaforge/etc/profile.d/conda.sh" ]; then
  source "$HOME/mambaforge/etc/profile.d/conda.sh"
elif [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniforge3/etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
  # last resort; sometimes works if conda is already initialized by cluster
  eval "$(conda shell.bash hook)"
else
  echo "ERROR: conda.sh not found and conda not in PATH"
  exit 1
fi

conda activate splitseq   # <-- change to your exact env name if different

# Sanity check (optional but helpful)
echo "Activated env: $CONDA_DEFAULT_ENV"
