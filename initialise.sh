#!/usr/bin/env bash
set -euo pipefail

#----------------------------------------
# initialise.sh
#
# Bootstrap script to initialise a new Ubuntu EC2
# node for the GPU pipeline.
#----------------------------------------

# 1) Update package lists and install system tools
sudo apt-get update -y
sudo apt-get install -y \
    git \
    wget \
    bzip2 \
    ca-certificates \
    build-essential \
    libssl-dev \
    libffi-dev \
    python3-dev \
    awscli \
    samtools \
    bwa \
    bcftools \
    vcftools \
    pigz \
    plink2 \
    unzip

# 2) Download and unpack the pre-built Conda environment
S3_BUCKET=djw-conda-artifacts-2025
CONDA_ARCHIVE=base_with_bio_20250624.tar.gz
TARGET_DIR="$HOME/conda_env"

aws s3 cp s3://$S3_BUCKET/$CONDA_ARCHIVE .
mkdir -p "$TARGET_DIR"
tar -xzf $CONDA_ARCHIVE -C "$TARGET_DIR"

# 3) Run conda-unpack to fix paths
"$TARGET_DIR/bin/conda-unpack"

# 4) Add the unpacked env to PATH
if ! grep -q 'conda_env/bin' ~/.bashrc; then
  echo "# Add local conda environment to PATH" >> ~/.bashrc
  echo "export PATH=\$HOME/conda_env/bin:\$PATH" >> ~/.bashrc
fi

# 5) Clone and prepare the GPU pipeline repo
PIPELINE_DIR="$HOME/gpu_pipeline"
if [ ! -d "$PIPELINE_DIR" ]; then
  git clone https://github.com/DJWalsh27/gpu-pipeline-djw.git "$PIPELINE_DIR"
fi

# Done
printf "\nInitialization complete. Please restart your shell or run 'source ~/.bashrc'.\n"
