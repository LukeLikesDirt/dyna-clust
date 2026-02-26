#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out

# Script name:  01.Prepare_eukaryome_ITS.sh
# Description:  Download and reformat EUKARYOME sequences in preparation for 
#               dnabarcoder predictions, chimera removal and taxanomic 
#               assignment.
# Note:         This script must be run from the source directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

# EUKARYOME PARAMETERS: Download URL (EUKARYOME EUK ITS v.2.0)
readonly DOWNLOAD_URL="https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_ITS_v2.0.zip"
readonly DOWNLOAD_FILE="./data/General_EUK_ITS_v2.0.zip"
readonly EXTRACTED_DIR="./data/"

# FILE PATHS
readonly INPUT_SEQS="./data/General_EUK_ITS_v2.0.fasta"
readonly OUTPUT_SEQS="./data/eukaryome_ITS_v2.0.fasta"
readonly OUTPUT_CLASSIFICATION="./data/eukaryome_ITS_v2.0.classification"

# Create data directory if it doesn't exist
mkdir -p ./data

# HELPER SCRIPT: R script for reformatting EUKARYOME headers
readonly REFORMAT_REFSEQS="./reformat_eukaryome.R"

# Activate conda environment
echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_env

# =============================================================================
# FILE DOWNLOAD: EUKARYOME
# =============================================================================

echo "=== DOWNLOADING EUKARYOME DATABASE ==="

echo "Downloading file from: $DOWNLOAD_URL"
if ! curl -o "$DOWNLOAD_FILE" "$DOWNLOAD_URL"; then
  echo "ERROR: Failed to download file from $DOWNLOAD_URL" >&2
  exit 1
fi

echo "Unzipping downloaded file..."
if ! 7z x "$DOWNLOAD_FILE" -o"$EXTRACTED_DIR" -y; then
  echo "ERROR: Failed to unzip $DOWNLOAD_FILE" >&2
  exit 1
fi

# Clean up zip file
rm -f "$DOWNLOAD_FILE"

echo "Download and extraction completed successfully!"
echo ""

# =============================================================================
# REFORMAT HEADERS
# =============================================================================

echo "=== REFORMATTING HEADERS ==="

# Check if R script exists
if [[ ! -f "$REFORMAT_REFSEQS" ]]; then
  echo "ERROR: R script not found: $REFORMAT_REFSEQS" >&2
  echo "Please ensure the script is located at: $REFORMAT_REFSEQS" >&2
  exit 1
fi

if [[ ! -f "$INPUT_SEQS" ]]; then
  echo "ERROR: EUKARYOME file not found: $INPUT_SEQS" >&2
  exit 1
fi

# Execute R script for reformatting
echo "Executing R script for header reformatting and merging..."
Rscript "$REFORMAT_REFSEQS" "$INPUT_SEQS" "$OUTPUT_SEQS" "$OUTPUT_CLASSIFICATION"

# Check if R script executed successfully
if [[ $? -ne 0 ]]; then
  echo "ERROR: R script execution failed!" >&2
  exit 1
fi

echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo ""
echo "Reformatted output saved to:"
echo " - Sequences: $OUTPUT_SEQS"
echo " - Classification: $OUTPUT_CLASSIFICATION"
echo ""