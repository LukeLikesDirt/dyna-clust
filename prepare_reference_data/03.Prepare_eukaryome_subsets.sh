#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=0-01:00:00
#SBATCH --partition=day
#SBATCH --output=slurm/%x.%j.out

# Script name:  03.Prepare_eukaryome_subsets.sh
# Description:  Generate ID subset files from the reformatted EUKARYOME
#               database for use in DNAbarcoder similarity predictions.
#               Produces one ID file per unique-sequence rank (STEP 1) and
#               one ID file per valid (target x parent) rank combination
#               (STEP 2).
# Note:         This script must be run from the source directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

# INPUT FILES (outputs from 01.Prepare_eukaryome_ITS.sh)
readonly INPUT_SEQS="./data/eukaryome_ITS_v2.0.fasta"
readonly INPUT_CLASSIFICATION="./data/eukaryome_ITS_v2.0.classification"

# OUTPUT DIRECTORY for ID files
readonly OUTPUT_DIR="./data"
mkdir -p "$OUTPUT_DIR"

# FILTER CONSTANTS
# min_subgroups: minimum unique child taxa a parent-taxon
readonly MIN_SUBGROUPS=10
# min_sequences: minimum sequences a parent-taxon must retain after max_proportion cap
readonly MIN_SEQUENCES=30
# max_sequences: maximum sequences per parent-taxon; excess is balanced-downsampled
readonly MAX_SEQUENCES=25000
# max_proportion: maximum fraction the proportion child taxon may represent per parent-taxon
readonly MAX_PROPORTION=0.5

# HELPER SCRIPT
readonly PREPARE_SUBSETS="./prepare_eukaryome_subsets.R"

# Activate conda environment
echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_env

# =============================================================================
# INPUT VALIDATION
# =============================================================================

echo "=== VALIDATING INPUTS ==="

if [[ ! -f "$PREPARE_SUBSETS" ]]; then
    echo "ERROR: R script not found: $PREPARE_SUBSETS" >&2
    echo "Please ensure the script is located at: $PREPARE_SUBSETS" >&2
    exit 1
fi

if [[ ! -f "$INPUT_SEQS" ]]; then
    echo "ERROR: FASTA file not found: $INPUT_SEQS" >&2
    echo "Please run 01.Prepare_eukaryome_ITS.sh first." >&2
    exit 1
fi

if [[ ! -f "$INPUT_CLASSIFICATION" ]]; then
    echo "ERROR: Classification file not found: $INPUT_CLASSIFICATION" >&2
    echo "Please run 01.Prepare_eukaryome_ITS.sh first." >&2
    exit 1
fi

echo "All input files found."
echo ""

# =============================================================================
# PREPARE EUKARYOME SUBSETS
# =============================================================================

echo "=== PREPARING EUKARYOME SUBSETS ==="
echo "Input FASTA          : $INPUT_SEQS"
echo "Input classification : $INPUT_CLASSIFICATION"
echo "Output directory     : $OUTPUT_DIR"
echo "min_subgroups        : $MIN_SUBGROUPS"
echo "min_sequences        : $MIN_SEQUENCES"
echo "max_sequences        : $MAX_SEQUENCES"
echo "max_proportion       : $MAX_PROPORTION"
echo ""

Rscript "$PREPARE_SUBSETS" \
    "$INPUT_SEQS" \
    "$INPUT_CLASSIFICATION" \
    "$OUTPUT_DIR" \
    "$MIN_SUBGROUPS" \
    "$MIN_SEQUENCES" \
    "$MAX_SEQUENCES" \
    "$MAX_PROPORTION"

if [[ $? -ne 0 ]]; then
    echo "ERROR: R script execution failed!" >&2
    exit 1
fi

echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo ""
echo "ID files written to: $OUTPUT_DIR"
echo "  STEP 1 unique-sequence IDs : <rank>_unique_id.txt"
echo "  STEP 2 prediction IDs      : <target>_pred_id_<parent>.txt"