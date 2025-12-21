#!/bin/bash

# Example workflow for dyna-clust
# This script demonstrates how to use dyna-clust with example data

set -e

echo "======================================"
echo "  dyna-clust Example Workflow"
echo "======================================"
echo ""

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# Set paths
ASV_FASTA="asvs.fasta"
REFERENCE="reference.fasta"
TAXONOMY="taxonomy.txt"
OUTPUT_DIR="output"

# Clean up previous output
if [ -d "$OUTPUT_DIR" ]; then
    echo "Removing previous output directory..."
    rm -rf "$OUTPUT_DIR"
fi

echo "Input files:"
echo "  - ASV FASTA: $ASV_FASTA"
echo "  - Reference: $REFERENCE"
echo "  - Taxonomy: $TAXONOMY"
echo ""

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    echo "Error: Rscript not found. Please install R."
    exit 1
fi

# Run dyna-clust
echo "Running dyna-clust pipeline..."
echo ""

../scripts/dyna_clust.R \
    --asv-fasta "$ASV_FASTA" \
    --reference "$REFERENCE" \
    --taxonomy "$TAXONOMY" \
    --output "$OUTPUT_DIR" \
    --threads 1

echo ""
echo "======================================"
echo "  Example workflow complete!"
echo "======================================"
echo ""
echo "Results are in: $OUTPUT_DIR"
echo ""
echo "Output files:"
echo "  1. Cutoffs: $OUTPUT_DIR/01_cutoffs/taxon_cutoffs.tsv"
echo "  2. Classifications: $OUTPUT_DIR/02_classification/asv_classifications.tsv"
echo "  3. OTU table: $OUTPUT_DIR/03_clustering/otu_table.tsv"
echo "  4. OTU representatives: $OUTPUT_DIR/03_clustering/otu_representatives.fasta"
echo ""
