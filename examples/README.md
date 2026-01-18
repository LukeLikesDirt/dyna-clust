# Example workflow for dyna-clust

This directory contains example data and scripts to demonstrate the dyna-clust workflow.

## Files

- `reference.fasta` - Small reference database with 20 sequences from 4 taxa
- `taxonomy.txt` - Taxonomy mapping for reference sequences
- `asvs.fasta` - Sample ASVs to cluster (10 sequences)
- `run_example.sh` - Shell script to run the complete workflow
- `example.R` - R script demonstrating the workflow

## Running the example

### Using the shell script:

```bash
./run_example.sh
```

### Using R:

```r
source("example.R")
```

## Expected output

The example will create an `output/` directory containing:

1. `01_cutoffs/` - Taxon-specific clustering cutoffs
2. `02_classification/` - ASV classifications
3. `03_clustering/` - Final OTU table and representative sequences

## Note

This example uses simulated data for demonstration purposes. The sequences are randomly generated and do not represent real biological data.
