# Testing and Validation

This document describes how to test the dyna-clust tool.

## Prerequisites

Before running tests, ensure you have:

1. **R** (>= 4.0.0) installed
2. Required R packages installed:
   ```r
   install.packages(c("dplyr", "readr", "stringr", "purrr", "tibble", "optparse"))
   ```
3. **BLAST+** (optional for full testing, will use fallback if not available)
4. **dnabarcoder** (optional for full testing, will use defaults if not available)

## Quick Test

Run the example workflow to verify basic functionality:

```bash
cd examples
./run_example.sh
```

Or using R:

```bash
cd examples
Rscript example.R
```

## Expected Behavior

The test should:

1. ✓ Parse input files (FASTA and taxonomy)
2. ✓ Compute cutoffs for each taxon (or use defaults if dnabarcoder not available)
3. ✓ Run BLAST classification (or skip if BLAST not available)
4. ✓ Cluster ASVs into OTUs using hierarchical clustering
5. ✓ Generate output files in three subdirectories

## Validation Checklist

- [ ] Script runs without errors
- [ ] Output directory is created with three subdirectories
- [ ] Cutoffs file contains entries for each taxon
- [ ] Classifications file contains all ASVs
- [ ] OTU table is generated
- [ ] Representative sequences FASTA is created

## Testing Individual Components

### Test cutoff computation

```r
source("R/compute_cutoffs.R")

cutoffs <- compute_cutoffs(
  fasta_file = "examples/reference.fasta",
  taxonomy_file = "examples/taxonomy.txt",
  output_dir = "/tmp/test_cutoffs"
)

# Should return a data frame with columns: taxon, cutoff, method, n_sequences
print(cutoffs)
```

### Test ASV classification

```r
source("R/compute_cutoffs.R")
source("R/classify_asvs.R")

# First get cutoffs
cutoffs <- compute_cutoffs(
  fasta_file = "examples/reference.fasta",
  taxonomy_file = "examples/taxonomy.txt",
  output_dir = "/tmp/test_cutoffs"
)

# Then classify
classifications <- classify_asvs(
  asv_fasta = "examples/asvs.fasta",
  reference_db = "examples/reference.fasta",
  cutoffs = cutoffs,
  output_dir = "/tmp/test_classification"
)

print(classifications)
```

### Test OTU clustering

```r
source("R/compute_cutoffs.R")
source("R/classify_asvs.R")
source("R/cluster_otus.R")

# Get cutoffs and classifications first
cutoffs <- compute_cutoffs(
  fasta_file = "examples/reference.fasta",
  taxonomy_file = "examples/taxonomy.txt",
  output_dir = "/tmp/test_cutoffs"
)

classifications <- classify_asvs(
  asv_fasta = "examples/asvs.fasta",
  reference_db = "examples/reference.fasta",
  cutoffs = cutoffs,
  output_dir = "/tmp/test_classification"
)

# Then cluster
otu_table <- cluster_otus(
  classifications = classifications,
  asv_fasta = "examples/asvs.fasta",
  output_dir = "/tmp/test_clustering"
)

print(otu_table)
```

## Known Limitations

1. If dnabarcoder is not available, the tool will use a default cutoff of 0.97 for all taxa
2. If BLAST is not available, the classification step will fail - ensure BLAST+ is installed
3. The sequence similarity calculation uses a simple identity metric - for production use, consider using more sophisticated alignment methods

## Troubleshooting

### Issue: "dnabarcoder not found"
**Solution**: Install dnabarcoder or the tool will automatically use default cutoffs (0.97)

### Issue: "BLAST output file not created"
**Solution**: Ensure BLAST+ is installed and in your PATH, or specify the full path with `--blast`

### Issue: "No sequences found in FASTA file"
**Solution**: Check that your FASTA file is properly formatted with headers starting with ">"

### Issue: "Taxonomy file not found"
**Solution**: Ensure the taxonomy file exists and is tab-separated with columns: seq_id and taxon
