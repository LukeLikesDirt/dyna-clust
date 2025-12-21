# dyna-clust

**Taxon-specific dynamic clustering of OTUs from metabarcoding data**

## Overview

`dyna-clust` is an R-based tool for performing dynamic clustering of Amplicon Sequence Variants (ASVs) into Operational Taxonomic Units (OTUs) using taxon-specific similarity thresholds. Unlike traditional OTU clustering approaches that use a single fixed threshold (e.g., 97% similarity), dyna-clust computes and applies different thresholds for different taxonomic groups, improving the biological accuracy of OTU definitions.

## Workflow

The dyna-clust pipeline consists of three main steps:

1. **Compute taxon-specific cutoffs**: Using [dnabarcoder](https://github.com/vuthuyduong/dnabarcoder), calculate optimal clustering thresholds for each taxonomic group in your reference database.

2. **Classify ASVs**: Use BLAST to align your ASV sequences against the reference database, then apply the taxon-specific cutoffs to confidently assign taxonomy.

3. **Cluster into OTUs**: Group ASVs into OTUs within each taxonomic group using the appropriate similarity threshold.

## Installation

### Prerequisites

- **R** (>= 4.0.0)
- **BLAST+** (>= 2.10.0) - [Download](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- **dnabarcoder** - [Installation instructions](https://github.com/vuthuyduong/dnabarcoder)

### Install R dependencies

```r
# Install required packages
install.packages(c("dplyr", "readr", "stringr", "purrr", "tibble", "optparse"))
```

### Install dyna-clust

```bash
# Clone the repository
git clone https://github.com/LukeLikesDirt/dyna-clust.git
cd dyna-clust

# Or install as R package
R CMD INSTALL .
```

## Usage

### Command-line interface

The easiest way to use dyna-clust is through the command-line script:

```bash
./scripts/dyna_clust.R \
  --asv-fasta my_asvs.fasta \
  --reference reference_db.fasta \
  --taxonomy reference_taxonomy.txt \
  --output output_directory \
  --threads 4
```

**Required arguments:**
- `--asv-fasta`: FASTA file containing your ASV sequences to cluster
- `--reference`: FASTA file with reference sequences for BLAST database
- `--taxonomy`: Tab-separated file mapping reference sequences to taxa (columns: seq_id, taxon)
- `--output`: Directory for output files

**Optional arguments:**
- `--dnabarcoder`: Path to dnabarcoder executable (default: "dnabarcoder")
- `--blast`: Path to BLAST executable (default: "blastn")
- `--threads`: Number of threads for BLAST (default: 1)
- `--method`: Clustering method - average, single, or complete (default: "average")
- `--skip-cutoffs`: Skip cutoff computation and use pre-computed cutoffs
- `--cutoffs-file`: Path to pre-computed cutoffs file (required if --skip-cutoffs)

### R interface

You can also use dyna-clust directly in R:

```r
# Load the package
library(dynaclust)

# Or source the functions
source("R/compute_cutoffs.R")
source("R/classify_asvs.R")
source("R/cluster_otus.R")
source("R/run_dyna_clust.R")

# Run the complete pipeline
results <- run_dyna_clust(
  asv_fasta = "my_asvs.fasta",
  reference_fasta = "reference_db.fasta",
  taxonomy_file = "reference_taxonomy.txt",
  output_dir = "dyna_clust_output",
  num_threads = 4
)

# Access results
head(results$cutoffs)
head(results$classifications)
head(results$otu_table)
```

### Step-by-step workflow

You can also run each step independently:

```r
# Step 1: Compute taxon-specific cutoffs
cutoffs <- compute_cutoffs(
  fasta_file = "reference_db.fasta",
  taxonomy_file = "reference_taxonomy.txt",
  output_dir = "cutoffs_output"
)

# Step 2: Classify ASVs using BLAST
classifications <- classify_asvs(
  asv_fasta = "my_asvs.fasta",
  reference_db = "reference_db.fasta",
  cutoffs = cutoffs,
  output_dir = "classification_output",
  num_threads = 4
)

# Step 3: Cluster into OTUs
otu_table <- cluster_otus(
  classifications = classifications,
  asv_fasta = "my_asvs.fasta",
  output_dir = "clustering_output"
)
```

## Input File Formats

### ASV FASTA file
Standard FASTA format with unique identifiers:
```
>ASV_001
ATCGATCGATCG...
>ASV_002
GCTAGCTAGCTA...
```

### Reference FASTA file
Standard FASTA format. Sequence IDs should be included in the taxonomy file:
```
>ref_seq_001
ATCGATCGATCG...
>ref_seq_002
GCTAGCTAGCTA...
```

### Taxonomy file
Tab-separated file with two columns (no header):
```
ref_seq_001    Bacteria
ref_seq_002    Fungi
ref_seq_003    Bacteria
ref_seq_004    Protist
```

## Output Files

The pipeline creates three subdirectories in the output folder:

### 1. `01_cutoffs/`
- `taxon_cutoffs.tsv`: Taxon-specific clustering thresholds
- Individual FASTA files for each taxon
- dnabarcoder output files

### 2. `02_classification/`
- `asv_classifications.tsv`: ASV taxonomy assignments with confidence scores
- `blast_results.txt`: Raw BLAST output

### 3. `03_clustering/`
- `otu_table.tsv`: Final OTU assignments for each ASV
- `otu_representatives.fasta`: Representative sequences for each OTU

## Example Workflow

See the `examples/` directory for a complete example workflow with test data.

```bash
# Run example
cd examples
./run_example.sh
```

## Citation

If you use dyna-clust in your research, please cite:

- This tool: (Publication pending)
- dnabarcoder: [Zhang et al. 2024](https://github.com/vuthuyduong/dnabarcoder)
- BLAST: Camacho C. et al. (2009) BLAST+: architecture and applications. BMC Bioinformatics 10:421.

## License

MIT License - see LICENSE file for details

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## Contact

For questions or issues, please open an issue on GitHub.

## References

- dnabarcoder: https://github.com/vuthuyduong/dnabarcoder
- BLAST+: https://blast.ncbi.nlm.nih.gov/

