# Installation Guide

This guide provides detailed instructions for installing dyna-clust and its dependencies.

## System Requirements

- Operating System: Linux, macOS, or Windows (with WSL)
- R version 4.0.0 or higher
- At least 4GB RAM (more for large datasets)
- Sufficient disk space for reference databases and output files

## Step 1: Install R

### Linux (Ubuntu/Debian)
```bash
sudo apt-get update
sudo apt-get install r-base r-base-dev
```

### macOS
```bash
brew install r
```

Or download from [CRAN](https://cran.r-project.org/bin/macosx/)

### Windows
Download and install from [CRAN](https://cran.r-project.org/bin/windows/base/)

## Step 2: Install R Packages

Open R and run:

```r
# Install required packages
install.packages(c("dplyr", "readr", "stringr", "purrr", "tibble", "optparse"))
```

## Step 3: Install BLAST+

BLAST+ is required for sequence classification.

### Linux (Ubuntu/Debian)
```bash
sudo apt-get update
sudo apt-get install ncbi-blast+
```

### macOS
```bash
brew install blast
```

### Manual Installation
Download from [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

Verify installation:
```bash
blastn -version
```

## Step 4: Install dnabarcoder

dnabarcoder is used to compute taxon-specific clustering cutoffs.

### Installation from source

1. Clone the repository:
   ```bash
   git clone https://github.com/vuthuyduong/dnabarcoder.git
   cd dnabarcoder
   ```

2. Follow the installation instructions in the dnabarcoder repository

3. Add dnabarcoder to your PATH or note its location

**Note**: dnabarcoder is optional. If not installed, dyna-clust will use default cutoff values (0.97 for all taxa).

## Step 5: Install dyna-clust

### Option A: Clone from GitHub (Recommended for development)

```bash
git clone https://github.com/LukeLikesDirt/dyna-clust.git
cd dyna-clust
```

Test the installation:
```bash
cd examples
./run_example.sh
```

### Option B: Install as R package

```bash
git clone https://github.com/LukeLikesDirt/dyna-clust.git
R CMD INSTALL dyna-clust
```

Then in R:
```r
library(dynaclust)
```

## Verification

Verify your installation by running the example workflow:

```bash
cd dyna-clust/examples
./run_example.sh
```

If successful, you should see output files in `examples/output/`.

## Troubleshooting

### R packages fail to install

Make sure you have the development tools installed:

**Linux**:
```bash
sudo apt-get install build-essential
```

**macOS**:
```bash
xcode-select --install
```

### BLAST not found

Make sure BLAST+ is in your PATH:
```bash
echo $PATH
which blastn
```

If not in PATH, you can specify the full path when running dyna-clust:
```bash
./scripts/dyna_clust.R --blast /path/to/blastn ...
```

### dnabarcoder not found

dnabarcoder is optional. The tool will use default cutoffs if it's not available. To use dnabarcoder, ensure it's in your PATH or specify the full path:
```bash
./scripts/dyna_clust.R --dnabarcoder /path/to/dnabarcoder ...
```

### Permission denied

Make sure the scripts are executable:
```bash
chmod +x scripts/dyna_clust.R
chmod +x examples/run_example.sh
```

## Getting Help

If you encounter issues:

1. Check the [TESTING.md](TESTING.md) guide
2. Review the [README.md](README.md) for usage examples
3. Open an issue on [GitHub](https://github.com/LukeLikesDirt/dyna-clust/issues)

