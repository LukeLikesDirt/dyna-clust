#!/usr/bin/env Rscript

#' Command-line interface for dyna-clust
#' 
#' This script provides a command-line interface to the dyna-clust package
#' for dynamic clustering of OTUs from metabarcoding data.

# Load required packages
suppressPackageStartupMessages({
  library(optparse)
})

# Define command-line options
option_list <- list(
  make_option(c("-a", "--asv-fasta"), type="character", default=NULL,
              help="Path to ASV FASTA file [required]", metavar="FILE"),
  make_option(c("-r", "--reference"), type="character", default=NULL,
              help="Path to reference FASTA file [required]", metavar="FILE"),
  make_option(c("-t", "--taxonomy"), type="character", default=NULL,
              help="Path to taxonomy file (tab-separated: seq_id, taxon) [required unless --skip-cutoffs]", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="dyna_clust_output",
              help="Output directory [default: %default]", metavar="DIR"),
  make_option(c("-d", "--dnabarcoder"), type="character", default="dnabarcoder",
              help="Path to dnabarcoder executable [default: %default]", metavar="PATH"),
  make_option(c("-b", "--blast"), type="character", default="blastn",
              help="Path to BLAST executable [default: %default]", metavar="PATH"),
  make_option(c("-n", "--threads"), type="integer", default=1,
              help="Number of threads for BLAST [default: %default]", metavar="N"),
  make_option(c("-m", "--method"), type="character", default="average",
              help="Clustering method: average, single, or complete [default: %default]", metavar="METHOD"),
  make_option(c("--skip-cutoffs"), action="store_true", default=FALSE,
              help="Skip cutoff computation and use pre-computed cutoffs"),
  make_option(c("-c", "--cutoffs-file"), type="character", default=NULL,
              help="Path to pre-computed cutoffs file (required if --skip-cutoffs)", metavar="FILE")
)

opt_parser <- OptionParser(
  usage = "usage: %prog [options]",
  option_list = option_list,
  description = "\ndyna-clust: Taxon-specific dynamic clustering of OTUs\n\nThis tool performs dynamic clustering of Amplicon Sequence Variants (ASVs) using:\n  1. Taxon-specific cutoffs computed with dnabarcoder\n  2. BLAST-based classification\n  3. Custom clustering within taxonomic groups"
)

opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$`asv-fasta`)) {
  print_help(opt_parser)
  stop("Missing required argument: --asv-fasta", call.=FALSE)
}
if (is.null(opt$reference)) {
  print_help(opt_parser)
  stop("Missing required argument: --reference", call.=FALSE)
}
if (!opt$`skip-cutoffs` && is.null(opt$taxonomy)) {
  print_help(opt_parser)
  stop("Missing required argument: --taxonomy (or use --skip-cutoffs with --cutoffs-file)", call.=FALSE)
}
if (opt$`skip-cutoffs` && is.null(opt$`cutoffs-file`)) {
  print_help(opt_parser)
  stop("When using --skip-cutoffs, --cutoffs-file must be provided", call.=FALSE)
}

# Source the R functions (if not in a package)
script_dir <- dirname(sys.frame(1)$ofile)
r_dir <- file.path(dirname(script_dir), "R")

if (dir.exists(r_dir)) {
  source(file.path(r_dir, "compute_cutoffs.R"))
  source(file.path(r_dir, "classify_asvs.R"))
  source(file.path(r_dir, "cluster_otus.R"))
  source(file.path(r_dir, "run_dyna_clust.R"))
} else {
  # Try loading as package
  tryCatch({
    library(dynaclust)
  }, error = function(e) {
    stop("Could not find dyna-clust R functions. Make sure the package is installed or R files are in the correct location.", call.=FALSE)
  })
}

# Run the pipeline
cat("\n")
cat("╔════════════════════════════════════════════════════════════════════╗\n")
cat("║                          DYNA-CLUST                                ║\n")
cat("║        Taxon-specific Dynamic Clustering of OTUs                   ║\n")
cat("╚════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

tryCatch({
  results <- run_dyna_clust(
    asv_fasta = opt$`asv-fasta`,
    reference_fasta = opt$reference,
    taxonomy_file = opt$taxonomy,
    output_dir = opt$output,
    dnabarcoder_path = opt$dnabarcoder,
    blast_path = opt$blast,
    num_threads = opt$threads,
    clustering_method = opt$method,
    skip_cutoff_computation = opt$`skip-cutoffs`,
    cutoffs_file = opt$`cutoffs-file`
  )
  
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════════════╗\n")
  cat("║                      Pipeline Successful!                          ║\n")
  cat("╚════════════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
}, error = function(e) {
  cat("\n")
  cat("╔════════════════════════════════════════════════════════════════════╗\n")
  cat("║                          ERROR                                     ║\n")
  cat("╚════════════════════════════════════════════════════════════════════╝\n")
  cat("\nError: ", e$message, "\n")
  quit(status = 1)
})
