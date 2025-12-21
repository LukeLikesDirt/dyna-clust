#!/usr/bin/env Rscript

# Example R script demonstrating the dyna-clust workflow

cat("\n")
cat("======================================\n")
cat("  dyna-clust Example Workflow (R)\n")
cat("======================================\n")
cat("\n")

# Get script directory
script_dir <- dirname(sys.frame(1)$ofile)
setwd(script_dir)

# Source the dyna-clust functions
cat("Loading dyna-clust functions...\n")
source("../R/compute_cutoffs.R")
source("../R/classify_asvs.R")
source("../R/cluster_otus.R")
source("../R/run_dyna_clust.R")

# Set up paths
asv_fasta <- "asvs.fasta"
reference <- "reference.fasta"
taxonomy <- "taxonomy.txt"
output_dir <- "output"

# Clean up previous output
if (dir.exists(output_dir)) {
  cat("Removing previous output directory...\n")
  unlink(output_dir, recursive = TRUE)
}

cat("\nInput files:\n")
cat("  - ASV FASTA:", asv_fasta, "\n")
cat("  - Reference:", reference, "\n")
cat("  - Taxonomy:", taxonomy, "\n")
cat("\n")

# Run the complete pipeline
cat("Running dyna-clust pipeline...\n")
cat("\n")

results <- run_dyna_clust(
  asv_fasta = asv_fasta,
  reference_fasta = reference,
  taxonomy_file = taxonomy,
  output_dir = output_dir,
  num_threads = 1
)

cat("\n")
cat("======================================\n")
cat("  Example workflow complete!\n")
cat("======================================\n")
cat("\n")

# Display summary
cat("Results summary:\n")
cat("  - Taxa with cutoffs:", nrow(results$cutoffs), "\n")
cat("  - ASVs classified:", nrow(results$classifications), "\n")
cat("  - OTUs created:", length(unique(results$otu_table$otu_id)), "\n")
cat("\n")

cat("Cutoffs by taxon:\n")
print(results$cutoffs)
cat("\n")

cat("OTU summary:\n")
print(table(results$otu_table$taxon))
cat("\n")

cat("Output files saved to:", output_dir, "\n")
cat("\n")
