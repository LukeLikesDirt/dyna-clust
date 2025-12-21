#' Run Complete Dynamic Clustering Workflow
#'
#' This is the main wrapper function that runs the complete dyna-clust workflow:
#' 1. Compute taxon-specific cutoffs using dnabarcoder
#' 2. Classify ASVs using BLAST and cutoffs
#' 3. Cluster ASVs into OTUs using taxon-specific thresholds
#'
#' @param asv_fasta Path to FASTA file containing ASV sequences to cluster
#' @param reference_fasta Path to reference FASTA file for BLAST and cutoff calculation
#' @param taxonomy_file Path to taxonomy file mapping reference sequences to taxa
#'   (tab-separated with columns: seq_id, taxon)
#' @param output_dir Main output directory for all results
#' @param dnabarcoder_path Path to dnabarcoder executable (default: "dnabarcoder")
#' @param blast_path Path to BLAST executable (default: "blastn")
#' @param num_threads Number of threads for BLAST (default: 1)
#' @param clustering_method Clustering method: "average", "single", or "complete" (default: "average")
#' @param skip_cutoff_computation Skip cutoff computation and use provided cutoffs file (default: FALSE)
#' @param cutoffs_file Path to pre-computed cutoffs file (used if skip_cutoff_computation is TRUE)
#'
#' @return A list containing:
#'   \itemize{
#'     \item cutoffs: Data frame of taxon-specific cutoffs
#'     \item classifications: Data frame of ASV classifications
#'     \item otu_table: Data frame of OTU assignments
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run complete workflow
#' results <- run_dyna_clust(
#'   asv_fasta = "my_asvs.fasta",
#'   reference_fasta = "reference_database.fasta",
#'   taxonomy_file = "reference_taxonomy.txt",
#'   output_dir = "dyna_clust_output",
#'   num_threads = 4
#' )
#' 
#' # Access results
#' head(results$otu_table)
#' }
run_dyna_clust <- function(asv_fasta,
                          reference_fasta,
                          taxonomy_file,
                          output_dir,
                          dnabarcoder_path = "dnabarcoder",
                          blast_path = "blastn",
                          num_threads = 1,
                          clustering_method = "average",
                          skip_cutoff_computation = FALSE,
                          cutoffs_file = NULL) {
  
  # Check inputs
  if (!file.exists(asv_fasta)) {
    stop("ASV FASTA file not found: ", asv_fasta)
  }
  if (!file.exists(reference_fasta)) {
    stop("Reference FASTA file not found: ", reference_fasta)
  }
  if (!skip_cutoff_computation && !file.exists(taxonomy_file)) {
    stop("Taxonomy file not found: ", taxonomy_file)
  }
  
  # Create main output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("===== DYNA-CLUST: Dynamic OTU Clustering Pipeline =====")
  message("Started at: ", Sys.time())
  
  # Step 1: Compute taxon-specific cutoffs
  if (skip_cutoff_computation) {
    if (is.null(cutoffs_file) || !file.exists(cutoffs_file)) {
      stop("When skipping cutoff computation, a valid cutoffs_file must be provided")
    }
    message("\n--- Step 1: Loading pre-computed cutoffs ---")
    cutoffs <- readr::read_tsv(cutoffs_file, col_types = readr::cols())
  } else {
    message("\n--- Step 1: Computing taxon-specific cutoffs with dnabarcoder ---")
    cutoffs_dir <- file.path(output_dir, "01_cutoffs")
    cutoffs <- compute_cutoffs(
      fasta_file = reference_fasta,
      taxonomy_file = taxonomy_file,
      output_dir = cutoffs_dir,
      dnabarcoder_path = dnabarcoder_path
    )
  }
  
  message(sprintf("Loaded %d taxon-specific cutoffs", nrow(cutoffs)))
  
  # Step 2: Classify ASVs using BLAST
  message("\n--- Step 2: Classifying ASVs using BLAST ---")
  blast_dir <- file.path(output_dir, "02_classification")
  classifications <- classify_asvs(
    asv_fasta = asv_fasta,
    reference_db = reference_fasta,
    cutoffs = cutoffs,
    output_dir = blast_dir,
    blast_path = blast_path,
    num_threads = num_threads
  )
  
  message(sprintf("Classified %d ASVs", nrow(classifications)))
  confident <- sum(classifications$classification_status == "confident")
  uncertain <- sum(classifications$classification_status == "uncertain")
  message(sprintf("  - Confident: %d (%.1f%%)", confident, 100 * confident / nrow(classifications)))
  message(sprintf("  - Uncertain: %d (%.1f%%)", uncertain, 100 * uncertain / nrow(classifications)))
  
  # Step 3: Cluster into OTUs
  message("\n--- Step 3: Clustering ASVs into OTUs ---")
  clustering_dir <- file.path(output_dir, "03_clustering")
  otu_table <- cluster_otus(
    classifications = classifications,
    asv_fasta = asv_fasta,
    output_dir = clustering_dir,
    clustering_method = clustering_method
  )
  
  # Summary
  message("\n===== Pipeline Complete =====")
  message("Finished at: ", Sys.time())
  message(sprintf("Total ASVs processed: %d", nrow(classifications)))
  message(sprintf("Total OTUs created: %d", length(unique(otu_table$otu_id))))
  message(sprintf("Results saved to: %s", output_dir))
  
  # Return results
  results <- list(
    cutoffs = cutoffs,
    classifications = classifications,
    otu_table = otu_table
  )
  
  return(invisible(results))
}
