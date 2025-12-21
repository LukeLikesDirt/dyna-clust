#' Cluster OTUs Using Taxon-Specific Cutoffs
#'
#' This function performs dynamic clustering of classified ASVs into OTUs using
#' taxon-specific similarity cutoffs. ASVs are grouped by their assigned taxon,
#' and clustering is performed within each group using the appropriate cutoff.
#'
#' @param classifications Data frame of ASV classifications (output from classify_asvs)
#' @param asv_fasta Path to FASTA file containing ASV sequences
#' @param output_dir Directory for output files
#' @param clustering_method Method for clustering: "hierarchical", "single", "complete", or "average" (default: "average")
#' @param min_cluster_size Minimum number of ASVs to form a cluster (default: 1)
#'
#' @return A data frame with OTU assignments containing:
#'   \itemize{
#'     \item asv_id: ASV identifier
#'     \item otu_id: Assigned OTU identifier
#'     \item taxon: Taxonomic group
#'     \item cluster_size: Number of ASVs in the OTU
#'     \item representative_seq: ID of the representative sequence for the OTU
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' otu_table <- cluster_otus(
#'   classifications = asv_classifications,
#'   asv_fasta = "asvs.fasta",
#'   output_dir = "otu_clustering"
#' )
#' }
cluster_otus <- function(classifications,
                        asv_fasta,
                        output_dir,
                        clustering_method = "average",
                        min_cluster_size = 1) {
  
  # Check inputs
  if (!file.exists(asv_fasta)) {
    stop("ASV FASTA file not found: ", asv_fasta)
  }
  if (!is.data.frame(classifications)) {
    stop("classifications must be a data frame")
  }
  if (!all(c("asv_id", "taxon", "cutoff_used") %in% names(classifications))) {
    stop("classifications must contain 'asv_id', 'taxon', and 'cutoff_used' columns")
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Read ASV sequences
  message("Reading ASV sequences...")
  asv_sequences <- read_fasta(asv_fasta)
  
  # Filter to only confidently classified ASVs
  confident_asvs <- classifications %>%
    dplyr::filter(classification_status == "confident")
  
  message(sprintf("Clustering %d confidently classified ASVs across %d taxa",
                 nrow(confident_asvs), 
                 length(unique(confident_asvs$taxon))))
  
  # Cluster by taxon
  otu_assignments_list <- list()
  otu_counter <- 1
  
  taxa <- unique(confident_asvs$taxon)
  
  for (taxon in taxa) {
    message(sprintf("Processing taxon: %s", taxon))
    
    # Get ASVs for this taxon
    taxon_asvs <- confident_asvs %>%
      dplyr::filter(taxon == !!taxon)
    
    # Get cutoff for this taxon
    cutoff <- unique(taxon_asvs$cutoff_used)[1]
    
    # Get sequences for these ASVs
    taxon_seq_ids <- taxon_asvs$asv_id
    taxon_sequences <- asv_sequences[names(asv_sequences) %in% taxon_seq_ids]
    
    if (length(taxon_sequences) == 0) {
      warning("No sequences found for taxon: ", taxon)
      next
    }
    
    # Compute pairwise similarities
    similarity_matrix <- compute_sequence_similarities(taxon_sequences)
    
    # Perform clustering
    clusters <- perform_clustering(
      similarity_matrix = similarity_matrix,
      cutoff = cutoff,
      method = clustering_method
    )
    
    # Assign OTU IDs
    for (cluster_id in unique(clusters)) {
      cluster_members <- names(clusters)[clusters == cluster_id]
      
      # Skip if cluster is too small
      if (length(cluster_members) < min_cluster_size) {
        next
      }
      
      # Select representative sequence (most abundant or first)
      representative <- cluster_members[1]
      
      otu_assignments_list[[otu_counter]] <- data.frame(
        asv_id = cluster_members,
        otu_id = sprintf("OTU_%04d", otu_counter),
        taxon = taxon,
        cluster_size = length(cluster_members),
        representative_seq = representative,
        stringsAsFactors = FALSE
      )
      
      otu_counter <- otu_counter + 1
    }
  }
  
  # Combine results
  if (length(otu_assignments_list) == 0) {
    warning("No OTUs created")
    return(data.frame(
      asv_id = character(),
      otu_id = character(),
      taxon = character(),
      cluster_size = integer(),
      representative_seq = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  otu_assignments <- dplyr::bind_rows(otu_assignments_list)
  
  # Save OTU table
  output_file <- file.path(output_dir, "otu_table.tsv")
  readr::write_tsv(otu_assignments, output_file)
  message(sprintf("Created %d OTUs from %d ASVs", 
                 length(unique(otu_assignments$otu_id)),
                 nrow(otu_assignments)))
  message("OTU table saved to: ", output_file)
  
  # Create OTU representative sequences FASTA
  create_otu_fasta(otu_assignments, asv_sequences, output_dir)
  
  return(otu_assignments)
}

#' Read FASTA file
#' @keywords internal
read_fasta <- function(fasta_file) {
  lines <- readLines(fasta_file)
  
  # Find header lines
  header_idx <- which(startsWith(lines, ">"))
  
  if (length(header_idx) == 0) {
    stop("No sequences found in FASTA file")
  }
  
  sequences <- list()
  
  for (i in seq_along(header_idx)) {
    # Extract sequence ID
    header <- lines[header_idx[i]]
    seq_id <- sub("^>([^ ]+).*", "\\1", header)
    
    # Extract sequence
    start <- header_idx[i] + 1
    end <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(lines)
    
    seq <- paste(lines[start:end], collapse = "")
    sequences[[seq_id]] <- seq
  }
  
  return(sequences)
}

#' Compute pairwise sequence similarities
#' @keywords internal
compute_sequence_similarities <- function(sequences) {
  n <- length(sequences)
  seq_names <- names(sequences)
  
  # Initialize similarity matrix
  sim_matrix <- matrix(0, nrow = n, ncol = n)
  rownames(sim_matrix) <- seq_names
  colnames(sim_matrix) <- seq_names
  
  # Diagonal is 1 (identity)
  diag(sim_matrix) <- 1
  
  # Compute pairwise similarities
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      sim <- compute_similarity(sequences[[i]], sequences[[j]])
      sim_matrix[i, j] <- sim
      sim_matrix[j, i] <- sim
    }
  }
  
  return(sim_matrix)
}

#' Compute similarity between two sequences
#' @keywords internal
compute_similarity <- function(seq1, seq2) {
  # Simple identity calculation
  seq1_chars <- strsplit(seq1, "")[[1]]
  seq2_chars <- strsplit(seq2, "")[[1]]
  
  # Handle different lengths
  min_len <- min(length(seq1_chars), length(seq2_chars))
  max_len <- max(length(seq1_chars), length(seq2_chars))
  
  if (max_len == 0) return(0)
  
  # Count matches
  matches <- sum(seq1_chars[1:min_len] == seq2_chars[1:min_len])
  
  # Calculate similarity (matches / max length)
  similarity <- matches / max_len
  
  return(similarity)
}

#' Perform hierarchical clustering
#' @keywords internal
perform_clustering <- function(similarity_matrix, cutoff, method = "average") {
  n <- nrow(similarity_matrix)
  
  if (n == 1) {
    # Single sequence
    clusters <- c(1)
    names(clusters) <- rownames(similarity_matrix)[1]
    return(clusters)
  }
  
  # Convert similarity to distance
  dist_matrix <- as.dist(1 - similarity_matrix)
  
  # Perform hierarchical clustering
  hc <- hclust(dist_matrix, method = method)
  
  # Cut tree at the specified cutoff
  # cutoff is similarity, so distance threshold is (1 - cutoff)
  height_cutoff <- 1 - cutoff
  clusters <- cutree(hc, h = height_cutoff)
  
  return(clusters)
}

#' Create FASTA file with OTU representative sequences
#' @keywords internal
create_otu_fasta <- function(otu_assignments, asv_sequences, output_dir) {
  # Get unique representative sequences
  representatives <- otu_assignments %>%
    dplyr::select(otu_id, representative_seq, taxon) %>%
    dplyr::distinct()
  
  # Create FASTA output
  fasta_lines <- c()
  
  for (i in seq_len(nrow(representatives))) {
    otu_id <- representatives$otu_id[i]
    rep_seq_id <- representatives$representative_seq[i]
    taxon <- representatives$taxon[i]
    
    if (rep_seq_id %in% names(asv_sequences)) {
      header <- sprintf(">%s taxon=%s representative=%s", otu_id, taxon, rep_seq_id)
      sequence <- asv_sequences[[rep_seq_id]]
      
      fasta_lines <- c(fasta_lines, header, sequence)
    }
  }
  
  # Write FASTA file
  output_file <- file.path(output_dir, "otu_representatives.fasta")
  writeLines(fasta_lines, output_file)
  message("OTU representative sequences saved to: ", output_file)
}
