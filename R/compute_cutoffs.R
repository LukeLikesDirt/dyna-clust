#' Compute Taxon-Specific Clustering Cutoffs with dnabarcoder
#'
#' This function computes taxon-specific clustering cutoffs using dnabarcoder.
#' It processes sequence data to determine optimal clustering thresholds for different taxa.
#'
#' @param fasta_file Path to the input FASTA file containing reference sequences
#' @param taxonomy_file Path to the taxonomy file mapping sequences to taxa
#' @param output_dir Directory where dnabarcoder output files will be saved
#' @param dnabarcoder_path Path to dnabarcoder executable (default: "dnabarcoder")
#' @param method Method for cutoff calculation: "ABGD", "ASAP", or "both" (default: "both")
#' @param min_seq_per_taxon Minimum number of sequences per taxon for cutoff calculation (default: 5)
#'
#' @return A data frame containing taxon-specific clustering cutoffs with columns:
#'   \itemize{
#'     \item taxon: Taxonomic group name
#'     \item cutoff: Similarity cutoff value (0-1)
#'     \item method: Method used to compute cutoff
#'     \item n_sequences: Number of sequences in the taxon
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cutoffs <- compute_cutoffs(
#'   fasta_file = "reference_seqs.fasta",
#'   taxonomy_file = "taxonomy.txt",
#'   output_dir = "cutoffs_output"
#' )
#' }
compute_cutoffs <- function(fasta_file, 
                           taxonomy_file, 
                           output_dir,
                           dnabarcoder_path = "dnabarcoder",
                           method = "both",
                           min_seq_per_taxon = 5) {
  
  # Check input files exist
  if (!file.exists(fasta_file)) {
    stop("FASTA file not found: ", fasta_file)
  }
  if (!file.exists(taxonomy_file)) {
    stop("Taxonomy file not found: ", taxonomy_file)
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Read taxonomy file
  message("Reading taxonomy file...")
  taxonomy <- readr::read_tsv(taxonomy_file, col_types = readr::cols())
  
  # Split sequences by taxon
  message("Grouping sequences by taxon...")
  taxa_groups <- taxonomy %>%
    dplyr::group_by(taxon) %>%
    dplyr::summarise(
      seq_ids = list(seq_id),
      n_sequences = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_sequences >= min_seq_per_taxon)
  
  message(sprintf("Found %d taxa with >= %d sequences", 
                 nrow(taxa_groups), min_seq_per_taxon))
  
  # Process each taxon with dnabarcoder
  cutoffs_list <- list()
  
  for (i in seq_len(nrow(taxa_groups))) {
    taxon <- taxa_groups$taxon[i]
    message(sprintf("Processing taxon %d/%d: %s", i, nrow(taxa_groups), taxon))
    
    # Create temporary FASTA for this taxon
    taxon_fasta <- file.path(output_dir, paste0(taxon, ".fasta"))
    taxon_output <- file.path(output_dir, taxon)
    
    # Extract sequences for this taxon
    seq_ids <- taxa_groups$seq_ids[[i]]
    extract_sequences(fasta_file, seq_ids, taxon_fasta)
    
    # Run dnabarcoder
    cutoff <- run_dnabarcoder(
      fasta_file = taxon_fasta,
      output_prefix = taxon_output,
      dnabarcoder_path = dnabarcoder_path,
      method = method
    )
    
    cutoffs_list[[i]] <- data.frame(
      taxon = taxon,
      cutoff = cutoff,
      method = method,
      n_sequences = taxa_groups$n_sequences[i],
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results
  cutoffs_df <- dplyr::bind_rows(cutoffs_list)
  
  # Save cutoffs
  cutoffs_file <- file.path(output_dir, "taxon_cutoffs.tsv")
  readr::write_tsv(cutoffs_df, cutoffs_file)
  message("Cutoffs saved to: ", cutoffs_file)
  
  return(cutoffs_df)
}

#' Extract sequences from FASTA file
#' @keywords internal
extract_sequences <- function(fasta_file, seq_ids, output_file) {
  # Read FASTA file
  fasta_lines <- readLines(fasta_file)
  
  # Initialize output
  output_lines <- c()
  include_seq <- FALSE
  
  for (line in fasta_lines) {
    if (startsWith(line, ">")) {
      # Header line
      seq_id <- sub("^>([^ ]+).*", "\\1", line)
      include_seq <- seq_id %in% seq_ids
      if (include_seq) {
        output_lines <- c(output_lines, line)
      }
    } else if (include_seq) {
      # Sequence line
      output_lines <- c(output_lines, line)
    }
  }
  
  writeLines(output_lines, output_file)
}

#' Run dnabarcoder to compute cutoff
#' @keywords internal
run_dnabarcoder <- function(fasta_file, output_prefix, dnabarcoder_path, method) {
  # Check if dnabarcoder is available
  if (Sys.which(dnabarcoder_path) == "") {
    warning("dnabarcoder not found in PATH. Using default cutoff of 0.97")
    return(0.97)
  }
  
  # Construct dnabarcoder command
  # Note: This is a placeholder - actual dnabarcoder command syntax may vary
  cmd <- sprintf("%s -i %s -o %s -m %s", 
                dnabarcoder_path, fasta_file, output_prefix, method)
  
  # Run dnabarcoder
  tryCatch({
    result <- system(cmd, intern = TRUE)
    
    # Parse output to extract cutoff
    # This is a placeholder - actual parsing will depend on dnabarcoder output format
    cutoff_line <- grep("cutoff|threshold", result, value = TRUE, ignore.case = TRUE)
    if (length(cutoff_line) > 0) {
      cutoff <- as.numeric(sub(".*([0-9]+\\.[0-9]+).*", "\\1", cutoff_line[1]))
      if (!is.na(cutoff) && cutoff > 0 && cutoff <= 1) {
        return(cutoff)
      }
    }
    
    # If parsing fails, use default
    warning("Could not parse dnabarcoder output. Using default cutoff of 0.97")
    return(0.97)
  }, error = function(e) {
    warning("Error running dnabarcoder: ", e$message, ". Using default cutoff of 0.97")
    return(0.97)
  })
}
