#' Classify ASVs Using BLAST and Taxon-Specific Cutoffs
#'
#' This function classifies Amplicon Sequence Variants (ASVs) by comparing them
#' to a reference database using BLAST, then applying taxon-specific similarity
#' cutoffs to assign taxonomy.
#'
#' @param asv_fasta Path to FASTA file containing ASV sequences to classify
#' @param reference_db Path to BLAST database or FASTA file to use as reference
#' @param cutoffs Data frame with taxon-specific cutoffs (output from compute_cutoffs)
#' @param output_dir Directory for output files
#' @param blast_path Path to BLAST executable (default: "blastn")
#' @param num_threads Number of threads for BLAST (default: 1)
#' @param max_target_seqs Maximum number of BLAST hits to consider (default: 100)
#' @param evalue E-value threshold for BLAST (default: 1e-10)
#'
#' @return A data frame with ASV classifications containing:
#'   \itemize{
#'     \item asv_id: ASV identifier
#'     \item taxon: Assigned taxonomic group
#'     \item similarity: Percent identity to best match
#'     \item reference_id: ID of best matching reference sequence
#'     \item cutoff_used: Cutoff threshold applied
#'     \item classification_status: "confident" or "uncertain"
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' classifications <- classify_asvs(
#'   asv_fasta = "asvs.fasta",
#'   reference_db = "reference.fasta",
#'   cutoffs = cutoffs_df,
#'   output_dir = "blast_output"
#' )
#' }
classify_asvs <- function(asv_fasta,
                         reference_db,
                         cutoffs,
                         output_dir,
                         blast_path = "blastn",
                         num_threads = 1,
                         max_target_seqs = 100,
                         evalue = 1e-10) {
  
  # Check inputs
  if (!file.exists(asv_fasta)) {
    stop("ASV FASTA file not found: ", asv_fasta)
  }
  if (!file.exists(reference_db) && !file.exists(paste0(reference_db, ".nhr"))) {
    stop("Reference database not found: ", reference_db)
  }
  if (!is.data.frame(cutoffs) || !all(c("taxon", "cutoff") %in% names(cutoffs))) {
    stop("cutoffs must be a data frame with 'taxon' and 'cutoff' columns")
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Check if BLAST database exists, if not create it
  if (!file.exists(paste0(reference_db, ".nhr"))) {
    message("Creating BLAST database...")
    create_blast_db(reference_db, blast_path)
  }
  
  # Run BLAST
  message("Running BLAST search...")
  blast_output <- file.path(output_dir, "blast_results.txt")
  run_blast(
    query = asv_fasta,
    database = reference_db,
    output = blast_output,
    blast_path = blast_path,
    num_threads = num_threads,
    max_target_seqs = max_target_seqs,
    evalue = evalue
  )
  
  # Parse BLAST results
  message("Parsing BLAST results...")
  blast_results <- parse_blast_results(blast_output)
  
  if (nrow(blast_results) == 0) {
    warning("No BLAST hits found")
    return(data.frame(
      asv_id = character(),
      taxon = character(),
      similarity = numeric(),
      reference_id = character(),
      cutoff_used = numeric(),
      classification_status = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Get reference taxonomy (assuming reference IDs contain taxon info)
  message("Assigning taxonomy based on cutoffs...")
  classifications <- apply_cutoffs(blast_results, cutoffs)
  
  # Save classifications
  output_file <- file.path(output_dir, "asv_classifications.tsv")
  readr::write_tsv(classifications, output_file)
  message("Classifications saved to: ", output_file)
  
  return(classifications)
}

#' Create BLAST database
#' @keywords internal
create_blast_db <- function(fasta_file, blast_path) {
  makeblastdb_path <- file.path(dirname(blast_path), "makeblastdb")
  if (Sys.which(makeblastdb_path) == "" && Sys.which("makeblastdb") != "") {
    makeblastdb_path <- "makeblastdb"
  }
  
  cmd <- sprintf("%s -in %s -dbtype nucl -parse_seqids", 
                makeblastdb_path, fasta_file)
  
  result <- system(cmd, intern = TRUE)
  if (length(result) > 0 && any(grepl("Error|error", result))) {
    stop("Error creating BLAST database: ", paste(result, collapse = "\n"))
  }
}

#' Run BLAST search
#' @keywords internal
run_blast <- function(query, database, output, blast_path, 
                     num_threads, max_target_seqs, evalue) {
  
  # Custom output format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
  outfmt <- "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
  
  cmd <- sprintf(
    "%s -query %s -db %s -out %s -outfmt '%s' -num_threads %d -max_target_seqs %d -evalue %e",
    blast_path, query, database, output, outfmt, num_threads, max_target_seqs, evalue
  )
  
  result <- system(cmd, intern = TRUE)
  
  if (!file.exists(output)) {
    stop("BLAST output file not created")
  }
}

#' Parse BLAST results
#' @keywords internal
parse_blast_results <- function(blast_output) {
  if (!file.exists(blast_output) || file.info(blast_output)$size == 0) {
    return(data.frame(
      qseqid = character(),
      sseqid = character(),
      pident = numeric(),
      length = integer(),
      evalue = numeric(),
      bitscore = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  blast_results <- readr::read_tsv(
    blast_output,
    col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", 
                 "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"),
    col_types = readr::cols(
      qseqid = readr::col_character(),
      sseqid = readr::col_character(),
      pident = readr::col_double(),
      length = readr::col_integer(),
      mismatch = readr::col_integer(),
      gapopen = readr::col_integer(),
      qstart = readr::col_integer(),
      qend = readr::col_integer(),
      sstart = readr::col_integer(),
      send = readr::col_integer(),
      evalue = readr::col_double(),
      bitscore = readr::col_double()
    )
  )
  
  return(blast_results)
}

#' Apply taxon-specific cutoffs to BLAST results
#' @keywords internal
apply_cutoffs <- function(blast_results, cutoffs) {
  # Get best hit for each ASV
  best_hits <- blast_results %>%
    dplyr::group_by(qseqid) %>%
    dplyr::slice_max(bitscore, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Extract taxon from reference ID (assuming format like "taxon|seqid" or similar)
  # This is a placeholder - actual extraction will depend on reference format
  best_hits <- best_hits %>%
    dplyr::mutate(
      ref_taxon = extract_taxon_from_id(sseqid)
    )
  
  # Join with cutoffs
  classifications <- best_hits %>%
    dplyr::left_join(cutoffs, by = c("ref_taxon" = "taxon")) %>%
    dplyr::mutate(
      similarity = pident / 100,
      cutoff_used = dplyr::coalesce(cutoff, 0.97),  # Default if no cutoff found
      classification_status = dplyr::if_else(
        similarity >= cutoff_used,
        "confident",
        "uncertain"
      ),
      taxon = dplyr::if_else(
        classification_status == "confident",
        ref_taxon,
        "Unclassified"
      )
    ) %>%
    dplyr::select(
      asv_id = qseqid,
      taxon,
      similarity,
      reference_id = sseqid,
      cutoff_used,
      classification_status
    )
  
  return(classifications)
}

#' Extract taxon from reference sequence ID
#' @keywords internal
extract_taxon_from_id <- function(seq_ids) {
  # Attempt multiple common formats
  # Format 1: taxon|seqid
  taxon <- stringr::str_extract(seq_ids, "^([^|]+)\\|")
  taxon <- stringr::str_remove(taxon, "\\|$")
  
  # Format 2: seqid_taxon or seqid-taxon
  if (all(is.na(taxon))) {
    taxon <- stringr::str_extract(seq_ids, "[_-]([A-Za-z]+)$")
    taxon <- stringr::str_remove_all(taxon, "[_-]")
  }
  
  # Format 3: just return first part before space
  if (all(is.na(taxon))) {
    taxon <- stringr::str_extract(seq_ids, "^[^ ]+")
  }
  
  # Default to "Unknown" if still NA
  taxon[is.na(taxon)] <- "Unknown"
  
  return(taxon)
}

