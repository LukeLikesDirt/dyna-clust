#!/usr/bin/env Rscript
# ==============================================================================
# FASTA Header Formatting for Classification Script
# Description: Format FASTA headers to create a classification file
#              and modify headers to match DNABARCODER requirements
# ==============================================================================

# Check for required packages
required_packages <- c("dplyr", "readr", "stringr", "furrr", "future")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop("ERROR: Missing required packages: ", paste(missing_packages, collapse = ", "), 
       "\nInstall them with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))")
}

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(furrr)
  library(future)
})

# Input and output file paths
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments provided
if (length(args) != 3) {
  stop("Usage: Rscript reformat_eukaryome.R <input_fasta> <output_fasta> <classification_output>")
}

# Parse arguments
fasta_in <- args[1]
fasta_out <- args[2]
classification_out <- args[3]

# Check if input file exists
if (!file.exists(fasta_in)) {
  stop("ERROR: Input file does not exist: ", fasta_in)
}

# Check if output directories exist, create if needed
output_dir <- dirname(fasta_out)
if (!dir.exists(output_dir) && output_dir != ".") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

classification_dir <- dirname(classification_out)
if (!dir.exists(classification_dir) && classification_dir != ".") {
  dir.create(classification_dir, recursive = TRUE, showWarnings = FALSE)
}

# Increase global size limit for large FASTA files in future
options(future.globals.maxSize = 2000 * 1024^2)  # 2 GB limit

# ==============================================================================
# OPTIMIZED PARALLEL FASTA READING
# ==============================================================================

cat("Reading FASTA file...\n")
fasta_lines <- readLines(fasta_in)

cat("Extracting headers and sequences...\n")
header_indices <- which(grepl("^>", fasta_lines))

# Set up parallel processing
# Use multicore on Linux/HPC (memory efficient via forking)
plan(multicore, workers = availableCores() - 1)

cat("Processing", length(header_indices), "sequences in parallel...\n")

# Parallel processing - creates list of tibbles, then combines
fasta_df <- future_map_dfr(seq_along(header_indices), function(i) {
  current_header_idx <- header_indices[i]
  
  # Determine sequence end
  if (i < length(header_indices)) {
    seq_end_idx <- header_indices[i + 1] - 1
  } else {
    seq_end_idx <- length(fasta_lines)
  }
  
  # Get header and sequence
  tibble(
    header = fasta_lines[current_header_idx],
    sequence = paste(fasta_lines[(current_header_idx + 1):seq_end_idx], collapse = "")
  )
}, .progress = TRUE, .options = furrr_options(seed = TRUE))

# ==============================================================================
# TAXONOMIC CLASSIFICATION PROCESSING
# ==============================================================================

cat("\nProcessing taxonomic classifications...\n")

# Process headers to create classification data
classification_df <- fasta_df %>%
  # Extract ID (everything before the first semicolon)
  mutate(id = str_extract(header, "(?<=^>)[^;]+")) %>%
  # Extract taxonomic information
  mutate(
    kingdom = str_extract(header, "k__([^;]+)") %>% str_remove("k__"),
    phylum = str_extract(header, "p__([^;]+)") %>% str_remove("p__"),
    class = str_extract(header, "c__([^;]+)") %>% str_remove("c__"),
    order = str_extract(header, "o__([^;]+)") %>% str_remove("o__"),
    family = str_extract(header, "f__([^;]+)") %>% str_remove("f__"),
    genus = str_extract(header, "g__([^;]+)") %>% str_remove("g__"),
    species_epithet = str_extract(header, "s__([^;]+)$") %>% str_remove("s__")
  ) %>%
  # Mutate "-" to "." for consistent use of delimiters
  mutate(
    across(c(kingdom, phylum, class, order, family, genus, species_epithet), 
           ~str_replace_all(.x, "-", "."))
  ) %>%
  # Replace "unclassified" with "unidentified" across all taxonomic ranks
  mutate(across(c(kingdom, phylum, class, order, family, genus, species_epithet), 
                ~ifelse(.x == "unclassified", "unidentified", .x))) %>%
  # Replace EUKARYOME taxon predictions, "Incertae sedis", and tentative IDs to "unidentified"
  mutate(
    # Tentative IDs (cf. & nr.), affinities (aff.) and invalidly published names (nov.inval.)
    species_epithet = ifelse(
      str_detect(species_epithet, "cf\\.|nr\\.|aff\\.|nov\\.inval\\."), 
      "unidentified", 
      species_epithet
    ),
    # EUKARYOME predictions and Incertae sedis
    phylum = ifelse(grepl(".phy", phylum), "unidentified", phylum),
    class = ifelse(grepl(".cl", class), "unidentified", class),
    order = ifelse(grepl(".ord", order), "unidentified", order),
    family = ifelse(grepl(".fam", family), "unidentified", family),
    genus = ifelse(grepl(".gen", genus), "unidentified", genus)
  ) %>%
  # Reformat genus names that contain kingdom name in parentheses
  mutate(
    genus = case_when(
      # If parentheses exist
      grepl("\\(", genus) ~ {
        prefix <- str_extract(genus, "^[^\\(]+")  # Before "("
        inside <- str_extract(genus, "(?<=\\()[^\\)]+")  # Inside "()"
        
        # If there is no prefix, parentheses are a typo and the genus instead of the 
        # kingdom is inside the parentheses, so remove the parenthesis and keep 
        # only genus. Otherwise use: prefix_kng_inside
        ifelse(is.na(prefix) | prefix == "", 
               inside,
               paste0(prefix, "_kng_", inside))
      },
      # No parentheses - keep as is
      TRUE ~ genus
    )
  ) %>%
  # Create proper species names (genus + species epithet for identified species)
  mutate(
    species = case_when(
      species_epithet == "unidentified" ~ "unidentified",
      is.na(species_epithet) | species_epithet == "" ~ "unidentified",
      genus == "unidentified" ~ "unidentified",
      !is.na(genus) & !is.na(species_epithet) & 
        species_epithet != "unidentified" & genus != "unidentified" ~ 
        paste(genus, species_epithet, sep = " "),
      TRUE ~ "unidentified"
    )
  ) %>%
  # Calculate taxonomic completeness score
  mutate(
    tax_score = 
      (!is.na(kingdom) & kingdom != "" & kingdom != "unidentified") +
      (!is.na(phylum) & phylum != "" & phylum != "unidentified") +
      (!is.na(class) & class != "" & class != "unidentified") +
      (!is.na(order) & order != "" & order != "unidentified") +
      (!is.na(family) & family != "" & family != "unidentified") +
      (!is.na(genus) & genus != "" & genus != "unidentified") +
      (!is.na(species) & species != "" & species != "unidentified")
  ) %>%
  # Select unique IDs (keeping highest taxonomic score)
  group_by(id) %>%
  arrange(desc(tax_score), id) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  # Replace any NA values with empty strings if needed
  mutate(across(everything(), ~ifelse(is.na(.x), "", .x)))

# Check for duplicated IDs and exit with a warning if found
if (any(duplicated(classification_df$id))) {
  dup_ids <- classification_df$id[duplicated(classification_df$id)]
  stop(paste("Error: Duplicate IDs found in classification data:", paste(dup_ids, collapse = ", ")))
}

# ==============================================================================
# WRITE OUTPUT FILES
# ==============================================================================

cat("Writing output files...\n")

# Create new FASTA file using header and sequence columns
fasta_content <- classification_df %>%
  mutate(fasta_header = paste0(">", id)) %>%
  select(fasta_header, sequence) %>%
  tidyr::pivot_longer(everything(), names_to = NULL, values_to = "line") %>%
  pull(line)

writeLines(fasta_content, fasta_out)

# Write the classification file
classification_df %>%
  select(id, kingdom, phylum, class, order, family, genus, species) %>%
  write_tsv(classification_out)

# ==============================================================================
# PRINT SUMMARY
# ==============================================================================

cat("\n")
cat("Processing complete!\n")
cat("Input FASTA:", fasta_in, "\n")
cat("Output FASTA:", fasta_out, "\n")
cat("Classification file:", classification_out, "\n")
cat("Number of sequences processed:", nrow(classification_df), "\n\n")
cat("   Phyla identified:", sum(classification_df$phylum != "unidentified"), "\n")
cat("   Phyla unidentified:", sum(classification_df$phylum == "unidentified"), "\n\n")
cat("   Classes identified:", sum(classification_df$class != "unidentified"), "\n")
cat("   Classes unidentified:", sum(classification_df$class == "unidentified"), "\n\n")
cat("   Orders identified:", sum(classification_df$order != "unidentified"), "\n")
cat("   Orders unidentified:", sum(classification_df$order == "unidentified"), "\n\n")
cat("   Families identified:", sum(classification_df$family != "unidentified"), "\n")
cat("   Families unidentified:", sum(classification_df$family == "unidentified"), "\n\n")
cat("   Genera identified:", sum(classification_df$genus != "unidentified"), "\n")
cat("   Genera unidentified:", sum(classification_df$genus == "unidentified"), "\n\n")
cat("   Species identified:", sum(classification_df$species != "unidentified"), "\n")
cat("   Species unidentified:", sum(classification_df$species == "unidentified"), "\n\n")