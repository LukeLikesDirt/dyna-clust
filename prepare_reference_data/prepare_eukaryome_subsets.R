#!/usr/bin/env Rscript
# ==============================================================================`
# Prepare datasets for prediction
# Description: Produce ID lists used to subset sequences for DNAbarcoder
#              similarity-based taxonomic predictions.
#
# STEP 1 - Unique sequences:
#   Remove duplicate sequences at each rank, keeping one representative per
#   taxonomic group based on abundance. Outputs one ID file per rank.
#
# STEP 2 - Prediction ID lists:
#   For each target rank (e.g. species), evaluate within every valid parent
#   rank above it (genus, family, ... kingdom), applying four filters per
#   parent-taxon chunk in this order:
#     1. min_subgroups  - checked on raw child-taxon count before any removal
#     2. Ambiguous taxa - remove unclassified / sp. / incertae sedis children
#     3. max_proportion - cap the dominant child taxon to <= max_proportion
#     4. min_sequences  - checked after capping; drop chunk if too few remain
#     5. max_sequences  - balanced round-robin downsample if still too many
#   Outputs one ID file per (target x parent) combination, e.g.
#     species_pred_id_gen.txt, species_pred_id_fam.txt, ...
#     genus_pred_id_fam.txt,   genus_pred_id_ord.txt,   ...
#
# Usage:
#   Rscript prepare_eukaryome_subsets.R <fasta_in> <classification_in> \
#       <output_dir> <min_subgroups> <min_sequences> <max_sequences> \
#       <max_proportion>
#
# All outputs are plain-text ID files (one sequence ID per line) so the
# global FASTA and classification table can be filtered at prediction time,
# keeping intermediate file sizes small.
# ==============================================================================`

# Check for required packages before attempting to load
required_packages <- c("readr", "Biostrings", "dplyr", "data.table")
missing_packages  <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "ERROR: Missing required packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them with: install.packages(c('",
    paste(missing_packages, collapse = "', '"), "'))"
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(Biostrings)
  library(dplyr)
  library(data.table)
})

# ==============================================================================`
# ARGUMENT PARSING -------------------------------------------------------------
# ==============================================================================`

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
  stop(
    "Usage: Rscript prepare_eukaryome_subsets.R ",
    "<fasta_in> <classification_in> <output_dir> ",
    "<min_subgroups> <min_sequences> <max_sequences> <max_proportion>"
  )
}

fasta_in          <- args[1]
classification_in <- args[2]
output_dir        <- args[3]
min_subgroups     <- as.integer(args[4])
min_sequences     <- as.integer(args[5])
max_sequences     <- as.integer(args[6])
max_proportion    <- as.numeric(args[7])

# ==============================================================================`
# INPUT AND CONSTANT VALIDATION ------------------------------------------------
# ==============================================================================`

if (!file.exists(fasta_in)) {
  stop("FASTA file not found: ", fasta_in)
}
if (!file.exists(classification_in)) {
  stop("Classification file not found: ", classification_in)
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}
if (!is.numeric(min_subgroups) || min_subgroups <= 0 || min_subgroups != floor(min_subgroups)) {
  stop("min_subgroups must be a positive integer (e.g. 10)")
}
if (!is.numeric(min_sequences) || min_sequences <= 0 || min_sequences != floor(min_sequences)) {
  stop("min_sequences must be a positive integer (e.g. 30)")
}
if (!is.numeric(max_sequences) || max_sequences <= 0 || max_sequences != floor(max_sequences)) {
  stop("max_sequences must be a positive integer (e.g. 25000)")
}
if (is.na(max_proportion) || max_proportion <= 0 || max_proportion >= 1) {
  stop("max_proportion must be between 0 and 1 (e.g. 0.5 for 50%)")
}

cat("=== PARAMETERS ===\n")
cat(sprintf("  fasta_in          : %s\n", fasta_in))
cat(sprintf("  classification_in : %s\n", classification_in))
cat(sprintf("  output_dir        : %s\n", output_dir))
cat(sprintf("  min_subgroups     : %d\n", min_subgroups))
cat(sprintf("  min_sequences     : %d\n", min_sequences))
cat(sprintf("  max_sequences     : %d\n", max_sequences))
cat(sprintf("  max_proportion    : %.2f\n", max_proportion))
cat("\n")

# ==============================================================================`
# RANK METADATA ----------------------------------------------------------------
# ==============================================================================`

# Ordered low-to-high; used to derive parent rank relationships.
rank_hierarchy <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

# Short suffixes used in output file names (e.g. species_pred_id_gen.txt).
rank_abbr <- c(
  kingdom = "kng", phylum = "phy", class = "cls", order = "ord",
  family  = "fam", genus  = "gen", species = "spe"
)

# All valid parent ranks for each target rank (every rank above it).
# kingdom has no parents so produces no prediction files.
parent_ranks_map <- list(
  species = c("genus", "family", "order", "class", "phylum", "kingdom"),
  genus   = c("family", "order", "class", "phylum", "kingdom"),
  family  = c("order", "class", "phylum", "kingdom"),
  order   = c("class", "phylum", "kingdom"),
  class   = c("phylum", "kingdom"),
  phylum  = c("kingdom"),
  kingdom = character(0)
)

# ==============================================================================`
# FUNCTIONS --------------------------------------------------------------------
# ==============================================================================`

# ------------------------------------------------------------------------------`
# filter_unique_sequences
# ------------------------------------------------------------------------------`
# Removes duplicate sequences at a given taxonomic rank by keeping the single
# most abundant representative within each duplicate group. Sequences sharing
# an identical string are grouped; within each group the taxon with the most
# sequences is kept and all others are removed. If no duplicates exist the
# input is returned unchanged.
#
# Args:
#   classification_df : data.frame with columns id + all seven rank columns
#   fasta_seqs        : DNAStringSet named by sequence ID
#   rank              : target rank at which to deduplicate (default "species")
#
# Returns: named list with two elements - <rank>_df and <rank>_seqs

filter_unique_sequences <- function(classification_df, fasta_seqs, rank = "species") {
  
  valid_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  if (!rank %in% valid_ranks) {
    stop(paste("Invalid rank. Must be one of:", paste(valid_ranks, collapse = ", ")))
  }
  
  # All ranks from kingdom down to (and including) the target rank
  rank_hier <- valid_ranks[1:which(valid_ranks == rank)]
  
  # Retain only rows with a valid classification at the target rank
  rank_df <- classification_df %>%
    filter(
      !is.na(.data[[rank]]),
      .data[[rank]] != "",
      .data[[rank]] != "unclassified",
      .data[[rank]] != "unidentified"
    )
  
  # Group by the rank immediately above the target (or by kingdom itself)
  group_var <- if (rank == "kingdom") "kingdom" else rank_hier[length(rank_hier) - 1]
  rank_df   <- rank_df %>% group_by(across(all_of(group_var)))
  
  rank_seqs <- fasta_seqs[names(fasta_seqs) %in% rank_df$id]
  
  # Identify sequences that appear more than once (exact string match)
  duplicated_seqs <- tibble(
    id       = names(rank_seqs),
    sequence = as.character(rank_seqs)
  ) %>%
    group_by(sequence) %>%
    filter(n() > 1) %>%
    mutate(group_id = paste0("group_", cur_group_id())) %>%
    ungroup() %>%
    select(-sequence)
  
  # No duplicates - return as-is
  if (nrow(duplicated_seqs) == 0) {
    message(paste("No duplicated sequences found at", rank, "level"))
    result <- list(rank_df, rank_seqs)
    names(result) <- c(paste0(rank, "_df"), paste0(rank, "_seqs"))
    return(result)
  }
  
  # Within each duplicate group, keep the taxon with the highest sequence
  # count; if tied, keep the first occurrence
  unique_seqs <- duplicated_seqs %>%
    left_join(
      classification_df %>% select(id, all_of(rank_hier)),
      by = "id"
    ) %>%
    group_by(across(c(group_id, all_of(rank_hier)))) %>%
    mutate(n_seqs = n()) %>%
    filter(n_seqs == max(n_seqs)) %>%
    arrange(desc(n_seqs)) %>%
    ungroup() %>%
    group_by(group_id) %>%
    slice(1)
  
  seqs_to_remove <- setdiff(duplicated_seqs$id, unique_seqs$id)
  
  rank_df   <- rank_df %>% filter(!id %in% seqs_to_remove)
  rank_seqs <- rank_seqs[names(rank_seqs) %in% rank_df$id]
  
  result <- list(rank_df, rank_seqs)
  names(result) <- c(paste0(rank, "_df"), paste0(rank, "_seqs"))
  return(result)
}

# ------------------------------------------------------------------------------'
# balanced_downsample
# ------------------------------------------------------------------------------'
# Reduces a data.frame to n rows while keeping the taxon distribution as
# even as possible. Works by round-robin: assigns a within-taxon index to
# every row (1st sequence per taxon, 2nd, 3rd, ...), sorts by that index,
# then takes the first n rows. Rare taxa (with only 1-2 sequences) are
# therefore always included before any taxon receives a third draw.
#
# Args:
#   df          : data.frame to downsample (must contain target_rank column)
#   target_rank : column name of the rank being balanced
#   n           : number of rows to return
#
# Returns: data.frame with n rows

balanced_downsample <- function(df, target_rank, n) {
  df %>%
    ungroup() %>%
    group_by(!!sym(target_rank)) %>%
    mutate(.draw_order = row_number()) %>%        # position within each taxon
    ungroup() %>%
    arrange(.draw_order, !!sym(target_rank)) %>%  # round-robin sort
    slice_head(n = n) %>%
    select(-.draw_order)
}

# ------------------------------------------------------------------------------`
# nested_prediction_filter
# ------------------------------------------------------------------------------`
# Filters a data.frame for one (target_rank x parent_rank) combination by
# iterating over every parent-taxon chunk and applying five sequential filters.
# Returns a combined data.frame of passing rows, or NULL if nothing survives.
#
# Filter order within each chunk:
#   1. min_subgroups  - raw unique child-taxon count before ambiguous removal;
#                       drop chunk immediately if below threshold
#   2. Ambiguous taxa - remove unclassified / unidentified / sp. / incertae sedis
#   3. max_proportion - randomly subsample the dominant child taxon so it
#                       represents no more than max_proportion of the chunk
#   4. min_sequences  - drop chunk if too few sequences survive step 3
#   5. max_sequences  - balanced round-robin downsample if chunk is still large
#
# Parent-level ambiguous filtering (step 0) is applied to the whole data.frame
# before splitting, so chunks with uncertain parent placement are excluded
# regardless of how well-resolved their child taxa are.
#
# Args:
#   df             : data.frame from filter_unique_sequences for target_rank
#   target_rank    : rank being filtered (e.g. "species")
#   parent_rank    : rank defining the chunks (e.g. "family")
#   max_proportion : max dominant-taxon fraction per chunk (default 0.5)
#   min_subgroups  : min raw child taxa per chunk (default 10)
#   min_sequences  : min sequences per chunk after step 3 (default 30)
#   max_sequences  : max sequences per chunk (default 25000)
#
# Returns: data.frame of retained rows across all passing chunks, or NULL

nested_prediction_filter <- function(
    df,
    target_rank,
    parent_rank,
    max_proportion = 0.5,
    min_subgroups  = 10,
    min_sequences  = 30,
    max_sequences  = 25000
) {
  
  # Step 0: remove rows whose parent-rank classification is ambiguous.
  # A species can be well-resolved while its family is "unidentified";
  # such rows cannot be meaningfully grouped by parent rank.
  df <- df %>%
    filter(
      !is.na(!!sym(parent_rank)) &
        !!sym(parent_rank) != "" &
        !tolower(!!sym(parent_rank)) %in% c("unidentified", "unclassified") &
        !grepl("incertae sedis", !!sym(parent_rank), ignore.case = TRUE) &
        !grepl("sp\\.", !!sym(parent_rank))
    )
  
  if (nrow(df) == 0) return(NULL)
  
  df_split      <- split(df, df[[parent_rank]])
  chunk_results <- list()
  
  for (parent_taxon in names(df_split)) {
    
    chunk_raw <- df_split[[parent_taxon]]
    
    # Filter 1: min_subgroups on raw child count (before ambiguous removal).
    # Rejects sparse parent groups early without further computation.
    if (length(unique(chunk_raw[[target_rank]])) < min_subgroups) next
    
    # Filter 2: remove ambiguous child taxa within this chunk
    chunk <- chunk_raw %>%
      filter(
        !is.na(!!sym(target_rank)) &
          !!sym(target_rank) != "" &
          !tolower(!!sym(target_rank)) %in% c("unidentified", "unclassified") &
          !grepl("incertae sedis", !!sym(target_rank), ignore.case = TRUE) &
          !grepl("sp\\.", !!sym(target_rank))
      )
    
    if (nrow(chunk) == 0) next
    
    # Filter 3: max_proportion - cap the dominant child taxon.
    # All sequences from non-dominant taxa are kept; the dominant taxon is
    # randomly subsampled so it reaches at most max_proportion of the chunk.
    # Formula: max_allowed = floor((p / (1-p)) * n_smaller), where p = max_proportion
    rank_counts <- chunk %>%
      group_by(!!sym(target_rank)) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n))
    
    if (nrow(rank_counts) == 1) {
      # Only one child taxon - proportion filter not applicable
      selected_ids <- chunk$id
    } else {
      largest_taxon    <- rank_counts[[1, target_rank]]
      largest_count    <- rank_counts$n[1]
      smaller_count    <- sum(rank_counts$n[-1])
      max_from_largest <- floor((max_proportion / (1 - max_proportion)) * smaller_count)
      
      non_dominant_ids <- chunk %>%
        filter(!!sym(target_rank) != largest_taxon) %>%
        pull(id)
      largest_ids <- chunk %>%
        filter(!!sym(target_rank) == largest_taxon) %>%
        pull(id)
      
      selected_ids <- if (largest_count <= max_from_largest) {
        c(non_dominant_ids, largest_ids)           # already within limit
      } else {
        c(non_dominant_ids, sample(largest_ids, max_from_largest))
      }
    }
    
    chunk <- chunk %>% filter(id %in% selected_ids)
    
    # Filter 4: min_sequences after proportion cap
    if (nrow(chunk) < min_sequences) next
    
    # Filter 5: max_sequences - balanced round-robin downsample
    if (nrow(chunk) > max_sequences) {
      chunk <- balanced_downsample(chunk, target_rank, max_sequences)
    }
    
    chunk_results[[parent_taxon]] <- chunk
  }
  
  if (length(chunk_results) == 0) return(NULL)
  bind_rows(chunk_results)
}

# ==============================================================================`
# READ INPUT FILES -------------------------------------------------------------
# ==============================================================================`

cat("Reading classification file...\n")
classification_df <- fread(classification_in) %>%
  select(id, kingdom, phylum, class, order, family, genus, species)

cat("Reading FASTA file...\n")
fasta_seqs <- readDNAStringSet(fasta_in)

cat(sprintf("  Sequences loaded: %d\n\n", length(fasta_seqs)))

# ==============================================================================`
# STEP 1: Select unique sequences ----------------------------------------------
# ==============================================================================`
# For each rank, remove duplicate sequences and save the retained IDs.
# These ID files are used to subset the global FASTA for BLAST-based
# taxonomic assignment. Saving IDs rather than FASTA keeps file sizes small.

cat("STEP 1: Selecting unique sequences...\n")

rank_dfs <- list()

for (rank in rank_hierarchy) {
  
  result           <- filter_unique_sequences(classification_df, fasta_seqs, rank = rank)
  rank_dfs[[rank]] <- result[[paste0(rank, "_df")]]
  
  out_path <- sprintf("data/%s_unique_id.txt", rank)
  writeLines(rank_dfs[[rank]]$id, out_path)
  
  cat(sprintf("  %-10s  %6d unique sequences  ->  %s\n",
              rank, nrow(rank_dfs[[rank]]), out_path))
}

# ==============================================================================`
# STEP 2: Prepare prediction ID lists -----------------------------------------
# ==============================================================================`
# For each valid (target x parent) combination, apply nested_prediction_filter()
# and write the retained IDs to a text file. One file is produced per
# combination; combinations where no chunks survive are skipped.
#
# Output naming: <target>_pred_id_<parent_abbr>.txt
#   e.g.  species_pred_id_gen.txt  - species filtered within each genus
#         species_pred_id_fam.txt  - species filtered within each family
#         class_pred_id_phy.txt    - class filtered within each phylum

cat("\nSTEP 2: Preparing prediction ID lists...\n")

for (target_rank in names(parent_ranks_map)) {
  
  valid_parents <- parent_ranks_map[[target_rank]]
  
  if (length(valid_parents) == 0) {
    cat(sprintf("  %-10s  no valid parent ranks - skipping\n", target_rank))
    next
  }
  
  df_target <- rank_dfs[[target_rank]]
  
  for (parent_rank in valid_parents) {
    
    out_path <- sprintf("data/%s_pred_id_%s.txt",
                        target_rank, rank_abbr[[parent_rank]])
    
    result <- nested_prediction_filter(
      df             = df_target,
      target_rank    = target_rank,
      parent_rank    = parent_rank,
      max_proportion = max_proportion,
      min_subgroups  = min_subgroups,
      min_sequences  = min_sequences,
      max_sequences  = max_sequences
    )
    
    if (is.null(result)) {
      cat(sprintf("  %-10s within %-10s ->  no groups passed filters\n",
                  target_rank, parent_rank))
      next
    }
    
    writeLines(result$id, out_path)
    
    cat(sprintf("  %-10s within %-10s ->  %6d IDs  |  %d unique %s  ->  %s\n",
                target_rank, parent_rank,
                nrow(result),
                length(unique(result[[target_rank]])),
                target_rank,
                out_path))
  }
}

cat("\nDone. All ID files saved to data/\n")