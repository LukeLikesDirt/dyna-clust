# Script name: 02.Check_infraspecific_annotations.R
# Description: Check and re-annotate species epithet to handle infraspecific 
#              annotations in the EUKARYOME reference database.
# Note:        This script must be run from the source directory.

# Required packages
require(data.table)
require(tidyverse)

# Read in the classification file
taxa <- fread("./data/eukaryome_ITS_v2.0.classification")

# Subset to species column and handle front-end issues
species <- taxa %>%
  select(id, species, genus) %>%  # Keep genus for reconstruction later
  mutate(
    
    # Handle "ë" formatting issue (assuming +1/2 is meant to be ë)
    species = str_replace_all(species, "\\+.1/2.", "e"),
    
    # Handle "thaliana.x.arenosa.x.arenosa" typo
    species = str_replace_all(species, "thaliana\\.x\\.arenosa\\.x\\.arenosa", "thaliana.x.arenosa"),
    
    # Handle "metallica.x.Begonia.sanguinea" typo
    species = str_replace_all(species, "metallica\\.x\\.Begonia\\.sanguinea", "metallica.x.sanguinea"),
    
    # Extract infraspecific annotation (if present)
    infraspecific = case_when(
      str_detect(species, "\\.subf\\.") ~ str_extract(species, "\\.subf\\.[^\\s]+"),
      str_detect(species, "\\.f\\.sp\\.") ~ str_extract(species, "\\.f\\.sp\\.[^\\s]+"),
      str_detect(species, "\\.f\\.") ~ str_extract(species, "\\.f\\.[^\\s]+"),
      str_detect(species, "\\.var\\.") ~ str_extract(species, "\\.var\\.[^\\s]+"),
      str_detect(species, "\\.subsp\\.") ~ str_extract(species, "\\.subsp\\.[^\\s]+"),
      TRUE ~ NA_character_
    ),
    
    # Remove leading period from infraspecific annotation
    infraspecific = str_remove(infraspecific, "^\\."),
    
    # Remove infraspecific annotation from species name
    species = case_when(
      !is.na(infraspecific) ~ str_remove(species, "\\.subf\\.[^\\s]+|\\.f\\.sp\\.[^\\s]+|\\.f\\.[^\\s]+|\\.var\\.[^\\s]+|\\.subsp\\.[^\\s]+"),
      TRUE ~ species
    ),
    
    # Clean up any trailing/leading whitespace
    species = str_trim(species),
    infraspecific = str_trim(infraspecific)
  )

# Show summary of infraspecific annotations found
cat("Infraspecific annotations summary:\n")
print(table(species$infraspecific, useNA = "ifany"))
cat("\n")

# Generate a unique species vector
unique_species <- unique(species$species)

# Split species into genus and species epithet
species_split <- tibble(species_name = unique_species) %>%
  filter(species_name != "unidentified") %>%
  mutate(
    genus = str_extract(species_name, "^\\S+"),  # First word
    species_epithet = str_remove(species_name, "^\\S+\\s*")  # Everything after first word
  )

# Pull all species with special characters in the species epithet
# (not equal to "_" or alphanumeric)
infraspecific_candidates <- species_split %>%
  filter(str_detect(species_epithet, "[^A-Za-z0-9_]")) %>%
  pull(species_name)

# Check the first 1000 entries
cat("Species with potential infraspecific annotations (first 1000):\n")
print(head(infraspecific_candidates, 1000))

# Check the second 1000 entries
cat("\nSpecies with potential infraspecific annotations (second 1000):\n")
print(head(infraspecific_candidates[-(1:1000)], 1000))

# Address any issues with manual corrections
species <- species %>%
  mutate(
    # Apply manual corrections
    species = case_when(
      species == "Oryza sativa.Indica" ~ "Oryza sativa",
      species == "Oryza sativa.Japonica" ~ "Oryza sativa",
      species == "Marchandiomyces nt59.1784" ~ "Marchandiomyces unidentified",
      TRUE ~ species
    ),
    # Update infraspecific annotations for corrected species
    infraspecific = case_when(
      species == "Oryza sativa" & str_detect(lag(species, default = ""), "sativa\\.Indica") ~ "subsp.Indica",
      species == "Oryza sativa" & str_detect(lag(species, default = ""), "sativa\\.Japonica") ~ "subsp.Japonica",
      TRUE ~ infraspecific
    )
  )

# How many species with potential infraspecific annotations?
cat("Total unique species:", length(unique_species) - sum(unique_species == "unidentified"), "\n")
cat("Species with potential infraspecific annotations:", length(infraspecific_candidates), "\n\n")

# Extract abbreviations between periods from species epithet only
unique_abbreviations <- species_split %>%
  filter(species_name %in% infraspecific_candidates) %>%
  pull(species_epithet) %>%
  str_extract("\\.[^.]+(\\.([^.]+))?\\.") %>% # Matches .x. or .x.y.
  .[!is.na(.)] %>%
  unique(.) %>%
  sort()

cat("Unique abbreviation patterns found:\n")
print(unique_abbreviations)
cat("\n")

# RECONSTRUCT CLASSIFICATION FILE WITH UPDATED TAXONOMY ------------------------

cat("Reconstructing classification file with updated taxonomy...\n")

# Join the updated species information back to the original taxa
taxa_updated <- taxa %>%
  select(-species) %>%  # Remove old species column
  left_join(
    species %>% select(id, species, infraspecific),
    by = "id"
  ) %>%
  # Reorder columns
  select(id, kingdom, phylum, class, order, family, genus, species, infraspecific)

# Write updated classification file
fwrite(taxa_updated, "./data/eukaryome_ITS_v2.0.classification", sep = "\t")
