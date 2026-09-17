# Chapter 1 Eco-Evo Phase 1: Soil microbial data setup
# Author: Jennah Brown
# Date: 9/17/2026
#
# Purpose: Prepare fungal ITS sequencing data for the Eco-Evo Phase 1 analysis.
#
# This workflow was adapted from microbial community analysis code developed  by Allison Simler-Williamson (ASW). The starting rarefied ITS phyloseq object was generated using the BRC-BIO ITS preprocessing pipeline.

library(phyloseq)
library(tidyverse)

#---------------
# Load and filter ITS sequencing data
#---------------

# Load the rarefied ITS phyloseq object containing BRC-BIO 2024-25 samples.
ps.all <- readRDS(
  "Data/phyloseq_rare_ITS_BRCBIO_2024-25.rds"
)

# Keep the Spring 2025 sequencing project used for the Eco-Evo analysis.
ps.spring <- subset_samples(
  ps.all,
  Project.ID == "BRCBIO Spring 2025"
)

# Remove taxa that are absent from the Spring 2025 subset.
ps.spring <- prune_taxa(taxa_sums(ps.spring) > 0, ps.spring)

# Keep the five soil collection sites used in Eco-Evo Phase 1.
ps.focal <- subset_samples(
  ps.spring,
  SoilID %in% c("BB", "HG", "NG", "OH", "SR")
)

# Extract sample metadata from the filtered phyloseq object.
ps.metadata <- data.frame(sample_data(ps.focal))


# Load BRTE biomass data for later matching with microbial metrics.
biomass <- read.csv(
  "Data/Eco-Evo Phase 1 Biomass Data Fall 2025.csv"
)

#---------------
# Calculate fungal alpha diversity
#---------------

# Calculate fungal alpha-diversity metrics for each sequencing sample. Extraction IDs are reformatted to match the sample metadata before joining.
rich <- estimate_richness(ps.focal) %>%
  rownames_to_column(var = "Extraction.ID") %>%
  mutate(
    Extraction.ID = as.character(str_sub(Extraction.ID, 2, 4))
  ) %>%
  left_join(ps.metadata, by = "Extraction.ID")

#---------------
# Calculate alpha diversity by replicate
#---------------

# Calculate average alpha diversity for each site, replicate, and soil source.
alpha.measured <- rich %>%
  group_by(SoilID, Replicate, Sample.Type) %>%
  summarize(
    alphadiv = mean(Chao1, na.rm = TRUE),
    .groups = "drop"
  )

#---------------
# Create alpha-diversity values for combined replicates
#---------------

# Phase 1 used A-B and B-C soil replicate combinations, while sequencing was done on individual replicates. Combined values are the average of the corresponding individual replicates.
combo.key <- tribble(
  ~Replicate, ~rep1, ~rep2,
  "A-B", "A", "B",
  "B-C", "B", "C"
)

alpha.combos <- combo.key %>%
  left_join(alpha.measured %>%
              select(SoilID, Sample.Type, rep1 = Replicate, alphadiv.1 = alphadiv),
            by = "rep1") %>%
  left_join(alpha.measured %>%
              select(SoilID, Sample.Type, rep2 = Replicate, alphadiv.2 = alphadiv),
            by = c("SoilID", "Sample.Type", "rep2")) %>%
  mutate(alphadiv = (alphadiv.1 + alphadiv.2) / 2) %>%
  select(SoilID, Replicate, Sample.Type, alphadiv)

# Combine measured and combo values.
alpha.df <- bind_rows(alpha.measured, alpha.combos)

# OH Interspace C and SR Rhizosphere C were not sequenced, so those C reps and SR rhizosphere B-C combinations remain NA.

#---------------
# Match alpha diversity to biomass data
#---------------

# Match soil source naming used in the biomass dataset.
alpha.df <- alpha.df %>%
  mutate(
    Rhizo.Int = recode(
      Sample.Type,
      "Interspace" = "INT",
      "Rhizosphere" = "RHIZO"
    )
  ) %>%
  select(SoilID, Replicate, Rhizo.Int, alphadiv)

# Add alpha diversity to the corresponding biomass observations.
biomass <- biomass %>%
  left_join(alpha.df, by = c("SoilID", "Replicate", "Rhizo.Int"))