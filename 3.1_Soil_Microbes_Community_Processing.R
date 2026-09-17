# Chapter 1 Eco-Evo Phase 1: Soil microbial community processing
# Author: Jennah Brown
# Date: 9/17/2026
#
# Purpose: Calculate Bray-Curtis community dissimilarity, NMDS ordination, and microbial community projection scores.
#
# This workflow was adapted from microbial community analysis code developed by Allison Simler-Williamson.

library(phyloseq)
library(tidyverse)
library(vegan)

# Run microbial data setup
source("3.0_Soil_Microbes_Data_Setup.R")

#---------------
# Bray-Curtis NMDS
#---------------

# Convert sequence counts to relative abundance within each sample.
ps.focal <- transform_sample_counts(ps.focal, function(otu) otu / sum(otu))

# Run NMDS using Bray-Curtis dissimilarity. The seed makes the ordination reproducible.
set.seed(1)
ord.nmds.bray <- ordinate(ps.focal, method = "NMDS", distance = "bray")

ord.nmds.bray

#---------------
# Bray-Curtis distances
#---------------

# Extract sample metadata for use with the Bray-Curtis distance calculations.
sampledf <- data.frame(sample_data(ps.focal))

# Calculate pairwise Bray-Curtis dissimilarities among samples.
ps_bray <- phyloseq::distance(ps.focal, method = "bray")

#---------------
# Microbial translocation distances
#---------------

# Calculate each sample's Bray-Curtis distance to the site x soil-source centroids.
dist.seed2soil <- usedist::dist_to_centroids(ps_bray, sampledf$SitePar) %>%
  dplyr::rename(Extraction.ID = Item) %>%
  mutate(Extraction.ID = as.character(Extraction.ID)) %>%
  
  # Separate centroid labels into comparison site and soil source.
  mutate(ParInt = str_sub(CentroidGroup, start = 3, end = 13),
         CompareSite = str_sub(CentroidGroup, start = 1, end = 2)) %>%
  select(-CentroidGroup) %>%
  pivot_wider(values_from = CentroidDistance, names_from = ParInt) %>%
  dplyr::rename(dist.PAR = Rhizosphere, dist.INT = Interspace) %>%
  
  # Add sample metadata and calculate mean distance to the rhizosphere centroid for each soil replicate and comparison site.
  left_join(sampledf, by = "Extraction.ID") %>%
  group_by(SoilID, Sample.Type, Replicate, CompareSite) %>%
  summarize(seed2soil.dist = mean(dist.PAR)) %>%
  
  # Rename variables to match the experimental datasets.
  dplyr::rename(SeedID = CompareSite) %>%
  mutate(Project.ID = "BRCBIO Spring 2025")

#---------------
# Soil-conditioning projection scores
#---------------

# Extract each sample's NMDS coordinates and pair them with the sample metadata.
bray_scores <- data.frame(scores(ord.nmds.bray, display = "sites"))
bray_scores$Extraction.ID <- as.character(row.names(bray_scores))
bray_scores <- bray_scores %>% left_join(sampledf)

#---------------
# Calculate projection scores within each site
#---------------

# Projection scores describe community change along the NMDS direction from the interspace centroid toward the rhizosphere centroid.
compute.projection <- function(df) {
  
  # Return NA if a site does not contain both interspace and rhizosphere samples.
  if (!all(c("Interspace", "Rhizosphere") %in% df$Sample.Type)) {
    df$Projection <- NA
    df$Projection.std <- NA
    return(df)
  }
  
  # Calculate the interspace and rhizosphere centroids in NMDS space.
  ref.centroid <- colMeans(df[df$Sample.Type == "Interspace", c("NMDS1", "NMDS2")])
  treat.centroid <- colMeans(df[df$Sample.Type == "Rhizosphere", c("NMDS1", "NMDS2")])
  
  # Define the direction from the interspace centroid toward the rhizosphere centroid.
  direction.vector <- treat.centroid - ref.centroid
  direction.unit <- direction.vector / sqrt(sum(direction.vector^2))
  
  # Center each sample relative to the interspace centroid.
  centered.scores <- sweep(df[, c("NMDS1", "NMDS2")], 2, ref.centroid)
  centered.scores <- as.matrix(centered.scores)
  
  # Project each sample onto the interspace-to-rhizosphere direction.
  projection.scores <- as.vector(centered.scores %*% direction.unit)
  df$Projection <- projection.scores
  
  # Standardize projection scores relative to variation among interspace samples.
  control.proj <- projection.scores[df$Sample.Type == "Interspace"]
  control.mean <- mean(control.proj, na.rm = TRUE)
  control.sd <- sd(control.proj, na.rm = TRUE)
  
  # Avoid dividing by zero if interspace samples have no variation.
  if (control.sd == 0) {
    df$Projection.std <- NA
  } else {
    df$Projection.std <- (projection.scores - control.mean) / control.sd
  }
  
  return(df)
}

#---------------
# Apply projection calculation
#---------------

# Calculate projection scores separately within each site, then use the median when multiple sequencing observations represent the same soil replicate.
soil.conditioning <- bray_scores %>%
  group_by(SoilID) %>%
  group_modify(~ compute.projection(.x)) %>%
  ungroup() %>%
  select(SoilID, Sample.Type, Replicate,
         Projection, Projection.std) %>%
  group_by(SoilID, Replicate, Sample.Type) %>%
  summarize(Projection = median(Projection),
            Projection.std = median(Projection.std))

#---------------
# Create projection values for combined replicates
#---------------

# Phase 1 used A-B and B-C soil replicate combinations. Calculate their projection values by averaging the corresponding individual replicates.
projection.combos <- combo.key %>%
  left_join(soil.conditioning %>%
              select(SoilID, Sample.Type, rep1 = Replicate,
                     Projection.1 = Projection),
            by = "rep1") %>%
  left_join(soil.conditioning %>%
              select(SoilID, Sample.Type, rep2 = Replicate,
                     Projection.2 = Projection),
            by = c("SoilID", "Sample.Type", "rep2")) %>%
  mutate(Projection = (Projection.1 + Projection.2) / 2) %>%
  select(SoilID, Replicate, Sample.Type, Projection)

# Combine measured and combo projection values and match soil-source naming used in the biomass data.
projection.df <- bind_rows(
  soil.conditioning %>%
    select(SoilID, Replicate, Sample.Type, Projection),
  projection.combos
) %>%
  mutate(Rhizo.Int = recode(
    Sample.Type,
    "Interspace" = "INT",
    "Rhizosphere" = "RHIZO"
  )) %>%
  select(SoilID, Replicate, Rhizo.Int, Projection)

#---------------
# Match projection scores to biomass data
#---------------

# Add projection scores to the corresponding biomass observations.
biomass <- biomass %>%
  left_join(
    projection.df,
    by = c("SoilID", "Replicate", "Rhizo.Int")
  )