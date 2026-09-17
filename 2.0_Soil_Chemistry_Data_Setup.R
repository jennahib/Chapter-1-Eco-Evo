# Chapter 1 Eco-Evo Phase 1: Soil chemistry data setup
# Author: Jennah Brown
# Date: 9/12/2026
# Purpose: Import, summarize, and join soil chemistry data with BRTE biomass data for downstream analysis.

library(tidyverse)

#---------------
# Import data
#---------------

som <- read_csv(
  "Data/BRCBIO 23-24 - Soil chemistry_ pH, N, and C_LOI - C.csv",
  show_col_types = FALSE
)

nitrogen <- read_csv(
  "Data/BRCBIO 23-24 - Soil chemistry_ pH, N, and C_LOI - nitrogen.csv",
  show_col_types = FALSE
)

ph <- read_csv(
  "Data/BRCBIO 23-24 - Soil chemistry_ pH, N, and C_LOI - pH.csv",
  show_col_types = FALSE
)

biomass <- read_csv(
  "Data/Eco-Evo Phase 1 Biomass Data Fall 2025.csv",
  show_col_types = FALSE
)

#---------------
# Clean soil source labels
#---------------

nitrogen <- nitrogen %>%
  mutate(
    Treatment = recode(
      Treatment,
      "Interspace" = "INT",
      "Live Rhizosphere" = "RHIZO",
      "LIve Rhizosphere" = "RHIZO"
    )
  )

ph <- ph %>%
  mutate(
    Treatment = recode(
      Treatment,
      "Interspace" = "INT",
      "Live Rhizosphere" = "RHIZO",
      "LIve Rhizosphere" = "RHIZO"
    )
  )

som <- som %>%
  mutate(
    Treatment = recode(
      Treatment,
      "Interspace" = "INT",
      "Live Rhizosphere" = "RHIZO"
    )
  )

#---------------
# Calculate average soil chemistry by site and soil source
#---------------

nitrogen.site <- nitrogen %>%
  filter(Treatment %in% c("INT", "RHIZO")) %>%
  group_by(SiteID, Treatment) %>%
  summarize(
    NO3 = mean(mean.NO3, na.rm = TRUE),
    NH4 = mean(mean.NH4, na.rm = TRUE),
    .groups = "drop"
  )

ph.site <- ph %>%
  filter(Treatment %in% c("INT", "RHIZO")) %>%
  group_by(SiteID, Treatment) %>%
  summarize(
    pH = mean(pH.avg, na.rm = TRUE),
    .groups = "drop"
  )

som.site <- som %>%
  filter(Treatment %in% c("INT", "RHIZO")) %>%
  group_by(SiteID, Treatment) %>%
  summarize(
    SOM = mean(per.SOM, na.rm = TRUE),
    .groups = "drop"
  )

#---------------
# Combine soil chemistry data
#---------------

soil.chemistry.site <- nitrogen.site %>%
  left_join(
    ph.site,
    by = c("SiteID", "Treatment")
  ) %>%
  left_join(
    som.site,
    by = c("SiteID", "Treatment")
  )

#---------------
# Put interspace and rhizosphere chemistry in separate columns
# and calculate total nitrogen
#---------------

soil.chemistry <- soil.chemistry.site %>%
  pivot_wider(
    names_from = Treatment,
    values_from = c(NO3, NH4, pH, SOM)
  ) %>%
  mutate(
    totalN_INT = NH4_INT * 0.778 + NO3_INT * 0.226,
    totalN_RHIZO = NH4_RHIZO * 0.778 + NO3_RHIZO * 0.226
  )
#---------------
# Add soil chemistry values to BRTE biomass data
#---------------

biomass.chem <- biomass %>%
  left_join(
    soil.chemistry,
    by = c("SoilID" = "SiteID")
  )