# Chapter 1 Eco-Evo Phase 1: BRTE data setup
# Author: Jennah Brown
#Date: 9/12/2026
# Purpose: Import and clean BRTE germination, survival, biomass, and height data for downstream analyses.

library(tidyverse)
library(lubridate)

#---------------
#Germination
#--------------
# Import Phase 1 Version 1 germination data
eco.evo.p1.germ <- read_csv(
  "Data/Eco-Evo Phase 1 Spring 2025 GH Germination Data.csv",
  show_col_types = FALSE
)

# Reshape repeated germination data from wide to long format
long <- eco.evo.p1.germ %>%
  pivot_longer(
    cols = matches("^(Date|Live|New.Dead)\\d+$"),
    names_to = c(".value", "tmp"),
    names_pattern = "^(Date|Live|New.Dead)(\\d+)$"
  ) %>%
  mutate(Date = mdy(as.character(Date))) %>%
  select(-tmp) %>%
  arrange(ConeID, Date)

# Calculate maximum germination observed for each cone-tainer
long.germ.max <- long %>%
  group_by(
    ConeID, SoilID, Replicate, Treatment,
    Rhizo.Int, SeedID, PopID, SpeciesID
  ) %>%
  summarise(
    max_germ = max(Live, na.rm = TRUE),
    .groups = "drop"
  )

# Retain cheatgrass observations for downstream analysis
germ.brte <- long.germ.max %>%
  filter(SpeciesID == "BRTE")
#-------------
#survival
#------------
# Import Phase 1 Version 1 survival data
survival <- read_csv(
  "Data/Eco-Evo Phase 1 Spring 2025 GH Survival Data.csv",
  show_col_types = FALSE
)

# Reshape repeated survival censuses from wide to long format
sur.long <- survival %>%
  pivot_longer(
    cols = matches("^(Date|Live)\\d+$"),
    names_to = c(".value", "tmp"),
    names_pattern = "^(Date|Live)(\\d+)$"
  ) %>%
  mutate(
    Date = mdy(as.character(Date)),
    Live = as.integer(Live)
  ) %>%
  select(-tmp) %>%
  arrange(ConeID, Date)

# Use the final recorded survival census from July 31, 2025
end_date <- mdy("7/31/2025")

sur.final <- sur.long %>%
  filter(Date == end_date) %>%
  group_by(
    ConeID, SoilID, Replicate, Treatment,
    Rhizo.Int, SeedID, PopID, SpeciesID
  ) %>%
  summarise(
    brte.survival = max(Live, na.rm = TRUE),
    .groups = "drop"
  )

# Retain cheatgrass observations for downstream analysis
survival.brte <- sur.final %>%
  filter(SpeciesID == "BRTE")
#-----------
#Biomass
#-----------
# Import Phase 1 Version 1 cheatgrass biomass data
biomass <- read_csv(
  "Data/Eco-Evo Phase 1 Biomass Data Fall 2025.csv",
  show_col_types = FALSE
)
#-----------
#Height
#----------
# Import Phase 1 Version 1 cheatgrass height data
height <- read_csv(
  "Data/Eco-Evo Phase 1 Spring 2025 GH Height Data BRTE.csv",
  show_col_types = FALSE
)

# Reshape repeated height measurements from wide to long format
height.long <- height %>%
  pivot_longer(
    cols = matches("Height|Date"),
    names_to = c(".value", "time"),
    names_pattern = "(Height|Date)(\\d)"
  ) %>%
  rename(
    height = Height,
    date = Date
  ) %>%
  mutate(
    time = as.integer(time),
    time_f = factor(time),
    date = as.Date(date, format = "%m/%d/%Y")
  )

# Remove dates without a height measurement and retain cheatgrass
height.long <- height.long %>%
  filter(!is.na(height), SpeciesID == "BRTE") %>%
  arrange(ConeID, time)