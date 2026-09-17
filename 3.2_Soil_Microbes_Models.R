# Chapter 1 Eco-Evo Phase 1: Soil microbial models
# Author: Jennah Brown
# Date: 9/17/2026
# Purpose: Fit final models testing relationships between fungal community metrics and BRTE performance.

library(brms)
library(bayesplot)

# Run microbial community processing
source("3.1_Soil_Microbes_Community_Processing.R")

#---------------
# Standardize microbial predictors
#---------------

# Standardize alpha diversity and projection scores before modeling.
biomass$alphadiv.std <- as.numeric(scale(biomass$alphadiv))
biomass$projection.std <- as.numeric(scale(biomass$Projection))

#---------------
# Microbial predictors of BRTE biomass
#---------------

# Test whether fungal alpha diversity and microbial community position are associated with BRTE biomass, and whether these relationships depend on microbial treatment and soil source.
bbio.microbe <- brm(
  dry.mass ~ Treatment * Rhizo.Int * alphadiv.std +
    Treatment * Rhizo.Int * projection.std,
  data = biomass,
  family = Gamma(link = "log")
)

summary(bbio.microbe)

mcmc_plot(bbio.microbe, pars = "b_")

#---------------
# Fungal richness across sites
#---------------

# Test whether observed fungal ASV richness differs among soil collection sites.
rich.model <- brm(
  Observed ~ SoilID,
  data = rich,
  family = gaussian()
)

summary(rich.model)

mcmc_plot(rich.model, pars = "b_")