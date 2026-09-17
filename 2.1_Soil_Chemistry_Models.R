# Chapter 1 Eco-Evo Phase 1: Soil chemistry models
# Author: Jennah Brown
# Date: 9/17/2026
# Purpose: Fit final model testing relationships between soil chemistry and BRTE biomass.

library(brms)
library(bayesplot)

# Run soil chemistry data setup
source("2.0_Soil_Chemistry_Data_Setup.R")

#---------------
# Soil chemistry model
#---------------

# Test whether the effect of microbial treatment on BRTE biomass varies with soil pH, total nitrogen, and soil organic matter.
brte.chem <- brm(
  dry.mass ~ Treatment * Rhizo.Int * scale(pH_INT) +
    Treatment * Rhizo.Int * scale(totalN_INT) +
    Treatment * Rhizo.Int * scale(SOM_INT) +
    (1 | SeedID),
  data = biomass.chem,
  family = Gamma(link = "log")
)

# Review model estimates and convergence
summary(brte.chem)
mcmc_plot(brte.chem, pars = "b_")