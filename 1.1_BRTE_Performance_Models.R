# Chapter 1 Eco-Evo Phase 1: BRTE performance models
# Author: Jennah Brown
# Date: 9/12/2026
# Purpose: Fit final models for BRTE germination, biomass, and height.

library(brms)
library(marginaleffects)
library(bayesplot)

# Run BRTE data setup script
source("1.0_BRTE_Performance_Data_Setup.R")

#---------
# Germination
#-----------
# Fit Bayesian model for maximum BRTE germination.Maximum germination is modeled as the number of germinated seeds out of 8 seeds planted per cone-tainer.
bmod1 <- brm(
  max_germ | trials(8) ~ Treatment * Rhizo.Int * PopID +
    SeedID +
    (1 | ConeID),
  data = germ.brte,
  family = binomial(link = "logit")
)

# Review posterior estimates and convergence
summary(bmod1)
mcmc_plot(bmod1)
#-------------
# Survival
#------------
# Final BRTE survival was very high across treatments, resulting in insufficient variation for a useful model. Survival data were therefore retained and summarized, but no final Bayesian survival model was used for interpretation.

#-----------
#Biomass
#-----------
# Fit Bayesian model for BRTE aboveground dry biomass
brtebio <- brm(
  dry.mass ~ Treatment * Rhizo.Int * PopID + SeedID,
  data = biomass,
  family = Gamma(link = "log")
)

# Review posterior estimates and convergence
summary(brtebio)
mcmc_plot(brtebio)

# Save biomass model for use in the soil chemistry analysis
dir.create("outputs/models", recursive = TRUE, showWarnings = FALSE)
saveRDS(brtebio, "outputs/models/brtebio.rds")
#-----------
#Height
#-----------
# Fit Bayesian model for BRTE height through time. ConeID is included as a random effect to account for repeated height measurements from the same cone-tainer.
height.simp <- brm(
  height ~ time + Treatment * Rhizo.Int + (1 | ConeID),
  data = height.long,
  family = gaussian()
)

# Review posterior estimates and convergence
summary(height.simp)
mcmc_plot(height.simp)