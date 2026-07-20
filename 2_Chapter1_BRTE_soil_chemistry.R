# Chapter 1, Phase 1 (version 1)
# Soil chemistry assessment across sites and soil sources
# Author: Jennah Brown
#Date: 7/2026

# Purpose:
#   1. Summarize soil nitrogen, pH, and organic matter.
#   2. Explore associations between soil chemistry and estimated
#      microbial effects on cheatgrass biomass.
#   3. Explore associations between soil chemistry and mean biomass.

library(tidyverse)
library(marginaleffects)

# -------------------------------------------------------------------------
# Set sites and soil sources
# -------------------------------------------------------------------------

site.levels <- c("BB", "HG", "NG", "OH", "SR")
soil.source.levels <- c("INT", "RHIZO")

clean_treatment <- function(x) {
  recode(
    x,
    "LIve Rhizosphere" = "Live Rhizosphere"
  )
}

recode.soil.source <- function(x) {
  recode(
    x,
    "Interspace" = "INT",
    "Live Rhizosphere" = "RHIZO",
    "LIve Rhizosphere" = "RHIZO"
  )
}

# ---------------------
# Import data
# ----------------------
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
# --------------------------------------------------
# Summarize soil chemistry by site and soil source
# -----------------------------------------------------

nitrogen.site <- nitrogen %>%
  mutate(
    Treatment = recode.soil.source(Treatment)
  ) %>%
  filter(
    SiteID %in% site.levels,
    Treatment %in% soil.source.levels
  ) %>%
  group_by(
    SiteID,
    Treatment
  ) %>%
  summarize(
    NO3 = mean(mean.NO3, na.rm = TRUE),
    NH4 = mean(mean.NH4, na.rm = TRUE),
    .groups = "drop"
  )

ph.site <- ph %>%
  mutate(
    Treatment = recode.soil.source(Treatment)
  ) %>%
  filter(
    SiteID %in% site.levels,
    Treatment %in% soil.source.levels
  ) %>%
  group_by(
    SiteID,
    Treatment
  ) %>%
  summarize(
    pH = mean(pH.avg, na.rm = TRUE),
    .groups = "drop"
  )

som.site <- som %>%
  mutate(
    Treatment = recode.soil.source(Treatment)
  ) %>%
  filter(
    SiteID %in% site.levels,
    Treatment %in% soil.source.levels
  ) %>%
  group_by(
    SiteID,
    Treatment
  ) %>%
  summarize(
    SOM = mean(per.SOM, na.rm = TRUE),
    .groups = "drop"
  )

soil.chemistry.site <- nitrogen.site %>%
  left_join(
    ph.site,
    by = c("SiteID", "Treatment")
  ) %>%
  left_join(
    som.site,
    by = c("SiteID", "Treatment")
  ) %>%
  mutate(
    SiteID = factor(
      SiteID,
      levels = site.levels
    ),
    Treatment = factor(
      Treatment,
      levels = soil.source.levels
    )
  )


# ---------------------------------------------
# Soil nitrate by site and soil source
# --------------------------------------------

nitrate.site.fig <- ggplot(
  nitrogen %>%
    mutate(
      Treatment = recode.soil.source(Treatment),
      SiteID = factor(
        SiteID,
        levels = site.levels
      )
    ) %>%
    filter(
      SiteID %in% site.levels,
      Treatment %in% soil.source.levels
    ),
  aes(
    x = SiteID,
    y = mean.NO3,
    fill = Treatment
  )
) +
  stat_summary(
    fun = mean,
    geom = "col",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  labs(
    x = "Site",
    y = expression(
      "Soil nitrate"
    ),
    fill = "Soil source"
  ) +
  scale_fill_discrete(
    labels = c(
      "INT" = "Interspace",
      "RHIZO" = "Rhizosphere"
    )
  ) +
  theme_classic()

nitrate.site.fig


# -------------------------------------------------------------------------
# Estimated microbial effects on cheatgrass biomass
# -------------------------------------------------------------------------

bbio1 <- readRDS("outputs/models/bbio1.rds")
microbe.effect <- comparisons(
  bbio1,
  variables = "Treatment",
  by = c("SoilID", "Rhizo.Int")
)

microbe.effect.site <- microbe.effect %>%
  select(
    SoilID,
    Rhizo.Int,
    estimate,
    conf.low,
    conf.high
  ) %>%
  mutate(
    SoilID = factor(
      SoilID,
      levels = site.levels
    ),
    Rhizo.Int = factor(
      Rhizo.Int,
      levels = soil.source.levels
    )
  )


# -------------------------------------------------------------------------
# Join microbial effects with soil chemistry
# -------------------------------------------------------------------------

microbe.effect.chemistry <- microbe.effect.site %>%
  left_join(
    soil.chemistry.site,
    by = c(
      "SoilID" = "SiteID",
      "Rhizo.Int" = "Treatment"
    )
  )

microbe.effect.chemistry.long <- microbe.effect.chemistry %>%
  pivot_longer(
    cols = c(NO3, NH4, pH, SOM),
    names_to = "Chemistry",
    values_to = "Value"
  ) %>%
  mutate(
    Chemistry = factor(
      Chemistry,
      levels = c("NO3", "NH4", "pH", "SOM"),
      labels = c(
        "Nitrate",
        "Ammonium",
        "pH",
        "Soil organic matter"
      )
    )
  )


# -----------------------------------------
# microbial effects and soil chemistry
# -----------------------------------------

microbe.effect.chemistry.fig <- ggplot(
  microbe.effect.chemistry.long,
  aes(
    x = Value,
    y = estimate,
    color = Rhizo.Int
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.7
  ) +
  facet_wrap(
    ~ Chemistry,
    scales = "free_x"
  ) +
  labs(
    x = "Mean soil chemistry value",
    y = paste0(
      "Estimated microbial effect on biomass\n",
      "(Live - Dead)"
    ),
    color = "Soil source"
  ) +
  scale_color_discrete(
    labels = c(
      "INT" = "Interspace",
      "RHIZO" = "Rhizosphere"
    )
  ) +
  theme_classic()

microbe.effect.chemistry.fig


# Overall correlations
microbe.effect.correlations <- microbe.effect.chemistry %>%
  summarize(
    across(
      c(NO3, NH4, pH, SOM),
      ~ cor(
        estimate,
        .x,
        use = "complete.obs"
      ),
      .names = "{.col}.r"
    )
  )

microbe.effect.correlations


# Correlations calculated separately by soil source
microbe.effect.correlations.by.source <-
  microbe.effect.chemistry %>%
  group_by(Rhizo.Int) %>%
  summarize(
    across(
      c(NO3, NH4, pH, SOM),
      ~ cor(
        estimate,
        .x,
        use = "complete.obs"
      ),
      .names = "{.col}.r"
    ),
    .groups = "drop"
  )

microbe.effect.correlations.by.source


# -------------------------------------------------------
# Mean cheatgrass biomass in live and autoclaved soils
# --------------------------------------------------------

biomass.site <- biomass %>%
  filter(
    SoilID %in% site.levels,
    Treatment %in% c("Live", "Dead"),
    Rhizo.Int %in% soil.source.levels
  ) %>%
  group_by(
    SoilID,
    Treatment,
    Rhizo.Int
  ) %>%
  summarize(
    biomass = mean(
      dry.mass,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    SoilID = factor(
      SoilID,
      levels = site.levels
    ),
    Treatment = factor(
      Treatment,
      levels = c("Dead", "Live")
    ),
    Rhizo.Int = factor(
      Rhizo.Int,
      levels = soil.source.levels
    )
  )


# -------------------------------------------------------------------------
# Join mean biomass with soil chemistry
# -------------------------------------------------------------------------

biomass.chemistry <- biomass.site %>%
  left_join(
    soil.chemistry.site,
    by = c(
      "SoilID" = "SiteID",
      "Rhizo.Int" = "Treatment"
    )
  )

biomass.chemistry.long <- biomass.chemistry %>%
  pivot_longer(
    cols = c(NO3, NH4, pH, SOM),
    names_to = "Chemistry",
    values_to = "Value"
  ) %>%
  mutate(
    Chemistry = factor(
      Chemistry,
      levels = c("NO3", "NH4", "pH", "SOM"),
      labels = c(
        "Nitrate",
        "Ammonium",
        "pH",
        "Soil organic matter"
      )
    )
  )


# -------------------------------------------------------------------------
# explore mean biomass and soil chemistry
# -------------------------------------------------------------------------

biomass.chemistry.fig <- ggplot(
  biomass.chemistry.long,
  aes(
    x = Value,
    y = biomass,
    color = Rhizo.Int
  )
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.7
  ) +
  facet_grid(
    Treatment ~ Chemistry,
    scales = "free_x"
  ) +
  labs(
    x = "Mean soil chemistry value",
    y = "Mean cheatgrass biomass (g)",
    color = "Soil source"
  ) +
  scale_color_discrete(
    labels = c(
      "INT" = "Interspace",
      "RHIZO" = "Rhizosphere"
    )
  ) +
  theme_classic()

biomass.chemistry.fig


# Correlations calculated separately for live and autoclaved soils
biomass.chemistry.correlations <- biomass.chemistry %>%
  group_by(Treatment) %>%
  summarize(
    across(
      c(NO3, NH4, pH, SOM),
      ~ cor(
        biomass,
        .x,
        use = "complete.obs"
      ),
      .names = "{.col}.r"
    ),
    .groups = "drop"
  )

biomass.chemistry.correlations