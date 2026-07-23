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
library(brms)
library(ggplot2)
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


# ------------------------------------------
# Clean biomass data and join soil chemistry

soil.chemistry <- soil.chemistry.site %>%
  pivot_wider(
    names_from = Treatment,
    values_from = c(NO3, NH4, pH, SOM)
  ) %>%
  mutate(
    totalN_INT   = NH4_INT * 0.778 + NO3_INT * 0.226,
    totalN_RHIZO = NH4_RHIZO * 0.778 + NO3_RHIZO * 0.226
  )

biomass <- biomass %>%
  mutate(
    Treatment = clean_treatment(Treatment),
    SoilID = factor(SoilID, levels = site.levels),
    Rhizo.Int = factor(Rhizo.Int, levels = soil.source.levels)
  ) %>%
  left_join(
    soil.chemistry,
    by = c("SoilID" = "SiteID")
  )
biomass.chem <- biomass

#built your bayes####
brte.chem <- brm(dry.mass ~ Treatment*Rhizo.Int*scale(pH_INT) +
               Treatment*Rhizo.Int*scale(totalN_INT) +
               Treatment*Rhizo.Int*scale(SOM_INT) +
               (1|SeedID),
             data= biomass.chem, family= Gamma(link="log"))
summary(brte.chem)
mcmc_plot(brte.chem, pars="b_")

plot_comparisons(brte.chem, variable=c("Treatment"), condition=c("Rhizo.Int", "pH_INT"))
plot_predictions(brte.chem, condition=c("Rhizo.Int"))
plot_comparisons(brte.chem, variable=c("Treatment"), condition=c("Rhizo.Int", list("pH_INT"=c(6, 6.7, 7.3))))

#####plots
microbe.effects <- comparisons(
  brtebio,
  variables = list(
    Treatment = c("Dead", "Live")
  ),
  by = c(
    "SoilID",
    "Rhizo.Int"
  ),
  newdata = biomass.chem
) %>%
  left_join(
    soil.chemistry %>%
      select(
        SiteID,
        pH_INT
      ),
    by = c(
      "SoilID" = "SiteID"
    )
  )
##
effect.ph <- microbe.effects %>%
  left_join(
    soil.chemistry %>% select(SiteID, pH_INT),
    by = c("SoilID" = "SiteID")
  )

# figure
chem.fig <- plot_comparisons(
  brte.chem,
  variables = "Treatment",
  condition = c("pH_INT", "Rhizo.Int"),
  newdata = biomass.chem
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_color_manual(
    values = c(
      "INT" = "#8C6D46",
      "RHIZO" = "#2E8B57"
    ),
    name = "Soil inoculum source",
    labels = c(
      "INT" = "Interspace",
      "RHIZO" = "Rhizosphere"
    )
  ) +
  scale_fill_manual(
    values = c(
      "INT" = "#8C6D46",
      "RHIZO" = "#2E8B57"
    ),
    name = "Soil inoculum source",
    labels = c(
      "INT" = "Interspace",
      "RHIZO" = "Rhizosphere"
    )
  ) +
  guides(
    color = guide_legend(nrow = 1),
    fill = guide_legend(nrow = 1)
  ) +
  labs(
    x = "Soil pH",
    y = "Effect of microbes on\ncheatgrass biomass"
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(
      size = 15,
      color = "#2F2A20"
    ),
    axis.title = element_text(
      size = 20,
      face = "bold",
      color = "#2F3A1F"
    ),
    axis.title.x = element_text(
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      margin = margin(r = 8)
    ),
    legend.title = element_text(
      size = 15,
      face = "bold",
      color = "#2F3A1F"
    ),
    legend.text = element_text(
      size = 14,
      color = "#2F2A20"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.border = element_rect(
      color = "#C7B99D",
      fill = NA,
      linewidth = 0.6
    ),
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    )
  )

chem.fig

ggsave(
  filename = "chem.fig.png",
  plot = chem.fig,
  width = 6,
  height = 5,
  dpi = 600,
  bg = "transparent"
)

#probe more into chemistry with alpha diversity
alpha.df <- alpha.df %>%
  mutate(
    SoilID = factor(SoilID, levels = site.levels),
    Rhizo.Int = factor(Rhizo.Int, levels = soil.source.levels)
  ) %>%
  left_join(
    soil.chemistry,
    by = c("SoilID" = "SiteID")
  )

alpha.chem2 <- brm(
  alphadiv ~
    Rhizo.Int *
    (scale(pH_INT) +
       scale(totalN_INT) +
       scale(SOM_INT)),
  data = alpha.df,
  family = gaussian()
)
summary(alpha.chem2)
cor.test(soil.chemistry$pH_INT, soil.chemistry$totalN_INT)

#figure for poster
library(patchwork)

alpha.pH.fig <- plot_predictions(
  alpha.chem2,
  condition = c("pH_INT", "Rhizo.Int")
) +
  labs(
    x = "Soil pH",
    y = "Predicted fungal richness",
    color = "Soil inoculum source",
    fill = "Soil inoculum source"
  ) +
  scale_color_manual(
    values = c(
      "INT" = "#8B7355",
      "RHIZO" = "#2E8B57"
    ),
    labels = c(
      "INT" = "Interspace",
      "RHIZO" = "Rhizosphere"
    )
  ) +
  scale_fill_manual(
    values = c(
      "INT" = "#8B7355",
      "RHIZO" = "#2E8B57"
    ),
    labels = c(
      "INT" = "Interspace",
      "RHIZO" = "Rhizosphere"
    )
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")

alpha.N.fig <- plot_predictions(
  alpha.chem2,
  condition = c("totalN_INT", "Rhizo.Int")
) +
  labs(
    x = "Total soil nitrogen (%)",
    y = NULL,
    color = "Soil inoculum source",
    fill = "Soil inoculum source"
  ) +
  scale_color_manual(
    values = c(
      "INT" = "#8B7355",
      "RHIZO" = "#2E8B57"
    ),
    labels = c(
      "INT" = "Interspace",
      "RHIZO" = "Rhizosphere"
    )
  ) +
  scale_fill_manual(
    values = c(
      "INT" = "#8B7355",
      "RHIZO" = "#2E8B57"
    ),
    labels = c(
      "INT" = "Interspace",
      "RHIZO" = "Rhizosphere"
    )
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

alpha.chem.fig <- alpha.pH.fig + alpha.N.fig +
  plot_layout(guides = "collect")

alpha.chem.fig #for poster

ggsave(
  filename = "alpha.chem.fig.png",
  plot = alpha.chem.fig,
  width = 8,
  height = 5,
  dpi = 600,
  bg = "transparent"
)