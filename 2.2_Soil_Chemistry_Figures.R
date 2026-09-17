# Chapter 1 Eco-Evo Phase 1: Soil chemistry figures
# Author: Jennah Brown
# Date: 9/17/2026
# Purpose: Generate the soil chemistry figure used in the final ESA poster.

library(tidyverse)
library(marginaleffects)

# Run soil chemistry model
source("2.1_Soil_Chemistry_Models.R")

#---------------
# Soil pH figure
#---------------

# Poster Figure D: estimated effect of live microbial inoculum on BRTE biomass across the soil pH gradient.
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