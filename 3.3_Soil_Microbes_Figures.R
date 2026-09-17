# Chapter 1 Eco-Evo Phase 1: Soil microbial figures
# Author: Jennah Brown
# Date: 9/17/2026
# Purpose: Generate figures for the Eco-Evo Phase 1 soil microbial analysis and ESA poster.

library(tidyverse)

# Run microbial models
source("3.2_Soil_Microbes_Models.R")

#---------------
# Figure A: Fungal richness by site and soil source
#---------------

# Compare observed fungal ASV richness between interspace and rhizosphere soils across sagebrush populations.
rich.fig <- rich %>%
  filter(Sample.Type %in% c("Interspace", "Rhizosphere")) %>%
  mutate(
    Sample.Type = factor(
      Sample.Type,
      levels = c("Interspace", "Rhizosphere")
    )
  ) %>%
  ggplot(
    aes(
      x = SoilID,
      y = Observed,
      fill = Sample.Type
    )
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    linewidth = 0.6,
    outlier.shape = 16,
    outlier.size = 2
  ) +
  scale_x_discrete(
    labels = c(
      "BB" = "Bogus Basin",
      "HG" = "Hulls Gulch",
      "NG" = "Nancy Gulch",
      "OH" = "Owyhee High",
      "SR" = "Snake River"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Interspace" = "#8C7350",
      "Rhizosphere" = "#5F8F55"
    ),
    name = "Soil inoculum source"
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      ncol = 1
    )
  ) +
  labs(
    x = "Sagebrush population",
    y = "Observed fungal richness"
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(
      size = 15,
      color = "#2F2A20"
    ),
    axis.text.x = element_text(
      size = 14,
      angle = 0
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
      size = 13,
      color = "#2F2A20"
    ),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.margin = margin(t = 5),
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
    ),
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 15,
      l = 15
    )
  )
rich.fig