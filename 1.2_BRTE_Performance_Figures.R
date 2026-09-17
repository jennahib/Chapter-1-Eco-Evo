# Chapter 1 Eco-Evo Phase 1: BRTE performance figures
# Author: Jennah Brown
# Date: 9/12/2026
# Purpose: Generate figures and model-based predictions for BRTE germination, biomass, and height.

library(tidyverse)
library(marginaleffects)

# Run data setup and model scripts
source("1.0_BRTE_Performance_Data_Setup.R")
source("1.1_BRTE_Performance_Models.R")

#------------
# Germination
#------------
# Generate predicted BRTE germination across microbial treatment, inoculum source, and soil population.
germ.fig <- plot_predictions(
  bmod1,
  condition = c("Treatment", "Rhizo.Int", "PopID")
) +
  labs(
    x = NULL,
    y = "Predicted germinated",
    color = "Soil inoculum"
  ) +
  scale_x_discrete(
    labels = c(
      "Dead" = "Sterile",
      "Live" = "Live"
    )
  ) +
  scale_color_manual(
    values = c(
      "INT" = "#7A8B3A",
      "RHIZO" = "#8B5A2B"
    ),
    labels = c(
      "INT" = "Interspace",
      "RHIZO" = "Rhizosphere"
    )
  ) +
  facet_wrap(~PopID, nrow = 1) +
  theme_classic(base_size = 18) +
  theme(
    strip.background = element_rect(fill = "#D8C3A5", color = NA),
    strip.text = element_text(face = "bold", size = 20, color = "#2F3A1F"),
    axis.text = element_text(size = 15, color = "#2F2A20"),
    axis.text.x = element_text(size = 15, angle = 20, hjust = 1),
    axis.title = element_text(size = 20, face = "bold", color = "#2F3A1F"),
    legend.title = element_text(size = 15, face = "bold", color = "#2F3A1F"),
    legend.text = element_text(size = 14, color = "#2F2A20"),
    legend.position = "right",
    panel.border = element_rect(color = "#C7B99D", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

germ.fig

#--------------
#Survival
#--------------
# Final BRTE survival was very high across treatments, resulting in insufficient variation for a useful model. Because no final survival model was retained, no model based survival figure was produced.

#-------------
#Biomass
#-------------
# Plot the estimated effect of live microbes on BRTE biomass across inoculum source and soil population.
microbe.fig <- plot_comparisons(
  brtebio,
  variables = list(Treatment = c("Dead", "Live")),
  condition = c("Rhizo.Int", "PopID")
) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.6
  ) +
  scale_x_discrete(
    labels = c(
      "INT" = "Interspace",
      "RHIZO" = "Rhizosphere"
    )
  ) +
  scale_color_manual(
    values = c(
      "BB" = "#8C6D46",
      "HG" = "#D98C3F",
      "NG" = "#D8C8A8",
      "OH" = "#6E9A63",
      "SR" = "#A66A3A"
    ),
    labels = c(
      "BB" = "Bogus Basin",
      "HG" = "Hulls Gulch",
      "NG" = "Nancy Gulch",
      "OH" = "Owyhee High",
      "SR" = "Snake River"
    )
  ) +
  labs(
    x = "Soil inoculum source",
    y = "Effect of live microbes on\ncheatgrass biomass"
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(size = 15, color = "#2F2A20"),
    axis.text.x = element_text(size = 15),
    axis.title = element_text(size = 20, face = "bold", color = "#2F3A1F"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.position = "none",
    panel.border = element_rect(color = "#C7B99D", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 40)
  )

microbe.fig

# Poster Figure C: BRTE biomass in sterile soils
# Compare predicted BRTE biomass between interspace and rhizosphere soils without live microbes.
brte.performance <- avg_predictions(
  brtebio,
  newdata = datagrid(
    Treatment = "Dead",
    Rhizo.Int = c("INT", "RHIZO"),
    grid_type = "balanced"
  ),
  by = "Rhizo.Int"
)

brte.performance.fig <- ggplot(
  brte.performance,
  aes(
    x = Rhizo.Int,
    y = estimate,
    color = Rhizo.Int
  )
) +
  geom_point(size = 4) +
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    width = 0.08,
    linewidth = 0.8
  ) +
  scale_x_discrete(
    labels = c(
      "INT" = "Interspace",
      "RHIZO" = "Rhizosphere"
    )
  ) +
  scale_color_manual(
    values = c(
      "INT" = "#8B5A2B",
      "RHIZO" = "#7A8B3A"
    ),
    guide = "none"
  ) +
  labs(
    x = "Soil inoculum source",
    y = "Predicted cheatgrass biomass (g)\nin sterile soils"
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(size = 15, color = "#2F2A20"),
    axis.title = element_text(size = 20, face = "bold", color = "#2F3A1F"),
    panel.border = element_rect(color = "#C7B99D", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

brte.performance.fig

# Additional biomass performance figure
# Predicted BRTE biomass across microbial treatment, inoculum source, and soil population.
biomass.fig <- plot_predictions(
  brtebio,
  condition = c("Treatment", "Rhizo.Int", "PopID")
) +
  labs(
    x = NULL,
    y = "Predicted cheatgrass biomass (g)",
    color = "Soil inoculum"
  ) +
  scale_x_discrete(
    labels = c(
      "Dead" = "Sterile",
      "Live" = "Live"
    )
  ) +
  scale_color_manual(
    values = c(
      "INT" = "#7A8B3A",
      "RHIZO" = "#8B5A2B"
    ),
    labels = c(
      "INT" = "Interspace",
      "RHIZO" = "Rhizosphere"
    )
  ) +
  facet_wrap(~PopID, nrow = 1) +
  theme_classic(base_size = 18) +
  theme(
    strip.background = element_rect(fill = "#D8C3A5", color = NA),
    strip.text = element_text(face = "bold", size = 20, color = "#2F3A1F"),
    axis.text = element_text(size = 15, color = "#2F2A20"),
    axis.text.x = element_text(size = 15, angle = 20, hjust = 1),
    axis.title = element_text(size = 20, face = "bold", color = "#2F3A1F"),
    legend.title = element_text(size = 15, face = "bold", color = "#2F3A1F"),
    legend.text = element_text(size = 14, color = "#2F2A20"),
    legend.position = "none",
    panel.border = element_rect(color = "#C7B99D", fill = NA, linewidth = 0.6),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

biomass.fig

#-------------
# Height
#------------
# Plot predicted BRTE height across microbial treatment and inoculum source.
height.fig <- plot_predictions(
  height.simp,
  condition = c("Treatment", "Rhizo.Int")
)

height.fig