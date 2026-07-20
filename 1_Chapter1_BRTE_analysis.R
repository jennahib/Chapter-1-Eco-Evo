#Chapter 1 Eco-Evo phase 1: Sagebrush microbial resistance to cheatgrass
#Author: Jennah Brown
#Date: 7/2026

# Purpose:
#   1. Analyze cheatgrass germination, survival, biomass, and height.
#   2. Estimate effects of live versus sterile inoculum.
#   3. Generate figures for the Chapter 1 analysis and ESA poster.

library(tidyverse)
library(brms)
library(marginaleffects)
library(bayesplot)
library(modelr)
library(pROC)
library(MLmetrics)
library(lubridate)
library(lme4)
library(rstan)
library(sjPlot)

# -----------------------------
# germination analysis
# -----------------------------
eco.evo.p1.germ <- read_csv(
  "Data/Eco-Evo Phase 1 SPring 2025 GH Germination Data.csv",
  show_col_types = FALSE
)
# ------ moving dates, lives, and new deads into long format -----
long <- eco.evo.p1.germ %>%
  pivot_longer(                                      ### pivot to longer format
    cols = matches("^(Date|Live|New.Dead)\\d+$"),    
    names_to = c(".value", "tmp"),                
    names_pattern = "^(Date|Live|New.Dead)(\\d+)$"  ##take those matches and split them
  ) %>%
  mutate(
    Date = mdy(as.character(Date))   # just change to character, date
  ) %>%
  select(-tmp) %>%               ## remove the temp label
  arrange(ConeID, Date)     ## group by coneID here

#get to max germ
long <- long %>%
  group_by(ConeID) %>%
  mutate(max_germ = max(Live, na.rm = TRUE)) %>%  ### this gets you to long, max germ without cutting dates to one cone identity
  ungroup()

#clean the data now to drop dates
long.germ.max <- long %>%
  group_by(
    ConeID,
    Soil.ID,
    Replicate,
    Treatment,
    Rhizo.Int,
    SeedID,
    PopID,
    SpeciesID
  ) %>%
  summarise(
    max_germ = max(Live, na.rm = TRUE),
    .groups = "drop"
  )

#subset for BRTE
germ.brte <- long.germ.max %>%
  filter(SpeciesID == "BRTE")

# germination Bayes:
bmod1 <- brm(
  max_germ | trials(8) ~ Treatment * Rhizo.Int * PopID +
    SeedID +
    (1 | ConeID),
  data = germ.brte,
  family = binomial(link = "logit")
)
summary(bmod1)

mcmc_plot(bmod1)

#figure
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
  facet_wrap(~PopID,nrow = 1)+
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
ggsave(
  "germination_poster_figure.png",
  germ.fig,
  width = 13,
  height = 5,
  dpi = 600,
  bg = "transparent"
)


##plot comparison figure
plot_comparisons(bmod1, variable=c("Treatment"), condition=c("Rhizo.Int", "PopID"))

# ---- probing bayes ----
mod2.bayes <- brm(max_germ~Treatment*Rhizo.Int + (1|Soil.ID), data = germ.brte, 
                  family = poisson(link = "log"),
                  chains = 4,
                  cores = 4,
                  iter = 2000,
                  warmup= 1000,
                  control= list(adapt_delta = 0.99) #can try 0.995 if this doesn't work
               )
pairs(mod2.bayes)

summary(mod2.bayes)
plot(mod2.bayes)
pp_check(mod2.bayes) 
mcmc_plot(mod2.bayes, pars = "^b_", prob = 0.95)
# + increases germ, - decreases
#no effect of treatment on germ, credible interval overlaps 0
#rhizo microbes: post around 0, does not support a negative effect on germ
#no evidence that rhizo effects depend on treatment
#i.e. rhizo microbes do not resist BRTE invasion

mod2.bayes <- brm(max_germ ~ Treatment*Rhizo.Int*PopID + (1+Rhizo.Int|Soil.ID), data = germ.brte, 
                  family = poisson(link = "log"),
                  chains = 4,
                  cores = 4,
                  iter = 2000,
                  warmup= 1000,
                  control= list(adapt_delta = 0.99)) #can try 0.995 if this doesn't work
summary(mod2.bayes)
plot(mod2.bayes)
pp_check(mod2.bayes) 
mcmc_plot(mod2.bayes, pars = "^b_", prob = 0.95)


# -----------------------------
# survival analysis
# -----------------------------
#survival cleaning: 
survival <- read_csv(
  "Data/Eco-Evo Phase 1 Spring 2025 GH Survival Data.csv",
  show_col_types = FALSE
)

#just going to use final survival date 7/31/25
sur.long <- survival %>%
  pivot_longer(
    cols = matches("^(Date|Live)\\d+$"),
    names_to = c(".value", "tmp"),
    names_pattern = "^(Date|Live)(\\d+)$"
  )%>%
  mutate(
    Date = mdy(as.character(Date)),
    Live = as.integer(Live)
  )%>%
  select(-tmp)%>%
  arrange(ConeID,Date)

end_date <- mdy("7/31/2025")  

sur.final <- sur.long %>%
  filter(Date == end_date) %>%
  group_by(
    ConeID, SoilID, Replicate, Treatment, Rhizo.Int, SeedID, PopID, SpeciesID
  ) %>%
  summarise(
    brte.survival = max(Live, na.rm = TRUE),
    .groups = "drop"
  )
#subset brte
survival.brte <- sur.final %>%
  filter(SpeciesID == "BRTE")

# ---- Survival Bayes: -----
reg_priors <- c(
  set_prior("normal(0,5)", class="b"),
  set_prior("normal(0,5)", class="Intercept")
)

surv1 <- brm(brte.survival ~ Treatment*Rhizo.Int*PopID + SeedID + (1|ConeID), 
            data= survival.brte, family= bernoulli(),
            prior=reg_priors)
summary(surv1) #this model is funky, survival is very high. Tried to set priors, but didn't help

mcmc_plot(surv1)
plot_predictions(surv1, condition=c("Treatment", "Rhizo.Int", "PopID"))
plot_comparisons(surv1, variable=c("Treatment"), condition=c("Rhizo.Int", "PopID"))

# -----------------------------
# biomass analysis
# -----------------------------
biomass <- read_csv(
  "Data/Eco-Evo Phase 1 Biomass Data Fall 2025.csv",
  show_col_types = FALSE
)

#biomass Bayes
bbio1 <- brm(dry.mass ~ Treatment*Rhizo.Int*PopID + SeedID + (1|ConeID), 
            data= biomass, family= Gamma(link="log"))
summary(bbio1)
mcmc_plot(bbio1)

#save model for soil chemistry
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/models", recursive = TRUE, showWarnings = FALSE)
saveRDS(
  bbio1,
  "outputs/models/bbio1.rds"
)

#use for ESA 2026 poster:
plot_comparisons(
  bbio1,
  variables = "Rhizo.Int",
  condition = c("Treatment", "PopID")
)

plot_comparisons(
  bbio1,
  variables = "Treatment",
  condition = c("Rhizo.Int", "PopID")
)

#### Microbial effect figure for poster
microbe.fig <- plot_comparisons(
  bbio1,
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
  scale_color_discrete(
    name = "Sagebrush population",
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
    y = "Effect of microbes on\ncheatgrass biomass"
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(
      size = 15,
      color = "#2F2A20"
    ),
    axis.text.x = element_text(
      size = 15
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
    legend.position = "right",
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
      r = 10,
      b = 10,
      l = 40
    )
  )
microbe.fig
ggsave(
  filename = "microbe_effect_biomass_poster.png",
  plot = microbe.fig,
  width = 12,
  height = 4,
  units = "in",
  dpi = 600,
  bg = "transparent"
)

#####plot predictions figure####
plot_predictions(bbio1, condition=c("Treatment", "Rhizo.Int", "PopID"))
biomass.fig <- plot_predictions(
  bbio1,
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
  facet_wrap(~PopID,nrow = 1)+
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
#save it:
ggsave(
  "biomass_poster_figure.png",
  biomass.fig,
  width = 13,
  height = 5,
  dpi = 600,
  bg = "transparent"
)


plot_comparisons(bbio1, variable=c("Treatment"), condition=c("Rhizo.Int", "PopID"))
plot_comparisons(bbio1, variable=c("Rhizo.Int"), condition=c("PopID"))
avg_predictions(bbio1, variable=c("Treatment", "Rhizo.Int", "PopID"))

#tx effect:
  plot_predictions(
    bbio1,
    condition = c("Treatment", "Rhizo.Int")
  )


# -----------------------------
# height analysis
# -----------------------------
height <- read_csv(
    "Data/Eco-Evo Phase 1 Spring 2025 GH Height Data BRTE.csv",
    show_col_types = FALSE
  )

# need to make long, and filter for just BRTE
height.long <- height %>%
  pivot_longer(
    cols = matches ("Height|Date"),
    names_to = c(".value", "time"),
    names_pattern = "(Height|Date)(\\d)"
  )

#clean
height.long <- height.long %>%
  rename(height = Height,
         date = Date) %>%
  mutate(
    time = as.integer(time),
    time_f = factor(time),
    date = as.Date(date, format = "%m/%d/%Y")
  )

#deal with missing cells (not all weeks plants were measured)
height.long <- height.long %>%
  filter(!is.na(height))

#filter to brte
height.long <- height.long %>%
  filter(SpeciesID == "BRTE")

#order obs withing conetainer ID
height.long <- height.long %>%
  arrange(ConeID, time)

#lagged height## 
height.long <- height.long %>%
  group_by(ConeID) %>%
  mutate(height_lag = lag(height)) %>%
  ungroup() #not using lagged height only final height now

#brte final height bayes
height.simp <- brm(height ~ time + Treatment * Rhizo.Int + (1|ConeID), data = height.long, family = gaussian())

summary(height.simp)
mcmc_plot(height.simp)
plot_predictions(height.simp, condition=c("Treatment", "Rhizo.Int", "SoilID"))
plot_comparisons(height.simp, variable=c("Treatment"), condition=c("Rhizo.Int", "SoilID"))
plot_comparisons(height.simp, variable=c("Rhizo.Int"), condition=c("SoilID"))