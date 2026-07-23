# Chapter 1 Phase 1: Site, interspace, and rhizosphere comparisons
# Author: Jennah Brown
# Date: July 2026
#
# Purpose:
#   1. Quantify fungal alpha diversity across sagebrush populations.
#   2. Characterize fungal community composition using Bray-Curtis dissimilarity and NMDS.
#   3. Calculate soil-conditioning projection scores.
#   4. Generate figures for the Chapter 1 analysis and ESA poster.
# Data input:
# This analysis begins with a rarefied ITS phyloseq object generated from the BRC-BIO 1_Phyloseq_object_prep_ITS script preprocessing pipeline. The object contains all BRCBIO Spring 2025 samples,

library(phyloseq)
library(tidyverse)
library(vegan)
library(DESeq2)
library(brms)

##difference here is that i filtered data in the original asw scripts when i should extract her data then filter in this script because its confusing now

#load data####
#this is phyloseq object created and filtered to Spring 2025 project with specific sites
ps.spring <- readRDS(
  "Data/phyloseq_rare_ITS_BRCBIO_Spring2025.rds"
)
ps.focal <- subset_samples(
  ps.spring,
  SoilID %in% c("BB", "HG", "NG", "OH", "SR")
)

#metadata
# Extract metadata from the filtered phyloseq object
ps.metadata <- data.frame(sample_data(ps.focal))

#brte biomass
biomass <- read.csv("Data/Eco-Evo Phase 1 Biomass Data Fall 2025.csv")

# Alpha diversity calculation ###########
rich <- estimate_richness(ps.focal) %>% rownames_to_column(var="Extraction.ID") %>% 
  mutate(Extraction.ID = as.character(str_sub(Extraction.ID, 2, 4))) %>% 
  left_join(ps.focal@sam_data)

#plot it
ggplot(rich, aes(x=SoilID, y=Observed)) + 
  geom_boxplot() +
  facet_wrap(~Project.ID)

## Eco-Evo phase 1 used combos for soil replicates (A-B, B-C) we need to calculate these combos to match the biomass data
# Average sequencing replicates
alpha.measured <- rich %>%
  group_by(SoilID, Replicate, Sample.Type) %>%
  summarize(alphadiv = mean(Chao1, na.rm = TRUE),
            .groups = "drop")

# Create combo replicates
combo.key <- tribble(
  ~Replicate, ~rep1, ~rep2,
  "A-B", "A", "B",
  "B-C", "B", "C"
)

alpha.combos <- combo.key %>%
  left_join(alpha.measured %>%
              select(SoilID, Sample.Type,
                     rep1 = Replicate,
                     alphadiv.1 = alphadiv),
            by = "rep1") %>%
  left_join(alpha.measured %>%
              select(SoilID, Sample.Type,
                     rep2 = Replicate,
                     alphadiv.2 = alphadiv),
            by = c("SoilID", "Sample.Type", "rep2")) %>%
  mutate(alphadiv = (alphadiv.1 + alphadiv.2) / 2) %>%
  select(SoilID, Replicate, Sample.Type, alphadiv)

# Combine measured and combo values
alpha.df <- bind_rows(alpha.measured, alpha.combos)

# Note: OH Interspace C and SR Rhizosphere C were not sequenced,
# so those C reps and SR rhizo B-C combos remain NA.

# Match biomass naming
alpha.df <- alpha.df %>%
  mutate(
    Rhizo.Int = recode(
      Sample.Type,
      "Interspace" = "INT",
      "Rhizosphere" = "RHIZO"
    )
  ) %>%
  select(SoilID, Replicate, Rhizo.Int, alphadiv)

# Add alpha diversity to biomass
biomass <- biomass %>%
  left_join(alpha.df,
            by = c("SoilID", "Replicate", "Rhizo.Int"))

## Ordination / NMDS to visualize communities #######
#Split up ps object into the different experiments:

ps.focal <- transform_sample_counts(ps.focal, function(otu) otu/sum(otu))
# NMDS using ordinate():
set.seed(1)
ord.nmds.bray <- ordinate(ps.focal, method="NMDS", distance="bray")
ord.nmds.bray

ordplot <- plot_ordination(ps.focal, ord.nmds.bray, 
                           color="SoilID", shape="Sample.Type") +
  theme_minimal() +
  geom_point(size =2)     

ordplot

## Rhizosphere v. Interspace samples across sites:
ordplot1 <- ordplot+ 
  stat_ellipse(aes(group = SitePar), type = "t", linetype = 2) +
  facet_wrap(~ Sample.Type)
ordplot1

## Rhizosphere v. Interspace comparisons within each site:
ordplot2 <- ordplot+ 
  stat_ellipse(aes(group = SitePar), type = "t", linetype = 2) +
  facet_wrap(~ SoilID, scales = "free")
ordplot2

#Distance matrix and calculation of microbial translocation distance & soil conditioning metrics: ######
sampledf <- data.frame(sample_data(ps.focal)) 

##Calculate microbial translocation distances: ######
ps_bray <- phyloseq::distance(ps.focal, method = "bray")

## Distances for spring 2025 samples:
dist.seed2soil <- usedist::dist_to_centroids(ps_bray, sampledf$SitePar) %>% 
  dplyr::rename(Extraction.ID = Item) %>% 
  mutate(Extraction.ID = as.character(Extraction.ID)) %>% 
  ## Identify centroids as being for rhizosphere or interspace samples:
  mutate(ParInt = str_sub(CentroidGroup, start=3, end=13),
         CompareSite = str_sub(CentroidGroup, start=1, end=2)) %>%
  select(-CentroidGroup) %>%
  
  pivot_wider(values_from = CentroidDistance, names_from = ParInt) %>%
  dplyr::rename(dist.PAR = Rhizosphere, dist.INT = Interspace) %>%
  left_join(sampledf, by="Extraction.ID")  %>%
  
  group_by(SoilID, Sample.Type, Replicate, CompareSite) %>%
  summarize(seed2soil.dist = mean(dist.PAR))  %>%
  dplyr::rename(SeedID = CompareSite) %>%
  
  mutate(Project.ID = "BRCBIO Spring 2025")

#Soil conditioning as directional change in community composition: #######
## NMDS axes for greenhouse experiment sequencing data:
bray_scores <- data.frame(scores(ord.nmds.bray, display="sites"))
bray_scores$Extraction.ID <- as.character(row.names(bray_scores))
bray_scores <- bray_scores %>% left_join(sampledf)

# Function to compute projection scores within each site####
compute.projection <- function(df) {
  if (!all(c("Interspace", "Rhizosphere") %in% df$Sample.Type)) {
    df$Projection <- NA
    df$Projection_std <- NA
    return(df)
  }
  # Get centroids for each site
  ref.centroid <- colMeans(df[df$Sample.Type == "Interspace", c("NMDS1", "NMDS2")])
  treat.centroid <- colMeans(df[df$Sample.Type == "Rhizosphere", c("NMDS1", "NMDS2")])

  # Calculate the direction unit vector
  direction.vector <- treat.centroid - ref.centroid
  direction.unit <- direction.vector / sqrt(sum(direction.vector^2))

  # Center to the centroid for reference
  centered.scores <- sweep(df[, c("NMDS1", "NMDS2")], 2, ref.centroid)
  centered.scores <- as.matrix(centered.scores)

  # Projection score calculation:
  projection.scores <- as.vector(centered.scores %*% direction.unit)
  df$Projection <- projection.scores

  ## Standardize using "baseline" variability in interspace samples:
  control.proj <- projection.scores[df$Sample.Type == "Interspace"]
  control.mean <- mean(control.proj, na.rm = TRUE)
  control.sd <- sd(control.proj, na.rm = TRUE)

  # If control SD is 0 (no variation), avoid divide by zero
  if (control.sd == 0) {
    df$Projection.std <- NA
  } else {
    df$Projection.std <- (projection.scores - control.mean) / control.sd }
  
  return(df)
}  

# Apply the function by site for Spring 2025 nmds ordination:
soil.conditioning <- bray_scores %>%
  group_by(SoilID) %>%
  group_modify(~ compute.projection(.x)) %>%
  ungroup() %>%
  select(SoilID, Sample.Type, Replicate,
         Projection, Projection.std) %>%
  ## Several observations per replicate, so calculate the median value:
  group_by(SoilID, Replicate, Sample.Type) %>%
  summarize(Projection = median(Projection),
            Projection.std = median(Projection.std))
#we also have to address the replicate combos in the projections to be able to match with the biomass data
projection.combos <- combo.key %>%
  left_join(
    soil.conditioning %>%
      select(
        SoilID,
        Sample.Type,
        rep1 = Replicate,
        Projection.1 = Projection
      ),
    by = "rep1"
  ) %>%
  left_join(
    soil.conditioning %>%
      select(
        SoilID,
        Sample.Type,
        rep2 = Replicate,
        Projection.2 = Projection
      ),
    by = c("SoilID", "Sample.Type", "rep2")
  ) %>%
  mutate(
    Projection = (Projection.1 + Projection.2) / 2
  ) %>%
  select(SoilID, Replicate, Sample.Type, Projection)

# Combine original projection values with estimated A-B and B-C projection values
projection.df <- bind_rows(
  soil.conditioning %>%
    select(SoilID, Replicate, Sample.Type, Projection),
  projection.combos
) %>%
  mutate(
    Rhizo.Int = recode(
      Sample.Type,
      "Interspace" = "INT",
      "Rhizosphere" = "RHIZO"
    )
  ) %>%
  select(SoilID, Replicate, Rhizo.Int, Projection)

biomass <- biomass %>%
  left_join(
    projection.df,
    by = c("SoilID", "Replicate", "Rhizo.Int")
  )

#next build your model:####
#scaling
biomass$alphadiv.std <- as.numeric(scale(biomass$alphadiv))
biomass$projection.std <- as.numeric(scale(biomass$Projection))

## Reformat bray-curtis dissimilarity matrix for generalized dissimilarity modeling: #####
## This just takes the pairwise comparisons and turns them to long form, pairing with each sample's metadata
dist.mat <- ps_bray %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  dplyr::rename(Extraction.ID = sample,
                CompareID = name) %>%
  left_join(sampledf[c("Extraction.ID", "SoilID", "Replicate", "Sample.Type")], 
            by="Extraction.ID") %>%
  dplyr::rename(Sample.x = Extraction.ID,
                Extraction.ID = CompareID) %>%
  left_join(sampledf[c("Extraction.ID", "SoilID", "Replicate", "Sample.Type")], 
            by="Extraction.ID") %>%
  dplyr::rename(Sample.y = Extraction.ID) 

#build model
#all together with scaling
bbio.microbe <- brm(
  dry.mass ~ Treatment * Rhizo.Int * alphadiv.std +
    Treatment * Rhizo.Int * Projection.std,
  data = biomass,
  family = Gamma(link = "log")
)

summary(bbio.microbe)
mcmc_plot(bbio.microbe, pars="b_")
cor(
  biomass$alphadiv.std,
  biomass$Projection,
  use = "complete.obs"
)

#figures
plot_predictions(
  bbio.microbe,
  condition = c("alphadiv.std", "Treatment", "Rhizo.Int")
)


#Probe richnesses
rich.model <- brm(
  Observed ~ SoilID,
  data = rich,
  family = gaussian()
)
summary(rich.model)

#figure for poster
# Figure A
rich.fig <- ggplot(
  rich,
  aes(x = SoilID, y = Observed, fill = SoilID)
) +
  geom_boxplot(
    width = 0.65,
    linewidth = 0.6,
    outlier.shape = 16,
    outlier.size = 2
  ) +
  scale_fill_manual(
    values = c(
      "BB" = "#8C6D46",
      "HG" = "#D98C3F",
      "NG" = "#D8C8A8",
      "OH" = "#6E9A63",
      "SR" = "#A66A3A"
    ),
    name = "Sagebrush population",
    labels = c(
      "BB" = "Bogus Basin",
      "HG" = "Hulls Gulch",
      "NG" = "Nancy Gulch",
      "OH" = "Owyhee High",
      "SR" = "Snake River"
    )
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
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
      size = 13,
      color = "#2F2A20"
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
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

ggsave(
  filename = "rich.fig.png",
  plot = rich.fig,
  
  width = 9,
  height = 5,
  dpi = 600,
  bg = "transparent"
)