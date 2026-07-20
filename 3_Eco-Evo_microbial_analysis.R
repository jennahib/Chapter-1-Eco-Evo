# Chapter 1 Phase 1: Site, interspace, and rhizosphere comparisons
# Author: Jennah Brown
# Date: July 2026
#
# Purpose:
#   1. Compare fungal communities among sites.
#   2. Compare fungal communities between interspace and rhizosphere soils.
#   3. Visualize community composition using NMDS.
#   4. Calculate site-specific interspace rhizosphere projection scores.
#   5. Summarize abundant fungal genera.
#   6. Prepare fungal community data for comparison with BRTE performance.
#
# Data input:
# This analysis begins with a rarefied ITS phyloseq object generated from the BRC-BIO 1_Phyloseq_object_prep_ITS script preprocessing pipeline. The object contains all BRCBIO Spring 2025 samples,

library(phyloseq)
library(tidyverse)
library(vegan)
library(DESeq2)

# Create output folders
dir.create(
  "outputs/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# -----------------------------
# Import and subset data
# -----------------------------

# Import rarefied Spring 2025 fungal community samples
ps.spring <- readRDS(
  "Data/phyloseq_rare_ITS_BRCBIO_Spring2025.rds"
)

# Filter for sites
ps.focal <- subset_samples(
  ps.spring,
  SoilID %in% c("BB", "HG", "NG", "OH", "SR")
)

#Remove taxa with zero abundance after subsetting to the five sites
ps.focal <- prune_taxa(
  taxa_sums(ps.focal) > 0,
  ps.focal
)

# Check sample counts by site and soil source
dplyr::count(
  as.data.frame(sample_data(ps.focal)),
  SoilID,
  Sample.Type
)
# -----------------------------
# Sequencing depth checks
# -----------------------------
# Confirm equal sequencing depth after rarefaction
summary(sample_sums(ps.focal))

# -----------------------------
# Bray-Curtis dissimilarity
# -----------------------------

# Convert counts to relative abundance
ps.rel <- transform_sample_counts(
  ps.focal,
  function(x) x / sum(x)
)

bray <- phyloseq::distance(
  ps.rel,
  method = "bray"
)

meta <- data.frame(sample_data(ps.rel))

# Ensure metadata order matches distance matrix
meta <- meta[
  labels(bray),
  ,
  drop = FALSE
]

# -----------------------------
# Just looking at PERMANOVA
# -----------------------------

set.seed(1)

# testing site, soil source, and their interaction
permanova.results <- adonis2(
  bray ~ SoilID * Sample.Type,
  data = meta,
  permutations = 999
)

permanova.results

# -----------------------------
# Multivariate dispersion checks
# -----------------------------

# Check dispersion among sites
site.dispersion <- betadisper(
  bray,
  group = factor(meta$SoilID)
)

anova(site.dispersion)

set.seed(1)

site.dispersion.test <- permutest(
  site.dispersion,
  permutations = 999
)

site.dispersion.test

# Check dispersion between interspace and rhizosphere soils
source.dispersion <- betadisper(
  bray,
  group = factor(meta$Sample.Type)
)

anova(source.dispersion)

set.seed(1)

source.dispersion.test <- permutest(
  source.dispersion,
  permutations = 999
)

source.dispersion.test

# Check dispersion among site-by-source groups
site.source.group <- interaction(
  meta$SoilID,
  meta$Sample.Type,
  drop = TRUE
)

site.source.dispersion <- betadisper(
  bray,
  group = site.source.group
)

anova(site.source.dispersion)

set.seed(1)

site.source.dispersion.test <- permutest(
  site.source.dispersion,
  permutations = 999
)

site.source.dispersion.test

# -----------------------------
# NMDS ordination
# -----------------------------

# Match population colors used in the biomass figure
population.colors <- c(
  "BB" = "#E76F51",
  "HG" = "#A3A500",
  "NG" = "#009E73",
  "OH" = "#00A6D6",
  "SR" = "#CC79E0"
)

set.seed(1)

nmds <- ordinate(
  ps.rel,
  method = "NMDS",
  distance = "bray",
  trymax = 100
)

nmds

nmds.fig <- plot_ordination(
  ps.rel,
  nmds,
  color = "SoilID",
  shape = "Sample.Type"
) +
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  stat_ellipse(
    aes(
      group = SoilID,
      color = SoilID
    ),
    linewidth = 0.7,
    linetype = 2,
    alpha = 0.6,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = population.colors
  ) +
  scale_shape_manual(
    name = "Soil inoculum source",
    values = c(
      "Interspace" = 16,
      "Rhizosphere" = 17
    ),
    breaks = c(
      "Interspace",
      "Rhizosphere"
    )
  ) +
  labs(
    x = "NMDS1",
    y = "NMDS2"
  ) +
  guides(
    color = "none",
    shape = guide_legend(
      override.aes = list(
        color = "#2F2A20",
        size = 4,
        alpha = 1
      )
    )
  ) +
  theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(
      size = 20,
      face = "bold",
      color = "#2F3A1F"
    ),
    axis.text = element_text(
      size = 15,
      color = "#2F2A20"
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
    )
  )

nmds.fig

ggsave(
  filename = "outputs/figures/nmds_poster.png",
  plot = nmds.fig,
  width = 7,
  height = 5.5,
  units = "in",
  dpi = 600,
  bg = "transparent"
)

# -----------------------------
# Projection scores
# -----------------------------

# use NMDS site scores and sample metadata
scores.df <- plot_ordination(
  ps.rel,
  nmds,
  justDF = TRUE
)

bray.scores <- scores.df %>%
  select(
    SoilID,
    Sample.Type,
    Replicate,
    Subsample,
    Extraction.ID,
    NMDS1,
    NMDS2
  )

compute_projection <- function(df) {
  
  if (!all(c("Interspace", "Rhizosphere") %in% df$Sample.Type)) {
    df$Projection <- NA_real_
    return(df)
  }
  
  # Calculate interspace and rhizosphere centroids
  inter.centroid <- colMeans(
    df[
      df$Sample.Type == "Interspace",
      c("NMDS1", "NMDS2")
    ],
    na.rm = TRUE
  )
  
  rhizo.centroid <- colMeans(
    df[
      df$Sample.Type == "Rhizosphere",
      c("NMDS1", "NMDS2")
    ],
    na.rm = TRUE
  )
  
  # Define the interspace-to-rhizosphere direction
  direction.vector <- rhizo.centroid - inter.centroid
  direction.length <- sqrt(sum(direction.vector^2))
  
  if (
    is.na(direction.length) ||
    direction.length == 0
  ) {
    df$Projection <- NA_real_
    return(df)
  }
  
  direction.unit <- direction.vector / direction.length
  
  # Center observations on the interspace centroid
  centered <- sweep(
    df[, c("NMDS1", "NMDS2")],
    MARGIN = 2,
    STATS = inter.centroid
  )
  
  centered <- as.matrix(centered)
  
  # Project observations onto the site-specific direction vector
  df$Projection <- as.vector(
    centered %*% direction.unit
  )
  
  df
}

projection.df <- bray.scores %>%
  group_by(SoilID) %>%
  group_modify(
    ~ compute_projection(.x)
  ) %>%
  ungroup()

# Summarize subsamples within each replicate
projection.rep <- projection.df %>%
  group_by(
    SoilID,
    Replicate,
    Sample.Type
  ) %>%
  summarise(
    Projection = if (
      all(is.na(Projection))
    ) {
      NA_real_
    } else {
      median(
        Projection,
        na.rm = TRUE
      )
    },
    .groups = "drop"
  )

# Plot site projection scores
projection.fig <- ggplot(
  projection.rep,
  aes(
    x = SoilID,
    y = Projection,
    fill = Sample.Type
  )
) +
  geom_boxplot(
    alpha = 0.8,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.5
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Site",
    y = "Projection score",
    fill = "Soil source"
  )

projection.fig

# -----------------------------
# Genus-level relative abundance
# -----------------------------

ps.genus <- tax_glom(
  ps.focal,
  taxrank = "Genus"
)

# Convert genus counts to relative abundance
ps.genus.rel <- transform_sample_counts(
  ps.genus,
  function(x) x / sum(x)
)

# Calculate average relative abundance by genus
genus.abund <- psmelt(ps.genus.rel) %>%
  filter(
    !is.na(Genus),
    Genus != ""
  ) %>%
  group_by(Genus) %>%
  summarise(
    MeanAbundance = mean(Abundance),
    .groups = "drop"
  ) %>%
  arrange(
    desc(MeanAbundance)
  )

# Identify the 10 most abundant genera
top10 <- genus.abund %>%
  slice_head(n = 10) %>%
  pull(Genus)

# Calculate average abundance by site
heat.df <- psmelt(ps.genus.rel) %>%
  filter(
    Genus %in% top10
  ) %>%
  group_by(
    SoilID,
    Genus
  ) %>%
  summarise(
    MeanAbundance = mean(Abundance),
    .groups = "drop"
  )

heatmap.fig <- ggplot(
  heat.df,
  aes(
    x = SoilID,
    y = Genus,
    fill = MeanAbundance
  )
) +
  geom_tile(
    color = "white"
  ) +
  scale_fill_viridis_c(
    name = "Relative\nabundance"
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Site",
    y = "Genus"
  )

heatmap.fig

# -----------------------------
# probe for BRTE +/- microbes
# -----------------------------

# placeholder sites classified as suppressive or permissive
ps.deseq <- subset_samples(
  ps.genus,
  SoilID %in% c("BB", "HG", "NG", "OH")
)

# Remove taxa absent after subsetting
ps.deseq <- prune_taxa(
  taxa_sums(ps.deseq) > 0,
  ps.deseq
)

# Assign site-level comparison groups
sample_data(ps.deseq)$Group <- ifelse(
  sample_data(ps.deseq)$SoilID %in% c("HG", "NG"),
  "Suppressive",
  "Permissive"
)

sample_data(ps.deseq)$Group <- factor(
  sample_data(ps.deseq)$Group,
  levels = c(
    "Permissive",
    "Suppressive"
  )
)

# Check sample counts in each group
table(
  sample_data(ps.deseq)$Group
) #perhaps not as relevant