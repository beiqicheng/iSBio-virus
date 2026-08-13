#!/usr/bin/env Rscript

# ============================================================
# 04_mantel_env_correlogram.R
# Fig. 2D — Mantel correlogram: vOTU / vFunction composition vs.
# individual environmental variables, plus pairwise environmental
# variable correlations (upper-triangle heatmap)
# ============================================================

library(tidyverse)
library(linkET)
library(vegan)
library(RColorBrewer)

# ----------------------------
# 1) Read data
# ----------------------------
env <- read.table("isbiometa.txt",
                  row.names = 1, header = TRUE, sep = "\t")

# count.txt has vOTUs as rows / samples as columns. Read it in that
# orientation first so we can filter by row (vOTU) prevalence, then
# transpose to samples-as-rows afterward (matching vfunc's orientation).
votu_raw <- read.table("count.txt",
                       row.names = 1, header = TRUE, sep = "\t")

# Core virome filter: keep only vOTUs present in >=10% of samples,
# matching the same prevalence filter used for the community composition
threshold <- ceiling(0.10 * ncol(votu_raw))
prevalence <- apply(votu_raw, 1, function(x) sum(x > 0))
votu_core <- votu_raw[prevalence >= threshold, ]
cat("vOTUs retained after 10% prevalence filter:", nrow(votu_core), "out of", nrow(votu_raw), "\n")

votu <- t(votu_core)

# eggnog3.txt is functional annotation output, not protein clusters —
# "vFunctions" is the accurate label used throughout (not "vPCs").
vf_raw <- read.table("eggnog.txt",
                     row.names = 1, header = TRUE, sep = "\t")
vfunc <- t(vf_raw)

# ----------------------------
# 2) Align samples
# ----------------------------
common_samples <- Reduce(intersect, list(rownames(env), rownames(votu), rownames(vfunc)))

env   <- env[common_samples, , drop = FALSE]
votu  <- votu[common_samples, , drop = FALSE]
vfunc <- vfunc[common_samples, , drop = FALSE]

# keep numeric environmental variables only
env <- env %>% dplyr::select(where(is.numeric))

# ----------------------------
# 3) Combine vOTUs + vFunctions
# ----------------------------
spec <- cbind(votu, vfunc)
spec_select <- list(
  vOTUs      = 1:ncol(votu),
  vFunctions = (ncol(votu) + 1):(ncol(votu) + ncol(vfunc))
)

# ----------------------------
# 4) Mantel test (one per environmental variable)
# ----------------------------
mantel_df <- mantel_test(spec, env, spec_select = spec_select) %>%
  mutate(
    rd = cut(r,
             breaks = c(-Inf, 0.2, 0.4, Inf),
             labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
    pd = cut(p,
             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
             labels = c("p <= 0.001", "0.001 < p <= 0.01", "0.01 < p <= 0.05", "p > 0.05"))
  ) %>%
  filter(p <= 0.05)

# ----------------------------
# 5) Correlation heatmap + Mantel links
# ----------------------------
p <- qcorrplot(correlate(env), type = "upper", diag = FALSE) +
  geom_square() +
  geom_couple(
    data = mantel_df,
    aes(colour = spec, size = rd),
    curvature = nice_curvature()
  ) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    limits = c(-1, 1)
  ) +
  scale_size_manual(values = c(0.5, 1, 1.5)) +
  scale_colour_manual(values = c(
    "vOTUs"      = "#D73027",
    "vFunctions" = "#1B9E77"
  )) +
  guides(
    size   = guide_legend(title = "Mantel's r"),
    colour = guide_legend(title = "Community type"),
    fill   = guide_colorbar(title = "Spearman's rho")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  )

print(p)
ggsave("Fig2D_mantel_correlogram.png", p, width = 8, height = 6, dpi = 300)