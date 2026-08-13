#!/usr/bin/env Rscript

# ============================================================
# 02_dbrda_drivers.R
# dbRDA visualization of environmental drivers — Forest vs Grassland
# ============================================================

library(phyloseq)
library(vegan)
library(ggplot2)

# ------------------------------
# 1. Load data & build core virome (same logic as 01_alpha_beta_diversity.R)
# ------------------------------
otu <- read.table("count.txt",
                  row.names = 1, header = TRUE, sep = "\t")
meta <- read.table("isbiometa.txt",
                   row.names = 1, header = TRUE, sep = "\t")

OTU <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
SAM <- sample_data(meta)
ps  <- phyloseq(OTU, SAM)

threshold <- ceiling(0.10 * nsamples(ps))
prevalence <- apply(otu_table(ps), 1, function(x) sum(x > 0))
ps_core <- prune_taxa(prevalence >= threshold, ps)
ps_rel  <- transform_sample_counts(ps_core, function(x) x / sum(x))

ps_forest <- subset_samples(ps_rel, Biome2 == "Forest")
ps_grass  <- subset_samples(ps_rel, Biome2 == "Grassland")

dist_forest <- phyloseq::distance(ps_forest, method = "bray")
dist_grass  <- phyloseq::distance(ps_grass, method = "bray")

forest_meta <- as(sample_data(ps_forest), "data.frame")
grass_meta  <- as(sample_data(ps_grass),  "data.frame")

# ------------------------------
# 2. Forest dbRDA
# ------------------------------
# sand and nitrogen dropped due to multicollinearity; soc listed once
dbrda_forest <- dbrda(dist_forest ~ ph + mat + soc + swc + clay + silt + cn + map + elevation,
                      data = forest_meta)
anova(dbrda_forest, permutations = 999)
anova(dbrda_forest, by = "term", permutations = 999)

forest_var <- round(100 * eigenvals(dbrda_forest) / sum(eigenvals(dbrda_forest)), 1)

sites_f <- as.data.frame(scores(dbrda_forest, display = "sites"))
sites_f$Biome <- "Forest"
env_f <- as.data.frame(scores(dbrda_forest, display = "bp"))
env_f$var <- rownames(env_f)

p_forest <- ggplot() +
  geom_point(data = sites_f, aes(x = dbRDA1, y = dbRDA2),
             color = "#228B22", size = 3, alpha = 0.7) +
  geom_segment(data = env_f, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 0.8, color = "black") +
  geom_text(data = env_f, aes(x = dbRDA1 * 1.15, y = dbRDA2 * 1.15, label = var), size = 4) +
  labs(title = "dbRDA of Viral Communities (Forest)",
       x = paste0("dbRDA1 (", forest_var[1], "%)"),
       y = paste0("dbRDA2 (", forest_var[2], "%)")) +
  theme_bw(base_size = 14)
p_forest
ggsave("dbRDA_forest.png", p_forest, width = 6, height = 5, dpi = 300)

# ------------------------------
# 3. Grassland dbRDA
# ------------------------------
dbrda_grass <- dbrda(dist_grass ~ ph + mat + soc + swc + clay + silt + cn + map + elevation,
                     data = grass_meta)
anova(dbrda_grass, permutations = 999)
anova(dbrda_grass, by = "term", permutations = 999)

grass_var <- round(100 * eigenvals(dbrda_grass) / sum(eigenvals(dbrda_grass)), 1)

sites_g <- as.data.frame(scores(dbrda_grass, display = "sites"))
sites_g$Biome <- "Grassland"
env_g <- as.data.frame(scores(dbrda_grass, display = "bp"))
env_g$var <- rownames(env_g)

p_grass <- ggplot() +
  geom_point(data = sites_g, aes(x = dbRDA1, y = dbRDA2),
             color = "#DAA520", size = 3, alpha = 0.7) +
  geom_segment(data = env_g, aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               arrow = arrow(length = unit(0.25, "cm")),
               linewidth = 0.8, color = "black") +
  geom_text(data = env_g, aes(x = dbRDA1 * 1.15, y = dbRDA2 * 1.15, label = var), size = 4) +
  labs(title = "dbRDA of Viral Communities (Grassland)",
       x = paste0("dbRDA1 (", grass_var[1], "%)"),
       y = paste0("dbRDA2 (", grass_var[2], "%)")) +
  theme_bw(base_size = 14)
p_grass
ggsave("dbRDA_grassland.png", p_grass, width = 6, height = 5, dpi = 300)