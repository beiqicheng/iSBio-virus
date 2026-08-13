#!/usr/bin/env Rscript

# ============================================================
# 01_alpha_beta_diversity.R
# Alpha diversity, soil-property correlations, and beta diversity
# (PCoA/NMDS/PERMANOVA) for total samples and forest/grassland subsets
# ============================================================

library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggcorrplot)
library(Hmisc)

# ------------------------------
# 1. Load OTU table & metadata
# ------------------------------
otu <- read.table("count.txt",
                  row.names = 1, header = TRUE, sep = "\t")

meta <- read.table("isbiometa.txt",
                   row.names = 1, header = TRUE, sep = "\t")

OTU <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
SAM <- sample_data(meta)

ps <- phyloseq(OTU, SAM)

# ------------------------------
# 2. Alpha diversity
# ------------------------------
alpha_div <- estimate_richness(ps, measures = c("Shannon", "Observed"))
alpha_div$SampleID <- rownames(alpha_div)
alpha_div2 <- merge(alpha_div, meta, by.x = "SampleID", by.y = "row.names")
write.table(alpha_div2, "alpha_diversity_results.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Plot Shannon diversity by Biome
p_alpha <- ggplot(alpha_div2, aes(x = Biome1, y = Shannon, fill = Biome1)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, alpha = 0.4) +
  theme_bw() + xlab("Biome") + ylab("Shannon diversity")

p_alpha2 <- ggplot(alpha_div2, aes(x = Biome1, y = Observed, fill = Biome1)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, alpha = 0.4) +
  theme_bw() + xlab("Biome") + ylab("Observed diversity")

ggsave("alpha_diversity_boxplot.png", p_alpha, width = 6, height = 4, dpi = 300)
ggsave("alpha_diversity_observed_boxplot.png", p_alpha2, width = 6, height = 4, dpi = 300)

# ------------------------------
# 3. All samples correlations
# ------------------------------
numeric_vars <- c("Shannon", "elevation", "mat", "map", "clay", "silt",
                  "sand", "ph", "soc", "nitrogen", "cn", "swc")

soil_subset <- alpha_div2[, numeric_vars]
soil_subset <- as.data.frame(lapply(soil_subset, function(x) as.numeric(as.character(x))))

cor_matrix <- cor(soil_subset, method = "spearman", use = "pairwise.complete.obs")

p_corr <- ggcorrplot(cor_matrix,
                     hc.order = TRUE,
                     type = "lower",
                     lab = TRUE,
                     title = "Spearman Correlation: Shannon vs Soil Properties",
                     colors = c("#6D9EC1", "white", "#E46726"))
print(p_corr)
ggsave("correlation_heatmap.png", p_corr, width = 7, height = 6, dpi = 300)

# ------------------------------
# 4. Forest and grassland samples correlations
# ------------------------------
comparison_data <- alpha_div2 %>% filter(Biome2 %in% c("Forest", "Grassland"))

p1 <- ggplot(comparison_data, aes(x = ph, y = Shannon, color = Biome2)) +
  geom_point(alpha = 0.4) + geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman") + facet_wrap(~Biome2) + theme_bw()

p2 <- ggplot(comparison_data, aes(x = silt, y = Shannon, color = Biome2)) +
  geom_point(alpha = 0.4) + geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman") + facet_wrap(~Biome2) + theme_bw()

p3 <- ggplot(comparison_data, aes(x = mat, y = Shannon, color = Biome2)) +
  geom_point(alpha = 0.4) + geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman") + facet_wrap(~Biome2) + theme_bw()

p4 <- ggplot(comparison_data, aes(x = swc, y = Shannon, color = Biome2)) +
  geom_point(alpha = 0.4) + geom_smooth(method = "lm", se = TRUE) +
  stat_cor(method = "spearman") + facet_wrap(~Biome2) + theme_bw()

p_corr_panels <- ggarrange(p1, p2, p3, p4, ncol = 1)
print(p_corr_panels)
ggsave("Shannon_soil_correlations.png", p_corr_panels, width = 6, height = 16, dpi = 300)

# ------------------------------
# 5. Total samples beta diversity (PCoA / NMDS)
# ------------------------------
# Filtering to a "core virome" of vOTUs present in >=10% of samples
threshold <- ceiling(0.10 * nsamples(ps))
prevalence <- apply(otu_table(ps), 1, function(x) sum(x > 0))
ps_core <- prune_taxa(prevalence >= threshold, ps)

# Distance matrix (Bray-Curtis)
ps_rel <- transform_sample_counts(ps_core, function(x) x / sum(x))
bray_dist <- phyloseq::distance(ps_rel, method = "bray")

# PCoA
ord_pcoa <- ordinate(ps_rel, method = "PCoA", distance = bray_dist)
p_pcoa <- plot_ordination(ps_rel, ord_pcoa, color = "Biome1") +
  geom_point(size = 2) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("PCoA_bray.png", p_pcoa, width = 6, height = 4, dpi = 300)

# NMDS
ord_nmds <- ordinate(ps_rel, method = "NMDS", distance = bray_dist, trymax = 100)
p_nmds <- plot_ordination(ps_rel, ord_nmds, color = "Biome1") +
  geom_point(size = 4) + theme_bw()
ggsave("NMDS_bray.png", p_nmds, width = 6, height = 4, dpi = 300)

# ------------------------------
# 6. PERMANOVA (adonis2) — total samples
# ------------------------------
meta_rel <- as(sample_data(ps_rel), "data.frame")
adonis_res <- adonis2(as.matrix(bray_dist) ~ Biome1, data = meta_rel)
write.table(adonis_res, file = "permanova_results.txt", sep = "\t", quote = FALSE)

# ------------------------------
# 7. Forest and grassland beta diversity (PCoA / NMDS) + PERMANOVA
# ------------------------------
ps_forest <- subset_samples(ps_rel, Biome2 == "Forest")
ps_grass  <- subset_samples(ps_rel, Biome2 == "Grassland")

dist_forest <- phyloseq::distance(ps_forest, method = "bray")
dist_grass  <- phyloseq::distance(ps_grass, method = "bray")

forest_meta <- as(sample_data(ps_forest), "data.frame")
grass_meta  <- as(sample_data(ps_grass),  "data.frame")

# NOTE ON PREDICTOR SETS: this is a reduced, biome-specific variable set

perm_forest <- adonis2(dist_forest ~ ph + mat + soc + swc,
                       data = forest_meta,
                       permutations = 999)

perm_grass <- adonis2(dist_grass ~ silt + sand + swc + ph,
                      data = grass_meta,
                      permutations = 999)

print(perm_forest)
print(perm_grass)