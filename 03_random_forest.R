#!/usr/bin/env Rscript

# ============================================================
# 03_random_forest.R
# Random forest variable importance for viral community structure
# (PCoA1), Total samples / Forest / Grassland
# ============================================================

library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(randomForest)
library(rfPermute)

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

ord_pcoa <- ordinate(ps_rel, method = "PCoA", distance = "bray")

# ------------------------------
# 2. RF importance function
# ------------------------------

run_rf_importance <- function(ps_obj, pcoa_res, title_name) {
  y_var <- pcoa_res$vectors[, 1]
  env_vars <- c("ph", "mat", "map", "elevation", "soc", "cn", "clay", "silt", "swc")
  x_data <- as.data.frame(sample_data(ps_obj))[, env_vars]
  x_data <- as.data.frame(lapply(x_data, function(x) as.numeric(as.character(x))))
  valid_idx <- complete.cases(x_data)
  x_clean <- x_data[valid_idx, ]
  y_clean <- y_var[valid_idx]
  rf_model <- rfPermute(y_clean ~ ., data = x_clean, ntree = 1000, importance = TRUE)
  imp_tab <- as.data.frame(importance(rf_model, scale = TRUE))
  imp_tab$Variable <- rownames(imp_tab)
  imp_tab$Dataset <- title_name
  return(list(model = rf_model, importance = imp_tab))
}

# ------------------------------
# 3. Run RF for total samples, forest, and grassland
# ------------------------------
rf_all <- run_rf_importance(ps_rel, ord_pcoa, "Total Samples")

ps_forest <- subset_samples(ps_rel, Biome2 == "Forest")
ord_pcoa_forest <- ordinate(ps_forest, method = "PCoA", distance = "bray")
rf_forest <- run_rf_importance(ps_forest, ord_pcoa_forest, "Forest")

ps_grass <- subset_samples(ps_rel, Biome2 == "Grassland")
ord_pcoa_grass <- ordinate(ps_grass, method = "PCoA", distance = "bray")
rf_grass <- run_rf_importance(ps_grass, ord_pcoa_grass, "Grassland")

all_imp <- rbind(rf_all$importance, rf_forest$importance, rf_grass$importance)

# ------------------------------
# 4. Plot
# ------------------------------

all_imp <- all_imp %>%
  group_by(Dataset) %>%
  mutate(Variable = factor(Variable, levels = Variable[order(`%IncMSE`)])) %>%
  ungroup()

p_rf <- ggplot(all_imp, aes(x = Variable, y = `%IncMSE`, fill = Dataset)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  facet_wrap(~Dataset, scales = "free_y") +
  theme_bw() +
  labs(x = "Soil/Environmental Factors",
       y = "Importance (% Increase in MSE)",
       title = "Variable Importance for Viral Community Structure (PCoA1)") +
  scale_fill_manual(values = c("Total Samples" = "grey30", "Forest" = "#228B22", "Grassland" = "#DAA520"))

print(p_rf)
ggsave("RF_Variable_Importance.png", p_rf, width = 10, height = 6, dpi = 300)