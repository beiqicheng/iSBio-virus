#!/usr/bin/env Rscript

# ============================================================
# 05_mantel_distance_decay.R
# Geographic and elevational distance-decay: Mantel & partial Mantel tests
# ============================================================

library(phyloseq)
library(vegan)
library(geosphere)
library(ggplot2)
library(scales)

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

forest_meta <- as(sample_data(ps_forest), "data.frame")
grass_meta  <- as(sample_data(ps_grass),  "data.frame")

# ------------------------------
# 2. Geographic + environmental Mantel suite
# ------------------------------
run_mantel_suite <- function(ps_obj, meta_df, env_vars, label) {
  comm_dist <- phyloseq::distance(ps_obj, method = "bray")
  
  # Geographic distance — great-circle distance in km, not Euclidean on
  # raw lat/long, since samples span multiple continents
  coords <- as.matrix(meta_df[, c("longitude", "latitude")])
  geo_dist <- as.dist(distm(coords, fun = distHaversine) / 1000)  # km
  
  # Environmental distance — Euclidean on z-scored soil/climate variables
  # (same variable set as dbRDA, so results are directly comparable)
  env_scaled <- scale(meta_df[, env_vars])
  env_dist <- dist(env_scaled, method = "euclidean")
  
  stopifnot(attr(comm_dist, "Size") == attr(geo_dist, "Size"),
            attr(comm_dist, "Size") == attr(env_dist, "Size"))
  
  mantel_geo <- mantel(comm_dist, geo_dist, method = "spearman", permutations = 999)
  mantel_env <- mantel(comm_dist, env_dist, method = "spearman", permutations = 999)
  pmantel_geo_given_env <- mantel.partial(comm_dist, geo_dist, env_dist,
                                          method = "spearman", permutations = 999)
  pmantel_env_given_geo <- mantel.partial(comm_dist, env_dist, geo_dist,
                                          method = "spearman", permutations = 999)
  mantel_geo_env <- mantel(geo_dist, env_dist, method = "spearman", permutations = 999)
  
  cat("\n===========================================\n", label, "\n===========================================\n")
  cat("Simple Mantel, community ~ geography:   r =", round(mantel_geo$statistic, 3), " p =", mantel_geo$signif, "\n")
  cat("Simple Mantel, community ~ environment: r =", round(mantel_env$statistic, 3), " p =", mantel_env$signif, "\n")
  cat("Partial Mantel, geography | environment: r =", round(pmantel_geo_given_env$statistic, 3), " p =", pmantel_geo_given_env$signif, "\n")
  cat("Partial Mantel, environment | geography: r =", round(pmantel_env_given_geo$statistic, 3), " p =", pmantel_env_given_geo$signif, "\n")
  cat("Geography-environment collinearity:      r =", round(mantel_geo_env$statistic, 3), " p =", mantel_geo_env$signif, "\n")
  
  list(label = label, simple_geo = mantel_geo, simple_env = mantel_env,
       partial_geo_given_env = pmantel_geo_given_env,
       partial_env_given_geo = pmantel_env_given_geo,
       geo_env_collinearity = mantel_geo_env,
       comm_dist = comm_dist, geo_dist = geo_dist, env_dist = env_dist)
}

# sand/nitrogen excluded for multicollinearity, matching dbRDA (02_dbrda_drivers.R)
env_vars <- c("ph", "mat", "map", "elevation", "soc", "swc", "clay", "silt", "cn")

mantel_forest <- run_mantel_suite(ps_forest, forest_meta, env_vars, "Forest")
mantel_grass  <- run_mantel_suite(ps_grass,  grass_meta,  env_vars, "Grassland")

# ------------------------------
# 3. Elevation-specific Mantel suite
# ------------------------------
run_elevation_mantel <- function(comm_dist, geo_dist, meta_df, other_env_vars, label) {
  elev_dist <- dist(meta_df$elevation, method = "euclidean")
  other_env_scaled <- scale(meta_df[, other_env_vars])
  other_env_dist <- dist(other_env_scaled, method = "euclidean")
  
  stopifnot(attr(comm_dist, "Size") == attr(elev_dist, "Size"))
  
  mantel_elev <- mantel(comm_dist, elev_dist, method = "spearman", permutations = 999)
  pmantel_elev_given_geo <- mantel.partial(comm_dist, elev_dist, geo_dist,
                                           method = "spearman", permutations = 999)
  pmantel_elev_given_env <- mantel.partial(comm_dist, elev_dist, other_env_dist,
                                           method = "spearman", permutations = 999)
  
  cat("\n===========================================\n", label, "- elevation effect\n===========================================\n")
  cat("Simple Mantel, community ~ elevation distance:   r =", round(mantel_elev$statistic, 3), " p =", mantel_elev$signif, "\n")
  cat("Partial Mantel, elevation | geography:           r =", round(pmantel_elev_given_geo$statistic, 3), " p =", pmantel_elev_given_geo$signif, "\n")
  cat("Partial Mantel, elevation | rest of environment: r =", round(pmantel_elev_given_env$statistic, 3), " p =", pmantel_elev_given_env$signif, "\n")
  
  list(label = label, elev_dist = elev_dist, comm_dist = comm_dist,
       simple = mantel_elev, partial_given_geo = pmantel_elev_given_geo,
       partial_given_env = pmantel_elev_given_env)
}

# Same env set as above, minus elevation itself
other_env_vars <- c("ph", "mat", "map", "soc", "swc", "clay", "silt", "cn")

elev_forest <- run_elevation_mantel(mantel_forest$comm_dist, mantel_forest$geo_dist,
                                    forest_meta, other_env_vars, "Forest")
elev_grass  <- run_elevation_mantel(mantel_grass$comm_dist, mantel_grass$geo_dist,
                                    grass_meta, other_env_vars, "Grassland")

# ------------------------------
# 4. Plots
# ------------------------------
make_decay_df <- function(mantel_res, biome_label) {
  data.frame(geo_km = as.vector(mantel_res$geo_dist), bray = as.vector(mantel_res$comm_dist), Biome = biome_label)
}
decay_df <- rbind(make_decay_df(mantel_forest, "Forest"), make_decay_df(mantel_grass, "Grassland"))

p_decay <- ggplot(decay_df, aes(x = geo_km, y = bray, color = Biome)) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_x_continuous(trans = "log1p", labels = scales::comma, breaks = c(0, 10, 100, 1000, 10000)) +
  scale_color_manual(values = c("Forest" = "#228B22", "Grassland" = "#DAA520")) +
  labs(x = "Geographic distance (km, log scale)", y = "Bray-Curtis dissimilarity",
       title = "Distance-decay of soil viral community similarity") +
  theme_bw(base_size = 14)
print(p_decay)
ggsave("distance_decay.png", p_decay, width = 7, height = 5, dpi = 300)

make_elev_df <- function(res, biome_label) {
  data.frame(elev_diff_m = as.vector(res$elev_dist), bray = as.vector(res$comm_dist), Biome = biome_label)
}
elev_df <- rbind(make_elev_df(elev_forest, "Forest"), make_elev_df(elev_grass, "Grassland"))

p_elev <- ggplot(elev_df, aes(x = elev_diff_m, y = bray, color = Biome)) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_x_continuous(labels = scales::comma) +
  scale_color_manual(values = c("Forest" = "#228B22", "Grassland" = "#DAA520")) +
  labs(x = "Elevation difference between sites (m)", y = "Bray-Curtis dissimilarity",
       title = "Elevation-decay of soil viral community similarity") +
  theme_bw(base_size = 14)
print(p_elev)
ggsave("elevation_decay.png", p_elev, width = 7, height = 5, dpi = 300)

# ------------------------------
# 5. Save numeric summaries
# ------------------------------
mantel_summary <- data.frame(
  Biome = c("Forest", "Grassland"),
  Mantel_r_geo = c(mantel_forest$simple_geo$statistic, mantel_grass$simple_geo$statistic),
  Mantel_p_geo = c(mantel_forest$simple_geo$signif, mantel_grass$simple_geo$signif),
  Mantel_r_env = c(mantel_forest$simple_env$statistic, mantel_grass$simple_env$statistic),
  Mantel_p_env = c(mantel_forest$simple_env$signif, mantel_grass$simple_env$signif),
  PartialMantel_r_geo_given_env = c(mantel_forest$partial_geo_given_env$statistic, mantel_grass$partial_geo_given_env$statistic),
  PartialMantel_p_geo_given_env = c(mantel_forest$partial_geo_given_env$signif, mantel_grass$partial_geo_given_env$signif),
  PartialMantel_r_env_given_geo = c(mantel_forest$partial_env_given_geo$statistic, mantel_grass$partial_env_given_geo$statistic),
  PartialMantel_p_env_given_geo = c(mantel_forest$partial_env_given_geo$signif, mantel_grass$partial_env_given_geo$signif),
  GeoEnv_collinearity_r = c(mantel_forest$geo_env_collinearity$statistic, mantel_grass$geo_env_collinearity$statistic),
  GeoEnv_collinearity_p = c(mantel_forest$geo_env_collinearity$signif, mantel_grass$geo_env_collinearity$signif)
)
write.table(mantel_summary, "mantel_results_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
print(mantel_summary)

elev_summary <- data.frame(
  Biome = c("Forest", "Grassland"),
  Mantel_r_elev = c(elev_forest$simple$statistic, elev_grass$simple$statistic),
  Mantel_p_elev = c(elev_forest$simple$signif, elev_grass$simple$signif),
  PartialMantel_r_elev_given_geo = c(elev_forest$partial_given_geo$statistic, elev_grass$partial_given_geo$statistic),
  PartialMantel_p_elev_given_geo = c(elev_forest$partial_given_geo$signif, elev_grass$partial_given_geo$signif),
  PartialMantel_r_elev_given_env = c(elev_forest$partial_given_env$statistic, elev_grass$partial_given_env$statistic),
  PartialMantel_p_elev_given_env = c(elev_forest$partial_given_env$signif, elev_grass$partial_given_env$signif)
)
write.table(elev_summary, "mantel_elevation_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)
print(elev_summary)