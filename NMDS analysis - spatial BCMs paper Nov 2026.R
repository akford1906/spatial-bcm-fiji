# Script for analysing spatial BCM data in Fiji

library(dplyr)
library(ggplot2)
library(gtsummary)
library(dplyr)
library(vegan)

# Import reef data as dat
dat <- read.csv(file.choose()) 

# Organise data for NMDS
dat_nmds <- dat %>%
  mutate(
    Seascape = case_when(
      Province %in% c("Lomaiviti", "Bua", "Ra", "Tailevu") ~ "VIR_Seascape",
      Province == "Lau" ~ "LAU_Seascape",
      TRUE ~ NA_character_),
    BS  = `Bare.Substrate`, # Rename benthic variables with spaces
    OI  = `Other.Inverts`,
    CCA = Crustose) %>%
  filter(!is.na(Seascape))

table(dat_nmds$Seascape) # Check site numbers
str(dat_nmds)
head(dat_nmds)

# VATU-I-RA: 2D and 3D NMDS of benthic community composition
# Prepare data

dat_vir <- dat_nmds %>%
  filter(Seascape == "VIR_Seascape") %>%
  mutate(BCM_pct = as.numeric(Cyanobacteria))

benthic_cols <- c(
  "HC", "TA", "MA", "CCA", "SC",
  "BS", "OI", "SD", "RB")

# Extract benthic community data
benthic_vir <- dat_vir[, benthic_cols, drop = FALSE]

# Ensure all benthic variables are numeric
benthic_vir[] <- lapply(benthic_vir, as.numeric)

# Standardise non-BCM benthic data
# Calculate total cover represented by the non-BCM benthic
# categories at each site
row_tot <- rowSums(benthic_vir, na.rm = TRUE)

# Retain only sites with:
#   1. complete benthic composition data
#   2. a BCM-cover value
#   3. non-zero total non-BCM benthic cover
keep <- complete.cases(benthic_vir, dat_vir$BCM_pct) &
  row_tot > 0
benthic_vir_keep <- benthic_vir[keep, , drop = FALSE]
dat_vir_keep <- dat_vir[keep, , drop = FALSE]

# Convert benthic cover to relative composition after excluding BCM.
# Now each row sums to 1.
benthic_rel_vir <- benthic_vir_keep / row_tot[keep]
nrow(benthic_rel_vir) # Number of sites included in the NMDS

# Run 2D NMDS
set.seed(123)
nmds_vir_2d <- vegan::metaMDS(
  benthic_rel_vir,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE)

# Run 3D NMDS
set.seed(123)
nmds_vir_3d <- vegan::metaMDS(
  benthic_rel_vir,
  distance = "bray",
  k = 3,
  trymax = 200,
  autotransform = FALSE)

# Compare ordination stress
stress_comparison <- data.frame(
  Dimensions = c("2-D", "3-D"),
  Stress = c(
    nmds_vir_2d$stress,
    nmds_vir_3d$stress))

stress_comparison

# Current result:
# 2-D stress = 0.204
# 3-D stress = 0.126
# Interpretation: Adding a third dimension substantially improves the representation
# of the Bray-Curtis dissimilarities among Vatu-i-Ra reef sites.

# Fit BCM cover to the ordinations
set.seed(123)
bcm_fit_vir_2d <- vegan::envfit(
  nmds_vir_2d,
  dat_vir_keep$BCM_pct,
  permutations = 999,
  choices = 1:2)

set.seed(123)
bcm_fit_vir_3d <- vegan::envfit(
  nmds_vir_3d,
  dat_vir_keep$BCM_pct,
  permutations = 999,
  choices = 1:3)

bcm_fit_vir_2d
bcm_fit_vir_3d

# Current results:
# 2-D:
# r2 = 0.041
# p  = 0.102

# 3-D:
# r2 = 0.098
# p  = 0.014

# In the 3-D solution, the BCM vector is oriented predominantly
# along NMDS3, indicating that an important component of the
# BCM-associated compositional gradient is poorly represented
# in the 2-D projection.

# Extract direction of BCM vector in the 3-D ordination
bcm_vector_3d <- vegan::scores(
  bcm_fit_vir_3d,
  display = "vectors",
  choices = 1:3)
bcm_vector_3d

bcm_fit_vir_3d$vectors$arrows

# Current direction:
# NMDS1 = 0.291
# NMDS2 = 0.303
# NMDS3 = 0.907
# The strong NMDS3 component explains why the BCM relationship
# becomes clearer in the 3-D solution.

# Fit benthic components to ordinations
set.seed(123)
benthic_fit_vir_2d <- vegan::envfit(
  nmds_vir_2d,
  benthic_rel_vir,
  permutations = 999,
  choices = 1:2)

set.seed(123)
benthic_fit_vir_3d <- vegan::envfit(
  nmds_vir_3d,
  benthic_rel_vir,
  permutations = 999,
  choices = 1:3)

benthic_fit_vir_2d
benthic_fit_vir_3d

# Extract envfit results for reporting
envfit_2d <- data.frame(
  Variable = rownames(benthic_fit_vir_2d$vectors$arrows),
  R2 = benthic_fit_vir_2d$vectors$r,
  P = benthic_fit_vir_2d$vectors$pvals)

envfit_3d <- data.frame(
  Variable = rownames(benthic_fit_vir_3d$vectors$arrows),
  R2 = benthic_fit_vir_3d$vectors$r,
  P = benthic_fit_vir_3d$vectors$pvals)

envfit_2d
envfit_3d

# Extract site scores from the 3D solution
site_scores_3d <- as.data.frame(
  vegan::scores(
    nmds_vir_3d,
    display = "sites",
    choices = 1:3))

# Add BCM cover for visualisation
site_scores_3d$BCM_pct <- dat_vir_keep$BCM_pct

# Plot the three pairwise views of the 3D solution
library(ggplot2)
library(patchwork)
library(ggrepel)

bcm_colours <- c(
  "#C7CEC4",
  "#7FA074",
  "#526B49",
  "#2B4920")

community_axes <- as.data.frame(
  vegan::scores(
    benthic_fit_vir_3d,
    display = "vectors",
    choices = 1:3))
community_axes$Category <- rownames(community_axes)

# for manuscript - NMDS1 versus NMDS3 with community shown
vir_1_3 <- ggplot(
  site_scores_3d,
  aes(x = NMDS1, y = NMDS3)) +
  geom_point(
    aes(colour = BCM_pct),
    size = 3,
    alpha = 0.8) +
  geom_segment(
    data = community_axes,
    aes(x = 0, y = 0, xend = NMDS1, yend = NMDS3),
    inherit.aes = FALSE,
    arrow = arrow(length = grid::unit(0.2, "cm")),
    colour = "black") +
  geom_text_repel(
    data = community_axes,
    aes(x = NMDS1, y = NMDS3, label = Category),
    inherit.aes = FALSE,
    size = 4) +
  scale_colour_gradientn(
    colours = bcm_colours,
    name = "BCM cover (%)") +
  coord_equal() +
  theme_classic(base_size = 14) +
  labs(
    x = "NMDS1",
    y = "NMDS3",
    title = "NMDS1 × NMDS3",
    subtitle = sprintf(
      "3-D solution; stress = %.3f",
      nmds_vir_3d$stress))

# NMDS1 versus NMDS2
vir_1_2 <- ggplot(
  site_scores_3d,
  aes(x = NMDS1, y = NMDS2)) +
  geom_point(
    aes(colour = BCM_pct),
    size = 3,
    alpha = 0.8) +
  scale_colour_gradientn(
    colours = bcm_colours,
    name = "BCM cover (%)") +
  coord_equal() +
  theme_classic(base_size = 14) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    title = "NMDS1 × NMDS2")

# NMDS1 versus NMDS3
vir_1_3_supp <- ggplot(
  site_scores_3d,
  aes(x = NMDS1, y = NMDS3)) +
  geom_point(
    aes(colour = BCM_pct),
    size = 3,
    alpha = 0.8) +
  scale_colour_gradientn(
    colours = bcm_colours,
    name = "BCM cover (%)") +
  coord_equal() +
  theme_classic(base_size = 14) +
  labs(
    x = "NMDS1",
    y = "NMDS3",
    title = "NMDS1 × NMDS3")

# NMDS2 versus NMDS3
vir_2_3 <- ggplot(
  site_scores_3d,
  aes(x = NMDS2, y = NMDS3)) +
  geom_point(
    aes(colour = BCM_pct),
    size = 3,
    alpha = 0.8) +
  scale_colour_gradientn(
    colours = bcm_colours,
    name = "BCM cover (%)") +
  coord_equal() +
  theme_classic(base_size = 14) +
  labs(
    x = "NMDS2",
    y = "NMDS3",
    title = "NMDS2 × NMDS3")

# Combined supplementary figure
figure_nmds_3d_vir <- vir_1_2 + vir_1_3_supp + vir_2_3
figure_nmds_3d_vir

### --------------------------------------------------------------
# LAU: 2D and 3D NMDS of benthic community composition
# Prepare Lau data

dat_lau <- dat_nmds %>%
  filter(Seascape == "LAU_Seascape") %>%
  mutate(BCM_pct = as.numeric(as.character(Cyanobacteria)))

benthic_cols <- c(
  "HC", "TA", "MA", "CCA", "SC",
  "BS", "OI", "SD", "RB")

# Ensure all benthic variables are numeric
dat_lau[benthic_cols] <- lapply(
  dat_lau[benthic_cols],
  function(x) as.numeric(as.character(x)))

# Extract benthic community data
benthic_lau <- dat_lau[, benthic_cols, drop = FALSE]

# Standardise non-BCM benthic composition
# Calculate total cover represented by the non-BCM benthic
# categories at each site
row_tot <- rowSums(benthic_lau, na.rm = TRUE)

# Retain only sites with:
#   1. complete benthic composition data
#   2. a BCM-cover value
#   3. non-zero total non-BCM benthic cover
keep <- complete.cases(benthic_lau, dat_lau$BCM_pct) &
  row_tot > 0
benthic_lau_keep <- benthic_lau[keep, , drop = FALSE]
dat_lau_keep <- dat_lau[keep, , drop = FALSE]

# Convert benthic cover to relative composition after excluding BCM.
# Each row sums to 1.
benthic_rel_lau <- benthic_lau_keep / row_tot[keep]
nrow(benthic_rel_lau) # Number of sites included in the NMDS

# 2D NMDS
set.seed(123)
nmds_lau_2d <- vegan::metaMDS(
  benthic_rel_lau,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE)

# 3D NMDS
set.seed(123)
nmds_lau_3d <- vegan::metaMDS(
  benthic_rel_lau,
  distance = "bray",
  k = 3,
  trymax = 200,
  autotransform = FALSE)

# 5. COMPARE ORDINATION STRESS
stress_comparison_lau <- data.frame(
  Dimensions = c("2-D", "3-D"),
  Stress = c(
    nmds_lau_2d$stress,
    nmds_lau_3d$stress))

stress_comparison_lau

# Current result:
# 2-D stress = 0.120
# 3-D stress = 0.086
# Interpretation: Adding a third dimension improves the representation
# of the Bray-Curtis dissimilarities among Lau reef sites.

# Fit BCM cover to each ordination
set.seed(123)
bcm_fit_lau_2d <- vegan::envfit(
  nmds_lau_2d,
  dat_lau_keep$BCM_pct,
  permutations = 999,
  choices = 1:2)

set.seed(123)
bcm_fit_lau_3d <- vegan::envfit(
  nmds_lau_3d,
  dat_lau_keep$BCM_pct,
  permutations = 999,
  choices = 1:3)

bcm_fit_lau_2d
bcm_fit_lau_3d

# Current results:
# 2-D:
# r2 = 0.094
# p  = 0.197

# 3-D:
# r2 = 0.140
# p  = 0.196

# BCM association remains non-significant in both solutions.
# Retain the 2-D solution for manuscript reporting

# Fit benthic components to each ordination
set.seed(123)
benthic_fit_lau_2d <- vegan::envfit(
  nmds_lau_2d,
  benthic_rel_lau,
  permutations = 999,
  choices = 1:2)

set.seed(123)
benthic_fit_lau_3d <- vegan::envfit(
  nmds_lau_3d,
  benthic_rel_lau,
  permutations = 999,
  choices = 1:3)

benthic_fit_lau_2d
benthic_fit_lau_3d

# Extract envfit results for reporting
envfit_lau_2d <- data.frame(
  Variable = rownames(benthic_fit_lau_2d$vectors$arrows),
  R2 = benthic_fit_lau_2d$vectors$r,
  P = benthic_fit_lau_2d$vectors$pvals)

envfit_lau_3d <- data.frame(
  Variable = rownames(benthic_fit_lau_3d$vectors$arrows),
  R2 = benthic_fit_lau_3d$vectors$r,
  P = benthic_fit_lau_3d$vectors$pvals)

envfit_lau_2d
envfit_lau_3d

# Extract site scores from the 2D solution
site_scores_2d <- as.data.frame(
  vegan::scores(
    nmds_lau_2d,
    display = "sites",
    choices = 1:2))

# Add BCM cover for visualisation
site_scores_2d$BCM_pct <- dat_lau_keep$BCM_pct

# Plot the 2D solution for Lau
bcm_colours <- c(
  "#EAF7FC", 
  "#9EBCD7",
  "#759EBC", 
  "#2A5B7F")

# NMDS1 versus NMDS2
community_axes <- as.data.frame(
  scores(benthic_fit_lau_2d, display = "vectors"))
community_axes$Category <- rownames(community_axes)
community_axes

lau_1_2 <- ggplot(
  site_scores_2d,
  aes(x = NMDS1, y = NMDS2)) +
  geom_point(
    aes(colour = BCM_pct),
    size = 3,
    alpha = 0.8) +
  geom_segment(
    data = community_axes,
    aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
    inherit.aes = FALSE,
    arrow = arrow(length = grid::unit(0.2, "cm")),
    colour = "black") +
  geom_text_repel(
    data = community_axes,
    aes(x = NMDS1, y = NMDS2, label = Category),
    inherit.aes = FALSE,
    size = 4) +
  scale_colour_gradientn(
    colours = bcm_colours,
    name = "BCM cover (%)") +
  coord_equal() +
  theme_classic(base_size = 14) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    title = "NMDS1 × NMDS2",
    subtitle = sprintf(
      "2-D solution; stress = %.3f",
      nmds_lau_2d$stress))

lau_1_2

#need to plot lau_1_2 and vir_1_3

# Test correlations
library(purrr)
library(tibble)

benthic_vars <- c("HC", "CCA", "TA", "MA", "SC",
                  "BS", "OI", "SD", "RB")

run_spearman <- function(data, seascape) {
  map_dfr(benthic_vars, function(var) {
    test <- cor.test(
      data$BCM_pct,
      data[[var]],
      method = "spearman",
      exact = FALSE)
    tibble(
      Seascape = seascape,
      Category = var,
      rho = unname(test$estimate),
      p = test$p.value)})}

spearman_results <- bind_rows(
  run_spearman(dat_lau_keep, "Lau"),
  run_spearman(dat_vir_keep, "Vatu-i-Ra")) %>%
  group_by(Seascape) %>%
  mutate(p_adj = p.adjust(p, method = "holm")) %>%
  ungroup()

print(spearman_results, n = Inf)

# To be more ecologically meaningful, look at the relationship with unconsolidated
# substrate (sand + rubble) and reef-builders (CCA + hard coral)

run_combined_spearman <- function(dat_nmds, seascape) {
  dat_nmds <- dat_nmds %>%
    mutate(
      Unconsolidated = SD + RB,
      Reef_builders = HC + CCA)
  vars <- c("Unconsolidated", "Reef_builders")
  map_dfr(vars, function(var) {
    test <- cor.test(
      dat_nmds$BCM_pct,
      dat_nmds[[var]],
      method = "spearman",
      exact = FALSE)
    tibble(
      Seascape = seascape,
      Category = var,
      rho = unname(test$estimate),
      p = test$p.value)})}

combined_results <- bind_rows(
  run_combined_spearman(dat_lau_keep, "Lau"),
  run_combined_spearman(dat_vir_keep, "Vatu-i-Ra")) %>%
  group_by(Seascape) %>%
  mutate(p_adj = p.adjust(p, method = "holm")) %>%
  ungroup()

print(combined_results, n = Inf)

# Plot the functional categories

dat_lau_plot <- dat_lau_keep %>%
  mutate(
    Reef_builders = HC + CCA,
    Unconsolidated = SD + RB)

lau_RB <- ggplot(dat_lau_plot,
                 aes(Reef_builders, BCM_pct)) +
  geom_point(size = 2, alpha = 0.6, colour = "#9EBCD7") +
  geom_smooth(method = "gam",
              formula = y ~ s(x, k = 2),
              colour = "black") +
  theme_classic() +
  labs(x = "Reef-building taxa (%)\n(Hard coral + crustose coralline algae)",
    y = "BCM cover (%)")

lau_US <- ggplot(dat_lau_plot,
                 aes(Unconsolidated, BCM_pct)) +
  geom_point(size = 2, alpha = 0.6, colour = "#9EBCD7") +
  geom_smooth(method = "gam",
              formula = y ~ s(x, k = 2),
              colour = "black") +
  theme_classic() +
  labs(x = "Unconsolidated substrate (%)\n(Sand + Rubble)",
    y = "BCM cover (%)")

dat_vir_plot <- dat_vir_keep %>%
  mutate(
    Reef_builders = HC + CCA,
    Unconsolidated = SD + RB)

vir_RB <- ggplot(dat_vir_plot,
                 aes(Reef_builders, BCM_pct)) +
  geom_point(size = 2, alpha = 0.6, colour = "#7FA074") +
  geom_smooth(method = "gam",
              formula = y ~ s(x, k = 2),
              colour = "black") +
  theme_classic() +
  labs(x = "Reef-building taxa (%)\n(Hard coral + crustose coralline algae)",
    y = "BCM cover (%)")

vir_US <- ggplot(dat_vir_plot,
                 aes(Unconsolidated, BCM_pct)) +
  geom_point(size = 2, alpha = 0.6, colour = "#7FA074") +
  geom_smooth(method = "gam",
              formula = y ~ s(x, k = 2),
              colour = "black") +
  theme_classic() +
  labs(x = "Unconsolidated substrate (%)\n(Sand + Rubble)",
    y = "BCM cover (%)")

# Make final combined figure 3
library(cowplot)

lau_RB <- lau_RB +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = "rho = -0.69\np < 0.001",
    hjust = 1.1, vjust = 1.1,
    size = 4)

lau_US <- lau_US +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = "rho = 0.49\np = 0.002",
    hjust = 1.1, vjust = 1.1,
    size = 4)

vir_RB <- vir_RB +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = "rho = -0.13\np = 0.376",
    hjust = 1.1, vjust = 1.1,
    size = 4)

vir_US <- vir_US +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = "rho = -0.04\np = 0.679",
    hjust = 1.1, vjust = 1.1,
    size = 4)

top_row <- plot_grid(
  lau_1_2,
  plot_grid(
    lau_RB, lau_US,
    ncol = 1,
    labels = c("B", "C"),
    label_size = 14),
  ncol = 2,
  labels = c("A", ""),
  rel_widths = c(2, 1))

bottom_row <- plot_grid(
  vir_1_3,
  plot_grid(
    vir_RB, vir_US,
    ncol = 1,
    labels = c("E", "F"),
    label_size = 14),
  ncol = 2,
  labels = c("D", ""),
  rel_widths = c(2, 1))

figure_combined <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 1))

figure_combined

png("figure_combined.png", width = 2500, height = 2500, res = 300)
print(figure_combined)  # renders the ggplot to the PNG device
dev.off()