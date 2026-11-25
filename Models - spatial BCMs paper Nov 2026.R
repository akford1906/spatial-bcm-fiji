# Script for analysing spatial BCM data in Fiji

library(dplyr)
library(ggplot2)
library(gtsummary)

#import mastersheet as dat
dat <- read.csv(file.choose()) 

dat <- dat %>%
  mutate(Seascape = case_when(
      Province %in% c("Lomaiviti", "Bua", "Ra", "Tailevu") ~ "VIR_Seascape",
      Province == "Lau" ~ "LAU_Seascape",
      TRUE ~ NA_character_),
    BS  = `Bare.Substrate`, # Rename benthic variables with spaces
    OI  = `Other.Inverts`,
    CCA = Crustose) %>%
  filter(!is.na(Seascape))
head(dat)

#import wave data and merge for sites
wave_data <- read.csv(file.choose()) 

wave_data_clean <- wave_data %>%
  distinct(Year, Sample_month, Site, .keep_all = TRUE)

dat_merged <- dat %>%
  left_join(wave_data_clean, by = c("Year", "Sample_month", "Site"))

#include new column for seascape (two - VatuiRa and Lau)
dat2 <- dat_merged %>%
  mutate(Seascape = case_when(
      Province %in% c("Lomaiviti", "Bua", "Ra", "Tailevu") ~ "VIR_Seascape",
      Province == "Lau" ~ "LAU_Seascape",
      TRUE ~ NA_character_))

# Ensure wave exposure is numeric
dat2$Mean_HS <- as.numeric(as.character(dat2$Mean_HS))

# Convert NaN values to NA
dat2$Mean_HS[is.nan(dat2$Mean_HS)] <- NA_real_

# Check conversion
str(dat2$Mean_HS)
summary(dat2$Mean_HS)

dat2 %>%
  group_by(Seascape) %>%
  summarise(
    N_sites = n(),
    Missing_wave = sum(is.na(Mean_HS)),
    Missing_herbivore = sum(is.na(Herbivore_biomass)),
    Missing_both = sum(is.na(Mean_HS) & is.na(Herbivore_biomass)),
    Complete_both = sum(complete.cases(Mean_HS, Herbivore_biomass)))

dat2 %>%
  select(
    Seascape,
    Depth, Mean_HS, Market_gravity, Total_biomass, Herbivore_biomass,
    Distance_river, Distance_town, Distance_village,
    Cyanobacteria, HC, MA, TA, RB, SD, SC,
    BS, OI, CCA,
    Island_Type, Season, `Reef.Type`, `Reef.Zone`) %>%
  tbl_summary(
    by = Seascape,
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}"))

#clear differences in predictors across seascapes - Table S1


######################################

#Differences by seascape

dat2 %>%
  group_by(Seascape) %>%
  summarise(
    N_sites = sum(!is.na(Cyanobacteria)),
    Mean_BCM = mean(Cyanobacteria, na.rm = TRUE),
    SE_BCM = sd(Cyanobacteria, na.rm = TRUE) / sqrt(N_sites))

dat2 %>%
  filter(Cyanobacteria > 0) %>%
  summarise(
    N_sites = n(),
    Mean_BCM = mean(Cyanobacteria),
    SE_BCM = sd(Cyanobacteria) / sqrt(N_sites))

dat2 %>%
  summarise(
    N_sites = sum(!is.na(Cyanobacteria)),
    N_above10 = sum(Cyanobacteria > 10, na.rm = TRUE),
    Pct_above10 = 100 * N_above10 / N_sites)

dat2 %>%
  filter(Cyanobacteria > 0) %>%
  summarise(
    N_present = n(),
    N_above10 = sum(Cyanobacteria >= 10),
    N_above20 = sum(Cyanobacteria >= 20),
    N_above30 = sum(Cyanobacteria >= 30),
    Pct_above10 = 100 * N_above10 / N_present,
    Pct_above20 = 100 * N_above20 / N_present,
    Pct_above30 = 100 * N_above30 / N_present)

dat2 %>%
  group_by(Seascape) %>%
  summarise(
    N_sites = sum(!is.na(Cyanobacteria)),
    N_above10 = sum(Cyanobacteria >= 10, na.rm = TRUE),
    Pct_above10 = 100 * N_above10 / N_sites)

dat2 %>%
  filter(Cyanobacteria > 0) %>%
  group_by(Seascape) %>%
  summarise(
    N_present = n(),
    N_above10 = sum(Cyanobacteria >= 10),
    N_above20 = sum(Cyanobacteria >= 20),
    N_above30 = sum(Cyanobacteria >= 30),
    Pct_above10 = 100 * N_above10 / N_present,
    Pct_above20 = 100 * N_above20 / N_present,
    Pct_above30 = 100 * N_above30 / N_present)

#Figure for differences across seascapes (Figure 2e)

seascape_colors <- c(
  "LAU_Seascape" = "#9EBCD7", 
  "VIR_Seascape" = "#7FA074")

figure2e <- ggplot(dat2, aes(x = Seascape, y = Cyanobacteria, fill = Seascape)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(
    aes(fill = Seascape),
    shape = 21,
    colour = "black",
    width = 0.1,
    alpha = 0.5,
    size = 2,
    show.legend = FALSE) +
  scale_x_discrete(
    limits = c("VIR_Seascape", "LAU_Seascape")) +
  scale_fill_manual(values = seascape_colors) +
  labs(x = "Seascape", y = "BCM Coverage (%)") +
  theme_classic() +
  theme(legend.position = "none")
figure2e

#check statistical differences across seascapes
library(glmmTMB)
library(DHARMa)
library(geosphere)

# Retain complete observations for the seascape models
dat_seascape <- dat2 %>%
  filter(
    complete.cases(
      Cyanobacteria,
      Seascape,
      longitude,
      latitude))

m_seascape <- glmmTMB(
  Cyanobacteria ~ Seascape,
  family = tweedie(link = "log"),
  data = dat_seascape)
summary(m_seascape)

sim <- simulateResiduals(m_seascape, n = 1000)
plot(sim)

# Check spatial correlation
# Geographic distances in metres; longitude comes first
distance_matrix <- geosphere::distm(
  as.matrix(dat_seascape[, c("longitude", "latitude")]),
  fun = geosphere::distGeo)

library(ncf)

set.seed(123)
correlog_res <- correlog(
  x = dat_seascape$longitude,
  y = dat_seascape$latitude,
  z = sim$scaledResiduals,
  increment = 10,       # 10-km distance classes
  resamp = 999,
  latlon = TRUE,
  quiet = TRUE)
plot(correlog_res)
print(correlog_res) #clear spatial correlation up to 40km

# incorporate random effect for spatial autocorrelation
# create random effect so that nearby sites to share similar unexplained 
# conditions, with that similarity declining as distance increases
library(sf)

# Project coordinates to UTM Zone 60S and express in kilometres
points <- sf::st_as_sf(
  dat_seascape,
  coords = c("longitude", "latitude"),
  crs = 4326)

coords_m <- sf::st_coordinates(
  sf::st_transform(points, crs = 32760))

coords_km <- coords_m / 1000 #change to km

dat_seascape$pos <- glmmTMB::numFactor(
  coords_km[, 1],
  coords_km[, 2])

dat_seascape$spatial_group <- factor(rep(1, nrow(dat_seascape)))

m_spatial <- glmmTMB(
  Cyanobacteria ~ Seascape + exp(pos + 0 | spatial_group),
  family = tweedie(link = "log"),
  data = dat_seascape)

# Compare spatial and non-spatial seascape models

# Model comparison
AIC(m_seascape, m_spatial)

# Model sample size

nobs(m_spatial)

# Convergence
m_spatial$fit$convergence
m_spatial$sdr$pdHess

# Check model assumptions
sim <- simulateResiduals(m_spatial, n = 1000, rotation = "estimated")
plot(sim) 

# Account for observations at shared coordinates
coord_id <- interaction(
  dat_seascape$longitude,
  dat_seascape$latitude,
  drop = TRUE)

sim_space <- recalculateResiduals(
  sim,
  group = coord_id,
  rotation = "estimated")

coord_unique <- aggregate(
  dat_seascape[, c("longitude", "latitude")],
  by = list(coord_id = coord_id),
  FUN = mean)

dist_seascape <- geosphere::distm(
  as.matrix(coord_unique[, c("longitude", "latitude")]),
  fun = geosphere::distGeo)

testSpatialAutocorrelation(
  sim_space,
  distMat = dist_seascape)

sim_seascape_final <- DHARMa::simulateResiduals(
  fittedModel = m_spatial,
  n = 1000,
  seed = 123,
  rotation = "estimated")

# Confirm model sample sizes
nobs(m_seascape)
nobs(m_spatial)

# Overall seascape spatial model
summary(m_spatial)$coefficients$cond

library(emmeans)
emmeans(m_spatial, ~ Seascape, type = "response")

# Statistical contrast
pairs(emmeans(m_spatial, ~ Seascape),
  type = "response")

# Test hypothesis-driven Predictors

par(mfrow = c(1,1))

dat_lau <- dat2 %>% 
  filter(Seascape == "LAU_Seascape") %>%
  mutate(BCM_pct = as.numeric(Cyanobacteria))

dat_ViR <- dat2 %>% 
  filter(Seascape == "VIR_Seascape") %>%
  mutate(BCM_pct = as.numeric(Cyanobacteria))

hist(dat_lau$Cyanobacteria)
shapiro.test(dat_lau$Cyanobacteria)
summary(dat_lau$Cyanobacteria)

hist(dat_ViR$Cyanobacteria)
shapiro.test(dat_ViR$Cyanobacteria)
summary(dat_ViR$Cyanobacteria)

# check for spatial autoclustering
library(ape)
library(ade4)

#check for lau
coords <- dat_lau[, c("longitude", "latitude")]
geo_dist <- as.matrix(dist(coords))
bcm_dist <- dist(dat_lau$BCM_pct)
mantel.rtest(as.dist(geo_dist),
             bcm_dist,
             nrepet = 999) #no issues

#check for vir
coords <- dat_ViR[, c("longitude", "latitude")]
geo_dist <- as.matrix(dist(coords))
bcm_dist <- dist(dat_ViR$BCM_pct)
mantel.rtest(as.dist(geo_dist),
             bcm_dist,
             nrepet = 999) #no issues

# Model for Lau
dat_lau$Reef_builders <- dat_lau$HC + dat_lau$CCA
dat_lau$distance_villagekm <- dat_lau$Distance_village/1000

# Retain complete observations for model fitting
vars <- c("Cyanobacteria", "Mean_HS", "Reef_builders",
  "distance_villagekm", "Herbivore_biomass",
  "longitude", "latitude")

d_lau <- dat_lau[complete.cases(dat_lau[, vars]),]

m <- lm(Cyanobacteria ~ Mean_HS +
    Reef_builders +
    distance_villagekm +
    Herbivore_biomass,
  data = d_lau)

summary(m)

# Lau model explained variance
summary(m)$r.squared
summary(m)$adj.r.squared

nobs(m)

par(mfrow = c(1,1))

shapiro.test(residuals(m))
library(lmtest)
bptest(m)
plot(cooks.distance(m))
abline(h = 4/nobs(m), col = "red")
library(performance)
library(see)
performance::check_collinearity(m)
check_model(m)

m1 <- lm(Cyanobacteria ~ Reef_builders,
  data = d_lau)

m2 <- lm(Cyanobacteria ~ Reef_builders + Mean_HS,
  data = d_lau)

m3 <- lm(Cyanobacteria ~ Reef_builders + Mean_HS + Distance_village,
  data = d_lau)

AIC(m1, m2, m3) #m3 is best

# Check spatial autocorrelation

set.seed(123)

cor_lau <- ncf::correlog(
  x = d_lau$longitude,
  y = d_lau$latitude,
  z = residuals(m),
  increment = 10,
  resamp = 999,
  latlon = TRUE,
  quiet = TRUE)

plot(cor_lau)

data.frame(
  distance_km = cor_lau$mean.of.class,
  Moran_I = cor_lau$correlation,
  p = cor_lau$p,
  p_Holm = p.adjust(cor_lau$p, method = "holm"))

p.adjust(cor_lau$p, method = "holm")

library(ggeffects)

pred_hs <- ggpredict(m, terms = "Mean_HS")
pred_rb <- ggpredict(m, terms = "Reef_builders")
pred_vil <- ggpredict(m, terms = "distance_villagekm")

p_hs <- ggplot(pred_hs,
               aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = "#9EBCD7",
    alpha = 0.2) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = "β = -11.16\np = 0.028",
    hjust = 1.1,
    vjust = 1.3,
    size = 4) +
  geom_line(
    colour = "#005f80",
    linewidth = 1 ) +
  theme_classic() +
  labs(
    x = "Mean wave exposure (Mean_HS)",
    y = "Predicted BCM cover (%)")

p_rb <- ggplot(pred_rb,
               aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = "#9EBCD7",
    alpha = 0.2  ) + 
  annotate(
    "text",
    x = Inf, y = Inf,
    label = "β = -0.302\np < 0.001",
    hjust = 1.1,
    vjust = 1.3,
    size = 4) +
  geom_line(
    colour = "#005f80",
    linewidth = 1) +
  theme_classic() +
  labs(
    x = "Reef-building taxa (%)",
    y = "Predicted BCM cover (%)")

p_vil <- ggplot(pred_vil,
                aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = "#9EBCD7",
    alpha = 0.2) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = "β = -0.3361\np = 0.012",
    hjust = 1.1,
    vjust = 1.3,
    size = 4) +
  geom_line(
    colour = "#005f80",
    linewidth = 1) +
  theme_classic() +
  labs(
    x = "Distance to village (km)",
    y = "Predicted BCM cover (%)")

library(cowplot)

plot_grid(p_hs, p_rb, p_vil,ncol = 3, labels = c("A","B","C"))

# for VIR
dat_ViR$Reef_builders <- dat_ViR$HC + dat_ViR$CCA
dat_ViR$distance_villagekm <- dat_ViR$Distance_village/1000

m_VIR <- glm(
  Cyanobacteria ~ Mean_HS + 
    Reef_builders + 
    distance_villagekm + 
    Herbivore_biomass,
  data = dat_ViR)
summary(m_VIR)
par(mfrow = c(2,2))
plot(m_VIR)
performance::check_collinearity(m_VIR)
check_model(m_VIR) #not good

# Retain complete observations for model fitting
vars <- c("BCM_pct", "Reef_builders", "Mean_HS",
  "distance_villagekm", "Herbivore_biomass",
  "longitude", "latitude")

datViR_model <- dat_ViR[
  complete.cases(dat_ViR[, vars]),]

nrow(datViR_model)

# Try Tweedie GLM (for right-skewed and zero data)
m_vir_tweedie <- glmmTMB(
  BCM_pct ~ Reef_builders +
    Mean_HS +
    distance_villagekm +
    Herbivore_biomass,
  family = tweedie(link = "log"),
  data = datViR_model)

summary(m_vir_tweedie)

#check residual spatial autocorrelation
nobs(m_vir_tweedie)

sim <- simulateResiduals(
  fittedModel = m_vir_tweedie,
  n = 1000)

plot(sim)
length(sim$scaledResiduals)

# aggregate residuals at shared horizontal locations
datViR_model$coord_id <- interaction(
  datViR_model$longitude,
  datViR_model$latitude,
  drop = TRUE)

sim_space <- recalculateResiduals(
  sim,
  group = datViR_model$coord_id)

coord_unique <- aggregate(
  datViR_model[, c("longitude", "latitude")],
  by = list(coord_id = datViR_model$coord_id),
  FUN = mean)

# GLOBAL MORAN'S I TEST using geographic distances
dist_vir <- geosphere::distm(
  as.matrix(coord_unique[, c("longitude", "latitude")]),
  fun = geosphere::distGeo)

testSpatialAutocorrelation(
  sim_space,
  distMat = dist_vir)

# CORRELOGRAM using the same aggregated DHARMa residuals
set.seed(123)

corr_vir <- ncf::correlog(
  x = coord_unique$longitude,
  y = coord_unique$latitude,
  z = sim_space$scaledResiduals,
  increment = 10,
  resamp = 999,
  latlon = TRUE,
  quiet = TRUE)

plot(corr_vir)
print(corr_vir)

vir_spatial_table <- data.frame(
  distance_km = corr_vir$mean.of.class,
  Moran_I = corr_vir$correlation,
  p = corr_vir$p,
  p_Holm = p.adjust(corr_vir$p, method = "holm"))

vir_spatial_table

# create random effect so that nearby sites to share similar unexplained 
# conditions, with that similarity declining as distance increases
library(sf)
points <- sf::st_as_sf(
  datViR_model,
  coords = c("longitude", "latitude"),
  crs = 4326)

coords_m <- sf::st_coordinates(
  sf::st_transform(points, crs = 32760))

coords_km <- coords_m / 1000

stopifnot(nrow(coords_km) == nrow(datViR_model))

datViR_model$pos <- glmmTMB::numFactor(
  coords_km[, 1],
  coords_km[, 2])

datViR_model$spatial_group <- factor(
  rep(1, nrow(datViR_model)))

# Replace the cluster intercept with exponential spatial correlation
m_vir_exp <- glmmTMB(
  BCM_pct ~ Reef_builders + Mean_HS +
    distance_villagekm + Herbivore_biomass +
    exp(pos + 0 | spatial_group), #add distance-dependent spatial correlation
  family = tweedie(link = "log"),
  data = datViR_model)

summary(m_vir_exp)
summary(m_vir_exp)$coefficients$cond

# Compare spatial and non-spatial Vatu-i-Ra models
AIC(m_vir_tweedie, m_vir_exp)

# Confirm consistent model sample sizes
nobs(m_vir_tweedie)
nobs(m_vir_exp)

# Check convergence
m_vir_exp$fit$convergence   
m_vir_exp$sdr$pdHess       

# Final Vatu-i-Ra spatial model diagnostics
sim_vir_exp <- DHARMa::simulateResiduals(
  fittedModel = m_vir_exp,
  n = 1000,
  seed = 123,
  rotation = "estimated")

plot(sim_vir_exp)

# Aggregate residuals at shared locations
sim_final_space <- DHARMa::recalculateResiduals(
  sim_vir_exp,
  group = datViR_model$coord_id,
  rotation = "estimated") #included because of spatial covariate

# Confirm dimensions match
stopifnot(length(sim_final_space$scaledResiduals) == nrow(dist_vir))

# Test residual spatial autocorrelation
DHARMa::testSpatialAutocorrelation(
  sim_final_space,
  distMat = dist_vir)

# Other model diagnostics
DHARMa::testUniformity(sim_vir_exp)
DHARMa::testDispersion(sim_vir_exp)
DHARMa::testZeroInflation(sim_vir_exp)

library(car)

vif_mod <- lm(
  BCM_pct ~ Mean_HS +
    distance_villagekm +
    Reef_builders +
    Herbivore_biomass,
  data = datViR_model)

vif(vif_mod)

pred_hs <- ggpredict(m_vir_exp, terms = "Mean_HS")
pred_rb <- ggpredict(m_vir_exp, terms = "Reef_builders")
pred_vil <- ggpredict(m_vir_exp, terms = "distance_villagekm")

p_hs_VIR <- ggplot(pred_hs,
                   aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = "#7FA074",
    alpha = 0.2) + 
  annotate(
    "text",
    x = Inf, y = Inf,
    label = "β = 1.426\np = 0.358",
    hjust = 1.1,
    vjust = 1.3,
    size = 4) +
  geom_line(
    colour = "#456B3B",
    linewidth = 1) +
  theme_classic() +
  labs(
    x = "Mean wave exposure (Mean_HS)",
    y = "Predicted BCM cover (%)")

p_rb_VIR <- ggplot(pred_rb,
                   aes(x = x, y = predicted)) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = "#7FA074",
    alpha = 0.2) + 
  annotate(
    "text",
    x = Inf, y = Inf,
    label = "β = -0.035\np < 0.001",
    hjust = 1.1,
    vjust = 1.3,
    size = 4) +
  geom_line(
    colour = "#456B3B",
    linewidth = 1) +
  theme_classic() +
  labs(
    x = "Reef-building taxa (%)",
    y = "Predicted BCM cover (%)")

p_vil_VIR <- ggplot(pred_vil,
                    aes(x = x, y = predicted)) +
  geom_line(linewidth = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high),
    fill = "#7FA074",
    alpha = 0.2) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = "β = 0.098\np = 0.002",
    hjust = 1.1,
    vjust = 1.3,
    size = 4) +
  geom_line(
    colour = "#456B3B",
    linewidth = 1) +
  theme_classic() +
  labs(
    x = "Distance to village (km)",
    y = "Predicted BCM cover (%)")

library(cowplot)

plot_grid(
  p_hs_VIR,
  p_rb_VIR,
  p_vil_VIR,
  ncol = 3,
  labels = c("A","B","C"))

# Final combined figure 4

plot_grid(
  p_hs,
  p_rb,
  p_vil,
  p_hs_VIR,
  p_rb_VIR,
  p_vil_VIR,
  ncol = 3,
  labels = c("A","B","C", "D", "E", "F"))

# Vatu-i-Ra marginal and conditional R2

formula(m_vir_exp)
nobs(m_vir_exp)
AIC(m_vir_exp)

r2_vir <- performance::r2_nakagawa(m_vir_exp)
print(r2_vir)
