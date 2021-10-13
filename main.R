# Code to simulate the antibody and vaccine efficacy component of the Imperial Malaria Transmission Model.

# Author: Alexandra Hogan
# Date: 13 October 2021
# Note that this version has the updated vaccine AB model parameters corresponding to White et al 2016 (https://doi.org/10.1016/S1473-3099(15)00239-X)

library(tidyverse)

source("functions.R")

tvec <- seq(0, 365*5, length.out=365*5) # time vector
reps <- 1000 # number of runs

# Fitted parameter values: from White et al (2015) Lancet ID
ab_mu <- 621  # median of the geometric means of observed antibody titres
ab_mu_boost <- 277 # median of the boosted antibody titre
t_boost <- 548 #time of booster dose in days after third vaccine dose
ab_sigma <- 0.35 # observational variance of antibody titre (log-normal)
ab_sigma_boost <- 0.35
d1_mu <- 45 # half-life of the short-lived component of antibody response 
d1_sigma <- 16 # standard deviation of HL of short-lived response
d2_mu <- 591 # half-life of long-lived component of antibody response
d2_sigma <- 245 # standard deviation of HL of long-lived response
rho_mu <- 2.378 # proportion of short-lived component following primary schedule
rho_sigma <- 1.008 
rho_mu_boost <- 1.034 # proportion of short-lived component following booster
rho_sigma_boost <- 1.027
alpha <- 0.74 # Shape parameter of dose-response curve
beta <- 99.2 # Scale parameter of dose-response curve
Vmax <- 0.93    # maximum efficacy against infection

# Initialise dataframes to store simulations of antibody titre and vaccine efficacy over time.
df_titre <- NULL
df_efficacy <- NULL

for (i in 1:reps) {
  titre <- antibody_titre(i = i,
                          tvec = tvec,
                          d1_mu = d1_mu,
                          d1_sigma = d1_sigma,
                          d2_mu = d2_mu,
                          d2_sigma = d2_sigma,
                          rho_mu = rho_mu,
                          rho_sigma = rho_sigma,
                          ab_mu = ab_mu,
                          ab_sigma= ab_sigma,
                          rho_mu_boost = rho_mu_boost,
                          rho_sigma_boost = rho_sigma_boost,
                          t_boost = t_boost,
                          ab_mu_boost = ab_mu_boost,
                          ab_sigma_boost = ab_sigma_boost)
  
  df_titre <- rbind(df_titre, titre)
  
  eff <- efficacy_profile(i = i,
                          alpha = alpha,
                          beta = beta,
                          Vmax = Vmax,
                          ab = titre$ab)
  
  df_efficacy <- rbind(df_efficacy, eff)
}

# summarise efficacies over runs
efficacy_median <- df_efficacy %>%
  group_by(time) %>%
  summarise(efficacy_med = median(efficacy),
            efficacy_upper = quantile(efficacy, 0.975),
            efficacy_lower = quantile(efficacy, 0.025))

write_csv(efficacy_median, "modelled_efficacy.csv")

# plot outputs
ggplot(data = efficacy_median, aes(x = time, y = efficacy_med*100)) +
  geom_line() +
  geom_ribbon(aes(ymin = efficacy_lower*100, ymax = efficacy_upper*100), alpha = 0.3, fill = "seagreen") +
  labs(x = "time (days)", y = "efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 365*5, by = 365)) +
  lims(y = c(0, 100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line())
  