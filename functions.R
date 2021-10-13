
# Functions from Bob Verity:

# draw from logitnormal distribution. Input arguments include the mean and standard deviation of the RAW normal variate that will be converted to the logit scale. Therefore the mean and standard deviation of the final random variable will not equal these values.
rlogitnorm <- function(n, mean_raw, sd_raw) {
  x <- exp(rnorm(n, mean_raw, sd_raw))
  x/(1+x)
}

# draw from lognormal distribution. Input arguments include the desired mean and standard deviation of the final lognormal random variable; the mean and standard deviation of the raw normal variate are calculated from these values.
rlnorm2 <- function(n, mean, sd) {
  meanlog <- log(mean^2/sqrt(mean^2+sd^2)) # mean of the lognormal distribution
  sdlog <- sqrt(log(1+(sd/mean)^2)) # standard deviation of the lognormal distribution
  rlnorm(n,meanlog,sdlog) # draw from the lognormal distribution
}

# draw single antibody titre from the model
antibody_titre <- function(i = 1,
                           tvec = 1:1825,
                           d1_mu = 45,
                           d1_sigma = 16,
                           d2_mu = 591,
                           d2_sigma = 245,
                           rho_mu = 2.37832,
                           rho_sigma = 1.00813,
                           ab_mu = 621,
                           ab_sigma = 0.35,
                           rho_mu_boost = 1.034,
                           rho_sigma_boost = 1.027,
                           t_boost = 548, 
                           ab_mu_boost = 277,
                           ab_sigma_boost = 0.35)
{
  # draw parameters from distributions
  rho <- rlogitnorm(1, mean_raw=rho_mu, sd_raw=rho_sigma)
  rho_boost <- rlogitnorm(1, mean_raw=rho_mu_boost, sd_raw=rho_sigma_boost)
  d1 <- rlnorm2(1, mean=d1_mu, sd=d1_sigma) # half-life of short-lived component of antibody response
  d2 <- rlnorm2(1, mean=d2_mu, sd=d2_sigma) # half-life of long-lived component of antibody response
  ab0 <- exp(rnorm(1, log(ab_mu)-ab_sigma^2/2, sd=ab_sigma))
  ab0_boost <- exp(rnorm(1, log(ab_mu_boost)-ab_sigma_boost^2/2, sd=ab_sigma_boost))
  
  # Simulate antibody titre over time
  r1 <- log(2)/d1
  r2 <- log(2)/d2
  tmax <- length(tvec)
  
  ab <- vector(mode = "numeric", length = tmax)
  tvec1 <- 1:t_boost
  ab[1:t_boost] <- ab0*(rho*exp(-r1*tvec1)+(1-rho)*exp(-r2*tvec1))
  tvec2 <- seq((t_boost + 1), tmax, by = 1)
  ab[(t_boost+1):tmax] <- ab0_boost*(rho_boost*exp(-r1*(tvec2-t_boost))+(1-rho_boost)*exp(-r2*(tvec2-t_boost)))
  
  ret <- data.frame(ab = ab, rep = i, time = tvec)
  return(ret)
}

# calculate efficacy profile from the antibody titre

efficacy_profile <- function(i=1, alpha=0.74, beta=99.2, Vmax=0.93, ab=1) {
  
  # Convert antibody titre ab to vaccine efficacy (Hill function)
  efficacy<-Vmax*(1 - 1/(1+(ab/beta)^alpha))
  
  ret <- data.frame(efficacy = efficacy, rep = i, time = 1:length(efficacy))
  return(ret)
}
