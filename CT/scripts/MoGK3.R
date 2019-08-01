# Load packages
library(tidyverse)
library(rstan)
library(plyr)
library(loo)
library(magrittr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("CopyTask2019/")

# Sampling parameters
set.seed(125)
n_cores = 3
n_chain = 3 
iterations = 10000

# Load df
d <- read_csv("data/ct.csv") %>% 
  dplyr::filter(IKI > 0 
         ,target == 1
#         ,subj %in% c(1:10)
         ) %>%
  group_by( subj, bigram, component ) %>%
  dplyr::summarise(IKI = mean(IKI),
                   N = n()) %>%
  ungroup() %>%
  mutate(component = factor(component, levels = c("Tapping", "Sentence", "HF", "LF", "Consonants"), ordered = TRUE)); d 

  # Data as list for stan input
components <- d %$% as.numeric(mapvalues(component, from = unique(component), to= 1:length(unique(component))))
dat <-  within( list(), {
    N <- nrow(d)
    y <- d$IKI # DV
    components <- components
    K <- max(components)
    subj <- d$subj
    S <- length(unique(d$subj))
    bigram <- as.integer(factor(d$bigram))
    B <- length(unique(d$bigram))
  } );str(dat)


# Initialise start values
start <- 
  function(chain_id = 1){
    list(beta = rep(mean(log(dat$y)),dat$K)
         , delta = rep(log(50), dat$K)
         , sigma = rep(sd(log(dat$y)), dat$K)
         , sigma_diff = rep(.01, dat$K)
         , sigma_diff2 = rep(.01, dat$K)
         , theta = matrix(rep(1/3, 3*dat$K), nrow = dat$K)
         , z_u = matrix(rep(0.1,dat$S*dat$K), nrow = dat$K)
         , w = rep(0.1, length(unique(d$bigram)))
         , L_u = matrix(rep(0, dat$K*dat$K), nrow = dat$K)
         , sigma_u = rep(1, dat$K)
         , sigma_w = 1
         , alpha = chain_id
    )
  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# --------------
# Stan models ##
# --------------
#---- 
# Mixture of three gaussians with unequal variance (constraint)
#---- 
# Load model
mog <- stan_model(file = "stanin/MoGK3.stan")

# Check model
#m <- sampling(mog, chain = 1, iter = 1, data = dat) 

# Parameters to omit in output
omit <- c("mu1", "mu2", "mu3", "L_u", "z_u", "sigma_e_comp", "sigmap_e_comp", "sigmap_ee_comp", "log_theta", "log_theta_comp")

# Fit model
m <- sampling(mog, 
              data = dat,
              init = start_ll,
              iter = iterations,
              warmup = iterations/2,
              chains = n_chain, 
              cores = n_cores, 
	      refresh = 250,
              save_warmup = FALSE, # Don't save the warmup
              include = FALSE, # Don't include the following parameters in the output
              pars = omit,
              thin = 1,
              seed = 81,
              control = list(max_treedepth = 16,
                             adapt_delta = 0.99)
)



# Save model
saveRDS(m, 
        file = "stanout/MoGK3.rda",
        compress = "xz")


# Extract and save posterior and log likelihood seperately
# Get log likelihood
log_lik <- extract_log_lik(m) 
saveRDS(log_lik, 
        file = "stanout/log_lik/MoGK3_loglik.rda",
        compress = "xz")


# Get parameter posterior
param <- c("beta1", "beta", "beta3", "delta", "gamma", "theta", "sigma", "sigma_diff1", "sigma_diff2", "sigmap_e", "sigmap_ee", "sigma_e")
samps <- as.data.frame(m, pars = param) 
saveRDS(samps, 
        file = "stanout/posterior/MoGK3_posterior.rda",
        compress = "xz")


# Get random effects
param <- c("u", "w", "sigma_u", "sigma_w")
re <- as.data.frame(m, pars = param) 
saveRDS(re, 
        file = "stanout/RE/MoGK3_re.rda",
        compress = "xz")


# Get posterior predicted values
param <- c("y_tilde")
y_tilde <- as.data.frame(m, pars = param) 
saveRDS(y_tilde, 
        file = "stanout/y_tilde/MoGK3_y_tilde.rda",
        compress = "xz")

