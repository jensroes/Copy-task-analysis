#Load packages
library(tidyverse)
library(rstan)
library(loo)
library(plyr)
library(magrittr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("CopyTask2019/")

# Sampling parameters
set.seed(125)
n_chain = 3 # number of cores/chains
n_cores = 3
iterations = 6000

d <- read_csv("data/ct.csv") %>% 
  dplyr::filter(IKI > 0 
                ,!is.na(IKI)
                ,target == 1
  ) %>%
  mutate(subj = as.numeric(factor(subj))) %>%
  group_by( subj, bigram, component ) %>%
  dplyr::summarise(IKI = mean(IKI),
                   N = n()) %>%
  ungroup() %>%
  mutate(component = factor(component, levels = c("Tapping", "Sentence", "HF", "LF", "Consonants"), ordered = TRUE)); d 

set.seed(123)
d %<>% filter(subj %in% sample(1:max(subj), 500)) %>%
  mutate(subj = as.numeric(factor(subj))) 

dat <-  within( list(), 
                {
                  N <- nrow(d)
                  y <- d$IKI # DV
                  subj <- d$subj
                  S <- max(d$subj)
                  bigram <- as.integer(factor(d$bigram))
                  B <- length(unique(d$bigram))
                } );str(dat)

# Initialise start values
start <- 
  function(chain_id = 1){
    list(beta = mean(log(d$IKI)), sigma = sd(log(d$IKI)),
         u = rep(0, max(d$subj)),
         w = rep(0, length(unique(d$bigram))),
         sigma_u = 1,
         sigma_w = 1,
         alpha = chain_id
    )
  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# --------------
# Stan models ##
# --------------
#---- 

#---- 
# Load model
lmm <- stan_model(file = "stanin/LMMsubjintercepts.stan")

# Check model
#m <- sampling(lmm, chain = 1, iter = 1, data = dat) 

# Fit model
m <- sampling(lmm, 
              data = dat,
              init = start_ll,
              iter = iterations,
              warmup = iterations/2,
              chains = n_chain, 
	      cores = n_cores,
	      refresh = 250,
              thin = 1,
              seed = 81,
              save_warmup = FALSE,
              include = FALSE,
              pars = c("mu"),
              control = list(max_treedepth = 16,
                             adapt_delta = 0.99)
)

# Save model
saveRDS(m, 
        file = "stanout/LMMsubjintercepts.rda",
        compress = "xz")


# Extract and save posterior and log likelihood seperately
# Get log likelihood
log_lik <- extract_log_lik(m) 
#loo(log_lik)
saveRDS(log_lik, 
        file = "stanout/log_lik/LMMsubjintercepts_loglik.rda",
        compress = "xz")

# Get parameter posterior
param <- c("beta", "sigma")
samps <- as.data.frame(m, pars = param) 
saveRDS(samps, 
        file = "stanout/posterior/LMMsubjintercepts_posterior.rda",
        compress = "xz")


# Get random effects
param <- c("u", "w", "sigma_u", "sigma_w")
re <- as.data.frame(m, pars = param) 
saveRDS(re, 
        file = "stanout/RE/LMMsubjintercepts_re.rda",
        compress = "xz")


# Get posterior predicted values
param <- c("y_tilde")
y_tilde <- as.data.frame(m, pars = param) 
saveRDS(y_tilde, 
        file = "stanout/y_tilde/LMMsubjintercepts_y_tilde.rda",
        compress = "xz")

