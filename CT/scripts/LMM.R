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
n_chain = 3
n_cores = 3
iterations = 10000

# Load df
d <- read_csv("data/ct.csv") %>% 
  filter(IKI > 0, 
         target == 1
) %>%
  group_by( subj, bigram, component ) %>%
  dplyr::summarise(IKI = mean(IKI),
                   N = n()) %>%
  ungroup(); d 


d$component <- d %$% factor(component, levels = c("Tapping", "Sentence", "HF", "LF", "Consonants"), ordered = TRUE)

# Data as list for stan input
components <- d %$% as.numeric(mapvalues(component, from = levels(component), to= 1:length(unique(component))))
dat <-  within( list(), 
                {
                  N <- nrow(d)
                  y <- d$IKI # DV
                  components <- components
                  K <- max(components)
                  subj <- d$subj
                  S <- max(d$subj)
                  bigram <- as.integer(factor(d$bigram))
                  B <- length(unique(d$bigram))
                } );str(dat)

# Initialise start values
start <- 
  function(chain_id = 1){
    list(beta = rep(mean(log(dat$y)),dat$K) 
         , sigma = sd(log(dat$y)),
         z_u = matrix(rep(0.1,dat$S*dat$K), nrow = dat$K),
         w = rep(0.1, length(unique(d$bigram))),
         L_u = matrix(rep(0, dat$K*dat$K), nrow = dat$K),
         sigma_u = rep(1, dat$K),
         sigma_w = 1,
         alpha = chain_id
    )
  }

start_ll <- lapply(1:n_chain, function(id) start(chain_id = id) )

# --------------
# Stan models ##
# --------------
# Load model
lmm <- stan_model(file = "stanin/LMM.stan")

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
              save_warmup = FALSE,
              include = FALSE,
              pars = c("mu", "L_u", "z_u"),
              seed = 81,
              control = list(max_treedepth = 16,
                             adapt_delta = 0.99)
)

#Save model
saveRDS(m,
        file="stanout/LMM.rda",
        compress="xz")


# Extract and save posterior and log likelihood seperately
# Get log likelihood
log_lik <- extract_log_lik(m) 
saveRDS(log_lik, 
        file = "stanout/log_lik/LMM_loglik.rda",
        compress = "xz")

#ll <- readRDS(file="stanout/log_lik/LMMsubjslopes_loglik.rda")
#loo.ll <- loo(log_lik)

# Get parameter posterior
param <- c("beta", "sigma")
samps <- as.data.frame(m, pars = param) 
saveRDS(samps, 
        file = "stanout/posterior/LMM_posterior.rda",
        compress = "xz")


# Get random effects
param <- c("u", "w", "sigma_u", "sigma_w")
re <- as.data.frame(m, pars = param) 
saveRDS(re, 
        file = "stanout/RE/LMM_re.rda",
        compress = "xz")


# Get posterior predicted values
param <- c("y_tilde")
y_tilde <- as.data.frame(m, pars = param) 
saveRDS(y_tilde, 
        file = "stanout/y_tilde/LMM_y_tilde.rda",
        compress = "xz")

