# Load packages
library(tidyverse)
library(rstan)
library(plyr)
library(magrittr)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Sampling parameters
set.seed(125)
n_chain = 3 
iterations = 4000

# Load df
d <- read_csv("data/ct.csv") %>% 
  dplyr::filter(IKI > 0 
#         ,IKI < 10000
         ,target == 1
         ,subj %in% c(1:10)
         ) %>%
  group_by( subj, bigram, component ) %>%
  dplyr::summarise(IKI = mean(IKI),
                   N = n()) %>%
  ungroup() %>%
  mutate(component = factor(component, levels = c("Tapping", "Sentence", "HF", "LF", "Consonants"), ordered = TRUE)); d 

d %>% 
  mutate(logIKI = log(IKI)) %>%
  gather(DV, IKI, -N,-subj,-bigram,-component) %>%
#  filter(IKI < 1500) %>%
  mutate(DV = ifelse(DV == "IKI", "IKI (in ms)",
                     ifelse(DV == "logIKI", "log IKI", NA ))) %>%
  ggplot(aes(x=IKI, fill = component)) +
#  geom_histogram(alpha = .5, bins = 60) +
#  geom_density(alpha = .5) +
  geom_histogram(aes(y=..density..), alpha = .5, bins = 60) + 
  geom_density(aes(y=..density..), alpha = .5, size =.1 ) +
  facet_wrap(c("component","DV"),  scales = "free", ncol = 2, strip.position = "bottom") + 
  xlab("") +
  scale_fill_hue("Components: ", l = 20, c = 70) +
  theme_light() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

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

# --------------
# Stan models ##
# --------------
#---- 
# Mixture of three gaussians with unequal variance (unconstraint)
#---- 
# Load model
mog <- stan_model(file = "stanin/MoGK3var.stan")

# Check model
m <- sampling(mog, chain = 1, iter = 1, data = dat) 

# Fit model
m <- sampling(mog, 
              data = dat,
              iter = iterations,
              init = start_ll,
              warmup = iterations/2,
              chains = n_chain, 
              refresh = 250,
              thin = 1,
              seed = 81,
              control = list(max_treedepth = 16,
                             adapt_delta = 0.99)
)


# Save model
saveRDS(m, 
        file = "stanout/log_lik/MoGK3var.rda",
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

