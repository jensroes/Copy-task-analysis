library(tidyverse)

# Save model
#m <- readRDS("CopyTask2019/stanout/MoGK3.rda")

# Get parameter posterior
param <- c("beta1", "beta", "beta3", "delta", "gamma", "theta", "sigma", "sigma_diff", "sigma_diff2", "sigmap_e", "sigmap_ee", "sigma_e")
samps <- as.data.frame(m, pars = param) 
saveRDS(samps, 
        file = "CopyTask2019/stanout/posterior/MoGK3_posterior.rda",
        compress = "xz")


# Get random effects
param <- c("u", "w", "sigma_u", "sigma_w")
re <- as.data.frame(m, pars = param) 
saveRDS(re, 
        file = "stanout/RE/MoGK3_re.rda",
        compress = "xz")

# Get posterior predicted values
param <- c("y_tilde")
readRDS("stanout/MoGK3.rda") %>% as.data.frame(pars = param) -> y_tilde

saveRDS(y_tilde, 
        file = "stanout/y_tilde/MoGK3_y_tilde.rda",
        compress = "xz")

