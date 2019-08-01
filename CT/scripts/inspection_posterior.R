library(tidyverse)
library(rstan)

dir("stanout")
m <- readRDS("stanout/MoGK2.rda")
names(m)

excl.param <- c("log_lik", "y_tilde", "w", "u")
print(m, include = FALSE, pars = excl.param, probs = c(.025, .975))
traceplot(m, include = FALSE, pars = excl.param)
