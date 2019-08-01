# This script generates a table for the LMM and LMM0 comparison
# I'll probably overwrite the table that is used here (because there's an error) so don't run this script again.

library(tidyverse)

read_csv("stanout/loo_results_mog.csv") %>%
  filter(model %in% c("loo_LMM_loglik", "loo_LMMsubjintercepts_loglik")) %>%
  mutate(elpd_diff = 0) -> elpd

elpd$elpd_diff[2] <- diff(elpd$elpd_loo)
elpd$se_elpd_diff[1] <- 0

write_csv(elpd, "stanout/loo_results_lmm.csv")

