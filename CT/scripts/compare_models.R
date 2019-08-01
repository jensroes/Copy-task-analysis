library(loo)
options(mc.cores = 3)
library(tidyverse)
library(magrittr)

path <- "CopyTask2019/stanout/log_lik/"
(files <- dir(path, pattern = "LMM|MoG"))
ms <- gsub(files, pattern = ".rda", replacement = "")

for(i in 1:length(files)){
  tmp <- readRDS(paste0(path,files[i]))
  assign(paste0("loo_",ms[i]), loo(tmp))
  print(ms[i]); rm(list = "tmp")
}

loos <- ls(pattern = "loo_")
mc <- do.call(compare, lapply(loos, as.name))

mc %<>% as.data.frame() %>%
  round(2) %>%
  mutate(model=row.names(.)) %>%
  select(model, elpd_diff:se_elpd_loo);mc

# Get se diff
mc$se_elpd_diff <- 0
for(i in 2:nrow(mc)){
  mc$se_elpd_diff[i] <- do.call(compare, lapply(mc$model[(i-1):i], as.name))[2]
}

mc %<>% select(model, elpd_diff, se_elpd_diff, elpd_loo, se_elpd_loo)
#do.call(plot, lapply(loos[4], as.name))
#eval(as.name(loos[4]))
file_out <- "CopyTask2019/stanout/loo_results_mog.csv"
write_csv(mc, file_out)
#read_csv(file_out)


