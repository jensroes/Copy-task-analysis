library(rethinking)
library(tidyverse)
library(rstan)
library(gridExtra)
source("functions/functions.R")

# Model
m <- readRDS(file="stanout/posterior/MoGK3_posterior.rda")

names(m)
m %>% gather(Parameter, value) %>%
  separate(Parameter, into = c("Parameter", "Component"), sep = "\\[") %>%
  mutate(Component = gsub(Component, pattern = "\\]", replacement = ""),
         value = exp(value)) %>%
  filter(Parameter  %in% c("beta1", "beta", "beta3")) %>%
  mutate(Parameter = ifelse(Parameter == "beta", "beta2", Parameter)) %>%
  filter(Component == 1)%>%
  ggplot(aes(x=value, color = Parameter)) +
  geom_density() +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10))
  
  
