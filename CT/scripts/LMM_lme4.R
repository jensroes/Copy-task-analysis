# Load packages
library(tidyverse)
library(plyr)
library(magrittr)
library(lme4)
library(rethinking)
library(wesanderson)
source("functions/functions.R")

# Load df
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

# Linear mixed model
lmm0 <- lmer(log(IKI) ~ 1 + (1|subj) + (1|bigram), d)
lmm <- lmer(log(IKI) ~ 0 + component + (1|subj) + (1|bigram), d)
summary(lmm)
anova(lmm0,lmm)

#cis <- confint(lmm, method="boot")
#cis %>% as.data.frame() %>%
#  mutate(Parameter = rownames(.)) %>%
#  as_tibble() %>%
#  dplyr::rename(l.ci = `2.5 %`,
#                u.ci = `97.5 %`) %>%
#  select(Parameter, l.ci, u.ci) -> cis_df
#write_csv(cis_df, "data/lmm_confint.csv")
cis <- read_csv("data/lmm_confint.csv")

lmm_summary <- summary(lmm)$coef %>% round(2) %>% 
  as.data.frame() %>%
  mutate(Parameter = rownames(.)) %>%
  as_tibble() %>%
  left_join(cis) %>%
  mutate(Parameter = gsub(pattern = "component", "", Parameter)) %>%
  dplyr::rename(SE = `Std. Error`,
         t = `t value`) %>%
  select(Parameter, Estimate, SE, t, l.ci, u.ci) %>%
  mutate(Estimate = exp(Estimate),
         SE = exp(SE),
         l.ci = exp(l.ci),
         u.ci = exp(u.ci)); lmm_summary

# Get posterior of Bayesian LMM
blmm <- readRDS("stanout/posterior/LMM_posterior.rda")
blmm_summary <- blmm %>% as_tibble() %>%
  gather(Parameter, value) %>%
  group_by(Parameter) %>%
  dplyr::summarise(mean = mean(value),
                   map = dmode(value),
                   SE = sd(value),
                   l.hpdi = HPDI(value, prob = .95)[1],
                   u.hpdi = HPDI(value, prob = .95)[2]
 #                  l.pi = PI(value, prob = .95 )[1]
  #                 ,u.pi = PI(value, prob = .95)[2]
                   ) %>%
  mutate(mean = exp(mean),
         map = exp(map),
         SE = exp(SE),
         l.hpdi = exp(l.hpdi),
         u.hpdi = exp(u.hpdi)
         #,l.pi = exp(l.pi)
         #,u.pi = exp(u.pi)
         ) %>%
  filter(Parameter != "sigma"); blmm_summary

# Compare model parameters
blmm_summary$Parameter <- lmm_summary$Parameter
blmm_summary$model <- "BLMM"
lmm_summary$model <- "LMM"

mc <- blmm_summary %>% select(-SE, -mean) %>%
  dplyr::rename(Estimate = map,
                l.ci = l.hpdi,
                u.ci = u.hpdi) %>%
  bind_rows(lmm_summary[,-c(3,4)])

ctc <- c("Tapping", "Sentence", "HF bigrams", "LF bigrams", "Consonants")
p.cis <- mc %>% 
  mutate(Parameter = ifelse(Parameter == "LF", "LF bigrams", 
                            ifelse(Parameter == "HF", "HF bigrams", Parameter))) %>%
  mutate(Parameter = factor(Parameter, levels = rev(ctc), ordered = T),
         model = factor(model, levels = c("LMM", "BLMM"), ordered = T)) %>%
  ggplot(aes(x=Parameter, linetype = model)) +
  theme_bw() +
  geom_point(aes(y=Estimate), size = 2.5, position = position_dodge(-.5)) +
  geom_errorbar(aes(ymin = l.ci, ymax = u.ci), width = 0, position = position_dodge(-.5), size = .5) +
  labs(y = bquote("IKI in msecs"),
       x = "Copy-task\ncomponents:") +
  coord_flip() +
  scale_linetype_manual("Model: ", values = c("dashed", "dotted") ) +
  theme(axis.title.y = element_text(angle = 0, hjust = 0),
        legend.key.width = unit(1.5,"cm"),
        legend.position = "right",
        legend.background = element_blank());p.cis

p.probs <- blmm %>% as_tibble() %>%
  gather(Parameter, value) %>%
  filter(Parameter != "sigma") %>%
  mutate(value = exp(value),
         Parameter = mapvalues(Parameter, from = unique(Parameter), to = ctc)) %>%
  mutate(Parameter = factor(Parameter, levels = ctc, ordered = T)) %>%
  ggplot(aes(x=value)) +
  theme_linedraw() +
  geom_density(fill = "grey40", alpha = .15) +
  facet_wrap(~Parameter, scales = "free", ncol = 2) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "IKI in msecs"); p.probs





