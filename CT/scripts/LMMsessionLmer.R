# Load packages
library(tidyverse)
library(plyr)
library(magrittr)
library(lme4)

# Test-retest analysis only on subjects that took part in two sessions

# Load df
d <- read_csv("data/ct.csv") %>% 
  dplyr::filter(IKI > 0 
                ,!is.na(IKI)
                ,target == 1
  ) %>%
  mutate(subj = as.numeric(factor(subj))) %>%
  group_by( subj, bigram, component, session ) %>%
  dplyr::summarise(IKI = mean(IKI),
                   N = n()) %>%
  ungroup() %>%
  mutate(component = factor(component, levels = c("Tapping", "Sentence", "HF", "LF", "Consonants"), ordered = TRUE))

d %>% dplyr::select(subj, session) %>% unique() %>% dplyr::count(subj) %>% arrange(desc(n)) %>% filter(n > 1) -> keep.subj

d %>% 
  filter(subj %in% keep.subj$subj) %>%
  mutate(session = as.integer(factor(session)),
         subj = as.integer(factor(subj))) %>%
  filter(session %in% 1:2) %>%
  mutate(session = session -1) -> d

(subj_left <- length(unique(d$subj))) # 239 subjects

contrasts(d$component) <- contr.treatment(5)
d$session <- factor(d$session)
contrasts(d$session) <- contr.sum(2)

m <- lmer(log(IKI) ~ 0 + component * session + (1|subj) + (1|bigram), d)
round(summary(m)$coef,2)
