library(tidyverse)
library(plyr)
ctc <- c("Tapping", "HF bigrams", "Sentence", "LF bigrams", "Consonants")
d.full <- read_csv("data/ct.csv") 

d.full %>% 
  select(subj, age) %>% unique() %>%
  filter(!is.na(age)) %>%
  select(subj) %>% unique()

d.full %>% 
  dplyr::filter(IKI > 0 
                ,!is.na(IKI)
                ,target == 1
  ) %>%
  mutate(subj = as.numeric(factor(subj))) %>%
  group_by( subj, component, age ) %>%
  dplyr::summarise(IKI = mean(IKI),
                   N = n()) %>%
  ungroup() %>%
  mutate(component = mapvalues(component, from  = rev(c("Tapping", "HF", "Sentence", "LF", "Consonants")), to = rev(ctc))) %>%
  mutate(component = factor(component, levels = rev(ctc), ordered = TRUE)) %>%
  filter(!is.na(age)) %>%
  ggplot(aes(y=IKI, x=age, color = component, linetype = component)) +
  scale_x_continuous(trans='log10', breaks = seq(10, 80, 10), limits = c(13, 83)) +
  scale_y_continuous(trans='log10', breaks = seq(0, 4000, 500), limits = c(50, 3500)) +
  geom_jitter(size = .45, width = .1, alpha = .25) +
  theme_minimal() +
  scale_color_hue(name = "Copy-task\ncomponent: ", c =80, l = 40) +
  scale_linetype_manual("Copy-task\ncomponent: ", values = c("dashed", "solid", "solid", "solid", "dashed") ) +
  labs(y = "IKI in msecs",
       x = "Age in years"
  ) +
  stat_smooth(aes(x = age, y = IKI), method = "lm", formula = y ~ poly(x,2),  se = F, size = 1) +
  theme(legend.key=element_blank(),
        legend.position= "right", 
        legend.justification = "top",
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.width=unit(1.25,"cm"),
        legend.key.height = unit(1,"cm"))


ggsave("plots/crosscomp_age_lexical.png", family="Times", bg = "transparent", width = 8, height = 6)

