library(tidyverse)
library(magrittr)

d <- read_csv("data/ct.csv") %>% filter(target == 1) %>%
  mutate(component = factor(component, levels = c("Tapping", "Sentence", "HF", "LF", "Consonants"), ordered = TRUE))

minIKI <- 50
midIKI <- 2000
maxIKI <- 6000 
d %>% 
  filter(IKI <= minIKI | IKI >= midIKI) %>%
  mutate(Extr = ifelse(IKI >= maxIKI, paste(">", maxIKI,"ms", sep = ""),
                       ifelse(IKI >= midIKI & IKI <= maxIKI, 
                              paste(">", midIKI, "ms", ", <", maxIKI, "ms", sep = ""),
                                    ifelse(IKI <= minIKI, paste("<", minIKI, "ms", sep = "") , NA)))) %>%
  ggplot(aes(y=IKI, x=1)) +
  geom_jitter(size =.01, alpha =.5) +
  facet_grid(Extr~., scales = "free" ) +
  theme_minimal() +
  ylab("IKI (in ms)") +
  theme(strip.text.y = element_text(face = "bold", size = 12, angle = 0),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())


ggsave("plots/thresholds.pdf", family="Times", bg = "transparent", width = 10, height = 7)


d %>% filter(IKI > 0, IKI < 10000) %>%
  group_by( subj, bigram, component, session, age ) %>%
  dplyr::summarise(IKI = mean(IKI),
                   N = n()) %>%
  ungroup() -> d.meanIKI

d.meanIKI %>%
  count(component)

d.meanIKI %>%
filter(IKI < 2500) %>% 
  ggplot(aes(x = IKI, 
             color = factor(component), 
             fill = factor(component ))) + 
  #  geom_histogram(aes(y=..density..), alpha = .3, bins = 60) + 
  geom_density(aes(y=..density..), alpha = .3, size =.5 ) +
  ylab(" ") + 
  xlab("IKI (in ms)") + 
  theme_minimal() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        legend.key=element_blank(),
        legend.position= "top", 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size =12)) +
  scale_color_hue(name = "Components:", c =70, l = 40) +
  scale_fill_hue(name = "Components:", c =70, l = 40)
#  scale_linetype_discrete(name = "Session:") 

ggsave("plots/crosscompN2000.pdf", family="Times", bg = "transparent", width = 10, height = 7)


# Test Rest
with(d.meanIKI, table(subj, session))

d.meanIKI %>% select(subj, session) %>% unique() %>% count(subj) %>% arrange(desc(n)) %>% filter(n > 1) -> keep.subj

d.meanIKI %>% 
  filter(subj %in% keep.subj$subj) %>%
  group_by(subj) %>%
  mutate(session = as.integer(factor(session))) %>%
  ungroup() %>%
  filter(session %in% 1:2) -> d.2session
#  select(subj, session) %>%
#  unique() %>%
#  count(subj) %>%
#  filter(n != n)

(subj_left <- length(unique(d.2session$subj)))


d.2session %>% 
#  filter(IKI > 30, IKI < 1500) %>% 
  ggplot(aes(x = IKI, 
             color = factor(session) 
             #fill = factor(session) 
             )) + 
  geom_histogram(aes(y=..density..), alpha = .3, bins = 60,  position="identity", fill="white", color = "grey" ) + 
  geom_density(aes(y=..density..), alpha = .3, size =.5 ) +
  ylab(" ") + 
  xlab("IKI (in ms)") + 
  theme_light() +
  facet_wrap(~component, scales = "free") +
  scale_color_hue(name = "Session: ", c =70, l = 40) +
  ggtitle(paste0("N participants = ", subj_left)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        legend.key=element_blank(),
        legend.position= "top", 
        legend.justification = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        text = element_text(size = 12),
        plot.title = element_text(margin=margin(b = 1, unit = "pt"), hjust = 1),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size =12)) 

ggsave("plots/crosscompsession.pdf", family="Times", bg = "transparent", width = 10, height = 7)



d.meanIKI %>% 
  mutate(logIKI = log(IKI)) %>%
  gather(DV, IKI, -N,-subj,-bigram,-component,-session,-age) %>%
  #  filter(IKI < 1500) %>%
  mutate(DV = ifelse(DV == "IKI", "IKI (in ms)",
                     ifelse(DV == "logIKI", "log IKI", NA ))) %>%
  ggplot(aes(x=IKI, fill = component)) +
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

ggsave("plots/crosscomptransformation.pdf", family="Times", bg = "transparent", width = 6, height = 10)


p.age <- d.meanIKI %>%
  group_by(subj, component, age) %>%
  dplyr::summarise(IKI = mean(IKI)) %>%
  ungroup() %>%
  filter(!is.na(age)) %>%
  ggplot(aes(y=IKI, x=age, color = component)) +
  scale_x_continuous(trans='log10', breaks = seq(10, 80, 10), limits = c(13, 83)) +
  scale_y_continuous(trans='log10', breaks = seq(0, 2500, 250), limits = c(50, 2500)) +
#  geom_smooth(method = "gam",
#              formula = y ~ s(x, bs = "ad")) +
  geom_jitter(size = .25, width = .1, alpha = .2) +
  theme_minimal() +
  scale_color_hue(name = "Components: ", c =70, l = 40) +
  ggtitle(paste0("N participants = ", length(unique(d.meanIKI[!is.na(d.meanIKI$age),]$subj))
                 )) +
  labs(y = "IKI (in ms)",
      x = "Age (in year)"
      ) +
  stat_smooth(aes(x = age, y = IKI), method = "lm", formula = y ~ poly(x,2),  se = F, size = 1) +
  theme(axis.line.y = element_blank(),
        legend.key=element_blank(),
        legend.position= "top", 
        legend.justification = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.width=unit(.75,"cm"),
        legend.key.height = unit(.75,"cm"),
        text = element_text(size = 14),
        plot.title = element_text(margin=margin(b = 1, unit = "pt"), hjust = 1, size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size =14)) ; p.age

ggsave("plots/crosscomp_age.pdf", family="Times", bg = "transparent", width = 10, height = 8)



p.age <- d.meanIKI %>%
  mutate(component = ifelse(component %in% c("Tapping", "Consonants"), "Non-lexical", "Lexical")) %>%
  group_by(subj, component, age) %>%
  summarise(IKI = mean(IKI)) %>%
  ungroup() %>%
  filter(!is.na(age)) %>%
  ggplot(aes(y=IKI, x=age, color = component)) +
  scale_x_continuous(trans='log10', breaks = seq(10, 80, 10), limits = c(13, 83)) +
  scale_y_continuous(trans='log10', breaks = seq(0, 2500, 250), limits = c(50, 2500)) +
  #  geom_smooth(method = "gam",
  #              formula = y ~ s(x, bs = "ad")) +
  geom_jitter(size = .25, width = .1, alpha = .2) +
  theme_minimal() +
  scale_color_hue(name = "Components: ", c =70, l = 40) +
  ggtitle(paste0("N participants = ", length(unique(d.meanIKI[!is.na(d.meanIKI$age),]$subj))
  )) +
  labs(y = "IKI (in ms)",
       x = "Age (in year)"
  ) +
  stat_smooth(aes(x = age, y = IKI), method = "lm", formula = y ~ poly(x,2),  se = F, size = 1) +
  theme(axis.line.y = element_blank(),
        legend.key=element_blank(),
        legend.position= "top", 
        legend.justification = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.width=unit(.75,"cm"),
        legend.key.height = unit(.75,"cm"),
        text = element_text(size = 14),
        plot.title = element_text(margin=margin(b = 1, unit = "pt"), hjust = 1, size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size =14)) ; p.age

ggsave("plots/crosscomp_age_lexical.pdf", family="Times", bg = "transparent", width = 10, height = 8)

  
