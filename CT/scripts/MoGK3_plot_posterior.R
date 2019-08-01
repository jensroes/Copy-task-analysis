library(rethinking)
library(tidyverse)
library(rstan)
library(gridExtra)
source("functions/functions.R")

# Data
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

# Model
m <- readRDS(file="stanout/posterior/MoGK2_posterior.rda")

#(param <- names(m)[1:6])

names(m)
m %>% gather(Parameter, value) %>%
  group_by(Parameter) %>%
  summarise_all(list(map = dmode, 
                     l.hpdi = ~HPDI(., prob = .95)[1], 
                     u.hpdi = ~HPDI(., prob = .95)[2])) %>%
  ungroup() %>%
  as.data.frame() %>%
  mutate(map = round(map, 2),
         l.hpdi = round(l.hpdi, 2),
         u.hpdi = round(u.hpdi, 2)) %>%
  filter(!(Parameter %in% c(paste0("gamma[",1:5,"]"),
                            paste0("delta[",1:5,"]"),
                            paste0("sigma_diff[",1:5,"]"),
                            paste0("sigma_diff2[",1:5,"]"),
                            paste0("sigma[",1:5,"]")))) %>%
  separate(Parameter, into = c("Parameter", "pmix"), sep = "[,]") %>%
  mutate(Parameter = gsub(pattern = "_", ".", Parameter)) %>%
  mutate(Parameter= gsub(pattern = "]", replacement = "", Parameter)) %>%
  mutate(Parameter= gsub(pattern = "\\[", replacement = "_", Parameter)) %>%
  mutate(pmix = ifelse(is.na(pmix), "", pmix)) %>%
  mutate(pmix = gsub(pattern = "]", replacement = "", pmix)) %>%
  separate(Parameter, into = c("Parameter", "Component"), sep = "_") %>%
  mutate(Parameter = paste0(Parameter, pmix)) %>%
  select(Parameter, Component, map) %>%
  mutate(Parameter = ifelse(Parameter == "beta", "beta2", 
                            ifelse(Parameter == "sigma.e", "sigma2",
                                   ifelse(Parameter == "sigmap.e", "sigma3",
                                          ifelse(Parameter == "sigmap.ee", "sigma1", Parameter))))) -> samps
  
for(p in unique(samps$Parameter)){
  for(i in unique(samps$Component)){
    tmp <- samps[samps$Parameter == p & samps$Component == i,]$map
  #  if(!(p %in% paste0("theta",1:3))){
  #    tmp <- exp(tmp)
  #  }
    assign(paste0(p, "_", i), tmp)
  }
}

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dlnorm(x, mu, sigma, log= FALSE)
}

xmax = 5000
binwidth = 20
alpha = .25  
theme_set(theme_minimal() +
            theme(axis.text.y = element_blank(),
                  plot.margin = margin(0.1, 0, 0, 0, "cm"),
                  legend.position = c(.9,.9),
                  legend.justification = c("right"),
                  legend.box.just = "right",
#                  legend.direction="horizontal",
                  legend.margin = margin(6, 6, 6, 6),
                  legend.title = element_text(face = "bold", size = 10, vjust = .5),
                  legend.text = element_text(size = 10),
                  legend.key.height =  unit(.75,"line"),
                  legend.key.width = unit(.75, "cm"),
                  panel.grid.minor = element_blank(),
                  plot.subtitle = element_text(size = 10),
                  plot.title = element_text(size = 14, face = "bold")) )

p.tap <- d %>% 
  filter(component == "Tapping") %>%
  mutate(x = IKI) %>%
  ggplot() +
  geom_histogram(aes(x, ..density.., fill = "data"), binwidth = binwidth, color = "white", alpha = alpha) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta1_1, sigma = sigma1_1, lam = theta1_1 ), aes(fill = "K1"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta2_1, sigma = sigma2_1, lam = theta2_1 ), aes(fill= "K2"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta3_1, sigma = sigma3_1, lam = theta3_1 ), aes(fill = "K3"), lwd = .25) +
  scale_fill_manual("Mixture\ncomponents:", 
                    breaks = c("K1", "K2", "K3", "data"),
                    values = c(
                      "data" = alpha("grey70", alpha),
                      "K1" = alpha("darkred", alpha),
                      "K2" = alpha("darkblue", alpha),
                      "K3" = alpha("forestgreen", alpha)),
                    labels = c(bquote(italic("K")[1]),
                               bquote(italic("K")[2]),
                               bquote(italic("K")[3]),
                               "data") 
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.1))) +
  scale_x_continuous(breaks = seq(0, 1000, 100), limits = c(0, 800)) +
  scale_colour_manual("Mixture components:", 
                      values = c("grey80", "darkred", "darkblue", "forestgreen"),
                      breaks = c("K1", "K2", "K3", "data"),
                      labels = c(bquote(italic("K")[1]),
                                 bquote(italic("K")[2]),
                                 bquote(italic("K")[3]),
                                 "data") 
                      ) +
  ggtitle("Tapping") +
  labs(y = "", x = "")


p.sent <- d %>% 
  filter(component == "Sentence") %>%
  mutate(x = IKI) %>%
  ggplot() +
  geom_histogram(aes(x, ..density.., fill = "data"), binwidth = binwidth, color = "white", alpha = alpha) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta1_2, sigma = sigma1_2, lam = theta1_2 ), aes(fill = "K1"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta2_2, sigma = sigma2_2, lam = theta2_2 ), aes(fill= "K2"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta3_2, sigma = sigma3_2, lam = theta3_2 ), aes(fill = "K3"), lwd = .25) +
  scale_fill_manual("Mixture\ncomponents:", 
                    breaks = c("K1", "K2", "K3", "data"),
                    values = c(
                      "data" = alpha("grey70", alpha),
                      "K1" = alpha("darkred", alpha),
                      "K2" = alpha("darkblue", alpha),
                      "K3" = alpha("forestgreen", alpha)),
                    labels = c(bquote(italic("K")[1]),
                               bquote(italic("K")[2]),
                               bquote(italic("K")[3]),
                               "data") 
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.1))) +
  scale_x_continuous(breaks = seq(0, 1000, 100), limits = c(0, 1000)) +
  scale_colour_manual("Mixture components:", 
                      values = c("grey80", "darkred", "darkblue", "forestgreen"),
                      breaks = c("K1", "K2", "K3", "data"),
                      labels = c(bquote(italic("K")[1]),
                                 bquote(italic("K")[2]),
                                 bquote(italic("K")[3]),
                                 "data") 
  ) +
  ggtitle("Sentence") +
  labs(y = "", x = ""); p.sent


p.hf <- d %>% 
  filter(component == "HF") %>%
  mutate(x = IKI) %>%
  ggplot() +
  geom_histogram(aes(x, ..density.., fill = "data"), binwidth = binwidth, color = "white", alpha = alpha) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta1_3, sigma = sigma1_3, lam = theta1_3 ), aes(fill = "K1"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta2_3, sigma = sigma2_3, lam = theta2_3 ), aes(fill= "K2"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta3_3, sigma = sigma3_3, lam = theta3_3 ), aes(fill = "K3"), lwd = .25) +
  scale_fill_manual("Mixture\ncomponents:", 
                    breaks = c("K1", "K2", "K3", "data"),
                    values = c(
                      "data" = alpha("grey70", alpha),
                      "K1" = alpha("darkred", alpha),
                      "K2" = alpha("darkblue", alpha),
                      "K3" = alpha("forestgreen", alpha)),
                    labels = c(bquote(italic("K")[1]),
                               bquote(italic("K")[2]),
                               bquote(italic("K")[3]),
                               "data") 
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.1))) +
  scale_x_continuous(breaks = seq(0, 1000, 100), limits = c(0, 1000)) +
  scale_colour_manual("Mixture components:", 
                      values = c("grey80", "darkred", "darkblue", "forestgreen"),
                      breaks = c("K1", "K2", "K3", "data"),
                      labels = c(bquote(italic("K")[1]),
                                 bquote(italic("K")[2]),
                                 bquote(italic("K")[3]),
                                 "data") 
  ) +
  ggtitle("High-frequency bigrams") +
  labs(y = "", x = ""); p.hf

p.lf <- d %>% 
  filter(component == "LF") %>%
  mutate(x = IKI) %>%
  ggplot() +
  geom_histogram(aes(x, ..density.., fill = "data"), binwidth = binwidth, color = "white", alpha = alpha) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta1_4, sigma = sigma1_4, lam = theta1_4 ), aes(fill = "K1"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta2_4, sigma = sigma2_4, lam = theta2_4 ), aes(fill= "K2"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta3_4, sigma = sigma3_4, lam = theta3_4 ), aes(fill = "K3"), lwd = .25) +
  scale_fill_manual("Mixture\ncomponents:", 
                    breaks = c("K1", "K2", "K3", "data"),
                    values = c(
                      "data" = alpha("grey70", alpha),
                      "K1" = alpha("darkred", alpha),
                      "K2" = alpha("darkblue", alpha),
                      "K3" = alpha("forestgreen", alpha)),
                    labels = c(bquote(italic("K")[1]),
                               bquote(italic("K")[2]),
                               bquote(italic("K")[3]),
                               "data") 
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.1))) +
  scale_x_continuous(breaks = seq(0, 1000, 100), limits = c(0, 1000)) +
  scale_colour_manual("Mixture components:", 
                      values = c("grey80", "darkred", "darkblue", "forestgreen"),
                      breaks = c("K1", "K2", "K3", "data"),
                      labels = c(bquote(italic("K")[1]),
                                 bquote(italic("K")[2]),
                                 bquote(italic("K")[3]),
                                 "data") 
  ) +
  ggtitle("Low-frequency bigrams") +
  labs(y = "", x = ""); p.lf

p.con <- d %>% 
  filter(component == "Consonants") %>%
  mutate(x = IKI) %>%
  ggplot() +
  geom_histogram(aes(x, ..density.., fill = "data"), binwidth = binwidth, color = "white", alpha = alpha) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta1_5, sigma = sigma1_5, lam = theta1_5 ), aes(fill = "K1"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta2_5, sigma = sigma2_5, lam = theta2_5 ), aes(fill= "K2"), lwd = .25) +
  geom_area(stat = "function", fun = plot_mix_comps, alpha = alpha, args = list(mu = beta3_5, sigma = sigma3_5, lam = theta3_5 ), aes(fill = "K3"), lwd = .25) +
  scale_fill_manual("Mixture\ncomponents:", 
                    breaks = c("K1", "K2", "K3", "data"),
                    values = c(
                      "data" = alpha("grey70", alpha),
                      "K1" = alpha("darkred", alpha),
                      "K2" = alpha("darkblue", alpha),
                      "K3" = alpha("forestgreen", alpha)),
                    labels = c(bquote(italic("K")[1]),
                               bquote(italic("K")[2]),
                               bquote(italic("K")[3]),
                               "data") 
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.1))) +
  scale_x_continuous(breaks = seq(0, 3000, 500), limits = c(0, 3000)) +
  scale_colour_manual("Mixture components:", 
                      values = c("grey80", "darkred", "darkblue", "forestgreen"),
                      breaks = c("K1", "K2", "K3", "data"),
                      labels = c(bquote(italic("K")[1]),
                                 bquote(italic("K")[2]),
                                 bquote(italic("K")[3]),
                                 "data") 
  ) +
  ggtitle("Consonants") +
  labs(y = "", x = ""); p.con


p <- grid.arrange(p.tap,p.sent,p.hf,p.lf,p.con) 



