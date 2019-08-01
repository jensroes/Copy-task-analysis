library(tidyverse)
library(magrittr)
iteratons <- 6e3
path <- "CopyTask2019/stanout/y_tilde/"
(files <- dir(path, pattern = "MoG"))

# Load df
d <- read_csv("CopyTask2019/data/ct.csv") %>% 
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

yt <- readRDS(paste0(path,files))
samps <- paste0("V", sample(1:iteratons, 500))
ytilde <- yt %>%
  as.matrix(pars = "y_tilde") %>%
  t() %>%
  as_tibble() %>%
  select(samps) %>%
  bind_cols(d,.) %>%
  gather(iteration, y_tilde, samps) 

ytilde %>%
  count(y_tilde < 10000)

ggplot(data = NULL, aes(x = y_tilde, 
                        group = interaction(iteration, component))) +
  theme_linedraw() +
  geom_density(data = ytilde, aes(linetype = "predicted", color = "predicted"), alpha = .05, size = .1) +
  geom_density(data = d, aes(x=IKI, group = component,linetype = "observed", color = "observed"),  alpha = .25) +
  scale_x_continuous(limits = c(0, 3000)) +
  labs(x = "IKI in msecs") +
#  labs(x = bquote(tilde("IKI"))) +
  facet_wrap(~component, scales = "free") +
  scale_linetype_manual("Copy-task\ncomponents:", 
                        values=c(predicted="dotted", observed="solid")) +
  scale_color_manual("Copy-task\ncomponents:", 
                        values=c(predicted="darkred", observed="black")) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(.85,.2))
  
  

