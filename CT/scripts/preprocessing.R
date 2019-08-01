library(tidyverse)
library(haven)
library(sjlabelled)
library(magrittr)
outputdf <- "ct.csv"
path_data <- "data"

(files <- dir(path_data, pattern = ".sav"))
(file <- paste(path_data, files[2], sep = "/"))
ppt_info <- paste(path_data, files[1], sep = "/")

ppt <- read_sav(ppt_info) %>% as_tibble() %>% 
  select(-Language, -Age_categories, -Handedness_groups, -Handedness_Score,-LogCreationDate, -Correctness, -Education,-ID) %>%
  mutate(participant = Partcipant) %>%
  select(-Partcipant) %>%
  mutate(participant = gsub('[[:punct:] ]+',' ',participant)) %>%
  mutate(participant = gsub('[[:digit:]]+', '', participant)) %>%
  mutate(participant = gsub(" ", "", participant, fixed = TRUE)) %>%
  mutate(subj = trimws(tolower(participant))) %>%
  mutate(logfile = gsub('[[:punct:] ]+',' ',logfile)) %>%
  mutate(logfile = gsub(" ", "", logfile, fixed = TRUE)) %>%
  mutate(logfile = trimws(tolower(logfile))) %>%
  dplyr::rename(sex=Gender__G,
         session = test_retest) %>%
  mutate(sex = ifelse(sex == "Vrouw", "f",
                      ifelse(sex == "Man", "m", NA))) %>%
  select(-participant, -subj)
ppt$session <- remove_all_labels(ppt$session)

str(ppt)
attr(ppt$session, "format.spss") <- NULL
attr(ppt$logfile, "format.spss") <- NULL
attr(ppt$logfile, "label") <- NULL
attr(ppt$logfile, "display_width") <- NULL
attr(ppt$Age, "format.spss") <- NULL
attr(ppt$Age, "label") <- NULL
attr(ppt$Age, "display_width") <- NULL


d <- read_sav(file) %>% as.tibble()
d %>% mutate(participant = gsub('[[:punct:] ]+',' ',participant)) %>%
      mutate(participant = gsub('[[:digit:]]+', '', participant)) %>%
      mutate(participant = gsub(" ", "", participant, fixed = TRUE)) %>%
     mutate(logfile = gsub('[[:punct:] ]+',' ',logfile)) %>%
     mutate(logfile = gsub(" ", "", logfile, fixed = TRUE)) %>%
     mutate(logfile = trimws(tolower(logfile))) %>%
     filter(!(logfile == "logfile")) %>%
     transmute(bigram = trimws(bigram),
            component = trimws(component),
            freq = trimws(bigr_frequency_class),
            trial = com_trial,
            subj = trimws(tolower(participant)),
            component = ifelse(component == "Zin", "Sentence",
                               ifelse(component %in% paste("Woorden", 1:3, sep = " "), "HF",
                                      ifelse(component == "Woorden 4", "LF",
                                             ifelse(component == "Snelheidstaak", "Tapping",
                                                    ifelse(component == "Medeklinkers", "Consonants", component))))),
            IKI = bigr_pausetime,
            Threshold = pauzegrens,
            target = bigr_is_targetted,
            logfile = logfile) %>% 
  filter(component %in% c("Tapping", "Consonants", "HF", "LF", "Sentence")) %>%
  filter(!(component == "HF" & freq %in% c("LF", "Indeterminate")),
         !(component == "LF" & freq %in% c("HF", "Indeterminate"))
  ) %>%
  select(-freq, -Threshold, -trial) -> d
str(d)
attr(d$bigram, "format.spss") <- NULL
attr(d$bigram, "display_width") <- NULL
attr(d$subj, "format.spss") <- NULL
attr(d$subj, "display_width") <- NULL
attr(d$IKI, "format.spss") <- NULL
attr(d$target, "format.spss") <- NULL
attr(d$logfile, "format.spss") <- NULL
attr(d$logfile, "display_width") <- NULL


#ppt$logfile = gsub(pattern = "ac14dutch2902201idfx", "ac14dutch29022016idfx",ppt$logfile)
#ppt$logfile = gsub(pattern = "admdutcht20180328idfx", "admjournalistdutch28032018idfx",ppt$logfile)

#removeL <- unique(ppt[which(!(ppt$logfile %in% d$logfile)),]$logfile)
#removeS <- unique(ppt[which(!(ppt$subj %in% d$subj)),]$subj)
#ppt %<>% filter(!(logfile %in% removeL))
#str(ppt)
#unique(d[which(!(d$logfile %in% ppt$logfile)),]$logfile)

files <- unique(d$logfile)
d$session <- NA
d$sex <- NA
d$age <- NA
for(f in files){
  if(f %in% unique(ppt$logfile)){
    d[d$logfile == f,]$session <- ppt[ppt$logfile == f,]$session[1]
    d[d$logfile == f,]$sex <- ppt[ppt$logfile == f,]$sex[1]
    d[d$logfile == f,]$age <- ppt[ppt$logfile == f,]$Age[1]
  }
}
d


d$subj <- as.integer(as.factor(d$subj))
length(unique(d$subj))
length(unique(d$bigram))

output <- paste(path_data, outputdf, sep = "/")
write_excel_csv(d, output)
