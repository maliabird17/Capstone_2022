## Loading Packages
# ----  
library(tidyverse)


# Loading Normalized Gene Data
# ----  

### ALL NORMALIZED DATA ###
all_organs <- read.csv("raw_data/Madalena_microarraysExpression2.csv", row.names = 1, sep = "\t") %>% 
  rename_all(tolower)  # change all column names to lowercase


### INDIVIDUAL ORGANS: NORMALIZED DATA ### 
gut <- read_tsv("raw_data/Madalena_microarraysExpression_gut.csv", row.names = 1, sep = "\t") %>% 
  rename_all(tolower) # change all column names to lowercase

muscle <- read_tsv("raw_data/Madalena_microarraysExpression_muscle.csv", row.names = 1, sep = "\t") %>% 
  rename_all(tolower) # change all column names to lowercase


testes <- read_tsv("raw_data/Madalena_microarraysExpression_testes.csv") %>% 
  rename_all(tolower) # change all column names to lowercase



# Sample Reference Table
# ---- 

### SAMPLE REFERENCE TABLE ###
labels <- colnames(all_organs) # extracting vector of labels

mega <- as.data.frame(labels) %>%
  # individual columns for genotype, age, and sample number
  separate(labels, c("genotype", "age", "replicate"), sep = "_", remove = FALSE) %>% 
  # remove sample number from sample types column
  mutate(sample = str_replace_all(labels, c("i" = "", "v" = ""))) %>%
  # extend organ values 
  mutate(organ = case_when(
    startsWith(replicate, "t") ~ "testes",
    startsWith(replicate, "m") ~ "muscle",
    startsWith(replicate, "g") ~ "gut")) %>%
  # replace roman numerals with numbers
  mutate(replicate = str_replace_all(replicate, "iii", "3")) %>%
  mutate(replicate = str_replace_all(replicate, "ii", "2")) %>%
  mutate(replicate = str_replace_all(replicate, "iv", "4")) %>%
  mutate(replicate = str_replace_all(replicate, "vi", "6")) %>%
  mutate(replicate = str_replace_all(replicate, "v", "5")) %>%
  mutate(replicate = str_replace_all(replicate, "i", "1")) %>%
  mutate(replicate = str_replace_all(replicate, c("g" = "", "t" = "", "m" = "")))











