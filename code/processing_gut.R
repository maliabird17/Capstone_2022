library(tidyverse)

gut <- read_tsv("raw_data/Madalena_microarraysExpression_gut.csv", 
         col_types = cols(...1 = col_character(), 
                          .default = col_double())) %>% # setting column types 
                                                        # (ascension number to character & rest to double)
  rename_all(tolower) %>% # change all column names to lowercase
  rename(ID = ...1)

muscle <- read_tsv("raw_data/Madalena_microarraysExpression_muscle.csv", 
         col_types = cols(...1 = col_character(), 
                          .default = col_double())) %>% # setting column types 
  # (ascension number to character & rest to double)
  rename_all(tolower) %>% # change all column names to lowercase
  rename(ID = ...1)

testes <- read_tsv("raw_data/Madalena_microarraysExpression_testes.csv", 
         col_types = cols(...1 = col_character(), 
                          .default = col_double())) %>% # setting column types 
  # (ascension number to character & rest to double)
  rename_all(tolower) %>% # change all column names to lowercase
  rename(ID = ...1)


colnames(gut)
colnames(muscle)
colnames(testes)

exp <- read_tsv("raw_data/Madalena_microarraysExpression2.csv")
colnames(exp)
