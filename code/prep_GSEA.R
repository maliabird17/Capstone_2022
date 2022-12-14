# PREPARING FOR GSEA - Fish
library(janitor)

df <- read.csv("processed_data/Gut_MUT-3-6-9_v_WT-3-6-36.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, Old_FC, end_padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(end_padj <= 0.05) %>%
  select(Human.gene.name, Old_FC)

write.table(df, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "gut_wt_mut_oldage_fc_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


df_muscle <- read.csv("processed_data/Muscle_MUT-3-6-9_v_WT-3-6-36.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, Old_FC, end_padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(end_padj <= 0.05) %>%
  select(Human.gene.name, Old_FC)

write.table(df_muscle, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "muscle_wt_mut_oldage_fc_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


df_testes <- read.csv("processed_data/Testes_MUT-3-6-9_v_WT-3-6-36.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, Old_FC, end_padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(end_padj <= 0.05) %>%
  select(Human.gene.name, Old_FC)

write.table(df_testes, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "testes_wt_mut_oldage_fc_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


#### AGEING SEPERATE GENOTYPES ###

# GUT #
df_wtgut <- read.csv("processed_data/Gut_WT-3-6-9-24-36.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC4, three_padj, four_padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(four_padj <= 0.05 & three_padj <= 0.05) %>%
  select(Human.gene.name, FC4)

write.table(df_wtgut, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "gut_wt_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_mutgut_young <- read.csv("processed_data/Gut_MUT-3-6-9.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC1, one_padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(one_padj <= 0.05) %>%
  select(Human.gene.name, FC1)

write.table(df_mutgut_young, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "gut_mut_young_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_mutgut_old <- read.csv("processed_data/Gut_MUT-3-6-9.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC2, two_padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(two_padj <= 0.05) %>%
  select(Human.gene.name, FC2)

write.table(df_mutgut_old, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "gut_mut_old_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)







# TESTES
df_wtmuscle <- read.csv("processed_data/Muscle_WT-3-6-9-24-36.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC3, three_padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(three_padj <= 0.05) %>%
  select(Human.gene.name, FC3)

write.table(df_wttestes, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "testes_wt_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_muttestes_young <- read.csv("processed_data/Testes_MUT-3-6-9.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC1, one_padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(one_padj <= 0.05) %>%
  select(Human.gene.name, FC1)

write.table(df_muttestes_young, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "testes_mut_young_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_muttestes_old <- read.csv("processed_data/Testes_MUT-3-6-9.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC2, two_padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(two_padj <= 0.05) %>%
  select(Human.gene.name, FC2)

write.table(df_muttestes_old, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "testes_mut_old_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# MUSCLE
df_wtmuscle <- read.csv("processed_data/Muscle_WT-3-6-9-24-36.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC4,three_padj, four_padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(four_padj <= 0.05 & three_padj <= 0.05) %>%
  select(Human.gene.name, FC4)

write.table(df_wtmuscle, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "muscle_wt_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_mutmuscle_young <- read.csv("processed_data/Muscle_MUT-3-6-9.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC1, one_padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(one_padj <= 0.05) %>%
  select(Human.gene.name, FC1)

write.table(df_mutmuscle_young, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "muscle_mut_young_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_mutmuscle_old <- read.csv("processed_data/Muscle_MUT-3-6-9.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC2, two_padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(two_padj <= 0.05) %>%
  select(Human.gene.name, FC2)


nrow(df_gut36 %>% filter(FC > 0))


###### WILDTYPES

WTposneg_gut <- c(nrow(df_gut36 %>% filter(FC > 0)),
                  - nrow(df_gut36 %>% filter(FC < 0)),
                  nrow(df_gut69 %>% filter(FC > 0)),
                  - nrow(df_gut69 %>% filter(FC < 0)),
                  nrow(df_gut924 %>% filter(FC > 0)),
                  - nrow(df_gut924 %>% filter(FC < 0)),
                  nrow(df_gut2436 %>% filter(FC > 0)),
                  - nrow(df_gut2436 %>% filter(FC < 0)))

WTposneg_testes <- c(nrow(df_testes36 %>% filter(FC > 0)),
                  - nrow(df_testes36 %>% filter(FC < 0)),
                  nrow(df_testes69 %>% filter(FC > 0)),
                  - nrow(df_testes69 %>% filter(FC < 0)),
                  nrow(df_testes924 %>% filter(FC > 0)),
                  - nrow(df_testes924 %>% filter(FC < 0)),
                  nrow(df_testes2436 %>% filter(FC > 0)),
                  - nrow(df_testes2436 %>% filter(FC < 0)))

WTposneg_muscle <- c(nrow(df_muscle36 %>% filter(FC > 0)),
                  - nrow(df_muscle36 %>% filter(FC < 0)),
                  nrow(df_muscle69 %>% filter(FC > 0)),
                  - nrow(df_muscle69 %>% filter(FC < 0)),
                  nrow(df_muscle924 %>% filter(FC > 0)),
                  - nrow(df_muscle924 %>% filter(FC < 0)),
                  nrow(df_muscle2436 %>% filter(FC > 0)),
                  - nrow(df_muscle2436 %>% filter(FC < 0)))



###### MUTANTS

MUTposneg_69gut <- c(nrow(df_mutgut_young %>% filter(df_mutgut_young[2] > 0)),
                      - nrow(df_mutgut_young %>% filter(df_mutgut_young[2] < 0)),
                      nrow(df_mutgut_old %>% filter(df_mutgut_old[2] > 0)),
                      - nrow(df_mutgut_old %>% filter(df_mutgut_old[2] < 0)))

MUTposneg_69testes <- c(nrow(df_muttestes_young %>% filter(df_muttestes_young[2] > 0)),
                   - nrow(df_muttestes_young %>% filter(df_muttestes_young[2] < 0)),
                   nrow(df_muttestes_old %>% filter(df_muttestes_old[2] > 0)),
                   - nrow(df_muttestes_old %>% filter(df_muttestes_old[2] < 0)))

MUTposneg_69muscle <- c(nrow(df_mutmuscle_young %>% filter(df_mutmuscle_young[2] > 0)),
                       - nrow(df_mutmuscle_young %>% filter(df_mutmuscle_young[2] < 0)),
                      nrow(df_mutmuscle_old %>% filter(df_mutmuscle_old[2] > 0)),
                      - nrow(df_mutmuscle_old %>% filter(df_mutmuscle_old[2] < 0)))


write.table(df_mutmuscle_old, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "muscle_mut_old_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




##### INDIVIDUAL GUT DEGS ####
df_gut36_mut <- read.csv("processed_data/Gut_MUT-3-6.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_gut36_mut, paste("/Users/malia/Desktop/Capstone_2022/processed_data/",
                                "gut_mut-3-6_sig", ".rnk", sep=""), sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

df_gut39_mut <- read.csv("processed_data/Gut_MUT-3-9.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_gut39_mut, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "gut_mut-3-9_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



df_gut36 <- read.csv("processed_data/Gut_WT-3-6.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_gut36, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "gut_wt-3-6_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_gut69 <- read.csv("processed_data/Gut_WT-3-9.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_gut69, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "gut_wt-3-9_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




df_gut924 <- read.csv("processed_data/Gut_WT-3-24.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_gut924, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "gut_wt-3-24_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_gut2436 <- read.csv("processed_data/Gut_WT-3-36.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_gut2436, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "gut_wt-3-36_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


##### INDIVIDUAL TESTES DEGS ####
df_testes36_mut <- read.csv("processed_data/Testes_MUT-3-6.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_testes36_mut, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "testes_mut-3-6_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_testes39_mut <- read.csv("processed_data/Testes_MUT-3-9.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_testes39_mut, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "testes_mut-3-9_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




df_testes36 <- read.csv("processed_data/Testes_WT-3-6.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_testes36, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "testes_wt-3-6_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_testes69 <- read.csv("processed_data/Testes_WT-3-9.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_testes69, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "testes_wt-3-9_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




df_testes924 <- read.csv("processed_data/Testes_WT-3-24.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_testes924, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "testes_wt-3-24_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_testes2436 <- read.csv("processed_data/Testes_WT-3-36.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_testes2436, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "testes_wt-3-36_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


##### INDIVIDUAL MUSCLE DEGS ####

df_muscle36_mut <- read.csv("processed_data/Muscle_MUT-3-6.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_muscle36_mut, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "muscle_mut-3-6_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_muscle39_mut <- read.csv("processed_data/Muscle_MUT-3-9.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_muscle39_mut, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "muscle_mut-3-9_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



df_muscle36 <- read.csv("processed_data/Muscle_WT-3-6.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_muscle36, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "muscle_wt-3-6_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_muscle69 <- read.csv("processed_data/Muscle_WT-3-9.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_muscle69, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "muscle_wt-3-9_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




df_muscle924 <- read.csv("processed_data/Muscle_WT-3-24.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_muscle924, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "muscle_wt-3-24_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

df_muscle2436 <- read.csv("processed_data/Muscle_WT-3-36.TXT", sep = "\t", na.strings=c("","NA")) %>% 
  select(Human.gene.name, FC, padj) %>%
  na.omit() %>%
  distinct() %>%
  filter(padj <= 0.05) %>%
  select(Human.gene.name, FC)

write.table(df_muscle2436, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", "muscle_wt-3-36_sig", ".rnk", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
