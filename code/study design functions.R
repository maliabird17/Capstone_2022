# Normalization & Annotation 
BiocManager::install("Biobase")
BiocManager::install("affy")
BiocManager::install("limma")

library(Biobase)
library(affy)
library(limma)

# Study Design 
# Function for defining groups depending on comparison
group_selection <- function(organ_, group1, time1, group2, time2) {
  subset1 <- mega %>% filter(organ == organ_ & genotype == group1 & age == time1)
  subset2 <- mega %>% filter(organ == organ_ & genotype == group2 & age == time2)
  subset <- rbind(subset1, subset2)
  
  return(subset)
}

# 9 MONTH COMPARISON #
testes_9m <- group_selection("testes", "mut", 9, "wt", 9)
groups <- testes_9m$genotype
### 

design <- model.matrix(~factor(groups)) # design experiment (set control vs treatment)
colnames(design) <- c(testes_9m$organ[1], "Mut vs WT") # rename columns

fit <- lmFit()


