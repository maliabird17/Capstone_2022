# Loading and installing packages 
# ----  
# Installing required packages from BiocManager
#BiocManager::install("limma")
#BiocManager::install("Biobase")
#BiocManager::install("Matrix")
#BiocManager::install("affy")
#BiocManager::install("marray")

# Loading packages
library(limma)
library(Biobase)
library(Matrix)
library(affy)
library(gplots)
library(marray)
library(dplyr)
library(multtest)
library(readxl)


###########
# ----  
#Producing heat map of the whole data

########  ########
# Comparison1: 3m WT vs 36m WT testes #
design.1 <- design[which(design[,"TimeGenoOrgan"]=="WT3T" | design[,"TimeGenoOrgan"]=="WT36Talt"),]
dataset.1 <- dataset[,match(design.1[,"Geo_accession"], colnames(dataset))]
design.1 <- as.data.frame(design.1)
CM <- design.1$TimeGenoOrgan
CM

# mine
wt_3_36 <- group_selection("testes", "wt", "3", "wt", "36")

data_wt_3_36 <- all_organs %>% select(all_of(wt_3_36$labels))
wt3_36_labels <- wt_3_36$sample

design_matrix <- model.matrix(~0+ wt3_36_labels)
colnames(design_matrix) <- unique(wt_3_36$sample)

fitting <- lmFit(data_wt_3_36, design_matrix)
contrast <- makeContrasts(wt_36_t - wt_3_t, levels = design_matrix)
fitting2 <- contrasts.fit(fitting, contrast)
fitting2 <- eBayes(fitting2)

fitting2$p.value


# p-val adjust
# Mult-test correction and determining p-values #
fitting_adjusted <- p.adjust(fitting2$p.value, method = "hochberg")
sig_table <- cbind(rownames(dataset), fitting2$coefficients, fitting2$p.value, fitting_adjusted, fitting2$lods, fitting2$t)
colnames(sig_table) <- c("ProbeID", "Log2R", "Rawpval", "AdjpvalBH", "Lods", "t-stats")
write.table(sig_table, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", NAME, "_Bstats.txt", sep=""), sep="\t", row.names=F)


# Annotating diff expressed genes #
ANOT <- as.data.frame(read_excel("/Users/malia/Desktop/Capstone_2022/raw_data/ZebGene-1_1-st-v1.na33.3.zv9.ANOT7.xlsx")) %>%
  select(transcript_cluster_id, Gene_symbol, transcript_ID, gene_assignement, GO_term, GO_description)

merged_df <- merge(ANOT, sig_table, by.x = "transcript_cluster_id", by.y ="ProbeID")
#NAMEd = "annotated_testes_WT_3_36"
#write.table(merged_df, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", NAMEd, ".txt", sep=""), sep="\t", row.names=F)

# Human (using ensembl)
human_zeb_map <- read.csv("/Users/malia/Desktop/Capstone_2022/raw_data/ensembl_human_zebra_map.txt", sep = ",", head = T) %>%
  dplyr::select(Gene.stable.ID, Human.gene.name)
final_df <- merge(merged_df, human_zeb_map, by.x = "gene_assignement", by.y = "Gene.stable.ID")
NAMEe = "annotated_testes_WT_3_36"
write.table(merged_df, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", NAMEd, ".txt", sep=""), sep="\t", row.names=F)
