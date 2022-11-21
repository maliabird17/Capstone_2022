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








fit_model <- function(organ_) {

    organ_keys <- setNames(as.list(c('_g', '_t', '_m')), c('gut', 'testis', 'muscle'))
    
    key <- organ_keys[organ_]
    organ_data <- all_organs %>% select(contains(paste(key)))
    organ_columns <- gsub("i", "", colnames(organ_data))
    colnames(organ_data) <- organ_columns
    
    organ_levels <- unique(colnames(organ_data))
    f <- factor(colnames(organ_data), levels = organ_levels)
    
    design <- model.matrix(~0 + f)
    colnames(design) <- organ_levels

    fit <- lmFit(organ_data, design)
    fit$levels <- organ_levels
    return(fit)
    
}



fitted <- fit_model('gut')


# Time series comparison (9.6.1 limma documentation)
timeseries_comparison <- function(fit_object, 
                                  baseline_wt,
                                  mid_wt,
                                  end_wt,
                                  baseline_mut,
                                  mid_mut,
                                  end_mut,
                                  write_data = FALSE,
                                  Name = "Empty") {
  
  organ_levels <- fit_object$levels 
  
  contrast <- makeContrasts(
    Diff_Mid = (mut_6_g - mut_3_g) - (wt_6_g - wt_3_g),
    Diff_End = (mut_9_g - mut_6_g) - (wt_36_g - wt_6_g),
    levels = organ_levels
  )
  
  fit2 <- contrasts.fit(fit_object, contrast)
  fit2 <- eBayes(fit2)
  table_pvals <- data.frame(fit2$p.value)
  
  adjust_table <- cbind(p.adjust(table_pvals$Diff_Mid, method = "BH"), 
                        p.adjust(table_pvals$Diff_End, method = "BH"))
  colnames(adjust_table) <- c("mid_padj", "end_padj")
  
  table_coefficients <- fit2$coefficients
  sig_table <- data.frame(cbind(table_coefficients, adjust_table))
  sig_table <- tibble::rownames_to_column(sig_table, "ProbeID")
  colnames(sig_table) <- c("ProbeID", "Mid_FC", "Old_FC", "mid_padj", "end_padj")
  
  ANOT <- as.data.frame(read_excel("/Users/malia/Desktop/Capstone_2022/raw_data/ZebGene-1_1-st-v1.na33.3.zv9.ANOT7.xlsx")) %>%
    select(transcript_cluster_id, Gene_symbol, transcript_ID, gene_assignement, GO_term, GO_description)
  
  merged_df <- merge(ANOT, sig_table, by.x = "transcript_cluster_id", by.y ="ProbeID")
  
  human_zeb_map <- read.csv("/Users/malia/Desktop/Capstone_2022/raw_data/ensembl_human_zebra_map.txt", sep = ",", head = T) %>%
    dplyr::select(Gene.stable.ID, Human.gene.name)
  final_df <- merge(merged_df, human_zeb_map, by.x = "gene_assignement", by.y = "Gene.stable.ID")
  
  if (write_data == T) {
    write.table(final_df, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", Name, ".txt", sep=""), sep="\t", row.names=F)
  }
  
  return(final_df)
}


gut_oldage_df <- timeseries_comparison(fitted,
                      'wt_3_g', 'wt_6_g', 'wt_36_g',
                      'mut_3_g', 'mut_6_g', 'mut_9_g',
                      TRUE,
                      "Gut_MUT-3-6-9_v_WT-3-6-36")



## GUT
######
gut_levels <- c("wt_3_g",  "wt_6_g", "wt_9_g", "wt_24_g", "wt_36_g", 
                "mut_3_g", "mut_6_g", "mut_9_g")
gut_data <- all_organs %>% select(contains('_g'))
gut_columns <- gsub("i", "", colnames(gut_data))

colnames(gut_data) <- gut_columns

gut_levels <- unique(colnames(gut_data))

f <- factor(colnames(gut_data), levels = gut_levels)
design <- model.matrix(~0 + f)
colnames(design) <- gut_levels

fit <- lmFit(gut_data, design)
######


### 3 -> 6 -> 9
gut_contrast_consistent <- makeContrasts(
  Diff_Mid = (wt_6_g - wt_3_g) - (mut_6_g - mut_3_g),
  Diff_M_Death = (wt_9_g - wt_6_g) - (mut_9_g - mut_6_g),
  levels = design
)

fit2_consistent_gut <- contrasts.fit(fit, gut_contrast_consistent)
fit2_consistent_gut <- eBayes(fit2_consistent_gut)
topTable(fit2_consistent_gut, adjust = "BH")


### 3 -> 6 -> 9/36
gut_contrast_oldage <- makeContrasts(
  Diff_Mid = (mut_6_g - mut_3_g) - (wt_6_g - wt_3_g),
  Diff_Old = (mut_9_g - mut_6_g) - (wt_36_g - wt_6_g),
  levels = design
)

fit2_oldage_gut <- contrasts.fit(fit, gut_contrast_oldage)
fit2_oldage_gut <- eBayes(fit2_oldage_gut)
table_pvals <- data.frame(fit2_oldage_gut$p.value)

adjust_table <- cbind(p.adjust(table_pvals$Diff_Mid, method = "BH"), 
                      p.adjust(table_pvals$Diff_Old, method = "BH"))
colnames(adjust_table) <- c("mid_padj", "old_padj")

table_gut_oa <- fit2_oldage_gut$coefficients
sig_gut_oldage_table <- data.frame(cbind(table_gut_oa, adjust_table))
sig_gut_oldage_table <- tibble::rownames_to_column(sig_gut_oldage_table, "ProbeID")
colnames(sig_gut_oldage_table) <- c("ProbeID", "Mid_FC", "Old_FC", "Mid_padj", "Old_padj")
#topTable(fit2_oldage_gut, adjust = "BH")


# Annotating diff expressed genes #
ANOT <- as.data.frame(read_excel("/Users/malia/Desktop/Capstone_2022/raw_data/ZebGene-1_1-st-v1.na33.3.zv9.ANOT7.xlsx")) %>%
  select(transcript_cluster_id, Gene_symbol, transcript_ID, gene_assignement, GO_term, GO_description)

merged_df <- merge(ANOT, sig_gut_oldage_table, by.x = "transcript_cluster_id", by.y ="ProbeID")
#NAMEd = "annotated_testes_WT_3_36"
#write.table(merged_df, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", NAMEd, ".txt", sep=""), sep="\t", row.names=F)

# Human (using ensembl)
human_zeb_map <- read.csv("/Users/malia/Desktop/Capstone_2022/raw_data/ensembl_human_zebra_map.txt", sep = ",", head = T) %>%
  dplyr::select(Gene.stable.ID, Human.gene.name)
final_df <- merge(merged_df, human_zeb_map, by.x = "gene_assignement", by.y = "Gene.stable.ID")



# ---- 


pairwise_DEGs <- function(organ, geno1, t1, geno2, t2) {
  # Creating name of file to export
  name <- paste("DEG", organ, paste(geno1, t1, sep =""), paste(geno2, t2, sep =""), sep = "_")
  
  # retrieve information of chosen groups
  ID_subset <- group_selection(organ, geno1, t1, geno2, t2)
  groups <- ID_subset$sample
  
  # Get microarray data of specified  
  microarray_data <- all_organs %>% select(all_of(ID_subset$labels))
  
  # Create design matrix (matrix defining which samples belong to each group)
  design_matrix <- model.matrix(~0 + groups)
  colnames(design_matrix) <- unique(groups)
  
  # Fitting the data to linear 
  
}

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
