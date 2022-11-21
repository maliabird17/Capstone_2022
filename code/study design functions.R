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



#### FITTING MODEL ####
fit_model <- function(organ_) {
  
  organ_keys <- setNames(as.list(c('_g', '_t', '_m')), c('gut', 'testis', 'muscle'))
  
  organ_data <- all_organs %>% select(contains(organ_keys$organ_))
  organ_columns <- gsub("i", "", colnames(organ_data))
  colnames(organ_data) <- organ_columns
  
  organ_levels <- unique(colnames(organ_data))
  f <- factor(colnames(organ_data), levels = organ_levels)
  
  design <- model.matrix(~0 + f)
  colnames(design) <- organ_levels
  
  fit <- lmFit(organ_data, design)
  return(fit)
}

### TIME SERIES COMPARISONS ###
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
  colnames(sig_table) <- c("ProbeID", "Mid_FC", "Old_FC", "Mid_padj", "Old_padj")
  
  sig_table <- sig_table %>% filter(Mid_padj <= 0.05 | End_padj <= 0.05)
  sig_table <- sig_table[order(sig_table$End_padj), ]
  
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

