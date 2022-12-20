source("/Users/malia/Desktop/Capstone_2022/code/function_hub.R")

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



############ MUSCLE ################
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
    Diff_Mid = (mut_6_m - mut_3_m) - (wt_6_m - wt_3_m),
    Diff_End = (mut_9_m - mut_6_m) - (wt_36_m - wt_6_m),
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


muscle_final_df <- timeseries_comparison(fitted,
                      'wt_3_m', 'wt_6_m', 'wt_36_m',
                      'mut_3_m', 'mut_6_m', 'mut_9_m',
                      TRUE,
                      "Muscle_MUT-3-6-9_v_WT-3-6-36")



############ TESTES ################
fitted <- fit_model('testis')


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
    Diff_Mid = (mut_6_t - mut_3_t) - (wt_6_t - wt_3_t),
    Diff_End = (mut_9_t - mut_6_t) - (wt_36_t - wt_6_t),
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


testes_final_df <- timeseries_comparison(fitted,
                      'wt_3_t', 'wt_6_t', 'wt_36_t',
                      'mut_3_t', 'mut_6_t', 'mut_9_t',
                      TRUE,
                      "Testes_MUT-3-6-9_v_WT-3-6-36")


############ GUT ################
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
  
  # inherit levels from the fit model 
  organ_levels <- fit_object$levels 
  
  # 
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

gut_final_df <- timeseries_comparison(fitted,
                                       'wt_3_g', 'wt_6_g', 'wt_36_g',
                                       'mut_3_g', 'mut_6_g', 'mut_9_g',
                                       TRUE,
                                       "Gut_MUT-3-6-9_v_WT-3-6-36")
