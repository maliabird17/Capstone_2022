##### DEG Tables RELATIVE TO 3 MONTH TIMEPOINT #####
fitted_muscle <- fit_model('muscle')
timeseries_comparison_test <- function(fit_object,
                                       genotype,
                                       organ,
                                       Time1,
                                       Time2,
                                       write_data = FALSE,
                                       Name = "Empty") {
  
  organ_levels <- fit_object$levels 
  
  contrast_string1 <- "G_T2_O - G_T1_O"
  contrast_string2 <- gsub("T1", Time1, contrast_string1)
  contrast_string3 <- gsub("T2", Time2, contrast_string2) 
  contrast_string4 <- gsub("G", genotype, contrast_string3)
  contrast_string5 <- gsub("O", organ, contrast_string4)
  
  
  contrast <- makeContrasts(contrasts=contrast_string5, levels=organ_levels)
  
  fit1 <- contrasts.fit(fit_object, contrast)
  fit2 <- eBayes(fit1)
  
  
  table_pvals <- data.frame(fit2$p.value)
  adjust_table <- p.adjust(table_pvals[,1], method = "BH")
  
  table_coefficients <- fit2$coefficients
  sig_table <- data.frame(cbind(table_coefficients, adjust_table))
  sig_table <- tibble::rownames_to_column(sig_table, "ProbeID")
  colnames(sig_table) <- c("ProbeID", "FC", "padj")
  
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
########## MUSCLE ##################################

muscle_wt_36 <- timeseries_comparison_test(fitted_muscle, "wt", "m",
                                        "3", "6",
                                        TRUE,
                                        "Muscle_WT-3-6")

muscle_wt_39 <- timeseries_comparison_test(fitted_muscle, "wt", "m",
                                        "3", "9",
                                        TRUE,
                                        "Muscle_WT-3-9")

muscle_wt_324 <- timeseries_comparison_test(fitted_muscle, "wt", "m",
                                         "3", "24",
                                         TRUE,
                                         "Muscle_WT-3-24")


muscle_wt_336 <- timeseries_comparison_test(fitted_muscle, "wt", "m",
                                         "3", "36",
                                         TRUE,
                                         "Muscle_WT-3-36")

muscle_mut_36 <- timeseries_comparison_test(fitted_muscle, "mut", "m",
                                         "3", "6",
                                         TRUE,
                                         "Muscle_MUT-3-6")

muscle_mut_39 <- timeseries_comparison_test(fitted_muscle, "mut", "m",
                                         "3", "9",
                                         TRUE,
                                         "Muscle_MUT-3-9")

