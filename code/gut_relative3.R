library(limma)
library(Biobase)
library(Matrix)
library(affy)
library(gplots)
library(marray)
library(dplyr)
library(multtest)
library(readxl)


##### DEG Tables RELATIVE TO 3 MONTH TIMEPOINT #####
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


fitted_gut <- fit_model('gut')


## Function to identify differencially expressed genes between two time points
timeseries_comparison_test <- function(fit_object,
                                       genotype,
                                       organ,
                                       Time1,
                                       Time2,
                                       write_data = FALSE,
                                       Name = "Empty") {
  
  #inherit levels from the fit model
  organ_levels <- fit_object$levels 
  
  # create string based on function's input parameters 
  contrast_string1 <- "G_T2_O - G_T1_O"
  contrast_string2 <- gsub("T1", Time1, contrast_string1)
  contrast_string3 <- gsub("T2", Time2, contrast_string2) 
  contrast_string4 <- gsub("G", genotype, contrast_string3)
  contrast_string5 <- gsub("O", organ, contrast_string4)
  
  
  # compare across selected levels
  contrast <- makeContrasts(contrasts=contrast_string5, levels=organ_levels)
  fit1 <- contrasts.fit(fit_object, contrast)
  fit2 <- eBayes(fit1)
  
  # adjust for multiple hypothesis testing and save p values
  table_pvals <- data.frame(fit2$p.value)
  adjust_table <- p.adjust(table_pvals[,1], method = "BH")
  
  # Select coefficients: fold change of each gene 
  table_coefficients <- fit2$coefficients
  
  # Create table with ProbeID, the fold change and adjusted p-value
  sig_table <- data.frame(cbind(table_coefficients, adjust_table))
  sig_table <- tibble::rownames_to_column(sig_table, "ProbeID")
  colnames(sig_table) <- c("ProbeID", "FC", "padj")
  
  # Annotate the probes with gene symbols and additional information on the gene 
  ANOT <- as.data.frame(read_excel("/Users/malia/Desktop/Capstone_2022/raw_data/ZebGene-1_1-st-v1.na33.3.zv9.ANOT7.xlsx")) %>%
    select(transcript_cluster_id, Gene_symbol, transcript_ID, gene_assignement, GO_term, GO_description)
  merged_df <- merge(ANOT, sig_table, by.x = "transcript_cluster_id", by.y ="ProbeID")
  
  # Map zebrafish genes to human ortholog 
  human_zeb_map <- read.csv("/Users/malia/Desktop/Capstone_2022/raw_data/ensembl_human_zebra_map.txt", sep = ",", head = T) %>%
    dplyr::select(Gene.stable.ID, Human.gene.name)
  final_df <- merge(merged_df, human_zeb_map, by.x = "gene_assignement", by.y = "Gene.stable.ID")
  
  # Save file of human and zebrafish gene information 
  if (write_data == T) {
    write.table(final_df, paste("/Users/malia/Desktop/Capstone_2022/processed_data/", Name, ".txt", sep=""), sep="\t", row.names=F)
  }
  
  return(final_df)
}

########## GUT #####################################


gut_wt_36 <- timeseries_comparison_test(fitted_gut, "wt", "g",
                                         "3", "6",
                                         TRUE,
                                         "Gut_WT-3-6")

gut_wt_39 <- timeseries_comparison_test(fitted_gut, "wt", "g",
                                        "3", "9",
                                        TRUE,
                                        "Gut_WT-3-9")

gut_wt_324 <- timeseries_comparison_test(fitted_gut, "wt", "g",
                                        "3", "24",
                                        TRUE,
                                        "Gut_WT-3-24")


gut_wt_336 <- timeseries_comparison_test(fitted_gut, "wt", "g",
                                          "3", "36",
                                          TRUE,
                                          "Gut_WT-3-36")

gut_mut_36 <- timeseries_comparison_test(fitted_gut, "mut", "g",
                                        "3", "6",
                                        TRUE,
                                        "Gut_MUT-3-6")

gut_mut_39 <- timeseries_comparison_test(fitted_gut, "mut", "g",
                                        "3", "9",
                                        TRUE,
                                        "Gut_MUT-3-9")

