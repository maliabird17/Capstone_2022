# WT AGEING 


#### GUT ###
fitted <- fit_model('gut') 
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


# Time series comparison (9.6.1 limma documentation)
timeseries_comparison <- function(fit_object,
                                  write_data = FALSE,
                                  Name = "Empty") {
  
  organ_levels <- fit_object$levels 
  
  contrast <- makeContrasts(
    "wt_6_g - wt_3_g",
    "wt_9_g - wt_6_g",
    "wt_24_g - wt_9_g",
    "wt_36_g - wt_24_g",
    levels = organ_levels
  )
  
  fit2 <- contrasts.fit(fitted, contrast)
  fit2 <- eBayes(fit2)
  table_pvals <- data.frame(fit2$p.value)
  top.table <- topTable(fit2, sort.by = "P", n = Inf)
  
  adjust_table <- cbind(p.adjust(table_pvals$wt_6_g...wt_3_g, method = "BH"), 
                        p.adjust(table_pvals$wt_9_g...wt_6_g, method = "BH"), 
                        p.adjust(table_pvals$wt_24_g...wt_9_g, method = "BH"), 
                        p.adjust(table_pvals$wt_36_g...wt_24_g, method = "BH"))
  colnames(adjust_table) <- c("one_padj", "two_padj", "three_padj", "four_padj")
  
  table_coefficients <- fit2$coefficients
  colnames(table_coefficients) <- c("FC1", "FC2", "FC3", "FC4")
  sig_table <- data.frame(cbind(table_coefficients, adjust_table))
  sig_table <- tibble::rownames_to_column(sig_table, "ProbeID")
  #colnames(sig_table) <- c("ProbeID", "Mid_FC", "Old_FC", "mid_padj", "end_padj")
  
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


# MUTANT Ageing 
gut_wt <- timeseries_comparison(fitted,
                                    TRUE,
                                    "Gut_WT-3-6-9-24-36")

timeseries_comparison <- function(fit_object,
                                  write_data = FALSE,
                                  Name = "Empty") {
  
  organ_levels <- fit_object$levels 
  
  contrast <- makeContrasts(
    "mut_6_g - mut_3_g",
    "mut_9_g - mut_6_g",
    levels = organ_levels
  )
  
  fit2 <- contrasts.fit(fitted, contrast)
  fit2 <- eBayes(fit2)
  table_pvals <- data.frame(fit2$p.value)
  
  adjust_table <- cbind(p.adjust(table_pvals$mut_6_g...mut_3_g, method = "BH"), 
                        p.adjust(table_pvals$mut_9_g...mut_6_g, method = "BH"))
  colnames(adjust_table) <- c("one_padj", "two_padj")
  
  table_coefficients <- fit2$coefficients
  colnames(table_coefficients) <- c("FC1", "FC2")
  sig_table <- data.frame(cbind(table_coefficients, adjust_table))
  sig_table <- tibble::rownames_to_column(sig_table, "ProbeID")
  #colnames(sig_table) <- c("ProbeID", "Mid_FC", "Old_FC", "mid_padj", "end_padj")
  
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

gut_mut <- timeseries_comparison(fitted,
                                TRUE,
                                "Gut_MUT-3-6-9")

#### TESTES #### 
fitted <- fit_model('testis')
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


# Time series comparison (9.6.1 limma documentation)
timeseries_comparison <- function(fit_object,
                                  write_data = FALSE,
                                  Name = "Empty") {
  
  organ_levels <- fit_object$levels 
  
  contrast <- makeContrasts(
    "wt_6_t - wt_3_t",
    "wt_9_t - wt_6_t",
    "wt_24_t - wt_9_t",
    "wt_36_t - wt_24_t",
    levels = organ_levels
  )
  
  fit2 <- contrasts.fit(fitted, contrast)
  fit2 <- eBayes(fit2)
  table_pvals <- data.frame(fit2$p.value)
  
  adjust_table <- cbind(p.adjust(table_pvals$wt_6_t...wt_3_t, method = "BH"), 
                        p.adjust(table_pvals$wt_9_t...wt_6_t, method = "BH"), 
                        p.adjust(table_pvals$wt_24_t...wt_9_t, method = "BH"), 
                        p.adjust(table_pvals$wt_36_t...wt_24_t, method = "BH"))
  colnames(adjust_table) <- c("one_padj", "two_padj", "three_padj", "four_padj")
  
  table_coefficients <- fit2$coefficients
  colnames(table_coefficients) <- c("FC1", "FC2", "FC3", "FC4")
  sig_table <- data.frame(cbind(table_coefficients, adjust_table))
  sig_table <- tibble::rownames_to_column(sig_table, "ProbeID")
  #colnames(sig_table) <- c("ProbeID", "Mid_FC", "Old_FC", "mid_padj", "end_padj")
  
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


# MUTANT Ageing 


testes_wt <- timeseries_comparison(fitted,
                                TRUE,
                                "Testes_WT-3-6-9-24-36")

timeseries_comparison <- function(fit_object,
                                  write_data = FALSE,
                                  Name = "Empty") {
  
  organ_levels <- fit_object$levels 
  
  contrast <- makeContrasts(
    "mut_6_t - mut_3_t",
    "mut_9_t - mut_6_t",
    levels = organ_levels
  )
  
  fit2 <- contrasts.fit(fitted, contrast)
  fit2 <- eBayes(fit2)
  table_pvals <- data.frame(fit2$p.value)
  
  adjust_table <- cbind(p.adjust(table_pvals$mut_6_t...mut_3_t, method = "BH"), 
                        p.adjust(table_pvals$mut_9_t...mut_6_t, method = "BH"))
  colnames(adjust_table) <- c("one_padj", "two_padj")
  
  table_coefficients <- fit2$coefficients
  colnames(table_coefficients) <- c("FC1", "FC2")
  sig_table <- data.frame(cbind(table_coefficients, adjust_table))
  sig_table <- tibble::rownames_to_column(sig_table, "ProbeID")
  #colnames(sig_table) <- c("ProbeID", "Mid_FC", "Old_FC", "mid_padj", "end_padj")
  
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

testes_mut <- timeseries_comparison(fitted,
                                 TRUE,
                                 "Testes_MUT-3-6-9")

#### MUSCLE ####
fitted <- fit_model('muscle')
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


# Time series comparison (9.6.1 limma documentation)
timeseries_comparison <- function(fit_object,
                                  write_data = FALSE,
                                  Name = "Empty") {
  
  organ_levels <- fit_object$levels 
  
  contrast <- makeContrasts(
    "wt_6_m - wt_3_m",
    "wt_9_m - wt_6_m",
    "wt_24_m - wt_9_m",
    "wt_36_m - wt_24_m",
    levels = organ_levels
  )
  
  fit2 <- contrasts.fit(fitted, contrast)
  fit2 <- eBayes(fit2)
  table_pvals <- data.frame(fit2$p.value)
  
  adjust_table <- cbind(p.adjust(table_pvals$wt_6_m...wt_3_m, method = "BH"), 
                        p.adjust(table_pvals$wt_9_m...wt_6_m, method = "BH"), 
                        p.adjust(table_pvals$wt_24_m...wt_9_m, method = "BH"), 
                        p.adjust(table_pvals$wt_36_m...wt_24_m, method = "BH"))
  colnames(adjust_table) <- c("one_padj", "two_padj", "three_padj", "four_padj")
  
  table_coefficients <- fit2$coefficients
  colnames(table_coefficients) <- c("FC1", "FC2", "FC3", "FC4")
  sig_table <- data.frame(cbind(table_coefficients, adjust_table))
  sig_table <- tibble::rownames_to_column(sig_table, "ProbeID")
  #colnames(sig_table) <- c("ProbeID", "Mid_FC", "Old_FC", "mid_padj", "end_padj")
  
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


# MUTANT Ageing 
muscle_wt <- timeseries_comparison(fitted,
                                TRUE,
                                "Muscle_WT-3-6-9-24-36")

timeseries_comparison <- function(fit_object,
                                  write_data = FALSE,
                                  Name = "Empty") {
  
  organ_levels <- fit_object$levels 
  
  contrast <- makeContrasts(
    "mut_6_m - mut_3_m",
    "mut_9_m - mut_6_m",
    levels = organ_levels
  )
  
  fit2 <- contrasts.fit(fitted, contrast)
  fit2 <- eBayes(fit2)
  table_pvals <- data.frame(fit2$p.value)
  
  adjust_table <- cbind(p.adjust(table_pvals$mut_6_m...mut_3_m, method = "BH"), 
                        p.adjust(table_pvals$mut_9_m...mut_6_m, method = "BH"))
  colnames(adjust_table) <- c("one_padj", "two_padj")
  
  table_coefficients <- fit2$coefficients
  colnames(table_coefficients) <- c("FC1", "FC2")
  sig_table <- data.frame(cbind(table_coefficients, adjust_table))
  sig_table <- tibble::rownames_to_column(sig_table, "ProbeID")
  #colnames(sig_table) <- c("ProbeID", "Mid_FC", "Old_FC", "mid_padj", "end_padj")
  
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

muscle_mut <- timeseries_comparison(fitted,
                                 TRUE,
                                 "Muscle_MUT-3-6-9")
