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

# Loading Data 
source("/Users/malia/Desktop/Capstone_2022/code/processing_all_data.R")

# Model Fitting -------------------------------------------

# Function to fit an organ's microarray data to model ## 

## INPUTS : ##
# organ_ : string corresponding to the organ you want to fit the model to (options: 'gut', 'muscle', or 'testis')
fit_model <- function(organ_) {
  
  # To select organ-specific data 
  organ_keys <- setNames(as.list(c('_g', '_t', '_m')), c('gut', 'testis', 'muscle'))
  key <- organ_keys[organ_]
  
  # Subset dataframe
  organ_data <- all_organs %>% select(contains(paste(key)))
  organ_columns <- gsub("i", "", colnames(organ_data)) # general across samples
  colnames(organ_data) <- organ_columns # name columns
  
  # specify levels to contrast (specific to each genotype, organ, and time)
  organ_levels <- unique(colnames(organ_data))
  f <- factor(colnames(organ_data), levels = organ_levels)
  
  # create matrix according to our levels
  design <- model.matrix(~0 + f)
  colnames(design) <- organ_levels
  
  # fit data to our model 
  fit <- lmFit(organ_data, design)
  fit$levels <- organ_levels
  return(fit)
  
}


## Function to identify differencially expressed genes between two time points ##
# Reference: Limma documentation, https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf 

## INPUTS : ##
# fit_object : output of fit_model function
# genotype : genotype being compared
# Time1 : first timepoint used in the comparison (should be the earlier time)
# Time2 : second timepoint used in the comparison (should be the later time)
# write_data : Boolean, by default FALSE to NOT write a txt file of differencially expressed genes to local directory specified
# Name : Name of data file to write, by default set to empty when no file is created (according to default above)
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

# Volcano Plots ------------------------------------------- 

# Volcano Plot function for interactive  visualization of DEGs yielded from the relative to 3 months comparison. 
# Modified from https://stevekm.github.io/2016/09/24/Plot.ly-Volcano-Plot.html

### INPUTS: ###
# final_df : output of timeseries_comparison_test 
# Name : string to be displayed as the plot title
# human : boolean, default set to true to display humanized gene names 
volcano_plot <- function(final_df, Name, human = TRUE) {
  
  # Mapping gene names 
  if (human == T) {
    gene <- as.character(final_df$Human.gene.name)
  } else {
    gene <- as.character(final_df$Gene_symbol)
  }
  p <- as.numeric(final_df$padj)
  FC <- as.numeric(final_df$FC)
  
  diff_df <- as.data.frame(cbind(gene, p, FC))
  diff_df <- diff_df %>% 
    dplyr::mutate_at(vars("p", "FC"), ~as.numeric(as.character(.))) 
  
  
  # add a grouping column; default value is "not significant"
  diff_df["group"] <- "NotSignificant"
  
  # for our plot, we want to highlight 
  # p-value < 0.05 (significance level)
  # Fold Change > 1.5
  
  # change the grouping for the entries with significance but not a large enough Fold change
  diff_df[which(diff_df['p'] < 0.05 & abs(diff_df['FC']) < 1.5 ),"group"] <- "Significant"
  
  # change the grouping for the entries a large enough Fold change but not a low enough p value
  diff_df[which(diff_df['p'] > 0.05 & abs(diff_df['FC']) > 1.5 ),"group"] <- "FoldChange"
  
  # change the grouping for the entries with both significance and large enough fold change
  diff_df[which(diff_df['p'] < 0.05 & abs(diff_df['FC']) > 1.5 ),"group"] <- "Significant&FoldChange"
  
  
  # Find and label the top peaks..
  top_peaks <- diff_df[with(diff_df, order(FC, p)),][1:5,]
  top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(-FC, p)),][1:5,])
  
  # Add gene labels for all of the top genes we found
  # here we are creating an empty list, and filling it with entries for each row in the dataframe
  # each list entry is another list with named items that will be used by Plot.ly
  a <- list()
  for (i in seq_len(nrow(top_peaks))) {
    m <- top_peaks[i, ]
    a[[i]] <- list(
      x = m[["FC"]],
      y = -log10(m[["p"]]),
      text = m[["gene"]],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0.5,
      ax = 20,
      ay = -40
    )
  }
  
  # make the Plot.ly plot
  plot_v <- plot_ly(data = diff_df, x = FC, y = -log10(p), text = gene, mode = "markers", color = diff_df$group) %>% 
    layout(title = Name) %>%
    layout(annotations = a)
  return(plot_v)
}



# Volcano Plot function for interactive  visualization of DEGs yielded from the alternative comparisons. 
# Modified from https://stevekm.github.io/2016/09/24/Plot.ly-Volcano-Plot.html

### INPUTS: ###
# final_df : output of timeseries_comparison 
# Name : string to be displayed as the plot title
# human : boolean, default set to true to display humanized gene names 
volcano_plot_alt <- function(final_df, Name, human = TRUE) {
  if (human == T) {
    gene <- as.character(final_df$Human.gene.name)
  } else {
    gene <- as.character(final_df$Gene_symbol)
  }
  p <- as.numeric(final_df$end_padj)
  FC <- as.numeric(final_df$Old_FC)
  
  diff_df <- as.data.frame(cbind(gene, p, FC))
  diff_df <- diff_df %>% 
    dplyr::mutate_at(vars("p", "FC"), ~as.numeric(as.character(.))) 
  
  
  # add a grouping column; default value is "not significant"
  diff_df["group"] <- "NotSignificant"
  
  # for our plot, we want to highlight 
  # p-value < 0.05 (significance level)
  # Fold Change > 1.5
  
  # change the grouping for the entries with significance but not a large enough Fold change
  diff_df[which(diff_df['p'] < 0.05 & abs(diff_df['FC']) < 1.5 ),"group"] <- "Significant"
  
  # change the grouping for the entries a large enough Fold change but not a low enough p value
  diff_df[which(diff_df['p'] > 0.05 & abs(diff_df['FC']) > 1.5 ),"group"] <- "FoldChange"
  
  # change the grouping for the entries with both significance and large enough fold change
  diff_df[which(diff_df['p'] < 0.05 & abs(diff_df['FC']) > 1.5 ),"group"] <- "Significant&FoldChange"
  
  
  # Find and label the top peaks..
  top_peaks <- diff_df[with(diff_df, order(FC, p)),][1:5,]
  top_peaks <- rbind(top_peaks, diff_df[with(diff_df, order(-FC, p)),][1:5,])
  
  # Add gene labels for all of the top genes we found
  # here we are creating an empty list, and filling it with entries for each row in the dataframe
  # each list entry is another list with named items that will be used by Plot.ly
  a <- list()
  for (i in seq_len(nrow(top_peaks))) {
    m <- top_peaks[i, ]
    a[[i]] <- list(
      x = m[["FC"]],
      y = -log10(m[["p"]]),
      text = m[["gene"]],
      xref = "x",
      yref = "y",
      showarrow = TRUE,
      arrowhead = 0.5,
      ax = 20,
      ay = -40
    )
  }
  
  # make the Plot.ly plot
  plot_v <- plot_ly(data = diff_df, x = FC, y = -log10(p), text = gene, mode = "markers", color = diff_df$group) %>% 
    layout(title = Name) %>%
    layout(annotations = a)
  return(plot_v)
}


