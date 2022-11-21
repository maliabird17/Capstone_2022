# VOLCANO PLOT #
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("plotly")
library(ggplot2)
library(gridExtra)
library(plotly)

volcano_plot <- function(final_df, human = TRUE) {
  #diff_df <- cbind(c(rownames(fitting2$p.value)), as.numeric(fitting2$F.p.value), as.numeric(fitting2$coefficients))
  if (human) {
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
    layout(title ="Volcano Plot") %>%
    layout(annotations = a)
  return(plot_v)
}


volcano_plot(gut_oldage_df)
