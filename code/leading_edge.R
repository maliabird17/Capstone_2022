# This function takes a list of file names and returns a merged data frame
# with only the rows that have values in the first column found in all input
# files, and with the sum of the values in the third column for each row
merge_files <- function(file_names) {
  
  # First, we read all the files and store the data in a list
  data_list <- lapply(file_names, read.table, sep="\t", header=TRUE)
  
  # Then, we take the data in the first column
  # of each file and compute the overlap
  overlap <- Reduce(intersect, lapply(data_list, "[[", 1))
  
  # Next, we create an empty data frame to store the merged data
  merged_data <- data.frame()
  
  # Then, we iterate over the overlapping values and compute the sum of the
  # values in the third column for each value
  for (value in overlap) {
    # We extract the rows from each data frame that have the current value
    # in the first column
    rows <- lapply(data_list, function(x) x[x[,1] == value, ])
    
    # We compute the sum of the values in the third column for each row
    sums <- sapply(rows, function(x) sum(x[,3]))
    
    # We add the sums to the merged data frame
    merged_data <- rbind(merged_data, sums)
  }
  
  data_frame <- cbind(overlap, merged_data)
  
  # Take the absolute sum of their rankings (impact on pathway) to rank them
  data_frame['Sum'] <- rowSums(data_frame[, 2:ncol(data_frame)])
  data_frame['Abs_sum'] <- abs(data_frame$Sum)
  data_frame <- data_frame[order(data_frame[,ncol(data_frame)], decreasing=TRUE),]
  
  # Finally, we add the overlapping values as the first column of the merged
  # data frame and return the result
  return(data_frame[1:10,])
}


mutant <- c("gsea/leading_edge/ranked_gene_list_MUT_gut-3-9.tsv",
                "gsea/leading_edge/ranked_gene_list_MUT_testes-3-9.tsv")

wt_gut <- c("/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_gut-3-36.tsv",
            "/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_gut-3-24.tsv")

wt_muscle <- c("/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_muscle-3-36.tsv",
               "/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_muscle-3-24.tsv")

wt <- c("/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_gut-3-24.tsv",
        "/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_muscle-3-24.tsv",
        "/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_testes-3-24.tsv")

testes <- c("gsea/leading_edge/ranked_gene_list_MUT_testes-3-9.tsv",
            "/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_testes-3-24.tsv")

gut <- c("/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_gut-3-36.tsv",
         "/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_gut-3-24.tsv",
         "gsea/leading_edge/ranked_gene_list_MUT_gut-3-9.tsv")

g <- read.csv("gsea/leading_edge/ranked_gene_list_MUT_gut-3-9.tsv", sep = "\t")


top_genes_MUT <- merge_files(mutant)
top_genes_WT <- merge_files(wt)
top_genes_WTgut <- merge_files(wt_gut)
top_genes_WTmuscle <- merge_files(wt_muscle)

top_genes_TESTES <- merge_files(testes)
top_genes_GUT <- merge_files(gut)

top_genes_MUT
