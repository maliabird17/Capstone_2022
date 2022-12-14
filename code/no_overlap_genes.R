merge_files_difference <- function(file_names) {
  
  # First, we read all the files and store the data in a list
  data_list <- lapply(file_names, read.table, sep="\t", header=TRUE)
  
  # Then, we make pairwise comparisons between the data in the first column
  # of each file and compute the overlap
  nooverlap <- Reduce(setdiff, lapply(data_list, "[[", 1))
  
  # Next, we create an empty data frame to store the merged data
  merged_data <- data.frame()
  
  # Then, we iterate over the overlapping values and compute the sum of the
  # values in the third column for each value
  for (value in nooverlap) {
    # We extract the rows from each data frame that have the current value
    # in the first column
    rows <- lapply(data_list, function(x) x[x[,1] == value, ])
    
    # We compute the sum of the values in the third column for each row
    sums <- sapply(rows, function(x) sum(x[,3]))
    
    # We add the sums to the merged data frame
    merged_data <- rbind(merged_data, sums)
  }
  
  data_frame <- cbind(nooverlap, merged_data)
  data_frame['Sum'] <- rowSums(data_frame[, 2:ncol(data_frame)])
  data_frame <- data_frame[order(data_frame[,ncol(data_frame)], decreasing=TRUE),]
  
  # Finally, we add the overlapping values as the first column of the merged
  # data frame and return the result
  return(data_frame[, 1])
}

justmut_genes_TESTES <- merge_files_difference(testes)
justmut_genes_GUT <- merge_files_difference(gut)
