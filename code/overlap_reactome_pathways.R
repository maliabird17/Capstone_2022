library(VennDiagram)

# Function for finding overlap between significantly deregulated 
findOverlap <- function(fileList) {
  # Read the first file and extract the first column
  firstFileData <- read.table(fileList[1], sep = "\t", header = T)%>%
    filter(NOM.p.val <= 0.05 & FDR.q.val <= 0.25)
  firstFileFirstColumn <- firstFileData[, 1]
  
  # Initialize a vector to store the overlaps
  overlaps <- c()
  
  # Loop through the rest of the files
  for (i in 2:length(fileList)) {
    # Read the current file and extract the first column
    currentFileData <- read.table(fileList[i], sep = "\t", header = T)%>%
      filter(NOM.p.val <= 0.05 & FDR.q.val <= 0.25)
    currentFileFirstColumn <- currentFileData[, 1]
    
    # Find the overlaps between the first column of the current file and the first file
    currentOverlaps <- intersect(firstFileFirstColumn, currentFileFirstColumn)
    
    # Add the overlaps to the list of overlaps
    overlaps <- union(overlaps, currentOverlaps)
  }
  
  # Return the list of overlaps
  return(overlaps)
}


# Function for finding overlap across groups (WT vs TERT or between organs)
shared_chars <- function(...) {
  # Convert the input vectors to a list
  vecs <- list(...)
  
  # Find the intersection of the characters in all vectors
  intersection <- Reduce(intersect, vecs)
  
  # Return the resulting character vector
  return(intersection)
}


# Plotting venn diagram for overlap between groups
plot_shared_chars <- function(...) {
  
  grid.newpage()
  
  # Create a venn diagram with the specified vectors
  v <- venn.diagram(..., filename = NULL)
  
  # Return the plot
  grid.draw(v)
  
}



#### MUTANT-ORGAN OVERLAP ####
overlaps_WTmuscle.up <- findOverlap(c("/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_WT_muscle-3-36_upreg.tsv",
                          "/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_WT_muscle-3-24_upreg.tsv" ))


overlaps_WTmuscle.down <- findOverlap(c("/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_WT_muscle-3-36_downreg.tsv",
                                   "/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_WT_muscle-3-24_downreg.tsv"))

overlaps_WTtestes.up <- findOverlap(c("/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_WT_testes-3-24_upreg.tsv",
                           "/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_WT_testes-3-24_upreg.tsv"))

overlaps_WTtestes.down <- findOverlap(c("/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_WT_testes-3-24_downreg.tsv",
                                        "/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_WT_testes-3-24_downreg.tsv"))

overlaps_WTgut.up <- findOverlap(c("/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_WT_gut-3-36_upreg.tsv",
                                   "/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_WT_gut-3-24_upreg.tsv"))
                                      
overlaps_WTgut.down <- findOverlap(c("/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_WT_gut-3-24_downreg.tsv",
                                     "/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_WT_gut-3-36_downreg.tsv"))

         
                     
#### WILDTYPE OVERLAP ####
overlaps_WT.up <- shared_chars(overlaps_WTmuscle.up, overlaps_WTtestes.up, overlaps_WTgut.up)   
plot_shared_chars(list(muscle = overlaps_WTmuscle.up, testes = overlaps_WTtestes.up, gut = overlaps_WTgut.up))

overlaps_WT.down <- shared_chars(overlaps_WTmuscle.down, overlaps_WTtestes.down, overlaps_WTgut.down)
plot_shared_chars(list(muscle = overlaps_WTmuscle.down, testes = overlaps_WTtestes.down, gut = overlaps_WTgut.down))
                                 
   


#### MUTANT-ORGAN OVERLAP ####
overlaps_MUTtestes.up <- findOverlap(c("/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_MUT_testes-3-9_upreg.tsv",
                                       "/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_MUT_testes-3-9_upreg.tsv"))

overlaps_MUTtestes.down <- findOverlap(c("/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_MUT_testes-3-9_downreg.tsv",
                                         "/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_MUT_testes-3-6_downreg.tsv"))

overlaps_MUTtestes.down_JUSTLATTER <- findOverlap(c("/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_MUT_testes-3-9_downreg.tsv",
                                                    "/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_MUT_testes-3-9_downreg.tsv"))

overlaps_MUTgut.up <- findOverlap(c("/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_MUT_gut-3-9_upreg.tsv",
                                    "/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_MUT_gut-3-9_upreg.tsv"))

overlaps_MUTgut.down <- findOverlap(c("/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_MUT_gut-3-9_downreg.tsv",
                                      "/Users/malia/Desktop/Capstone_2022/gsea/gsea_results/gsea_reactome_MUT_gut-3-9_downreg.tsv"))


#### MUTANT OVERLAP ####
overlaps_MUT.up <- shared_chars(overlaps_MUTtestes.up, overlaps_MUTgut.up)   
plot_shared_chars(list(testes = overlaps_MUTtestes.up, gut = overlaps_MUTgut.up))

overlaps_MUT.down <- shared_chars(overlaps_MUTtestes.down_JUSTLATTER, overlaps_MUTgut.down)
plot_shared_chars(list(testes = overlaps_MUTtestes.down_JUSTLATTER, gut = overlaps_MUTgut.down))




#### ORGAN OVERLAP ####
testes.up <- shared_chars(overlaps_MUTtestes.up, overlaps_WTtestes.up)   
plot_shared_chars(list(MUT = overlaps_MUTtestes.up, WT = overlaps_WTtestes.up))

testes.down <- shared_chars(overlaps_MUTtestes.down_JUSTLATTER, overlaps_WTtestes.down)
plot_shared_chars(list(MUT = overlaps_MUTtestes.down_JUSTLATTER, WT = overlaps_WTtestes.down))

gut.up <- shared_chars(overlaps_MUTgut.up, overlaps_WTgut.up)   
plot_shared_chars(list(MUT = overlaps_MUTgut.up, WT = overlaps_WTgut.up))

gut.down <- shared_chars(overlaps_MUTgut.down, overlaps_WTgut.down)
plot_shared_chars(list(MUT = overlaps_MUTgut.down, WT = overlaps_WTgut.down))

                                                               