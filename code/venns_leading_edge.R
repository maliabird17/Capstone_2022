# Muscle: AE5B83
# Gut: 679A7C,  
# Testes: 5267B4,  tert 505B81  wt 6BA5DB


my_venn <- function(...) {
  library(ggvenn)
  
  files <- list(...)
  
  # select first column of gene names
  columns <- lapply(files, function(file) {
    df <- read.delim(file, header = T)
    return(df[,1])
  })
  
  # name according to comparison
  input <- list(Gut = columns[[1]], Testes = columns[[2]])
  
  grid.newpage()
  # Create a venn diagram with the extracted columns
  v <- ggvenn(input, fill_color = c("#679A7C", "#5267B4"))
  grid.draw(v)

}




my_venn("gsea/leading_edge/ranked_gene_list_MUT_gut-3-9.tsv",
            "gsea/leading_edge/ranked_gene_list_MUT_testes-3-9.tsv")

my_venn("/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_gut-3-36.tsv",
            "/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_gut-3-24.tsv")

my_venn("/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_muscle-3-36.tsv",
               "/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_muscle-3-24.tsv")

my_venn("/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_gut-3-24.tsv",
        "/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_muscle-3-24.tsv",
        "/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_testes-3-24.tsv")

my_venn("/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_testes-3-24.tsv",
        "gsea/leading_edge/ranked_gene_list_MUT_testes-3-9.tsv")

my_venn( "/Users/malia/Desktop/Capstone_2022/gsea/leading_edge/ranked_gene_list_WT_gut-3-24.tsv",
         "gsea/leading_edge/ranked_gene_list_MUT_gut-3-9.tsv")
