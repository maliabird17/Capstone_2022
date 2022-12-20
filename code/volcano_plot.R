# VOLCANO PLOT #
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("plotly")
library(ggplot2)
library(gridExtra)
library(plotly)
source("/Users/malia/Desktop/Capstone_2022/code/function_hub.R")



######## Interactive Volcano Plots to explore DEGs ########

###### GUT ######
source("/Users/malia/Desktop/Capstone_2022/code/gut_relative3.R")

volcano_plot(gut_wt_36, "Volcano Plot: WT (3/6)")
volcano_plot(gut_wt_39, "Volcano Plot: WT (3/9)")
volcano_plot(gut_wt_324, "Volcano Plot: WT (3/24)")
volcano_plot(gut_wt_336, "Volcano Plot: WT (3/36)")

volcano_plot(gut_mut_36, "Volcano Plot: TERT (3/6)")
volcano_plot(gut_mut_39, "Volcano Plot: TERT (3/9)")


###### TESTES ######
source("/Users/malia/Desktop/Capstone_2022/code/testes_relative3.R")

volcano_plot3(testes_wt_36, "Volcano Plot: WT (3/6)")
volcano_plot3(testes_wt_39, "Volcano Plot: WT (3/9)")
volcano_plot3(testes_wt_324, "Volcano Plot: WT (3/24)")
volcano_plot3(testes_wt_336, "Volcano Plot: WT (3/36)")

volcano_plot3(testes_mut_36, "Volcano Plot: TERT (3/6)")
volcano_plot3(testes_mut_39, "Volcano Plot: TERT (3/9)")


###### MUSCLE ######
source("/Users/malia/Desktop/Capstone_2022/code/muscle_relative3.R")

volcano_plot3(muscle_wt_36, "Volcano Plot: WT (3/6)")
volcano_plot3(muscle_wt_39, "Volcano Plot: WT (3/9)")
volcano_plot3(muscle_wt_324, "Volcano Plot: WT (3/24)")
volcano_plot3(muscle_wt_336, "Volcano Plot: WT (3/36)")

volcano_plot3(testes_mut_36, "Volcano Plot: TERT (3/6)")
volcano_plot3(testes_mut_39, "Volcano Plot: TERT (3/9)")



###### ALTERNATIVE comparisons #####

##### Ageing across 3, 6, and 9 month timepoints comparison ####
source("/Users/malia/Desktop/Capstone_2022/code/alt-DEG_ageing-onset-369-3636.R")

volcano_plot_alt(gut_final_df, "Volcano Plot: TERT -/- (3/9) vs WT (3/36)")
volcano_plot_alt(testes_final_df, "Volcano Plot: TERT -/- (3/9) vs WT (3/36)")
volcano_plot_alt(muscle_final_df, "Volcano Plot: TERT -/- (3/9) vs WT (3/36)")

#### Comparison including onsets of ageing for each genotype ####
source("/Users/malia/Desktop/Capstone_2022/code/alt-DEG_369.R")

volcano_plot_alt(gfinal_369, "Volcano Plot: TERT -/- (3/9) vs WT (3/9)")
volcano_plot_alt(mfinal_369, "Volcano Plot: TERT -/- (3/9) vs WT (3/9)")
volcano_plot_alt(tfinal_369, "Volcano Plot: TERT -/- (3/9) vs WT (3/9)")





