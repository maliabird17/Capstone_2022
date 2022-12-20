

### code 

The R scripts for my original code in included in this folder. To replicate my results, all scripts can be run on R studio if pathnames are changed to your project location. The files are listed below in the same order as the pipeline seen in the analysis overview figure. 

- processing_all_data : loading all microarray data, basic data cleaning, and creating a matrix to store and reference all sample information. 

- function_hub : file of main functions used to identify DEGs. Organized into sections of the analysis pipeline: Model fitting (fit_model and timeseries_comparison_test), Volcano Plots (volcano_plot and volcano_plot_alt), 
 
- gut_relative3 : differencially expressed gene comparison in the gut, using the 3 month time point as the comparison point. 

- muscle_relative3 : differencially expressed gene comparison in the muscle, using the 3 month time point as the comparison point. 

- testes_relative3 : differencially expressed gene comparison in the testes, using the 3 month time point as the comparison point. 

- alt-DEG_369 : alternative analysis of differencially expressed genes (not used in results) which compares ageing across genotypes, using the 3, 6, and 9 month timepoints in both genotypes and contrasts that longterm progression. 

- alt-DEG_ageing-onset-369-36363 :  alternative analysis of differencially expressed genes (not used in results) which compares ageing across genotypes, using the 3, 6, and 9 month timepoints for the TERT mutant and the 3, 6, and 36 timepoints in the WT, representing their respective onsets of ageing.

- volcano_plot : plotting interactive volcano plots for DEGs analyses listed above. 
