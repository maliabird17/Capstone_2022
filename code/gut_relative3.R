source("/Users/malia/Desktop/Capstone_2022/code/function_hub.R")

########## GUT ################# 
##### DEG Tables RELATIVE TO 3 MONTH TIMEPOINT in GUT #####
fitted_gut <- fit_model('gut')



# Make contrasts between specific time points #
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

