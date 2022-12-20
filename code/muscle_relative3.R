source("/Users/malia/Desktop/Capstone_2022/code/function_hub.R")

########## MUSCLE ##################################
##### DEG Tables RELATIVE TO 3 MONTH TIMEPOINT #####
fitted_muscle <- fit_model('muscle')


# Make contrasts between specific time points #
muscle_wt_36 <- timeseries_comparison_test(fitted_muscle, "wt", "m",
                                        "3", "6",
                                        TRUE,
                                        "Muscle_WT-3-6")

muscle_wt_39 <- timeseries_comparison_test(fitted_muscle, "wt", "m",
                                        "3", "9",
                                        TRUE,
                                        "Muscle_WT-3-9")

muscle_wt_324 <- timeseries_comparison_test(fitted_muscle, "wt", "m",
                                         "3", "24",
                                         TRUE,
                                         "Muscle_WT-3-24")


muscle_wt_336 <- timeseries_comparison_test(fitted_muscle, "wt", "m",
                                         "3", "36",
                                         TRUE,
                                         "Muscle_WT-3-36")

muscle_mut_36 <- timeseries_comparison_test(fitted_muscle, "mut", "m",
                                         "3", "6",
                                         TRUE,
                                         "Muscle_MUT-3-6")

muscle_mut_39 <- timeseries_comparison_test(fitted_muscle, "mut", "m",
                                         "3", "9",
                                         TRUE,
                                         "Muscle_MUT-3-9")

