source("/Users/malia/Desktop/Capstone_2022/code/function_hub.R")

########## TESTES ##################################
##### DEG Tables RELATIVE TO 3 MONTH TIMEPOINT in TESTES #####
fitted_testes <- fit_model('testis')


# Make contrasts between specific time points #
testes_wt_36 <- timeseries_comparison_test(fitted_testes, "wt",
                                           "t",
                                           "3",
                                           "6",
                                          TRUE, "Testes_WT-3-6")

testes_wt_39 <- timeseries_comparison_test(fitted_testes, "wt",
                                            "t",
                                            "3",
                                            "9",
                                             TRUE,
                                            "Testes_WT-3-9")

testes_wt_324 <- timeseries_comparison_test(fitted_testes, "wt", "t",
                                          "3", "24",
                                          TRUE,
                                          "Testes_WT-3-24")


testes_wt_336 <- timeseries_comparison_test(fitted_testes, "wt", "t",
                                          "3", "36",
                                          TRUE,
                                          "Testes_WT-3-36")

testes_mut_36 <- timeseries_comparison_test(fitted_testes, "mut", "t",
                                         "3", "6",
                                        TRUE,
                                         "Testes_MUT-3-6")

testes_mut_39 <- timeseries_comparison_test(fitted_testes, "mut", "t",
                                         "3", "9",
                                         TRUE,
                                         "Testes_MUT-3-9")

