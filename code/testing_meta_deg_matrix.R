# MEGA DESIGN MATRIX

organ_info <- mega %>% 
  mutate(labels = str_replace_all(labels, "i", "")) %>%
  mutate(labels = str_replace_all(labels, "v", "")) %>%
  select(labels, genotype, age, organ) 

#organ_data['labels'] <- gsub("i", "", organ_data['labels'])

organ_columns <- gsub("i", "", colnames(all_organs))
colnames(all_organs) <- organ_columns

organ_levels_test <- c(unique(organ_info$labels))
f <- factor(colnames(all_organs), levels = organ_levels_test)

design <- model.matrix(~0 + f)
colnames(design) <- organ_levels_test


fit <- lmFit(organ_data, design)
fit$levels <- organ_levels
return(fit)





key <- organ_keys['gut' ]
organ_data <- all_organs %>% select(contains(paste(key)))
organ_columns <- gsub("i", "", colnames(organ_data))
colnames(organ_data) <- organ_columns

organ_levels <- unique(colnames(organ_data))
f <- factor(colnames(organ_data), levels = organ_levels)

design <- model.matrix(~0 + f)
colnames(design) <- organ_levels