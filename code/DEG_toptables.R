#### TOPTABLE ###
fit_model <- function(organ_) {
  
  organ_keys <- setNames(as.list(c('_g', '_t', '_m')), c('gut', 'testis', 'muscle'))
  
  key <- organ_keys[organ_]
  organ_data <- all_organs %>% select(contains(paste(key)))
  organ_columns <- gsub("i", "", colnames(organ_data))
  colnames(organ_data) <- organ_columns
  
  organ_levels <- unique(colnames(organ_data))
  f <- factor(colnames(organ_data), levels = organ_levels)
  
  design <- model.matrix(~0 + f)
  colnames(design) <- organ_levels
  
  fit <- lmFit(organ_data, design, plot = T)
  fit$levels <- organ_levels
  return(fit)
  
}
fitted <- fit_model('gut') 



# Time series comparison (9.6.1 limma documentation)
timeseries_comparison <- function(fit_object) {
  
  organ_levels <- fit_object$levels 
  
  contrast <- makeContrasts(
    first = "wt_36_g - wt_3_g",
    levels = organ_levels
  )
  
  fit2 <- contrasts.fit(fitted, contrast)
  fit2 <- eBayes(fit2)
  top.table <- topTable(fit2, n = Inf)
  
  return(top.table)
}


wttable_gut <- timeseries_comparison(fitted)
length(which(wttable_gut$adj.P.Val < 0.05))
length(which(wttable_gut$adj.P.Val < 0.05 & wttable_gut$logFC > 0))
length(which(wttable_gut$adj.P.Val < 0.05 & wttable_gut$logFC < 0))
