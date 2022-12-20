### PCA ###
library(Biobase)
library(GEOquery)
library(limma) 
library(gplots)


### INDIVIDUAL ORGANS: NORMALIZED DATA ### 
gut_read <- read.csv("raw_data/Madalena_microarraysExpression_gut.csv", row.names = 1, sep = "\t") %>% 
  rename_all(tolower) # change all column names to lowercase

muscle_read <- read.csv("raw_data/Madalena_microarraysExpression_muscle.csv", row.names = 1, sep = "\t") %>% 
  rename_all(tolower) # change all column names to lowercase


testes_read <- read.csv("raw_data/Madalena_microarraysExpression_testes.csv", row.names = 1, sep = "\t") %>% 
  rename_all(tolower) # change all column names to lowercase



gut <- as.matrix(gut_read)
col.order <- c("wt_3_gi","wt_3_gii","wt_3_giii",
               "wt_6_gi","wt_6_gii","wt_6_giii",
               "wt_9_gi","wt_9_gii","wt_9_giii",
               "wt_24_gi","wt_24_gii","wt_24_giii",
               "wt_36_gi","wt_36_gii","wt_36_giii",
               "mut_3_gi","mut_3_gii","mut_3_giii",
               "mut_6_gi","mut_6_gii","mut_6_giii",
               "mut_9_gi","mut_9_gii","mut_9_giii")
gut <- gut[ , col.order]

sum(is.na(gut))
w <- which(apply(is.na(gut), 1, sum) > 0 )
gut <- gut[-w, ]
sum(is.na(gut))

boxplot(as.data.frame(gut))
logdata <- log2(gut)
boxplot(as.data.frame(logdata))
probemeans <- apply(logdata, 1, mean)
probesd <- apply(logdata, 1, sd)
plot(probemeans, probesd)

q25 <- quantile(probemeans, 0.25, na.rm = T)
whichtosave <- which(probemeans > q25)
q25logdata <- logdata[whichtosave,]

mydata <- q25logdata[apply(q25logdata, 1, IQR) > 0.25, ]
tdata <- aperm(mydata)
pca <- prcomp(tdata, scale=T)
plot(pca$x, col=c("#009999", "#0000FF", "#999999"))
text(pca$x, rownames(pca$x), cex=0.5)




conditions <- colnames(gut) 
conditions <- str_replace_all(conditions, c("i" = "", "v" = "", "_g" = ""))
plot(pca$x, type = "p",  col = c("#9ACD32", "#9ACD32", "#9ACD32",
                                 "#8FBC8F", "#8FBC8F", "#8FBC8F",
                                 "#008080", "#008080", "#008080",
                                 "#6495ED", "#6495ED", "#6495ED",
                                 "#6A5ACD", "#6A5ACD", "#6A5ACD",
                                 "#A0522D", "#A0522D", "#A0522D",
                                 "#FF8C00", "#FF8C00", "#FF8C00",
                                 "#DC143C", "#DC143C", "#DC143C"))

text(pca$x, labels=conditions, col = c("#9ACD32", "#9ACD32", "#9ACD32",
                                       "#8FBC8F", "#8FBC8F", "#8FBC8F",
                                       "#008080", "#008080", "#008080",
                                       "#6495ED", "#6495ED", "#6495ED",
                                       "#6A5ACD", "#6A5ACD", "#6A5ACD",
                                       "#A0522D", "#A0522D", "#A0522D",
                                       "#FF8C00", "#FF8C00", "#FF8C00",
                                       "#DC143C", "#DC143C", "#DC143C"),  
     cex=1, font = 2, pos = 3)

heatmap(q25logdata, col=greenred(100))

pearsonCorr <- as.dist(1 - cor(mydata))
hC <- hclust(pearsonCorr)
plot(hC, labels = conditions)





### TESTES 

testes <- as.matrix(testes_read)
col.order <- c("wt_3_ti","wt_3_tii","wt_3_tiii",
               "wt_6_ti","wt_6_tii","wt_6_tiii",
               "wt_9_ti","wt_9_tii","wt_9_tiii",
               "wt_24_ti","wt_24_tii","wt_24_tiii",
               "wt_36_ti","wt_36_tii","wt_36_tiii",
               "mut_3_ti","mut_3_tii","mut_3_tiii",
               "mut_6_ti","mut_6_tii","mut_6_tiii",
               "mut_9_ti","mut_9_tii","mut_9_tiii")
testes <- testes[ , col.order]


sum(is.na(testes))
w <- which(apply(is.na(testes), 1, sum) > 0 )
testes <- testes[-w, ]
sum(is.na(testes))

boxplot(as.data.frame(testes))
logdata <- log2(testes)
boxplot(as.data.frame(logdata))
probemeans <- apply(logdata, 1, mean)
probesd <- apply(logdata, 1, sd)
plot(probemeans, probesd)

q25 <- quantile(probemeans, 0.25, na.rm = T)
whichtosave <- which(probemeans > q25)
q25logdata <- logdata[whichtosave,]

mydata <- q25logdata[apply(q25logdata, 1, IQR) > 0.25, ]
tdata <- aperm(mydata)
pca <- prcomp(tdata, scale=T)
#plot(pca$x, col=c("#009999", "#0000FF", "#999999"))
#text(pca$x, rownames(pca$x), cex=0.5)




conditions <- colnames(testes) 
conditions <- str_replace_all(conditions, c("i" = "", "v" = "", "_t" = ""))
plot(pca$x, type = "p",  col = c("#9ACD32", "#9ACD32", "#9ACD32",
                                 "#8FBC8F", "#8FBC8F", "#8FBC8F",
                                 "#008080", "#008080", "#008080",
                                 "#6495ED", "#6495ED", "#6495ED",
                                 "#6A5ACD", "#6A5ACD", "#6A5ACD",
                                 "#A0522D", "#A0522D", "#A0522D",
                                 "#FF8C00", "#FF8C00", "#FF8C00",
                                 "#DC143C", "#DC143C", "#DC143C"))
text(pca$x, labels=conditions, col = c("#9ACD32", "#9ACD32", "#9ACD32",
                                       "#8FBC8F", "#8FBC8F", "#8FBC8F",
                                       "#008080", "#008080", "#008080",
                                       "#6495ED", "#6495ED", "#6495ED",
                                       "#6A5ACD", "#6A5ACD", "#6A5ACD",
                                       "#A0522D", "#A0522D", "#A0522D",
                                       "#FF8C00", "#FF8C00", "#FF8C00",
                                       "#DC143C", "#DC143C", "#DC143C"),  
     cex=1, font = 2, pos = 3)

heatmap(q25logdata, col=greenred(100))

pearsonCorr <- as.dist(1 - cor(mydata))
hC <- hclust(pearsonCorr)
plot(hC, labels = conditions)


### MUSCLE 

muscle <- as.matrix(muscle_read)
col.order <- c("wt_3_mi","wt_3_mii","wt_3_miii",
               "wt_6_mi","wt_6_mii","wt_6_miii",
               "wt_9_mi","wt_9_mii","wt_9_miii",
               "wt_24_mi","wt_24_mii","wt_24_miii",
               "wt_36_mi","wt_36_mii","wt_36_miii",
               "mut_3_mi","mut_3_mii","mut_3_miii",
               "mut_6_mi","mut_6_mii","mut_6_miii",
               "mut_9_mi","mut_9_mii","mut_9_miii")

muscle <- muscle[ , col.order]


sum(is.na(muscle))
w <- which(apply(is.na(muscle), 1, sum) > 0 )
muscle <- muscle[-w, ]
sum(is.na(muscle))

boxplot(as.data.frame(muscle))
logdata <- log2(muscle)
boxplot(as.data.frame(logdata))
probemeans <- apply(logdata, 1, mean)
probesd <- apply(logdata, 1, sd)
plot(probemeans, probesd)

q25 <- quantile(probemeans, 0.25, na.rm = T)
whichtosave <- which(probemeans > q25)
q25logdata <- logdata[whichtosave,]

mydata <- q25logdata[apply(q25logdata, 1, IQR) > 0.25, ]
tdata <- aperm(mydata)
pca <- prcomp(tdata, scale=T)
#plot(pca$x, col=c("#009999", "#0000FF", "#999999"))
#text(pca$x, rownames(pca$x), cex=0.5)




conditions <- colnames(muscle) 
conditions <- str_replace_all(conditions, c("i" = "", "v" = "", "_m" = ""))
plot(pca$x, type = "p",  col = c("#9ACD32", "#9ACD32", "#9ACD32",
                                 "#8FBC8F", "#8FBC8F", "#8FBC8F",
                                 "#008080", "#008080", "#008080",
                                 "#6495ED", "#6495ED", "#6495ED",
                                 "#6A5ACD", "#6A5ACD", "#6A5ACD",
                                 "#A0522D", "#A0522D", "#A0522D",
                                 "#FF8C00", "#FF8C00", "#FF8C00",
                                 "#DC143C", "#DC143C", "#DC143C"))
text(pca$x, labels=conditions, col = c("#9ACD32", "#9ACD32", "#9ACD32",
                                       "#8FBC8F", "#8FBC8F", "#8FBC8F",
                                       "#008080", "#008080", "#008080",
                                       "#6495ED", "#6495ED", "#6495ED",
                                       "#6A5ACD", "#6A5ACD", "#6A5ACD",
                                       "#A0522D", "#A0522D", "#A0522D",
                                       "#FF8C00", "#FF8C00", "#FF8C00",
                                       "#DC143C", "#DC143C", "#DC143C"),  
     cex=1, font = 2, pos = 3)

heatmap(q25logdata, col=greenred(100))

pearsonCorr <- as.dist(1 - cor(mydata))
hC <- hclust(pearsonCorr)
plot(hC, labels = conditions)

