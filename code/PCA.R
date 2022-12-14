# PCA
library("affy")


gut <- as.matrix(gut)

gutdata <- mega %>% 
  filter(organ == 'gut')
  
pdata <- gutdata %>% select(sample, genotype, age)

rownames(pdata) <- gutdata$labels
rownames(pdata)

metadata <- data.frame(labelDescription =
                         c("Sample", "Genotype", "Age"), 
                       row.names = c("sample", "genotype", "age"))

phenodata <- new("AnnotatedDataFrame", data = pdata, varMetadata = metadata)
phenodata

exampleset <- ExpressionSet(assayData = gut, 
                            phenoData = phenodata)

tlog <- aperm(logdata)
pca <- prcomp(tlog, scale=T)


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

#mydata <- q25logdata[apply(q25logdata, 1, IQR) > 0.25, ]
tdata <- aperm(q25logdata)
pca <- prcomp(tdata, scale=T)
plot(pca$x, col=c("#009999", "#0000FF", "#999999"))
text(pca$x, rownames(pca$x), cex=0.5)




conditions <- colnames(gut) 
conditions <- str_replace_all(conditions, c("i" = "", "v" = ""))
plot(pca$x, type = "n")
text(pca$x, labels=conditions,  cex=0.5)

heatmap(mydata, col=greenred(100))
