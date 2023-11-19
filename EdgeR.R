rm(list = ls())

setwd("C:/Users/demik/Desktop/Project Bioinformatics II")


library(edgeR)
library(gplots)
library(caTools)
library(openxlsx)


annot<-read.xlsx("C:/Users/demik/Desktop/Project Bioinformatics II/ens_94_human.xlsx")
head(annot)

#load raw count data
data<-read.xlsx("C:/Users/demik/Desktop/Project Bioinformatics II/htseq_merged.xlsx")
colnames(data)
data[1:10,1:4]
rownames(data)<-data$gene_id

colnames(data)
names(data) = gsub(pattern = "H1_htseq.txt", replacement = "H1", x = names(data))
names(data) = gsub(pattern = "H2_htseq.txt", replacement = "H2", x = names(data))
names(data) = gsub(pattern = "H3_htseq.txt", replacement = "H3", x = names(data))
names(data) = gsub(pattern = "S1_htseq.txt", replacement = "S1", x = names(data))
names(data) = gsub(pattern = "S2_htseq.txt", replacement = "S2", x = names(data))
names(data) = gsub(pattern = "S3_htseq.txt", replacement = "S3", x = names(data))
names(data)


#remove first 5 lines(alignment info) and 1st column(gene_id)
df<-data[-c(1:5),-1]
df[1:10,1:4]
head(df)
colnames(df)
rownames(df)
write.xlsx(df,"raw_counts.xlsx", colNames=TRUE,rowNames=TRUE)

#load metadata file
covariates<-read.xlsx("metadata.xlsx")
covariates$sample


##Edit covariates object
rownames(covariates)<-covariates$sample
covariates$group<-as.factor(covariates$group)

#check and double-check if rownames(covariates) match colnames(df)
all(rownames(covariates) == colnames(df))

genomic_idx <- match(rownames(covariates), colnames(df))
genomic_idx

df_ord  <- df[,genomic_idx]

all(rownames(covariates) == colnames(df_ord))






######## Differential Gene Expression ########
df1 <- DGEList( counts=df_ord,group=covariates$group) 

#perform filtering before TMM normalization (calcNormFactors)
df1 <- calcNormFactors(df1)

## keep genes that have more than 1 CPM for at least 4 samples (10% samples)
keep <- rowSums( cpm( df1 ) > 1 ) >= 4
df2 <- df1[keep, ]
dim(df1$counts)
dim(df2$counts)


df2 <- DGEList( counts=df2, group=covariates$group)
df2 <- calcNormFactors(df2)

#export CPM counts
logcpm <- as.data.frame(cpm(df2, log=TRUE))
head(logcpm)
head(annot)

logcpm$Gene.stable.ID.version <-rownames(logcpm)
logcpm_m <-merge(logcpm,annot,by ="Gene.stable.ID.version")
colnames(logcpm_m)
head(logcpm_m)


cpm <- as.data.frame(cpm(df2, log=FALSE))
cpm$Gene.stable.ID.version<-rownames(cpm)
cpm_m<-merge(cpm,annot,by="Gene.stable.ID.version")
colnames(cpm_m)


##Export logCPM values for all 
write.xlsx(logcpm, file="logCPM.xlsx", colNames=TRUE)
##Export CPM values for all 
write.xlsx(cpm, file="CPM.xlsx", colNames=TRUE)

colnames(logcpm)
logcpm<-logcpm[,-7]
cpm<-cpm[,-7]


logcpm<-logcpm[,-7]
cpm<-cpm[,-7]

pdf("mds.pdf")
plotMDS(logcpm)
dev.off()

##### PCA ######

library(FactoMineR)
library(factoextra)
library(dplyr)

pca<-PCA(t(cpm),scale.unit = T,graph = F)

jpeg("PCA.jpeg", res = 600, height = 15, width = 15, units="cm")
fviz_pca_ind(pca,palette = c("green","orange"),axes = c(1,2),col.ind = covariates$group,legend.title="Project",mean.point=F)
dev.off()
jpeg("PCA1.jpeg", res = 600, height = 15, width = 15, units="cm")
fviz_pca_ind(pca,palette = c("pink","cyan"),axes = c(1,2),col.ind = covariates$sex,legend.title="Project",mean.point=F)
dev.off()




############# Healthy  vs Sle  ###############

design <- model.matrix(~0+covariates$group)
df3 <- estimateDisp( df2, design, verbose=TRUE)
design

fit <- glmQLFit(df3, design)
lrt <- glmQLFTest(fit, contrast = c(-1,1)) 
print(lrt)
de<-as.data.frame(lrt$table)
lrt$table$fdr <- p.adjust(lrt$table$PValue, method="BH")
top <- lrt$table[lrt$table$PValue < 0.01,]
top2 <- top[top$logFC >0.58,]

#merge with annotation
head(annot)
head(top)
top$Gene.stable.ID.version<-rownames(top)
top_m<-merge(top,annot,by="Gene.stable.ID.version")
head(top_m)

all<-lrt$table 
all$Gene.stable.ID.version <-rownames(all)
all_m<-merge(all,annot,by="Gene.stable.ID.version")


#save

write.xlsx(top_m, file="Healthy_vs_sle_pval.xlsx", colΝames=TRUE, rowΝames=TRUE)
write.xlsx(all_m, file="Healty_vs_Sle_all.xlsx", colΝames=TRUE)


#### Box Plot######

jpeg("gene_a.jpeg")
boxplot(as.numeric(logcpm["ENSG00000049192.14", ]) ~ covariates$group, main="gene_a", xlab="group",ylab="logcpm")
dev.off()