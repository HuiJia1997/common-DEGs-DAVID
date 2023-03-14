setwd("C:/Users/Janice/Desktop/SN/ARNI")

#setp 1 data prepration
rm(list = ls())
library(GEOquery)
library(limma)
library(Biobase)
library(dplyr)
library(data.table)
exp<-fread("GSE207283.txt")
exp<-data.frame(ID=exp$`Gene ID`,Symbol=exp$`Gene Symbol`,
                ARNI1=exp$`ARNI_A_1 Read Count`,
                ARNI2=exp$`ARNI_A_2 Read Count`,
                ARNI3=exp$`ARNI_A_3 Read Count`,
                control1=exp$`WKY_A_1 Read Count`,
                control2=exp$`WKY_A_2 Read Count`,
                control3=exp$`WKY_A_3 Read Count`)
exp$Symbol<-gsub("'","",exp$Symbol)
exp<-na.omit(exp)
EXP<-exp[,-c(1,2)]
rm<-rowMeans(EXP)
exp$rm<-rm
exp<- arrange(exp,Symbol,desc(rm))
exp<- exp[!duplicated(exp$Symbol),]
rownames(exp)<-exp$Symbol
exp<-exp[,-c(1,2,9)]

x<-as.matrix(exp)
f<-data.frame(symbol=rownames(x))
rownames(f)<-f$symbol
p<- read.csv("group.csv")
save(x,f,p,file = "DataPrep.Rdata")



#step 2 expressionset
rm(list = ls())
setwd("C:/Users/Janice/Desktop/SN/ARNI")
getwd()
load("ARNIDataPrep.Rdata")
rownames(p)<-p$sample
p$sample<-NULL
class(x)
class(p)
class(f)
boxplot(x[5,]~p[,"group"],main=f[5,"symbol"])
library(limma)
library(Biobase)
library(dplyr)
eSet<-ExpressionSet(assayData = x,
                    phenoData = AnnotatedDataFrame(p),
                    featureData = AnnotatedDataFrame(f))
dim(eSet)
plotDensities(eSet,legend = F)

exprs(eSet)<-log2(exprs(eSet)+1)
boxplot(exprs(eSet))
plotDensities(eSet,legend = F)

abline(v=1)
keep<-rowMeans(exprs(eSet))>1  #过滤低表达基因，原因见PMID: 26737772
eSet<-eSet[keep,]
plotDensities(eSet,legend = F)

exprs(eSet)<-normalizeBetweenArrays(exprs(eSet))
plotDensities(eSet,legend = F)
boxplot(exprs(eSet), col = c( "#DE480E","#DE480E","#DE480E","#70C2BB","#70C2BB","#70C2BB"))

plotMDS(eSet,labels = pData(eSet)[,"group"],gene.selection = "common", col=c(rep("#70C2BB",3), rep("#DE480E",3)))
save(eSet,file="expressionset.Rdata")

#step3 DEG analysis
rm(list = ls())
load("expressionset.Rdata")
P<-pData(eSet)
X=exprs(eSet)
F=fData(eSet)
design<-model.matrix(~0+group,data = P)
head(design)
colSums(design)
cm<-makeContrasts(status=groupARNI-groupcontrol,levels = design)
cm
fit<-lmFit(eSet,design)
head(fit)
fit2<-contrasts.fit(fit,contrasts = cm)
head(fit2$coefficients)
fit2<-eBayes(fit2)
head(fit2)
tempOutput=topTable(fit2,coef = 1,n=Inf,adjust="BH")


rowname<-rownames(X)
which(rowname=="Mettl24")
boxplot(X[8478,]~P[,"group"],main=F[8478,"symbol"])


nrDEG=na.omit(tempOutput)
write.csv(nrDEG,"ALLA.csv")
library(data.table)
nrDEG<-fread("ALLA.csv")
nrDEG$V1<-NULL
p=0.05
foldChange=1
diffSig = nrDEG[(nrDEG$P.Value < p & (nrDEG$logFC>foldChange |nrDEG$logFC<(-foldChange))),]
write.csv(diffSig,"SigA.csv")
diffUp = diffSig[diffSig$logFC>0,]
write.csv(diffUp,"UPA.csv")
diffDown = diffSig[diffSig$logFC<0,]
write.csv(diffDown,"DOWNA.csv")


#volcano plot
rm(list = ls())
setwd("C:/Users/Janice/Desktop/SN/ARNI")
library(ggVolcano)
# load the data
library(data.table)
data<-fread("ALLA.csv")
data<-as.data.frame(data)

# use the function -- add_regulate to add a regulate column 
# to the DEG result data. 
data <- add_regulate(data, log2FC_name = "logFC",
                     fdr_name = "P.Value",log2FC = 1,
                     fdr = 0.05)
# plot
ggvolcano(data, 
          x = "log2FoldChange", y = "padj",
          fills = c("#8EE3C8","#e0e0e0","#603162"),
          colors = c("#8EE3C8","#e0e0e0","#603162"),
          label = "symbol",add_label = F,
          legend_position="DR",y_lab="-log10Pvalue",
          x_lab="log2FC")





