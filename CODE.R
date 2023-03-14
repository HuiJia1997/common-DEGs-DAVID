setwd("C:/Users/Janice/Desktop/SN/common")
getwd()

#common DEGs
rm(list = ls())
library(readr)
library(data.table)
SigA<-fread("SigA.csv")
SigD<-fread("SigD.csv")
SigA$V1<-NULL
SigD$V1<-NULL
Sig<-merge(SigA,SigD,by="symbol")
write.csv(Sig,"Sig.csv")


#venn diagram
library("VennDiagram")                                                       
grid.newpage()                                                                
draw.pairwise.venn(area1 = 996, 
                   area2 = 1581,
                   cross.area = 116,
                   fill = c("#F07928","#B0D8E1"),
                   category=c("ARNI","Cardiotoxic"),
                   lty = "blank",                                            
                   alpha = 0.5,                                              
                   scaled = F) 
?draw.quad.venn


#volcano plot


#ARNIHeatmap
rm(list = ls())
setwd("C:/Users/Janice/Desktop/SN/common")
library(pheatmap)
library(ggplot2)
library(limma)
library(readr)
library(data.table)
load("ARNIDataPrep.Rdata")
ARNIp<-p
ARNIx<-x
Sig<-fread("Sig.csv")

ARNIx = log2(ARNIx+1)
boxplot(ARNIx)
ARNIx<-normalizeBetweenArrays(ARNIx)
ARNIx<-ARNIx[match(Sig$symbol,rownames(ARNIx)),]
rownames(ARNIp)<-ARNIp$sample
annotation_col<-data.frame(group=ARNIp)
colnames(ARNIx)
rownames(annotation_col)
sum(is.na(ARNIx))
hm<- ARNIx[apply(ARNIx, 1, function(x) sd(x)!=0),]
rm(f,p,x,Sig,ARNIp,ARNIx)

p <- pheatmap(hm,scale = "row",angle_col = 45,
              legend = T,fontsize = 7,
              color =  colorRampPalette(c("#88C63C", "#263165"))(100),
              treeheight_col = 20,
              treeheight_row = 20,
              cluster_cols=T,
              clustering_distance_rows="manhattan",
              clustering_method="ward.D2",
              annotation_col = annotation_col,
              show_rownames=F,annotation_legend=F)
?pheatmap
?colorRampPalette


#Doxorubicin Heatmap
rm(list = ls())
setwd("C:/Users/Janice/Desktop/SN/common")
library(pheatmap)
library(ggplot2)
library(limma)
library(readr)
library(data.table)
load("DoxorubicinDataPrep.Rdata")
DBCp<-p
DBCx<-x
Sig<-fread("Sig.csv")


DBCx = log2(DBCx+1)
boxplot(DBCx)
DBCx<-normalizeBetweenArrays(DBCx)
DBCx<-DBCx[match(Sig$symbol,rownames(DBCx)),]
rownames(DBCp)<-DBCp$sample
annotation_col<-data.frame(group=DBCp)
colnames(DBCx)
rownames(annotation_col)
sum(is.na(DBCx))
hm<-DBCx[apply(DBCx, 1, function(x) sd(x)!=0),]
rm(f,p,x,Sig,DBCp,DBCx)

p <- pheatmap(hm,scale = "row",angle_col = 45,
              legend = T,fontsize = 7,
              color =  colorRampPalette(c("#EC8C32", "#AAB6E0"))(100),
              treeheight_col = 20,
              treeheight_row = 20,
              cluster_cols=T,
              clustering_distance_rows="manhattan",
              clustering_method="ward.D2",
              annotation_col = annotation_col,
              show_rownames=F,annotation_legend=F)



#KEGG barplot
rm(list = ls())
library(readxl)
library(ggplot2)
library(viridis)
library(tidyverse)
library(data.table)
library(RColorBrewer)
setwd("C:/Users/Janice/Desktop/SN/common")
ed<-fread("edo.csv")
ed<-na.omit(ed)
p<-ggplot(ed,aes(x=Description,y=Count,fill=pvalue))+ 
          geom_bar(stat = "identity",width = 0.5)+
          coord_flip()+
          theme_bw()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
          labs(title = "KEGG Pathways")
p

#GO dotplot
rm(list = ls())
library(ggplot2)
library(tidyverse)
library(data.table)
library(viridis)
setwd("C:/Users/Janice/Desktop/SN/common")
BP<-fread("GO_BP.txt")
CC<-fread("GO_CC.txt")
MF<-fread("GO_MF.txt")
go<-rbind(BP[1:10,],CC[1:10,],MF[1:10,])
go$Category[which(go$Category=="GOTERM_BP_DIRECT")]="BP"
go$Category[which(go$Category=="GOTERM_CC_DIRECT")]="CC"
go$Category[which(go$Category=="GOTERM_MF_DIRECT")]="MF"
go$Category<-as.factor(go$Category)
str(go)
term<-data.frame(str_split_fixed(go$Term,"~",2))
go<-data.frame(category=go$Category,ID=term$X1,
               term=term$X2,pvalue=go$PValue,
               genes=go$Genes,count=go$Count)
rm(term,BP,CC,MF)
y <- ggplot(go, aes(term,pvalue,size=count,color=-1*log10(pvalue))) + 
     geom_point()+
     coord_flip()+
     scale_color_viridis_c(option = "C")+
     labs(color=expression(-log[10](pvalue)),size="Count",  
           x="Terms",y="pvalue",title="Gene Ontology Terms")+
     theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
     facet_grid(category~.,scales = "free_y",space =  "free_y")+
     theme_light()
y

