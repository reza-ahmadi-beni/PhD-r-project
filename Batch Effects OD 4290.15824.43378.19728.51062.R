setwd("C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Data")
library(data.table)
library(limma)
library(ggplot2)
library(sva)
data1 <- read.delim("GSE4290_series_matrix.txt", comment.char="!")
max(data1[,-1])
min(data1[,-1])
data1[,-1] <- log2(1+data1[,-1])
data1 <- data1[,c(1:8,10:11,15,17,19,22:24,27:28,30:33,35,38:41,44:48,50:54,56,60:61,63,65:181)]
data2 <- read.delim("GSE15824_series_matrix.txt", comment.char="!")
max(data2[,-1])
min(data2[,-1])
data2[,-1] <- log2(data2[,-1])
data2 <- data2[,1:28]
data3 <- read.delim("GSE43378_series_matrix.txt", comment.char="!")
max(data3[,-1])
min(data3[,-1])
data4 <- read.delim("GSE19728_series_matrix.txt", comment.char="!")
max(data4[,-1])
min(data4[,-1])
data4 <- data4[,c(1,5:19)]
data4[,-1] <- log2(data4[,-1])
data5 <- read.delim("GSE51062_series_matrix.txt", comment.char="!")
max(data5[,-1])
min(data5[,-1])
data5[,-1] <- log2(data5[,-1])
annot1 <- fread("GPL570.annot", skip="!platform_table_begin", data.table=F)
annot1 <- annot1[,c("ID", "Gene symbol", "Gene ID")]
colnames(annot1)[2] <- "Symbol"
rownames(annot1) <- annot1$ID
data1$Entrez <- annot1[as.character(data1$ID_REF), "Symbol"]
data2$Entrez <- annot1[as.character(data2$ID_REF), "Symbol"]
data3$Entrez <- annot1[as.character(data3$ID_REF), "Symbol"]
data4$Entrez <- annot1[as.character(data4$ID_REF), "Symbol"]
data5$Entrez <- annot1[as.character(data5$ID_REF), "Symbol"]
data1 <- data1[,-1]
data2 <- data2[,-1]
data3 <- data3[,-1]
data4 <- data4[,-1]
data5 <- data5[,-1]
data1 <- aggregate(.~Entrez, data1, mean)
data2 <- aggregate(.~Entrez, data2, mean)
data3 <- aggregate(.~Entrez, data3, mean)
data4 <- aggregate(.~Entrez, data4, mean)
data5 <- aggregate(.~Entrez, data5, mean)
####
listdata <- list(data1, data2, data3, data4, data5)
all <- Reduce(function(...) merge(..., by = "Entrez"), listdata)
all <- subset(all, Entrez != "")
#all <- merge(data1, data2, by="Entrez")

rownames(all) <- all$Entrez
all <- subset(all, select=-Entrez)
pdf("C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/all_boxplot.pdf", width=64)
boxplot(all)
dev.off()
#write.table(all, "allGSE4290.15824logaritmic.txt", row.names=T, sep="\t", quote = F)


gr <- read.delim("Batch.4290.15824.43378.19728.51062 - Copy.txt")
#gr1<-gr1[,-3]
Group<-gr$Grade

#all1=normalizeQuantiles(all)
#all1 is optional

####
batch <- factor(c(rep(1, ncol(data1)-1), rep(2, ncol(data2)-1), rep(3, ncol(data3)-1), rep(4, ncol(data4)-1), rep(5, ncol(data5)-1)))
allm <- all - rowMeans(all)
pc <- prcomp(allm)
pcr <- data.frame(pc$r[,1:3], batch, Group)
pdf("C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/all_ggplot.pdf")
ggplot(pcr, aes(PC1,PC2, color=batch, shape=Group))+geom_point()+theme_bw()
#+scale_shape_manual(values=seq(0,15))
###or +scale_shape_manual(values=LETTERS[1:7]) 
#or
#+scale_color_manual(values = c("Low"="#FF0000","High"="#0000FF", "BFP"="#E69F00", "NT"="#009E73"))
#or cbp1
#+scale_color_manual(values = cbp1)
#cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#or library(viridis)
#+scale_color_viridis(discrete = TRUE, option = "D")
dev.off()

memory.limit(size = 5000)
allc <- ComBat(as.matrix(all),batch)
pdf("C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/allc_boxplot.pdf", width=64)
boxplot(allc)
dev.off()

allm <- allc - rowMeans(allc)
pc <- prcomp(allm)
pcr <- data.frame(pc$r[,1:3], batch, Group)
pdf("C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/pcr_ggplot.pdf")
ggplot(pcr, aes(PC1,PC2, color=batch, shape=Group))+geom_point()+theme_bw()
#+scale_shape_manual(values=seq(0,15))
dev.off()

pdf("C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/pcr_ggplot1.pdf")
ggplot(pcr, aes(PC1,PC2, color=Group))+geom_point()+theme_bw()+scale_color_manual(values = c("Low"="#FF0000","High"="#0000FF"))
#+scale_shape_manual(values=seq(0,15))
dev.off()

gset <- new("ExpressionSet", exprs=as.matrix(allc))
gr1 <- gr
gr1 <- factor(gr1$Grade)
gset$description <- gr1
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(gr1)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(High-Low, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT$Gene.Symbol<- rownames(tT)
tT <- subset(tT, select=c("Gene.Symbol", "adj.P.Val", "logFC"))
dim(tT)
write.table(tT, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/Batch.4290.15824.43378.19728.51062.High-Low.txt", row.names=T, sep="\t", quote = F)

library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
pdf("C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/CorHeatmap.pdf", width = 30, height = 30)
pheatmap(cor(allc), labels_row = gr1, labels_col = gr1, color = greenred(256), border_color = NA)
dev.off()


tT$log10p=-log10(tT$adj.P.Val)
library(ggplot2)
#install.packages("ggrepel")
library(ggrepel)
y=ggplot(data = tT,aes(x = logFC,y = log10p )) + geom_point(size=1) +theme_bw(base_rect_size = 1,base_line_size = NA)  
y
up=subset(tT,logFC>=(1)&log10p> 1)
up$Genesymbol<-rownames(up)
down=subset(tT,logFC<=(-1)&log10p> 1)
down$Genesymbol<-rownames(down)
y1=y+annotate("point",up$logFC,up$log10p,col="green",size=1)
y1
plot=y1+annotate("point",down$logFC,down$log10p,col="red",size=1)
plot1=plot+ ylab("-log10 (adj.P.Value)")+ xlab("log2(fold change)") + theme(text =element_text(size = 16) )+geom_vline(xintercept = 1,lty=5)+geom_vline(xintercept = -1,lty=5)
#+geom_hline(yintercept = 0,lty=5)
pdf("C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/volgse.4290.15824.43378.19728.51062.pdf")
plot1
dev.off()
#plot2=plot+ ylab("-log10 (adj.P.Value)")+ xlab("log2(fold change)") + theme(text =element_text(size = 16) )+geom_vline(xintercept = 1.3,lty=5)+geom_vline(xintercept = -1.3,lty=5)
#+geom_hline(yintercept = 1.3,lty=5)
#pdf("C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/volcanocirc7.pdf")
#plot2
#dev.off()

gbm.up <- subset(tT, logFC > 2 & adj.P.Val < 0.001)
dim(gbm.up)
gbm.up.genes <- unique(gbm.up$Gene.Symbol)

write.table(gbm.up.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.UP_LogFC2.txt", row.names=F, col.names = F, sep="\t", quote = F)

gbm.down <- subset(tT, logFC < -2 & adj.P.Val < 0.001)
dim(gbm.down)
gbm.down.genes <- unique(gbm.down$Gene.Symbol)

write.table(gbm.down.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.DOWN_LogFC2.txt", row.names=F, col.names = F, sep="\t", quote = F)

gbm.up <- subset(tT, logFC > 1.9 & adj.P.Val < 0.001)
dim(gbm.up)
gbm.up.genes <- unique(gbm.up$Gene.Symbol)

write.table(gbm.up.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.UP_LogFC1.9.txt", row.names=F, col.names = F, sep="\t", quote = F)

gbm.down <- subset(tT, logFC < -1.9 & adj.P.Val < 0.001)
dim(gbm.down)
gbm.down.genes <- unique(gbm.down$Gene.Symbol)

write.table(gbm.down.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.DOWN_LogFC1.9.txt", row.names=F, col.names = F, sep="\t", quote = F)


gbm.up <- subset(tT, logFC > 1.8 & adj.P.Val < 0.001)
dim(gbm.up)
gbm.up.genes <- unique(gbm.up$Gene.Symbol)

write.table(gbm.up.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.UP_LogFC1.8.txt", row.names=F, col.names = F, sep="\t", quote = F)

gbm.down <- subset(tT, logFC < -1.8 & adj.P.Val < 0.001)
dim(gbm.down)
gbm.down.genes <- unique(gbm.down$Gene.Symbol)

write.table(gbm.down.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.DOWN_LogFC1.8.txt", row.names=F, col.names = F, sep="\t", quote = F)


gbm.up <- subset(tT, logFC > 1.7 & adj.P.Val < 0.001)
dim(gbm.up)
gbm.up.genes <- unique(gbm.up$Gene.Symbol)

write.table(gbm.up.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.UP_LogFC1.7.txt", row.names=F, col.names = F, sep="\t", quote = F)

gbm.down <- subset(tT, logFC < -1.7 & adj.P.Val < 0.001)
dim(gbm.down)
gbm.down.genes <- unique(gbm.down$Gene.Symbol)

write.table(gbm.down.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.DOWN_LogFC1.7.txt", row.names=F, col.names = F, sep="\t", quote = F)


gbm.up <- subset(tT, logFC > 1.6 & adj.P.Val < 0.001)
dim(gbm.up)
gbm.up.genes <- unique(gbm.up$Gene.Symbol)

write.table(gbm.up.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.UP_LogFC1.6.txt", row.names=F, col.names = F, sep="\t", quote = F)

gbm.down <- subset(tT, logFC < -1.6 & adj.P.Val < 0.001)
dim(gbm.down)
gbm.down.genes <- unique(gbm.down$Gene.Symbol)

write.table(gbm.down.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.DOWN_LogFC1.6.txt", row.names=F, col.names = F, sep="\t", quote = F)


gbm.up <- subset(tT, logFC > 1.5 & adj.P.Val < 0.001)
dim(gbm.up)
gbm.up.genes <- unique(gbm.up$Gene.Symbol)

write.table(gbm.up.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.UP_LogFC1.5.txt", row.names=F, col.names = F, sep="\t", quote = F)

gbm.down <- subset(tT, logFC < -1.5 & adj.P.Val < 0.001)
dim(gbm.down)
gbm.down.genes <- unique(gbm.down$Gene.Symbol)

write.table(gbm.down.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.DOWN_LogFC1.5.txt", row.names=F, col.names = F, sep="\t", quote = F)

gbm.up <- subset(tT, logFC > 1 & adj.P.Val < 0.001)
dim(gbm.up)
gbm.up.genes <- unique(gbm.up$Gene.Symbol)

write.table(gbm.up.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.UP_LogFC1.txt", row.names=F, col.names = F, sep="\t", quote = F)

gbm.down <- subset(tT, logFC < -1 & adj.P.Val < 0.001)
dim(gbm.down)
gbm.down.genes <- unique(gbm.down$Gene.Symbol)

write.table(gbm.down.genes, "C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/4290.15824.43378.19728.51062.DOWN_LogFC1.txt", row.names=F, col.names = F, sep="\t", quote = F)

###heatmap###
mdata5<-allc
#tT2=subset(tT,tT$adj.P.Val<0.001 & abs(tT$logFC)>=1.6079)
#rownames(tT2)
#colnames(tT2)
mdata5<-mdata5[c("LTF",                                   "IGFBP2"                                
                 ,"TIMP1",                                  "LOC101928916///NNMT"                                
                 ,"SERPINH1",                                  "COL4A1"                                
                 ,"CHI3L1",                                "EMP3"                                 
                 ,"IL13RA2",                                 "VMP1"                                
                 ,"SPOCD1",                               "SERPINA3"                                
                 ,"COL4A2",                                 "KIF20A"
                 ,"PTX3",                                    "FOXM1"                                  
                 ,"CA3",                                   "IGFBP3"                                 
                ,"FCGBP",                                  "TACC3"                                 
                 ,"FSTL5",                                 "NTSR2"                                
                 ,"SFRP2",                                  "ETNPPL"                                 
                 ,"KLRC2///KLRC1",                         "CSMD3"                                  
                 ,"SNAP91",                                "CNTN3"                                 
                 ,"HMP19",                                  "SPHKAP"                                   
                 ,"MGAT4C",                                 "SMOC1"                                  
                 ,"PTPRT",                                 "SH3GL2"                                 
                 ,"ZNF488",                                "KLRC3"                             
                 ,"FAM133A",                                 "INA"
                  ,"PRLHR",                                 "SELL"), ]
mdata5<-mdata5[c("LTF",                            "EMP3"                                
            ,"TIMP1",                                  "COL4A2"                                
            ,"CSMD3",                                  "COL4A1"                                
            ,"FAM133A",                                "PDPN"                                 
            ,"VMP1",                                   "LOC101930489///MIR4435-2HG///LINC00152"
            ,"TNFRSF12A",                              "CDC20"                                 
            ,"RASL10A",                                "SERPINH1"                              
            ,"CHI3L1",                                 "SHOX2"                                 
            ,"COL5A2",                                 "IGFBP3"                                
            ,"SERPINA3",                               "VEGFA"                                 
            ,"NCAPG",                                  "SPHKAP"                                
            ,"NARR///RAB34",                           "KIF20A"                                
            ,"LOC101928916///NNMT",                    "TACC3"                                 
            ,"SERPINE1",                               "CRNDE"                                 
            ,"ADM",                                    "GPX8"                                  
            ,"CA3",                                    "GALNT13"                               
            ,"FOXM1",                                  "PLP2"                                  
            ,"CCEPR",                                  "HPSE2"                                 
            ,"SELL",                                   "FCGBP"                                 
            ,"MGAT4C",                                 "FSTL5"                                 
            ,"FBXO17",                                 "SLITRK5"                               
            ,"SFRP2",                                  "HJURP"                                 
            ,"MIR155///MIR155HG",                      "COL5A1"                                
            ,"FABP5",                                  "CSDC2"                                 
            ,"PRLHR",                                  "PTX3"                                  
            ,"PLEKHA4",                                "KLRC3"                                 
            ,"NDC80",                                  "LTF"                                   
            ,"ETNPPL",                                 "RRM2"                                  
            ,"SPOCD1",                                 "SMOC1"                                 
            ,"METTL7B",                                "LINC00836"                             
            ,"AKR1C3",                                 "C8orf4"                                
            ,"PTPRT",                                  "COL3A1"                                
            ,"COL1A2",                                 "MELK"                                  
            ,"CHI3L2",                                 "DACH2"                                 
            ,"CTHRC1",                                 "PCOLCE"                                
            ,"HOXA5",                                  "LOX"                                   
            ,"PRR36",                                  "CEP55"                                 
            ,"HAR1A",                                  "CNTN3"                                 
            ,"INA",                                    "LOXL1"                                 
            ,"HMP19",                                  "KIF21B"                                
            ,"IBSP",                                   "RARRES2"                               
            ,"IL13RA2",                                "ACTL6B"                                
            ,"KLRC2///KLRC1",                          "SH3GL2"                                
            ,"DPP10",                                  "SNAP91"                                
            ,"ARHGAP36",                               "NTSR2"                                 
            ,"COL6A3",                                 "SLITRK1"                               
            ,"FLJ16779",                               "VSTM2A"                                
            ,"FREM3",                                  "ZNF488"                                
            ,"LOC286178",                              "SAA2///SAA1"                           
            ,"GPR17",                                  "GJB6"), ]
###pheatmap
#library(Biobase)
#library(GEOquery)
#library(parallel)
#library(base)
#library(reshape)
#colnames(mdata5)<-gr1$Accession	
#pheatmap(mdata5)
##colnames(mdata5)<-gr$Accession
#pheatmap(mdata5,color = c("red","black","green"))
#cols=colorRampPalette(c("red","black","green"))
#cols(30)
plot1<-pheatmap(mdata5,color =cols(10),cellwidth = 30,border_color = NA ,cellheight = 12,fontsize_row = 11,fontsize_col = 20,cluster_rows = T,cluster_cols = T,clustering_distance_cols = "correlation", clustering_distance_rows = "correlation" )
pdf("C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/pheatmap1.pdf",width = 130,height = 50)
plot1
dev.off()


plot1<-pheatmap(mdata5,color =cols(10),cellwidth = 30,cellheight = 12,border_color = NA,fontsize_row = 11,fontsize_col = 20,cluster_rows = T,cluster_cols = T,clustering_distance_cols ="euclidean" )
pdf("C:/Users/REZAHRA/Downloads/Data/GSEanalysis/Results/pheatmap2.pdf",width = 130,height = 50)
plot1
dev.off()

#colnames(mdata5)<-as.factor(gr1$Accession)
