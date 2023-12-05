# Code to download monocyte / macrophage microarray data from Xue et al. 2014 (GSE47189), normalise, collapse biological replicates and probes mapping to same gene (both by mean values), and use the resulting dataset to perform GSVA to identify conditions that most closely recapitulate state of monocytes from IBD.
library('GSVA')
library('edgeR')
library('affy')
library('GEOquery')
library('limma')
library('umap')
library('illuminaHumanv4.db')
data <- getGEO("GSE47189", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(data) > 1) idx <- grep("GPL6947", attr(data, "names")) else idx <- 1
data <- data[[idx]]
fvarLabels(data) <- make.names(fvarLabels(data))
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX00000000000000",
"00000000000000000000000000000000000000000000000000",
"00000000000000000000000000000000000000000000000000",
"00000000000000000000000000000000000000000000000000",
"00000000000000000000000000000000000000000000000000",
"00000000000000000000000000000000000000000000000000",
"00000000000000000000000000000000000000000000000000",
"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
sml <- strsplit(gsms, split="")[[1]]
sel <- which(sml != "X")
sml <- sml[sel]
data <- data[ ,sel]
ex <- exprs(data)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
(qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(data) <- log2(ex) }
exprs(data) <- normalizeBetweenArrays(exprs(data))
pData(data)$title=gsub("\\[.*\\]","",as.character(pData(data)$title))
data$group=paste(data$title,data$`activation stimuli:ch1`) #this comprises 67 unique combinations of activation stimuli and durations of exposure (including monocytes)
values=exprs(data)
colnames(values)=pData(data)$group
values2=aggregate(t(values), by = list(colnames(values)), FUN='median')
values3=as.matrix(t(values2[,-1]))
colnames(values3)=values2[,1]
genes=mget(rownames(values3),envir=illuminaHumanv4SYMBOL,ifnotfound=NA)
rownames(values3)=genes
values4=aggregate(values3, by=list(rownames(values3)), FUN = "max")
values5=noquote(as.matrix(values4[,-1]))
rownames(values5)=values4[,1]
macs.eset=new('ExpressionSet',exprs=as.matrix(unclass(values5)))
genelists=read.csv(file="ibd_genesets.csv",header=T,na.strings=c("")) #this contains lists of significantly upregulated and downregulated genes in IBD monocytes (differential expression analysis using limma from micrarray data from PMID 27015630).Microarray data are available in ArrayExpress, accession number E-MTAB-3554.
lists=lapply(genelists,function(col)col[!is.na(col)])
macs.gsva=gsva(macs.eset,lists,max.sz=2000)
res=exprs(macs.gsva)
res_dn=res[1,]*-1
res_up=res[2,]
res_sum=res_up+res_dn
res_sum=res_sum[order(-res_sum)]

# Co-expression of ETS2 across these 67 different macrophage activation conditions
library(Hmisc)
CorMatrix2 <- rcorr(t(as.matrix(unclass(values5))), type = "pearson")
ETS2_r <- CorMatrix2$r[, "ETS2"]
ETS2_P <- CorMatrix2$P[, "ETS2"]
ETS2_table2 <- cbind(ETS2_r, ETS2_P)
ETS2_table2 <- as.data.frame(ETS2_table2)
ETS2_table2$p.adj=p.adjust(ETS2_table2$ETS2_P,method="BH")
ETS2_table2=ETS2_table2[order(-ETS2_table2$ETS2_r),]
ETS2_table2$rank=c(length(ETS2_table2$ETS2_r):1)
oxphos=read.csv(file="oxphos.csv",header=T)
ox.phos=ETS2_table[rownames(ETS2_table) %in% oxphos$GOBP_OXIDATIVE_PHOSPHORYLATION,]
plot(ETS2_table$rank,ETS2_table$ETS2_r,pch=19,col="gray",cex=0.5,ylim=c(-1,1),xlab="Pearson r",ylab="")
abline(h=0.314572,lty=2,col="black") #Pearson r values corresponding to FDR P 0.05
abline(h=-0.31448,lty=2,col="black")
points(x=ox.phos$rank,y=ox.phos$ETS2_r,col="cornflowerblue",pch=19,cex=0.8)
points(ETS2_table[c("HIF1A","PFKFB3"),]$rank,ETS2_table[c("HIF1A","PFKFB3"),]$ETS2_r,col="red",pch=19,cex=0.8)

