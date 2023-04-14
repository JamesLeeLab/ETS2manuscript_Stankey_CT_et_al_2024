# Analysis of TPP macrophages treated with a MEK inhibitor
suppressMessages(library(limma))
suppressMessages(library(affy))
suppressMessages(library(edgeR))
suppressMessages(library(fgsea))
suppressMessages(library(data.table))
suppressMessages(library(plyr))
suppressMessages(library(ggplot2))
suppressMessages(library(gplots))
counts=data.matrix(read.table(file="raw_counts.txt",sep='\t',row.names=1,header=T))
names=paste0(rep(c("s1","s2","s3"),each=3),c("_100","_500","_ctrl"))
colnames(counts)=names
data2=as.matrix(counts)
eset=new('ExpressionSet',exprs=data2)
eset.expressed=rowSums(cpm(eset)>0.5)>=9
eset.exp=eset[eset.expressed,]
dim(eset.exp)
cpm.exp=cpm(eset)[eset.expressed,]
seq2=apply(counts,2,sum)
seq2
dge=DGEList(counts=eset.exp)
dge=calcNormFactors(dge)
dge$samples$donor=rep(c("s1","s2","s3"),each=3)
dge$samples$drug=rep(c("100nM","500nM","ctrl"),3)
design=model.matrix(~0+factor(dge$samples$drug)+factor(dge$samples$donor))
colnames(design)
colnames(design)=c("MEK_100","MEK_500","ctrl","s2","s3")
contrasts=makeContrasts(MEK_100nM=MEK_100-ctrl,MEK_500nM=MEK_500-ctrl, levels=colnames(design))
voom=voom(dge,design,plot=T)
fit=lmFit(voom,design)
fit2=contrasts.fit(fit,contrasts)
fit3=eBayes(fit2)
summary(decideTests(fit3, p=0.05, adjust.method = 'BH'))
MEK_100nM.res=topTable(fit3,number=12801,p.value=1,lfc=0,adjust.method="BH",coef=1)
MEK_500nM.res=topTable(fit3,number=12801,p.value=1,lfc=0,adjust.method="BH",coef=2)
#write.csv(MEK_100ng.res,file="PD_0.1_MEK.res.csv")
#write.csv(MEK_500ng.res,file="PD_0.5_MEK.res.csv")

# GSEA of ETS2-regulated macrophage gene sets
genelists=read.csv(file="macrophage_genesets_fgsea.csv",header=T,na.strings=c(""))
lists=lapply(genelists,function(col)col[!is.na(col)])
pd100=read.csv(file="PD_0.1_MEK.resSYMBOL.rank.csv",header=F,row.names=1)
pd100_num=unlist(pd100)
names(pd100_num)=rownames(pd100)
pd500=read.csv(file="PD_0.5_MEK.resSYMBOL.rank.csv",header=F,row.names=1)
pd500_num=unlist(pd500)
names(pd500_num)=rownames(pd500)
fgsea_pd100=fgsea(pathways=lists, stats=pd100_num, eps=0.0, minSize=15, maxSize=500)
fgsea_pd100[order(pval), ]
fgsea_pd500=fgsea(pathways=lists, stats=pd500_num, eps=0.0, minSize=15, maxSize=500)
fgsea_pd500[order(pval), ]
#fwrite(fgsea_pd500,"fgsea_PD_0.5_mac.pathways.csv")
#fwrite(fgsea_pd100,"fgsea_PD_0.1_mac.pathways.csv")

# Dot plot of GSEA results
data=read.csv(file="fgsea_PD_0.5_mac.pathways.csv",header=T)
data=data[-1,]
data$pathway = factor(data$pathway, levels=data[order(data$pval), "pathway"])
mapping=aes_string(x="NES",y="pathway",color="NES")
low="blue"
high="red"
mid="white"
limits=c(-3,3)
data$log10=-log10(data$pval)
ggp=ggplot(data, aes(x=NES, y=pathway, size=log10,color=NES))
ggp+ geom_point()+scale_size_area(limits=c(0.1,15),max_size = 10)+scale_color_gradient2(low = low, high = high, mid=mid, aesthetics = 'color', limits = limits)+ theme(axis.ticks.y=element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + xlim(-3,3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Heatmap of relative gene expression of chr21q22-regulated gene sets  
hm=read.csv(file="chr21_genesets_log2FC_MEKi.csv",header=T,row.names=1)
colours=ifelse(hm$chr21.geneset=="up","red","blue")
hmcols=colorRampPalette(c('blue','white','red'))(199)
breaks=seq(-8,8,length=200)
heatmap.2(t(hm[,3:5]), col=hmcols, breaks=breaks, Rowv=F, Colv=F, dendrogram="none", scale="none",trace="none",density.info="none",xlab="",key.title="", labCol=rownames(hm),labRow="",key=T,ColSideColors = colours)