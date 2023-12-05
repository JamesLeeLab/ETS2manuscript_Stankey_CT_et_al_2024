# R analysis of colonic biopsies treated with MEKi (figure 5)
suppressMessages(library(affy))
suppressMessages(library(edgeR))
suppressMessages(library(fgsea))
suppressMessages(library(GSVA))
suppressMessages(library(preprocessCore))
suppressMessages(library(limma))
counts=data.matrix(read.table(file="raw_counts.txt",sep='\t',row.names=1,header=T))
eset=new('ExpressionSet',exprs=counts)
eset.expressed=rowSums(cpm(eset)>1)>=27
eset.exp=eset[eset.expressed,]
cpm=cpm(counts)
cpm.filt=cpm[eset.expressed,]
cpm.filt.ensg=new("ExpressionSet",exprs=cpm.filt)
cpm.filt.ensg.qn=normalize.quantiles(exprs(cpm.filt.ensg))
rownames(cpm.filt.ensg.qn)=featureNames(cpm.filt.ensg)
colnames(cpm.filt.ensg.qn)=colnames(counts)
log.norm=log2(cpm.filt.ensg.qn)
norm=new("ExpressionSet",exprs=log.norm)
genelists=read.csv(file="bMIS_score_ETS2_ENSG.csv",header=T)#mapped Ensembl IDs for genes in bMIS scores and ETS2-regulated genes (chr21q22 deletion)
lists=lapply(genelists,function(col)col[!is.na(col)])
gsva.norm=gsva(norm,lists,verbose=F)
#write.csv(gsva.norm,file="bMIS_final.res.csv")
dge=DGEList(counts=eset.exp)
dge=calcNormFactors(dge)
dge$samples$treatment=rep(c("d","m","i"),9)
dge$samples$donor=rep(c(1:9),each=3)
design=model.matrix(~0+factor(dge$samples$treatment)+factor(dge$samples$donor))
colnames(design)
colnames(design)=c("d","i","m","s2","s3","s4","s5","s6","s7","s8","s9")
contrasts=makeContrasts(ifx=i-d,mek=m-d, levels=colnames(design))
voom=voom(dge,design,plot=T)
fit=lmFit(voom,design)
fit2=contrasts.fit(fit,contrasts)
fit3=eBayes(fit2)
MEKi.res=topTable(fit3,number=8902,p.value=1,lfc=0,adjust.method="BH",coef=2)
biopsy_MEK=MEKi.res[,"t"]
names(biopsy_MEK)=rownames(MEKi.res)
biopsy_MEK=sort(biopsy_MEK, decreasing=T)
ets2lists=read.csv(file="ets2_genesetsENSG_MASTER.csv",header=T,na.strings=c(""))
ets2.lists=lapply(ets2lists,function(col)col[!is.na(col)])
fgsea_biopsy_MEK=fgsea(pathways=ets2.lists, stats=biopsy_MEK, eps=0.0, minSize=15, maxSize=1000)
fgsea_biopsy_MEK
plotEnrichment(ets2.lists[["CHR21_DN"]],biopsy_MEK)
plotEnrichment(ets2.lists[["ETS2g1_DN"]],biopsy_MEK)
plotEnrichment(ets2.lists[["ETS2g2_DN"]],biopsy_MEK)






