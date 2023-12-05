# R analysis of ETS2 or revETS2 overexpression (250ng or 500ng) in M0 macrophages (figure 3)
library(limma)
library(edgeR)
library(DESeq2)
library(fgsea)
library(data.table)
library(ggplot2)
counts=data.matrix(read.table(file="raw_counts.txt",sep='\t',row.names=1,header=T))
dose=rep(c("d250ng","d500ng"),16)
rna=rep(rep(c("ets2","rev"),each=2),8)
donor=rep(c(1,2,3,4,5,6,7,8),each=4)
names=paste0(donor,"_",rna,"_",dose)
sample=paste0(rna,"_",dose)
pData=cbind(names,sample,donor,dose,rna)
pData=as.data.frame(pData)
rownames(pData)=names
counts.matrix=as.matrix(counts)
colnames(counts.matrix)=names
counts.eset=new("ExpressionSet",exprs=counts.matrix)
pData(counts.eset)=pData
macs.expressed=rowSums(cpm(counts.matrix)>0.5)>=32
counts.exp=counts.matrix[macs.expressed,]
dge=DGEList(counts=counts.exp)
dge=calcNormFactors(dge)
dge$samples$donor=donor
dge$samples$rna=rna
dge$samples$dose=dose
dge$samples$sample=sample
split.design=model.matrix(~0+factor(dge$samples$sample)+factor(dge$samples$donor))
colnames(split.design)=c("ets2_250","ets2_500","rev_250","rev_500","s2","s3","s4","s5","s6","s7","s8")
split.contrasts=makeContrasts(ets2_250=ets2_250-rev_250,ets2_500=ets2_500-rev_500,levels=colnames(split.design))
split.v=voom(dge,split.design,plot=T)
split.fit=lmFit(split.v,split.design)
split.fit2=contrasts.fit(split.fit,split.contrasts)
split.fit3=eBayes(split.fit2)
summary(decideTests(split.fit3))
design=model.matrix(~0+factor(dge$samples$rna)+factor(dge$samples$donor)+factor(dge$samples$dose))
colnames(design)=c('ets2','rev','s2','s3','s4','s5','s6','s7','s8','d500ng')
contrasts=makeContrasts(ETS2vsREV=ets2-rev,levels=colnames(design))
v=voom(dge,design,plot=T)
fit=lmFit(v,design)
fit2=contrasts.fit(fit,contrasts)
fit3=eBayes(fit2)
summary(decideTests(fit3))
res.table.1=topTable(fit3,adjust.method="BH",p=1,number=length(counts.exp[,1]),lfc=0)
res.table.2=topTable(split.fit3,adjust.method="BH",p=1,number=length(counts.exp[,1]),lfc=0,coef=1)
res.table.3=topTable(split.fit3,adjust.method="BH",p=1,number=length(counts.exp[,1]),lfc=0,coef=2)
#write.csv(res.table.1,file="ETS2_overexpression_DEGs_250and500.csv")
#write.csv(res.table.2,file="ETS2_250ng_overexpression_DEGs.csv")
#write.csv(res.table.3,file="ETS2_500ng_overexpression_DEGs.csv")

# fGSEA of selected macrophage pathways, including GIMATS inflammatory macrophage signature (figure 3)
genelists=read.csv(file="macrophage_genesets_fgsea.csv",header=T) # gene sets with genes listed as Gene Symbols
lists=lapply(genelists,function(col)col[!is.na(col)])
oe_500ng=read.csv(file="overexpression_500ng_rank.csv",header=F,row.names=1) # list of genes ranked by t statistic from limma analysis
oe_500ng_num=unlist(oe_500ng)
names(oe_500ng_num)=rownames(oe_500ng)
fgsea_500ng=fgsea(pathways=lists, stats=oe_500ng_num, eps=0.0, minSize=15, maxSize=1000)
fgsea_500ng[order(pval), ]
#fwrite(fgsea_500ng,"fgsea_overexpression_500ng_mac.pathways.csv")
oe_250ng=read.csv(file="overexpression_250ng_rank.csv",header=F,row.names=1) # list of genes ranked by t statistic from limma analysis
oe_250ng_num=unlist(oe_250ng)
names(oe_250ng_num)=rownames(oe_250ng)
fgsea_250ng=fgsea(pathways=lists, stats=oe_250ng_num, eps=0.0, minSize=15, maxSize=1000)
fgsea_250ng[order(pval), ]
#fwrite(fgsea_250ng,"fgsea_overexpression_250ng_mac.pathways.csv")

# Standard GSEA plot:
plotEnrichment(lists[["GIMATS"]],oe_500ng_num)

# GSEA dotplot:
fgsea_500ng$pathway = factor(fgsea_500ng$pathway, levels=fgsea_500ng[order(fgsea_500ng$order), "pathway"])
fgsea_250ng$pathway = factor(fgsea_250ng$pathway, levels=fgsea_250ng[order(fgsea_250ng$order), "pathway"])
fgsea_500ng$neg.log10.P = -log10(fgsea_500ng$pval)
fgsea_250ng$neg.log10.P = -log10(fgsea_250ng$pval)
mapping=aes_string(x="NES",y="pathway",color="NES")
low="midnightblue"
high="darkred"
mid="white"
limits=c(-3,3)
ggp=ggplot(fgsea_500ng, aes(x=NES, y=pathway, size=log10.P, color=NES))
ggp+ geom_point()+scale_size_area(limits=c(0.1,5.2),max_size = 10)+scale_color_gradient2(low = low, high = high, mid=mid, aesthetics = 'color', limits = limits)+ theme(axis.ticks.y=element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + xlim(-2,2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggp2=ggplot(fgsea_250ng, aes(x=NES, y=pathway, size=log10.P., color=NES))
ggp2+ geom_point()+scale_size_area(limits=c(0.1,5.2),max_size = 10)+scale_color_gradient2(low = low, high = high, mid=mid, aesthetics = 'color', limits = limits)+ theme(axis.ticks.y=element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + xlim(-2,2) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

# fGSEA of disease signatures from macrophages obtained from patients with the respective diseases and segment plot (figure 3) - sources of disease signatures provided in Supplementary Data of our manuscript 
genelists=read.csv(file="disease_genesets.csv",header=T,na.strings=c("")) # gene sets with genes listed as Gene Symbols
lists=lapply(genelists,function(col)col[!is.na(col)])
fgsea_500ng.disease=fgsea(pathways=lists, stats=oe_500ng_num, eps=0.0, minSize=15, maxSize=1000)
fgsea_500ng.disease$neg.log.P=-log10(fgsea_500ng.disease$padj)
fgsea_500ng.disease$type = c(rep("chr21q22_disease",4),"other autoimmune disease","bacterial infection","bacterial infection","viral infection","viral infection","atherosclerosis","tumour associated","tumour associated") # annotate each disease based on category
p = ggplot(fgsea_500ng.disease,aes(fct_reorder(pathway,order),neg.log10.P)) + 
  geom_bar(aes(y=7,fill=type),stat="identity",width=1,colour="white",alpha=0.1) + 
  geom_bar(stat="identity",width=1,aes(fill=type),colour="white") + 
  coord_polar() +
  scale_fill_manual(values=c("chr21q22 disease" = "firebrick2",                                   
                             "bacterial infection" = "orange",
                             "viral infection"="purple",
                             "atherosclerosis"="springgreen4",
                             "tumour associated"="dodgerblue2")) +  
   scale_y_continuous(limits = c(-1,7))+  
  geom_hline(yintercept=1.32, linetype="dashed",color="black",size=0.5) +  #equivalent to FDR 0.05
  geom_hline(yintercept=yticks, linetype=1,color="white",size=0.5) +
  geom_label(aes(label=round(NES,2),fill=type),size=2,color="white",show.legend = FALSE)+
  theme_minimal() +                                                                     
  theme(legend.position = "top",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
