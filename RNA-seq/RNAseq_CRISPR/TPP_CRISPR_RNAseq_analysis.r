# R analysis of ETS2 RNAseq data from primary human macrophages. Expression filter applied to all samples, then analyses on selected sample groups. 
# Used for figure 2
# NOTE: the same (edited) script was used to analyse the raw counts from the "Chr21q22 disease datasets" i.e. IBD macrophages (GSE123141), primary sclerosing cholangitis liver (GSE159676) and ankylosing spondylitis synovium (GSE41038) and then perform fGSEA analysis (supplementary figure 7)
suppressMessages(library(limma))
suppressMessages(library(affy))
suppressMessages(library(edgeR))
suppressMessages(library(biomaRt))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(plyr))
suppressMessages(library(fgsea))
counts=data.matrix(read.table(file="raw_counts.txt",sep='\t',row.names=1,header=T))
all.names=c('1_g1','1_g4','2_chr21','2_g1','2_g4','2_NTC','3_chr21','3_g1','3_g4','3_NTC','4_chr21','4_NTC','5_chr21','5_g1','5_g4','5_NTC','6_g1','6_g4','6_NTC','7_chr21','7_g1','7_g4','7_NTC','8_g1','8_g1_ROX','8_g4','8_g4_ROX','8_NTC','9_g1','9_g1_ROX','9_g4','9_g4_ROX','9_NTC','10_g1','10_g1_ROX','10_g4','10_g4_ROX','10_NTC')
colnames(counts)=all.names
data2=as.matrix(counts)
eset=new('ExpressionSet',exprs=data2)
eset.expressed=rowSums(cpm(eset)>0.5)>=38
eset.exp=eset[eset.expressed,]
reads=apply(counts,2,sum)
reads
ko=eset.exp[,c(1:2,4:6,8:10,12,14:19,21:24,26,28:29,31,33:34,36,38)]
ko_dge=DGEList(counts=ko)
ko_dge=calcNormFactors(ko_dge)
samples=c('g1','g2','g1','g2','NTC','g1','g2','NTC','NTC','g1','g2','NTC','g1','g2','NTC','g1','g2','NTC','g1','g2','NTC','g1','g2','NTC','g1','g2','NTC')
ko_dge$samples$target=samples
ko_dge$samples$donor=c(1,1,2,2,2,3,3,3,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10)
design=model.matrix(~0+factor(ko_dge$samples$target)+factor(ko_dge$samples$donor))
colnames(design)
colnames(design)=c('g1','g2','NTC','s2','s3','s4','s5','s6','s7','s8','s9','s10')
contrasts=makeContrasts(g1=g1-NTC,g2=g2-NTC, levels=colnames(design))
voom=voom(ko_dge,design,plot=T)
fit=lmFit(voom,design)
fit2=contrasts.fit(fit,contrasts)
fit3=eBayes(fit2)
summary(decideTests(fit3, p=0.05, adjust.method = 'BH'))
g1.res=topTable(fit3,number=12144,p.value=1,lfc=0,adjust.method="BH",coef=1)
g2.res=topTable(fit3,number=12144,p.value=1,lfc=0,adjust.method="BH",coef=2)
#write.csv(g1.res,file="ko.g1.res_final.csv")
#write.csv(g2.res,file="ko.g2.res_final.csv")
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensg = rownames(g1.res)
G_list = getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"),values=ensg,mart= mart)
g1.res.SYMBOL=g1.res[rownames(g1.res) %in% G_list$ensembl_gene_id,]
g1.res.SYMBOL$ensembl_gene_id=rownames(g1.res.SYMBOL)
g1.res_plot=join(g1.res.SYMBOL,G_list)
EnhancedVolcano(g1.res_plot,x="logFC",y="P.Value",lab=g1.res_plot$hgnc_symbol,pCutoff = 0.0066,legendLabels = F,title="",subtitle="", xlim=c(-4,4),ylim=c(0,15),labSize = 5,drawConnectors = T,max.overlaps=Inf, maxoverlapsConnectors=Inf, typeConnectors = 'open', FCcutoff=1,endsConnectors="last")

# Pathway analysis using fGSEA (figure 2)
genelists=read.csv(file="GOBP_genesets.csv",header=T) # csv of GO biological pathways (converted from gmt downloaded from MSigDB)
lists=lapply(genelists,function(col)col[!is.na(col)])
g1.fgsea=read.csv(file="ko.g1.res_finalSYMBOL_rank.csv",header=F,row.names=1)
g1.fgsea_num=unlist(g1.fgsea)
names(g1.fgsea_num)=rownames(g1.fgsea)
g1.fgsea_gobp=fgsea(pathways=lists, stats=g1.fgsea_num, eps=0.0, minSize=15, maxSize=1000)

# chr21q22 enhancer deletion analysis (figure 2)
chr21=eset.exp[,c(3,6,7,10,11,12,13,16,20,23)]
chr21_dge=DGEList(counts=chr21)
chr21_dge=calcNormFactors(chr21_dge)
chr21.samples=rep(c('chr21','NTC'),5)
chr21_dge$samples$target=chr21.samples
chr21_dge$samples$donor=rep(c("s2","s3","s4","s5","s7"),each=2)
chr21.design=model.matrix(~0+factor(chr21_dge$samples$target)+factor(chr21_dge$samples$donor))
colnames(chr21.design)=c('chr21','NTC','s3','s4','s5','s6')
chr21.contrasts=makeContrasts(chr21=chr21-NTC,levels=colnames(chr21.design))
voom=voom(chr21_dge,chr21.design,plot=T)
chr21.fit=lmFit(voom,chr21.design)
chr21.fit2=contrasts.fit(chr21.fit,chr21.contrasts)
chr21.fit3=eBayes(chr21.fit2)
summary(decideTests(chr21.fit3, p=0.05, adjust.method = 'BH'))
chr21.res=topTable(chr21.fit3,number=12144,p.value=1,lfc=0,adjust.method="BH")
#write.csv(chr21.res,file="chr21.res_final.csv")

# roxadustat-treated TPP analysis (figure 4)
rox=eset.exp[,c(24:38)]
rox_dge=DGEList(counts=rox)
rox_dge=calcNormFactors(rox_dge)
rox.samples=rep(c("g1","g1rox","g2","g2rox","ntc"),3)
rox_dge$samples$target=rox.samples
rox_dge$samples$donor=rep(c("s8","s9","s10"),each=5)
rox.design=model.matrix(~0+factor(rox_dge$samples$target)+factor(rox_dge$samples$donor))
colnames(rox.design)=c('g1','g1rox','g2','g2rox','ntc','s43',"s44")
rox.contrasts=makeContrasts(NTC_vs_g1=ntc-g1,NTC_vs_g1rox=ntc-g1rox,levels=colnames(rox.design))
voom=voom(rox_dge,rox.design,plot=T)
rox.fit=lmFit(voom,rox.design)
rox.fit2=contrasts.fit(rox.fit,rox.contrasts)
rox.fit3=eBayes(rox.fit2)
ntc.g1_rox.expt=topTable(rox.fit3,number=12144,p.value=1,lfc=0,adjust.method="BH",coef=1)
ntc.g1rox_rox.expt=topTable(rox.fit3,number=12144,p.value=1,lfc=0,adjust.method="BH",coef=2)
#write.csv(ntc.g1_rox.expt,file="ntc_vs_g1_rox.expt.csv")
#write.csv(ntc.g1rox_rox.expt,file="ntc_vs_g1rox_rox.expt.csv")
genelists=read.csv(file="macrophage_genesets_glycolysis_fgsea.csv",header=T)
lists=lapply(genelists,function(col)col[!is.na(col)])
ntc_g1=read.csv(file="ntc_vs_g1_rox.expt.res_final.rank.csv",header=F,row.names=1) # list of genes ranked by t statistic from limma analysis
ntc.g1=unlist(ntc_g1)
names(ntc.g1)=rownames(ntc_g1)
ntc_g1rox=read.csv(file="ntc_vs_g1rox_rox.expt.res_final.rank.csv",header=F,row.names=1) # list of genes ranked by t statistic from limma analysis
ntc.g1rox=unlist(ntc_g1rox)
names(ntc.g1rox)=rownames(ntc_g1rox)
fgsea_ntc.g1=fgsea(pathways=lists, stats=ntc.g1, eps=0.0, minSize=15, maxSize=1000)
fgsea_ntc.g1rox=fgsea(pathways=lists, stats=ntc.g1rox, eps=0.0, minSize=15, maxSize=1000)
fgsea_ntc.g1$neg.log.P=-log10(fgsea_ntc.g1$pval)
fgsea_ntc.g1rox$neg.log.P=-log10(fgsea_ntc.g1rox$pval)
lollipop=data.frame(x=fgsea_ntc.g1$pathway,g1=fgsea_ntc.g1$NES,g1rox=fgsea_ntc.g1rox$NES,order=c(11,9,8,7,10,6,5,4,3,1,2),g1_p=fgsea_ntc.g1$neg.log.P,g1rox_p=fgsea_ntc.g1rox$neg.log.P)
lollipop$x=factor(lollipop$x,levels=lollipop[order(lollipop$order),"x"])
ggplot(lollipop)+
  geom_segment( aes(x=x, xend=x, y=g1, yend=g1rox), color="grey") +
  geom_point( aes(x=x, y=g1), color=rgb(1,0.2,0.6,0.5), size=4 ) +
  geom_point( aes(x=x, y=g1rox), color=rgb(0,0.4,1,0.5), size=4 ) +
  coord_flip()+ theme_bw()+
  theme(
    legend.position = "none",
  ) +
  xlab("") +
  ylab("NES")