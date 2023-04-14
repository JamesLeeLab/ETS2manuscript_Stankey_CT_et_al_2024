# MPRA analysis of chr21q22 region containing 99% credible set SNPs in TPP macrophages
# Note, this region was included as part of a larger MPRA experiment involving other disease-associated loci that are not relevant to the chr21q22 manuscript. Results for oligos relating to the chr21q22 credible set SNPs were extracted from the larger dataset.
library(limma)
library(plotrix)
library(edgeR)
library(affy)
library(zoo)
library(les)
library(preprocessCore)
library(factoextra)
library(ggplot2)
library(gplots)
library(plyr)
library(magrittr)
library(data.table)
library(QuASAR)
library(GenomicRanges)
library(regioneR)
library(karyoploteR)
library(reshape2)

# preprocessing of raw barcode count data
raw=read.csv(file="all_barcodes_TPP_MPRA.csv",header=T,row.names=1,check.names=F)
pheno=raw[,1:10]
dim(raw)
tpp.matrix=as.matrix(raw[,11:22])
tpp.eset=new('ExpressionSet',exprs=tpp.matrix)
cpm.tpp=cpm(tpp.eset)
median.rna=apply(cpm.tpp[,1:8],1,median)
median.dna=apply(cpm.tpp[,9:12],1,median)
data=cbind(pheno,cpm.tpp)
data$median.rna=median.rna
data$median.dna=median.dna
filt=apply(data[,11:22],1,min)
data_filt_cpm1=data[!filt<1,]
tpp_element.level=aggregate(x=data_filt_cpm1[,11:24],by=list(unique.values=data_filt_cpm1$elementID),FUN=sum)
tpp_elements=tpp_element.level[,2:13]
rownames(tpp_elements)=tpp_element.level[,1]
tpp_elements.qn=normalize.quantiles(as.matrix(tpp_elements),copy=F)
rownames(tpp_elements.qn)=rownames(tpp_elements)
colnames(tpp_elements.qn)=colnames(tpp_elements)
dd <- dist(scale(t(tpp_elements.qn)), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
plot(hc) # dendrogram of RNA and DNA barcode count data for all samples

# PCA of samples and heatmap of correlation between donors
pca=prcomp(t(tpp_elements.qn))
status=c(rep('tpp',8),rep('dna',4))
pca$status=status
pairs(pca$x[,1:3],pch=20,cex=1.5,col=as.factor(pca$status))
pca.res=get_eigenvalue(pca)
fviz_eig(pca) # Scree plot of percentage of explained variance by PCs
cormat=round(cor(activity[,1:8]),2)
get_upper_tri <- function(cormat){
cormat[lower.tri(cormat)]<- NA
return(cormat)
}
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggheatmap=ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
geom_tile(color = "white")+
scale_fill_gradient2(low = "blue", high = "red", mid = "white",
midpoint = 0.5, limit = c(0,1), space = "Lab",
name="Pearson Correlation")+
theme_minimal()+
theme(axis.text.x = element_text(angle = 0,
size = 12),axis.text.y=element_text(angle=0,size=12))+
coord_fixed()
ggheatmap +
geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
theme(
axis.title.x = element_blank(),
axis.title.y = element_blank(),
panel.grid.major = element_blank(),
panel.border = element_blank(),
panel.background = element_blank(),
axis.ticks = element_blank(),
legend.justification = c(1, 0),
legend.position = c(0.6, 0.7),
legend.direction = "horizontal")+
guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
title.position = "top", title.hjust = 0.5))

# Enhancer tiling analysis 
dna=tpp_elements.qn[,grep("v",colnames(tpp_elements.qn))]
median.dna=apply(dna,1,median)
tpp=tpp_elements.qn[,grep("TPP",colnames(tpp_elements.qn))]
activity=apply(tpp,2,function(x)  ((x/median.dna)))
pvals=rep(NA,nrow(dna))
for(i in 1:nrow(activity)) pvals[i]=t.test(x=activity[i,],mu=1,alternative="greater")$p.value
activity=cbind(activity,pvals)
#write.csv(activity,file="element_activity.csv")
ets2=read.csv(file="ETS2_all.csv",header=T) # includes all enhancer tiles and credible SNP barcodes
activity_ets2=activity[rownames(activity)%in%ets2$ID,]
activity_ets2.1=cbind(apply(activity_ets2[,1:8],1,median))
log.activity=log2(activity_ets2.1[,1])
activity_ets2.1=cbind(activity_ets2.1,log.activity,activity_ets2[,9])
activity_ets2.1=as.data.frame(activity_ets2.1)
activity_ets2.1$elementID=rownames(activity_ets2.1)
colnames(activity_ets2.1)=c("FC","log.activity","p.value","elementID")
ets2_enh_pheno=pheno[pheno$elementID %in% rownames(activity_ets2.1),]
ets2_bp=as.data.frame(cbind(unique(ets2_enh_pheno$start),unique(ets2_enh_pheno$end),unique(ets2_enh_pheno$elementID)))
colnames(ets2_bp)=c("start","end","elementID")
ets2_enh=join(activity_ets2.1,ets2_bp)
les.input=cbind(ets2_enh$end,ets2_enh$p.value,21)
les.input=data.frame(les.input,stringsAsFactors = F)
les.input[] <- lapply(les.input, type.convert)
colnames(les.input)=c("pos","pval","chr")
les=Les(les.input$pos,les.input$pval)
win=300
les2=estimate(les,win)
threshold(les2) #to determine suitable theta threshold
res = les::regions(les2, limit=0.348781, verbose = TRUE)
region <- res["regions"]
region
borders <- c(region$start, region$end)
n=nrow(ets2_enh)
datalist = vector("list", length = n)
for (i in 1:n) {
  res=NULL
  res$i <- cbind(seq(ets2_enh[i,5],ets2_enh[i,6]),ets2_enh[i,"log.activity"])
  datalist[[i]]=res$i
}
big_data = do.call(rbind, datalist)
colnames(big_data)=c("coord","score")
df=aggregate(score~coord,data=big_data,median,na.rm=T)
win300slide10=rollapply(df,width=300,by=10,FUN=mean)
subset=win300slide10[win300slide10[,1]<40467100 & win300slide10[,1]>40465100,]
plot(x=subset[,1],y=subset[,2],ylim=c(-0.18,0.18),pch=19,cex=0.6,xlab="Chromosome 21",ylab="Expression Modulating Activity (log2)")
rect(region$start,-0.19,region$end,0.19,col='mistyrose',border=NA)#numbers from les for 300
points(x=subset[,1],y=subset[,2],pch=19,cex=0.6)
abline(h=0,lty=2,lwd=0.5)
ablineclip(v=40466570,y1=-0.15,y2=-0.1,lwd=3,col='red') #rs2836882 (highligted in red based on QuASAR-MPRA analysis results - code below)
ablineclip(v=40466468,y1=-0.15,y2=-0.1,lwd=3,col='black') #rs9808561
ablineclip(v=40466299,y1=-0.15,y2=-0.1,lwd=3,col='black') #rs2836881
ablineclip(v=40466744,y1=-0.15,y2=-0.1,lwd=1,col='black') #rs2836883
ablineclip(v=40465901,y1=-0.15,y2=-0.1,lwd=1,col='black') #rs4817987
ablineclip(v=40465534,y1=-0.15,y2=-0.1,lwd=1,col='black') #rs2836878
ablineclip(v=40465512,y1=-0.15,y2=-0.1,lwd=1,col='black') #rs4817986
raw2=raw[,-c(1:4,9,11:18)]
raw.filt=raw2[!filt<1,]
tpp_raw_aggregate=aggregate(x=raw.filt[,6:9],by=list(unique.values=raw.filt$elementID),FUN=sum)
pheno_unique=pheno[!duplicated(pheno$elementID),]
colnames(tpp_raw_aggregate)=c("elementID","V1","V2","V3","V4")
dna_coverage=merge(pheno_unique,tpp_raw_aggregate,by="elementID")
ets2=read.csv(file="ETS2_all.csv",header=T) #includes all enhancer tiles and credible SNP barcodes
ets2_elements=ets2[1:116,1]
ets2_cov=dna_coverage[dna_coverage$elementID %in% ets2_elements,]
ets2_cov$median=round(apply(ets2_cov[,11:14],1,median))
n=length(ets2_cov$median)
dat=vector("list",length=n)
for (i in 1:n) {
  res2=NULL
  res2$i <- rep(ets2_cov[i,]$localcoords,each=ets2_cov[i,]$median)
  dat[[i]]=res2$i
}
cov_list=unlist(dat)
cov=toGRanges(cov_list,genome="hg19")
snp.region=toGRanges("chr21:40465100-40467100", genome="hg19")
kp <- plotKaryotype(genome="hg19",plot.type=5, zoom = snp.region)
kpPlotCoverage(kp, data=cov,r0=0.0,r1=0.4,ymax=16000)
#kpPlotRegions(kp, data=cov, r0=0.48,r1=0,layer.margin = 0.01)
kpAxis(kp, ymin=0,numticks = 4, ymax=16000, r0=0.00, r1=0.4,tick.pos = c(0,5000,10000,15000))

# Expression-modulating variant analysis
tpp1=tpp_elements.qn[,grep("TPP",colnames(tpp_elements.qn))]
dna=tpp_elements.qn[,grep("v",colnames(tpp_elements.qn))]
median.dna=apply(dna,1,median)
tpp2=cbind(tpp1,median.dna)
tpp_snps=tpp2[-grep("Enhancer",rownames(tpp2)),]
ref=as.data.frame(tpp_snps[grep("ref",rownames(tpp_snps)),])
ref$name=rownames(ref)
out <- strsplit(as.character(ref$name),'_re')
ref=data.frame(ref,do.call(rbind,out))
ref$X2='ref'
new_names=gsub("X","allele2_",colnames(ref))
colnames(ref)=new_names
alt=as.data.frame(tpp_snps[grep("alt",rownames(tpp_snps)),])
alt$name=rownames(alt)
out <- strsplit(as.character(alt$name),'_al')
alt=data.frame(alt,do.call(rbind,out))
alt$X2='alt'
new_names=gsub("X","allele1_",colnames(alt))
colnames(alt)=new_names
mult=as.data.frame(tpp_snps[grep("mult",rownames(tpp_snps)),])
mult$name=rownames(mult)
out <- strsplit(as.character(mult$name),'_mu')
mult=data.frame(mult,do.call(rbind,out))
mult$X2='ref'
new_names=gsub("X","allele1_",colnames(mult))
colnames(mult)=new_names
alt_mult=rbind(alt,mult)
colnames(alt_mult)[11]="elementID"
colnames(alt_mult)[10]="alt_name"
colnames(alt_mult)[9]="alt_median.DNA"
colnames(alt_mult)[12]="alt_mult"
colnames(ref)[11]="elementID"
colnames(ref)[10]="ref_name"
colnames(ref)[9]="ref_median.DNA"
colnames(ref)[12]="ref"
tpp_quasar=join(ref,alt_mult,type='full')
#write.csv(tpp_quasar,file="PRE_8_sample_TPP_filt.csv")#this is for longhand removal of missing barcodes, but can also be done using MPRAeasieR.
tpp=read.csv(file="POST_8_sample_TPP_filt.csv", header=T,row.names=1)
tpp$DNA_A=tpp$alt_median.DNA
tpp$DNA_R=tpp$ref_median.DNA
tpp$DNA_prop=tpp$DNA_R/(tpp$DNA_R+tpp$DNA_A)
tpp$logit_prop=log2(tpp$DNA_R/tpp$DNA_A)
tpp=tpp[c((ncol(tpp)-4):ncol(tpp),1:(ncol(tpp)-5))]
head(tpp)
tpp_16.res=fitQuasarMpra(tpp$allele2_16TPP,tpp$allele1_16TPP,prop=tpp$DNA_prop)
tpp_17.res=fitQuasarMpra(tpp$allele2_17TPP,tpp$allele1_17TPP,prop=tpp$DNA_prop)
tpp_18.res=fitQuasarMpra(tpp$allele2_18TPP,tpp$allele1_18TPP,prop=tpp$DNA_prop)
tpp_22.res=fitQuasarMpra(tpp$allele2_22TPP,tpp$allele1_22TPP,prop=tpp$DNA_prop)
tpp_23.res=fitQuasarMpra(tpp$allele2_23TPP,tpp$allele1_23TPP,prop=tpp$DNA_prop)
tpp_25.res=fitQuasarMpra(tpp$allele2_25TPP,tpp$allele1_25TPP,prop=tpp$DNA_prop)
tpp_35.res=fitQuasarMpra(tpp$allele2_35TPP,tpp$allele1_35TPP,prop=tpp$DNA_prop)
tpp_29.res=fitQuasarMpra(tpp$allele2_29TPP,tpp$allele1_29TPP,prop=tpp$DNA_prop)
tpp_16.res$w=1/(tpp_16.res$betas_se^2)
tpp_17.res$w=1/(tpp_17.res$betas_se^2)
tpp_18.res$w=1/(tpp_18.res$betas_se^2)
tpp_22.res$w=1/(tpp_22.res$betas_se^2)
tpp_23.res$w=1/(tpp_23.res$betas_se^2)
tpp_25.res$w=1/(tpp_25.res$betas_se^2)
tpp_35.res$w=1/(tpp_35.res$betas_se^2)
tpp_29.res$w=1/(tpp_29.res$betas_se^2)
tpp_16.res$beta_0=log(tpp$DNA_prop/(1-tpp$DNA_prop))
tpp_17.res$beta_0=log(tpp$DNA_prop/(1-tpp$DNA_prop))
tpp_18.res$beta_0=log(tpp$DNA_prop/(1-tpp$DNA_prop))
tpp_22.res$beta_0=log(tpp$DNA_prop/(1-tpp$DNA_prop))
tpp_23.res$beta_0=log(tpp$DNA_prop/(1-tpp$DNA_prop))
tpp_25.res$beta_0=log(tpp$DNA_prop/(1-tpp$DNA_prop))
tpp_35.res$beta_0=log(tpp$DNA_prop/(1-tpp$DNA_prop))
tpp_29.res$beta_0=log(tpp$DNA_prop/(1-tpp$DNA_prop))
tpps=ls()[ls() %>% grep("^tpp_[0-9]",.)]
processQuasar <- function(x){
  rep <- gsub("tpp_([0-9]+).*","\\2",x)
  treat <- gsub("tpp_([0-9]+).*","\\1",x)
  tmp <- get(x) %>% data.table
  tmp[,rep:=rep]
  tmp[,c('rep','assay','treatment'):=list(rep,rownames(tpp),treat)]
  tmp
}
alltpps <- lapply(tpps,processQuasar) %>% rbindlist
plot(alltpps$betas_z,(alltpps$betas.beta.binom-alltpps$beta_0)/alltpps$betas_se,main="Stim",xlab="QuASAR Z",ylab="Computed Z") #sense check
alltpps[,beta_adj:=betas.beta.binom-beta_0]
meta.tpp <- alltpps[,list(sb=sum(beta_adj * w),sw=sum(w)),by='assay']
meta.tpp[,meta_beta:=sb/sw]
meta.tpp[,Z:=meta_beta/(1/sqrt(sw))]
FDR_THRESH <- 0.05
meta.tpp[,meta.p:=pnorm(abs(Z),lower.tail = FALSE) * 2]
meta.tpp[,meta.p.adj:=p.adjust(meta.p,method="fdr")]
meta.tpp.sig.count <- (meta.tpp$meta.p.adj<FDR_THRESH) %>% sum
tpp.res_quasar=cbind(meta.tpp$meta.p,meta.tpp$meta.p.adj,meta.tpp$meta_beta)
rownames(tpp.res_quasar)=rownames(tpp)
colnames(tpp.res_quasar)=c('meta.p','fdr.p','meta.beta')
tpp.res_quasar=tpp.res_quasar[rownames(tpp.res_quasar) %in% ets2$ID,]