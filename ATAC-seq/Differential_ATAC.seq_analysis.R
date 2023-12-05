#This is the code for genome-wide paired analysis of differential ATAC-seq data (peak-based) between ETS2-edited ("ETS2") and unedited ("NTC") inflammatory macrophages, used for figure S2 
#Adapted from https://github.com/reskejak/ATAC-seq. Raw data available from EGA. 
library(stringr)
library(GenomicRanges)
library(edgeR)
library(statmod)
bams.atac=dir(pattern="*.bam")#assumes you have only have relevant BAMs and index files in your working directory - otherwise need to select specific files for analysis
bams.atac=bams.atac[!str_detect(bams.atac,pattern="bai")]
ets2.d1.peaks=read.table("d1_ETS2_S5_L001_macs2_peaks.narrowPeak", sep="\t")[,1:3]
ets2.d2.peaks=read.table("d2_ETS2_S9_L001_macs2_peaks.narrowPeak", sep="\t")[,1:3]
ets2.d3.peaks=read.table("d3_ETS2_S11_L001_macs2_peaks.narrowPeak", sep="\t")[,1:3]
ntc.d1.peaks=read.table("d1_NTC_S4_L001_macs2_peaks.narrowPeak", sep="\t")[,1:3]
ntc.d2.peaks=read.table("d2_NTC_S8_L001_macs2_peaks.narrowPeak", sep="\t")[,1:3]
ntc.d3.peaks=read.table("d3_NTC_S10_L001_macs2_peaks.narrowPeak", sep="\t")[,1:3]
colnames(ets2.d1.peaks)=c("chrom","start","end")
colnames(ets2.d2.peaks)=c("chrom","start","end")
colnames(ets2.d3.peaks)=c("chrom","start","end")
colnames(ntc.d1.peaks)=c("chrom","start","end")
colnames(ntc.d2.peaks)=c("chrom","start","end")
colnames(ntc.d3.peaks)=c("chrom","start","end")
ets2.d1.peaks=GRanges(ets2.d1.peaks)
ets2.d2.peaks=GRanges(ets2.d2.peaks)
ets2.d3.peaks=GRanges(ets2.d3.peaks)
ntc.d3.peaks=GRanges(ntc.d3.peaks)
ntc.d2.peaks=GRanges(ntc.d2.peaks)
ntc.d1.peaks=GRanges(ntc.d1.peaks)
first.2.ets2=intersect(ets2.d1.peaks,ets2.d2.peaks)
ets2.peaks=intersect(first.2.ets2,ets2.d3.peaks)
first.2.ntc=intersect(ntc.d1.peaks,ntc.d2.peaks)
ntc.peaks=intersect(first.2.ntc,ntc.d3.peaks)
all.peaks=intersect(ets2.peaks,ntc.peaks)
blacklist=read.table("hg19-blacklist.v2.bed",sep="\t") #1-based
colnames(blacklist)=c("chrom","start","end")
blacklist=GRanges(blacklist)
standard.chr=paste0("chr",c(1:22,"X","Y"))
param=readParam(max.frag=1000,pe="both",discard=blacklist,restrict=standard.chr)
peak.counts=regionCounts(bams.atac,regions=all.peaks,param=param)
peak.abundances=aveLogCPM(asDGEList(peak.counts))
peak.counts.filt=peak.counts[peak.abundances>-1,] #equivalent to cpm > 0.5
binned=windowCounts(bams.atac,bin=T,width=10000,param=param)
macs.tmm=peak.counts.filt
macs.tmm=normFactors(binned, se.out=macs.tmm) #TMM normalisation
macs=asDGEList(macs.tmm)
colnames(macs$counts)=c("ETS2.d1","NTC.d1","ETS2.d2","NTC.d2","ETS2.d3","NTC.d3")
rownames(macs$samples)=c("ETS2.d1","NTC.d1","ETS2.d2","NTC.d2","ETS2.d3","NTC.d3")
macs$samples$group=rep(c("ets2","ntc"),3)
macs$samples$donor=rep(c("d1","d2","d3"),each=2)
macs.design=model.matrix(~0+group+donor, data=macs$samples)
colnames(macs.design)=c("ets2","ntc","d2","d3")
macs.fit=estimateDisp(macs,macs.design)
macs.fit2=glmQLFit(macs.fit,macs.design,robust=T)
macs.res=glmQLFTest(macs.fit2, contrast=makeContrasts(ets2-ntc,levels=macs.design))
rowData(macs.tmm)=cbind(rowData(macs.tmm),macs.res$table)
merged.macs=mergeWindows(rowRanges(macs.tmm),tol=500L, max.width=5000L)
tab.best=getBestTest(merged.macs$ids, macs.res$table)
final.merged.atac=GRanges(cbind(as.data.frame(merged.macs$regions), macs.res$table[tab.best$rep.test, -4], tab.best[,-c(7:8)]))
final.merged.atac@elementMetadata=cbind(final.merged.atac@elementMetadata,tab.best[,-1])
final.merged.atac=final.merged.atac[order(final.merged.atac@elementMetadata$FDR),]
final.merged.atac$sig="n.s."
plot(final.merged.atac$logFC,-log10(final.merged.atac$PValue),ylim=c(0,7),xlim=c(-1,1.2),pch=19,col="grey",cex=1.2)
