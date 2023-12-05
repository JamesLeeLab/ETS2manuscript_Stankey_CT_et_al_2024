# In our manuscript, both paired and unpaired differential ChIP-seq analysis are performed (depending on whether the comparator sample was derived from the same donor or not). 

# 1. Unpaired analysis of H3K27ac ChIP-seq data to assess the effect of rs2836882 genotype on enhancer activity at the extended chr21q22 locus
# This analysis was performed using MEDIPs - the software tool that performed best for differential analysis of sharp peak data (such as H3K27ac ChIP-seq) according to a recent benchmarking study (Eder and Grebien. Genome Biology 2022). 
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
maj1.h3="15_MAJ_TPP_H3K27ac_S4_L004_R1_001_trimmed_hg19.bam" #major hom rep1
maj1.in="15_MAJ_TPP_input_S8_L004_R1_001_trimmed_hg19.bam" #major hom rep1 input
min1.h3="19_MIN_TPP_H3K27ac_S3_L004_R1_001_trimmed_hg19.bam" #minor hom rep1
min1.in="19_MIN_TPP_input_S7_L004_R1_001_trimmed_hg19.bam" #minor hom rep1 input
maj2.h3="24_MAJ_TPP_H3K27ac_S12_L005_R1_001_trimmed_hg19.bam" #major hom rep2
maj2.in="24_MAJ_TPP_input_S16_L005_R1_001_trimmed_hg19.bam" #major hom rep2 input
min2.h3="30_MIN_TPP_H3K27ac_S11_L005_R1_001_trimmed_hg19.bam" #minor hom rep2
min2.in="30_MIN_TPP_input_S15_L005_R1_001_trimmed_hg19.bam" #minor hom rep2 input
BSgenome="BSgenome.Hsapiens.UCSC.hg19"
uniq=0
extend=300
shift=0
chr.select="chr21"
roi.wide=data.frame(chr="chr21",start=40150000,end=40710000,name="chr21q22.wide")
maj1.h3.roi=MEDIPS.createROIset(file=maj1.h3,BSgenome=BSgenome,extend=extend,shift=shift,uniq=uniq,ROI=roi.wide,chr.select=chr.select,bn=112)#this produces 112 x 5kb bins across the extended chr21q22 locus
maj2.h3.roi=MEDIPS.createROIset(file=maj2.h3,BSgenome=BSgenome,extend=extend,shift=shift,uniq=uniq,ROI=roi.wide,chr.select=chr.select,bn=112)
min2.h3.roi=MEDIPS.createROIset(file=min2.h3,BSgenome=BSgenome,extend=extend,shift=shift,uniq=uniq,ROI=roi.wide,chr.select=chr.select,bn=112)
min1.h3.roi=MEDIPS.createROIset(file=min1.h3,BSgenome=BSgenome,extend=extend,shift=shift,uniq=uniq,ROI=roi.wide,chr.select=chr.select,bn=112)
min1.in.roi=MEDIPS.createROIset(file=min1.in,BSgenome=BSgenome,extend=extend,shift=shift,uniq=uniq,ROI=roi.wide,chr.select=chr.select,bn=112)
min2.in.roi=MEDIPS.createROIset(file=min2.in,BSgenome=BSgenome,extend=extend,shift=shift,uniq=uniq,ROI=roi.wide,chr.select=chr.select,bn=112)
maj2.in.roi=MEDIPS.createROIset(file=maj2.in,BSgenome=BSgenome,extend=extend,shift=shift,uniq=uniq,ROI=roi.wide,chr.select=chr.select,bn=112)
maj1.in.roi=MEDIPS.createROIset(file=maj1.in,BSgenome=BSgenome,extend=extend,shift=shift,uniq=uniq,ROI=roi.wide,chr.select=chr.select,bn=112)
maj.h3.roi.wide=c(maj1.h3.roi,maj2.h3.roi)
min.h3.roi.wide=c(min1.h3.roi,min2.h3.roi)
min.in.roi.wide=c(min1.in.roi,min2.in.roi)
maj.in.roi.wide=c(maj1.in.roi,maj2.in.roi)
resultTable.roi.wide = MEDIPS.meth(MSet1 = maj.h3.roi.wide, MSet2 = min.h3.roi.wide, ISet1 = maj.in.roi.wide, ISet2 = min.in.roi.wide, chr = 'chr21', p.adj='BH', diff.method='edgeR', CNV=FALSE, MeDIP=FALSE, diffnorm='quantile',minRowSum=250)
sig = MEDIPS.selectSig(results=resultTable.roi.wide, p.value=0.05, adj=TRUE, ratio=NULL, bg.counts=NULL, CNV=FALSE)
resultTable.roi.wide$neg.log.p=-log10(resultTable.roi.wide$edgeR.adj.p.value)

# 2. Genome-wide paired analysis of differential H3K27ac ChIP-seq data (peak-based) between ETS2-edited and unedited inflammatory macrophages - code adapted from https://github.com/reskejak/ATAC-seq. Raw data available from EGA. 
library(GenomicRanges)
library(edgeR)
library(statmod)
library(stringr)
bams.h3k27.ko=dir(pattern="*.bam")#assumes you have only have relevant BAMs and index files in your working directory - otherwise need to select specific files for analysis 
bams.h3k27.ko=bams.h3k27.ko[!str_detect(bams.h3k27.ko,pattern="bai")]
d106.ets2.peaks=read.table("H3K27ac_TPP_NCI106ETS2_KO_hg19_peaks.narrowPeak", sep="\t")[,1:3]
d107.ets2.peaks=read.table("H3K27ac_TPP_NCI107ETS2_KO_hg19_peaks.narrowPeak", sep="\t")[,1:3]
d109.ets2.peaks=read.table("H3K27ac_TPP_NCI109ETS2_KO_hg19_peaks.narrowPeak", sep="\t")[,1:3]
d106.ntc.peaks=read.table("H3K27ac_TPP_NCI106_NTC_hg19_peaks.narrowPeak", sep="\t")[,1:3]
d107.ntc.peaks=read.table("H3K27ac_TPP_NCI107_NTC_hg19_peaks.narrowPeak", sep="\t")[,1:3]
d109.ntc.peaks=read.table("H3K27ac_TPP_NCI109_NTC_hg19_peaks.narrowPeak", sep="\t")[,1:3]
colnames(d107.ets2.peaks)=c("chrom","start","end")
colnames(d106.ets2.peaks)=c("chrom","start","end")
colnames(d109.ets2.peaks)=c("chrom","start","end")
colnames(d109.ntc.peaks)=c("chrom","start","end")
colnames(d107.ntc.peaks)=c("chrom","start","end")
colnames(d106.ntc.peaks)=c("chrom","start","end")
d107.ets2.peaks.ko=GRanges(d107.ets2.peaks)
d106.ets2.peaks.ko=GRanges(d106.ets2.peaks)
d109.ets2.peaks.ko=GRanges(d109.ets2.peaks)
d106.ntc.peaks.ko=GRanges(d106.ntc.peaks)
d109.ntc.peaks.ko=GRanges(d109.ntc.peaks)
d107.ntc.peaks.ko=GRanges(d107.ntc.peaks)
ntc.1.peaks.ko=intersect(d107.ntc.peaks.ko,d106.ntc.peaks.ko)
ntc.peaks.ko=intersect(ntc.1.peaks.ko,d109.ntc.peaks.ko)
ets2.1.peaks.ko=intersect(d107.ets2.peaks.ko,d106.ets2.peaks.ko)
ets2.peaks.ko=intersect(ets2.1.peaks.ko,d109.ets2.peaks.ko)
all.peaks.ko=intersect(ntc.peaks.ko,ets2.peaks.ko)
blacklist=read.table("hg19-blacklist.v2.bed",sep="\t") #1-based
colnames(blacklist)=c("chrom","start","end")
blacklist=GRanges(blacklist)
standard.chr=paste0("chr",c(1:22,"X","Y"))
param=readParam(max.frag=1000,pe="both",discard=blacklist,restrict=standard.chr)
peak.counts.ko=regionCounts(bams.h3k27.ko,regions=all.peaks.ko,param=param)
binned.ko=windowCounts(bams.h3k27.ko,bin=T,width=10000,param=param)
macs.tmm.ko=peak.counts.ko
macs.tmm.ko=normFactors(binned.ko, se.out=macs.tmm.ko)#TMM normalisation
macs.ko=asDGEList(macs.tmm.ko)
rownames(macs.ko$samples)=c("d106.ets2","d106.ntc","d107.ets2","d107.ntc","d109.ets2","d109.ntc")
colnames(macs.ko$counts)=c("d106.ets2","d106.ntc","d107.ets2","d107.ntc","d109.ets2","d109.ntc")
macs.ko$samples$group=rep(c("ets2","ntc"),3)
macs.ko$samples$donor=rep(c("d1","d2","d3"),each=2)
design.ko=model.matrix(~0+group+donor, data=macs.ko$samples)
colnames(design.ko)=c("ets2","ntc","d2","d3")
macs.ko.fit=estimateDisp(macs.ko,design.ko)
macs.ko.fit2=glmQLFit(macs.ko.fit,design.ko,robust=T)
macs.ko.res=glmQLFTest(macs.ko.fit2, contrast=makeContrasts(ets2-ntc,levels=design.ko))
rowData(macs.tmm.ko)=cbind(rowData(macs.tmm.ko),macs.ko.res$table)
merged.macs.ko=mergeWindows(rowRanges(macs.tmm.ko),tol=500L, max.width=5000L)
tab.best.ko=getBestTest(merged.macs.ko$ids, macs.ko.res$table)
final.merged.macs.ko=GRanges(cbind(as.data.frame(merged.macs.ko$regions), macs.ko.res$table[tab.best.ko$rep.test, -4], tab.best.ko[,-c(7:8)]))
final.merged.macs.ko@elementMetadata=cbind(final.merged.macs.ko@elementMetadata,tab.best.ko[,-1])
final.merged.macs.ko=final.merged.macs.ko[order(final.merged.macs.ko@elementMetadata$FDR),]
final.merged.macs.ko$sig="n.s."
final.merged.macs.ko$sig[final.merged.macs.ko$FDR<0.1]="s"
write.csv(final.merged.macs,file="H3K27ac_-1_KO_NTC_MACS.peak.csv")
sig=final.merged.macs[final.merged.macs$FDR<0.1,]
plot(final.merged.macs$logFC,-log10(final.merged.macs$PValue),ylim=c(0,7),xlim=c(-1,1.2), pch=19,col="grey",cex=1.2)
points(sig$logFC,-log10(sig$PValue),pch=19,col="red",cex=1.2)


# 3. Paired analysis of differential H3K27ac ChIP-seq signal at the disease-associated chr21q22 ETS2 enhancer between ETS2-edited and unedited inflammatory macrophages - code adapted from https://github.com/reskejak/ATAC-seq. Raw data available from EGA. 
# Note, this peak-independent analysis specifically focuses on the enhancer locus at chr21:40,465,000-40,470,000 (hg19)
library(GenomicRanges)
library(edgeR)
library(statmod)
library(stringr)
bams.h3k27.ko=dir(pattern="*.bam")#assumes you have only have relevant BAMs and index files in your working directory - otherwise need to select specific files for analysis 
bams.h3k27.ko=bams.h3k27.ko[!str_detect(bams.h3k27.ko,pattern="bai")]
chr21q22=toGRanges("chr21:40,465,000-40,470,000",genome="hg19")
blacklist=read.table("hg19-blacklist.v2.bed",sep="\t") #1-based
colnames(blacklist)=c("chrom","start","end")
blacklist=GRanges(blacklist)
standard.chr=paste0("chr",c(1:22,"X","Y"))
param=readParam(max.frag=1000,pe="both",discard=blacklist,restrict=standard.chr)
peak.counts.chr21q22=regionCounts(bams.h3k27.ko,regions=chr21q22,param=param)
binned.ko=windowCounts(bams.h3k27.ko,bin=T,width=10000,param=param)
macs.tmm.chr21q22=peak.counts.chr21q22
macs.tmm.chr21q22=normFactors(binned.ko, se.out=macs.tmm.chr21q22)#TMM normalisation
macs.chr21q22=asDGEList(macs.tmm.chr21q22)
colnames(macs.chr21q22$counts)=c("d106.ets2","d106.ntc","d107.ets2","d107.ntc","d109.ets2","d109.ntc")
rownames(macs.chr21q22$samples)=c("d106.ets2","d106.ntc","d107.ets2","d107.ntc","d109.ets2","d109.ntc")
macs.chr21q22$samples$group=rep(c("ets2","ntc"),3)
macs.chr21q22$samples$donor=rep(c("d106","d107","d109"),each=2)
design.chr21q22=model.matrix(~0+group+donor, data=macs.chr21q22$samples)
colnames(design.chr21q22)=c("ets2","ntc","d107","d109")
macs.chr21q22.fit=estimateDisp(macs.chr21q22,design.chr21q22)
macs.chr21q22.fit2=glmQLFit(macs.chr21q22.fit,design.chr21q22,robust=T)
macs.chr21q22.res=glmQLFTest(macs.chr21q22.fit2, contrast=makeContrasts(ets2-ntc,levels=design.chr21q22))
cpm(macs.chr21q22.res$fitted.values, lib.size=peak.counts.chr21q22$totals)
macs.chr21q22.res$table


# Paired analysis of differential H3K27ac ChIP-seq signal at the disease-associated chr21q22 ETS2 enhancer between resting (M0) macrophages overexpressing ETS2 or control mRNA - code adapted from https://github.com/reskejak/ATAC-seq. Raw data available from EGA. 
# Note, this peak-independent analysis specifically focuses on the enhancer locus at chr21:40,465,000-40,470,000 (hg19)
library(GenomicRanges)
library(edgeR)
library(statmod)
library(stringr)
bams.h3k27.oe=dir(pattern="*.bam")#assumes you have only have relevant BAMs and index files in your working directory - otherwise need to select specific files for analysis 
bams.h3k27.oe=bams.h3k27.oe[!str_detect(bams.h3k27.oe,pattern="bai")]
chr21q22=toGRanges("chr21:40,465,000-40,470,000",genome="hg19")
blacklist=read.table("hg19-blacklist.v2.bed",sep="\t") #1-based
colnames(blacklist)=c("chrom","start","end")
blacklist=GRanges(blacklist)
standard.chr=paste0("chr",c(1:22,"X","Y"))
param=readParam(max.frag=1000,pe="both",discard=blacklist,restrict=standard.chr)
peak.counts.oe=regionCounts(bams.h3k27.oe,regions=chr21q22,param=param)
binned.oe=windowCounts(bams.h3k27.oe,bin=T,width=10000,param=param)
macs.tmm.oe=peak.counts.oe
macs.tmm.oe=normFactors(binned.oe, se.out=macs.tmm.oe)#TMM normalisation
macs.oe=asDGEList(macs.tmm.oe)
rownames(macs.oe$samples)=c("d1.ets2","d1.rev","d2.ets2","d2.rev","d3.ets2","d3.rev")
colnames(macs.oe$counts)=c("d1.ets2","d1.rev","d2.ets2","d2.rev","d3.ets2","d3.rev")
macs.oe$samples$group=rep(c("ets2","rev"),3)
macs.oe$samples$donor=rep(c("d1","d2","d3"),each=2)
design.oe=model.matrix(~0+group+donor, data=macs.oe$samples)
colnames(design.oe)=c("ets2","rev","d2","d3")
macs.oe.fit=estimateDisp(macs.oe,design.oe)
macs.oe.fit2=glmQLFit(macs.oe.fit,design.oe,robust=T)
macs.oe.res=glmQLFTest(macs.oe.fit2, contrast=makeContrasts(ets2-rev,levels=design.oe))
cpm(macs.oe.res$fitted.values, lib.size=peak.counts.oe$totals)
macs.oe.res$table