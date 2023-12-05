## Analysis of TPP macrophages ETS2 CUT&RUN data, with hg19, used for figure 4 and supplementary figures 6 and 8
## libraries sequenced on a NovaSeq S2 with 2x100 bp paired-end reads

#FASTQC on initial fastq files:
module load FastQC/0.11.8-Java-1.8
folder="path/to/fastq/files"

for a in *.fastq.gz;
do 
        fastqc $a \
        -o $folder;
done

#to run MultiQC on these fastQC files:
cd $folder
module load MultiQC/0.9-foss-2016b-Python-2.7.12
module load matplotlib/1.5.3-foss-2016b-Python-2.7.12
multiqc .

#trimming:
module load TrimGalore/0.6.0

a=(*.fastq.gz)
for ((i=0; i<${#a[@]}; i+=2));
do 
        trim_galore --phred33 -q 24 --illumina \
        --length 25 --paired --stringency 6 \
        "${a[i]}" "${a[i+1]}" \
        --fastqc -o trimmed;
done

#alignment with bowtie2 using the same parameters as Skene & Henikoff, 2017 (Elife):
cd trimmed
module load Bowtie2/2.4.4-GCC-11.2.0

for ((i=0; i<${#a[@]}; i+=2));
do 
        bowtie2 -x /Bowtie2Index/genome -1 "${a[i]}" -2 "${a[i+1]}" \
        --local --very-sensitive-local --no-mixed --no-discordant \
        --phred33 -I 10 -X 700 -S ${a[i]%.fq.gz}_hg19_bowtie2.sam;
done

#convert SAM to BAM
module load GCC/10.2.0
module load SAMtools/1.11-GCC-10.2.0

for a in *.sam; do samtools view -bS - > bowtie2bam_hg19/${a%.sam}.bam; done

#sort
cd bowtie2bam_hg19
for a in *.bam; do samtools sort -o ${a%.bam}_sorted.bam $a; done

##Analysis on individual biological replicates:
#merging technical replicates (same libraries re-sequenced):
samtools merge CRETS2_Ther1_bowtie2hg19_STmerged.bam BOU4950A4_S45_L001_R1_001_val_1_hg19_bowtie2_sorted.bam BOU4950A41_S56_L002_R1_001_val_1_hg19_bowtie2_sorted.bam 
samtools merge CRETS2_Ther2_bowtie2hg19_STmerged.bam BOU4950A12_S53_L001_R1_001_val_1_hg19_bowtie2_sorted.bam BOU4950A43_S58_L002_R1_001_val_1_hg19_bowtie2_sorted.bam
samtools merge CRIgG_Ther1_bowtie2hg19_STmerged.bam BOU4950A1_S42_L001_R1_001_val_1_hg19_bowtie2_sorted.bam BOU4950A25_S45_L002_R1_001_val_1_hg19_bowtie2_sorted.bam
samtools merge CRIgG_Ther2_bowtie2hg19_STmerged.bam BOU4950A9_S50_L001_R1_001_val_1_hg19_bowtie2_sorted.bam BOU4950A27_S47_L002_R1_001_val_1_hg19_bowtie2_sorted.bam

#resort
for a in *STmerged.bam; do samtools sort -o ${a%.bam}_resorted.bam $a; done

#indexing:
for a in *resorted.bam; do samtools index $a; done

#Mark duplicates and unmapped reads with Picard:
ml load  picard/2.1.1-Java-1.8.0_92

for a in *resorted.bam;
do 
        java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates INPUT=$a \
        OUTPUT=${a%.bam}_markDup.bam METRICS_FILE=${a%.bam}_markDup.txt \
        REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=./ ;
done

#loop to remove unmapped reads 
#Note: duplicates are marked but not removed for samples, as advised on https://nf-co.re/cutandrun/3.1/usage 

for a in *markDup.bam;
do 
        samtools view $a -b -h -F 4 -F 256 -F 2048 -q 15 \
        -o ${a%sorted.bam}rmunm.bam;
done

#then loop to sort:
for a in *rmunm.bam; do samtools sort -o ${a%.bam}_resorted.bam $a; done

#then loop to index:
for a in *resorted.bam; do echo $a; cat $a | samtools index $a; done

#Create BW files:
ml deepTools/3.3.1-foss-2020a-Python-3.8.2
mkdir bigwig

for a in *resorted.bam;
do 
        bamCoverage -b $a \
        -o bigwig/${a%.bam}.bw -of bigwig \
        --binSize 10 --normalizeUsing CPM;
done

#call peaks with MACS2:
module load MACS2/2.2.5-foss-2018b-Python-3.6.6

macs2 callpeak -t CRETS2_Ther1_bowtie2hg19_STmerged_resorted_markDup.bamrmunm_resorted.bam -c CRIgG_Ther1_bowtie2hg19_STmerged_resorted_markDup.bamrmunm_resorted.bam -g hs --qvalue 0.05 -f BAMPE -n CRETS2_Ther1_bowtie2hg19_STmerged_rmunm_q05 --keep-dup all -B --nomodel --outdir PEAKS05 
macs2 callpeak -t CRETS2_Ther2_bowtie2hg19_STmerged_resorted_markDup.bamrmunm_resorted.bam -c CRIgG_Ther2_bowtie2hg19_STmerged_resorted_markDup.bamrmunm_resorted.bam -g hs --qvalue 0.05 -f BAMPE -n CRETS2_Ther2_bowtie2hg19_STmerged_rmunm_q05 --keep-dup all -B --nomodel --outdir PEAKS05


##analysis merging technical AND biological replicates, to improve sensitivity for peak detection:
samtools merge CRETS2_Ther_bowtie2hg19_STmerged.bam BOU4950A4_S45_L001_R1_001_val_1_hg19_bowtie2_sorted.bam BOU4950A41_S56_L002_R1_001_val_1_hg19_bowtie2_sorted.bam BOU4950A12_S53_L001_R1_001_val_1_hg19_bowtie2_sorted.bam BOU4950A43_S58_L002_R1_001_val_1_hg19_bowtie2_sorted.bam
samtools merge CRIgG_Ther_bowtie2hg19_STmerged.bam BOU4950A1_S42_L001_R1_001_val_1_hg19_bowtie2_sorted.bam BOU4950A25_S45_L002_R1_001_val_1_hg19_bowtie2_sorted.bam BOU4950A9_S50_L001_R1_001_val_1_hg19_bowtie2_sorted.bam BOU4950A27_S47_L002_R1_001_val_1_hg19_bowtie2_sorted.bam

#resort
for a in *STmerged.bam; do samtools sort -o ${a%.bam}_resorted.bam $a; done

#indexing:
for a in *resorted.bam; do samtools index $a; done

#Mark duplicates and unmapped reads with Picard:
for a in *markDup.bam;
do 
        samtools view $a -b -h -F 4 -F 256 -F 2048 -q 15 \
        -o ${a%sorted.bam}rmunm.bam;
done

#then loop to sort:
for a in *rmunm.bam; do samtools sort -o ${a%.bam}_resorted.bam $a; done

#then loop to index:
for a in *resorted.bam; do echo $a; cat $a | samtools index $a; done

#Create BW files:
ml deepTools/3.3.1-foss-2020a-Python-3.8.2
mkdir bigwig

for a in *resorted.bam;
do 
        bamCoverage -b $a \
        -o bigwig/${a%.bam}.bw -of bigwig \
        --binSize 10 --normalizeUsing CPM;
done

#call peaks with MACS2:
module load MACS2/2.2.5-foss-2018b-Python-3.6.6
macs2 callpeak -t CRETS2_Ther_bowtie2hg19_STmerged_resorted_markDup.bamrmunm_resorted.bam -c CRIgG_Ther_bowtie2hg19_STmerged_resorted_markDup.bamrmunm_resorted.bam -g hs --qvalue 0.05 -f BAMPE -n CRETS2_Ther_bowtie2hg19_STmerged_rmunm_05 --keep-dup all -B --nomodel --outdir PEAKS05 

## Use R to find idr using two biological replicates of ETS2 C&R:
R
library(idr)
file1 <- "CRETS2_Ther1_bowtie2hg19_STmerged_rmunm_q05_peaks.narrowPeak"
file2 <- "CRETS2_Ther2_bowtie2hg19_STmerged_rmunm_q05_peaks.narrowPeak"
df1 <- read_delim(file1, delim="\t", col_names=FALSE)
df2 <- read_delim(file2, delim="\t", col_names=FALSE)
peak1 <- GRanges(df1$X1, IRanges(df1$X2, df1$X3), score=df1$X7)
peak2 <- GRanges(df2$X1, IRanges(df2$X2, df2$X3), score=df2$X7)
peak1 <- keepStandardChromosomes(peak1, pruning.mode="coarse")
peak2 <- keepStandardChromosomes(peak2, pruning.mode="coarse")
fo <- findOverlaps(peak1, peak2)
write.csv(fo,file="fo_dups_kept.csv")
y1 <- peak1$score[fo[,1]]
y2 <- peak2$score[fo[,2]]
dat <- cbind(log10(y1), log10(y2))
res <- est.IDR(dat, mu=3, sigma=1, rho=.9, p=.5)
df <- data.frame(rep1=dat[,1],rep2=dat[,2], rank1=rank(-dat[,1]),rank2=rank(-dat[,2]), idr=res$idr)
write.csv(df,file="idr_of_overlapping_peaks_dups_KEPT.csv")

# heatmap of C&R idr data overlaid on chromatin marks from TPP macrophages
require(EnrichedHeatmap)
require(rtracklayer)
require(circlize)
require(data.table)
require(magick)
require(Cairo)
h3k27ac=read.table("h3k27ac_coords.txt",sep="\t",header=F)
h3k27ac=GRanges(h3k27ac$V1,IRanges(h3k27ac$V2,h3k27ac$V3))
h3k27ac <- keepStandardChromosomes(h3k27ac, pruning.mode="coarse")
h3k4me3=read.table("h3k4me3_coords.txt",sep="\t",header=F) # these coordinates are based on analysis of raw data (using ChIP-seq code provided) from a published study of H3K4me3 ChIP-seq in TPP macrophages (GSE47188)
h3k4me3=GRanges(h3k4me3$V1,IRanges(h3k4me3$V2,h3k4me3$V3))
h3k4me3 <- keepStandardChromosomes(h3k4me3, pruning.mode="coarse")
atac=read.table("atac_peaks.txt",sep="\t",header=F)
atac=GRanges(atac$V1,IRanges(atac$V2,atac$V3))
atac <- keepStandardChromosomes(atac, pruning.mode="coarse")
ExtendSize <- 2000
atac.extended  <- resize(atac, fix = "center", width = ExtendSize*2)
h3k27ac.extended  <- resize(h3k27ac, fix = "center", width = ExtendSize*2)
h3k4me3.extended  <- resize(h3k4me3, fix = "center", width = ExtendSize*2)
BigWig <- rtracklayer::import("CRETS2_Ther_hg19_STmerged_rmunm.bw", format = "BigWig", selection = BigWigSelection(atac.extended))
#can repeat for histone marks
normMatrix <- normalizeToMatrix(signal = BigWig, target = resize(atac, fix = "center", width = 1), background = 0, keep = c(0, 0.99), target_ratio = 0, mean_mode = "w0", value_column = "score", extend = ExtendSize)
col_fun = circlize::colorRamp2(quantile(normMatrix, c(0, .99)), c("darkblue", "darkgoldenrod1"))
EH <- EnrichedHeatmap( mat = normMatrix, pos_line = FALSE, border = T, col = col_fun, column_title = "ETS2 CUT&RUN", column_title_gp = gpar(fontsize = 15, fontfamily = "sans"), use_raster = TRUE, raster_quality = 10, raster_device = "CairoPNG", rect_gp = gpar(col = "transparent"), heatmap_legend_param = list(legend_direction = "horizontal", title = "normalized counts"), top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = "black", lty = 1, lwd=2),col="black")))
x11()
draw(EH, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", padding = unit(c(4, 4, 4, 4), "mm"))
dev.off()

# Use R to find where ETS2 binds relative to known functional and spatial genomic annotations
R
library(GenomicRanges)
prom=read.table("promoters_from_coords.txt",sep = "\t",header=F)
#-1000bp to +100bp relative to TSS (hg19)
ets2=read.table("ets2_peaks.txt",sep="\t",header=F)
#idr peaks from ETS2 C&R
h3k27=read.table("h3k27ac_from_coords.txt",sep="\t",header=F)
#H3K27ac peaks from ChIP-seq in H3k27ac_TPP_hg19_minor_fromMergedBAM_peaks.narrowpeak
body=read.table("genebody_from_coords.txt",sep="\t",header=F)
#gene body coordinates (+100 from TSS to end of gene, hg19)
prom=GRanges(prom$V1,IRanges(prom$V2,prom$V3))
ets2=GRanges(ets2$V1,IRanges(ets2$V2,ets2$V3))
h3k27=GRanges(h3k27$V1,IRanges(h3k27$V2,h3k27$V3))
body=GRanges(body$V1,IRanges(body$V2,body$V3))
ets2_body=findOverlaps(ets2,body)
ets2_prom=findOverlaps(ets2,prom1)
ets2_h3k27=findOverlaps(ets2,h3k27)

# to find where ETS2 binds relative to known functional and spatial genomic annotations 
R
library(GenomicRanges)
prom=read.table("promoters_from_coords.txt",sep = "\t",header=F)#-1000bp to +100bp relative to TSS (hg19)
ets2=read.table("ets2_peaks.txt",sep="\t",header=F)#idr peaks from ETS2 C&R
h3k27=read.table("h3k27ac_from_coords.txt",sep="\t",header=F)#H3K27ac peaks from ChIP-seq in H3k27ac_TPP_hg19_minor_fromMergedBAM_peaks.narrowpeak
body=read.table("genebody_from_coords.txt",sep="\t",header=F)#gene body coordinates (+100 from TSS to end of gene, hg19)
prom=GRanges(prom$V1,IRanges(prom$V2,prom$V3))
ets2=GRanges(ets2$V1,IRanges(ets2$V2,ets2$V3))
h3k27=GRanges(h3k27$V1,IRanges(h3k27$V2,h3k27$V3))
body=GRanges(body$V1,IRanges(body$V2,body$V3))
ets2_body=findOverlaps(ets2,body)
ets2_prom=findOverlaps(ets2,prom1)
ets2_h3k27=findOverlaps(ets2,h3k27)