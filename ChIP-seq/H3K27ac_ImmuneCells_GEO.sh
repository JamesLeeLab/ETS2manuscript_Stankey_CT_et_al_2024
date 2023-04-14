##Analysing H3K27ac ChIPseq data from various immune cell types  
#SRA files downloaded from GEO:
#Neutrophils: encode data GSE96014
#Tregs: datasets GSM1893223, GSM1893225, GSM1893236, GSM1893237
#Monocytes: datasets GSM1003559 and GSM1003475
#B cells, CD4 T cells, CD8 T cells, NK cells (from Hnisz et al. 2013): 
#datasets GSM1027287, GSM1027304, GSM772905, GSM772904, GSM1102781, GSM1102806, GSM1027288, GSM1027305

##Convert .sra to .fq files:
ml SRA-Toolkit/2.9.4-centos_linux64
for a in *.sra;
do 
    fastq-dump --gzip --split-files $a -O /fastq; 
done

##Follow first the ChIPseq_pipeline to process the data, until the bam files have been sorted, but don't run PICARD yet.
##IMPORTANT NOTE: for neutrophils, bwa alignment (single-end) was performed with bwa MEM since reads>70bp

#merging technical replicates (rerun of same libraries):
samtools merge H3K27_Neutro_S1_hg19_STmerged.bam H3K27_Neutro_S1_SRR5331498_1_trimmed_hg19_sorted.bam H3K27_Neutro_S1_SRR5331499_1_trimmed_hg19_sorted.bam
samtools merge H3K27_Neutro_S2_hg19_STmerged.bam H3K27_Neutro_S2_SRR5331500_1_trimmed_hg19_sorted.bam H3K27_Neutro_S2_SRR5331501_1_trimmed_hg19_sorted.bam H3K27_Neutro_S2_SRR5331502_1_trimmed_hg19_sorted.bam H3K27_Neutro_S2_SRR5331503_1_trimmed_hg19_sorted.bam
samtools merge H3K27_Neutro_S3_hg19_STmerged.bam H3K27_Neutro_S3_SRR5331504_1_trimmed_hg19_sorted.bam H3K27_Neutro_S3_SRR5331505_1_trimmed_hg19_sorted.bam
samtools merge H3K27_Monocytes-CD14_hg19_STmerged.bam SRR568409_GSM1003559_Broad_ChipSeq_Monocytes-CD14_RO01746_H3K27ac_trimmed_hg19_sorted.bam SRR568410_GSM1003559_Broad_ChipSeq_Monocytes-CD14_RO01746_H3K27ac_trimmed_hg19_sorted.bam
samtools merge Input_Monocytes-CD14_hg19_STmerged.bam SRR568246_GSM1003475_Broad_ChipSeq_Monocytes-CD14_RO01746_Control_trimmed_hg19_sorted.bam SRR568247_GSM1003475_Broad_ChipSeq_Monocytes-CD14_RO01746_Control_trimmed_hg19_sorted.bam 
samtools merge Input_CD19_Primary_Cells_hg19_STmerged.bam SRR609624_ChIP-Seq_Input_of_CD19_Primary_Cells_trimmed_hg19_sorted.bam SRR609625_ChIP-Seq_Input_of_CD19_Primary_Cells_trimmed_hg19_sorted.bam 
samtools merge Input_CD56_Primary_Cells_hg19_STmerged.bam SRR609626_ChIP-Seq_Input_of_CD56_Primary_Cells_trimmed_hg19_sorted.bam SRR609627_ChIP-Seq_Input_of_CD56_Primary_Cells_trimmed_hg19_sorted.bam
samtools merge H3K27_CD19_Primary_Cells_hg19_STmerged.bam SRR609630_Histone_H3K27ac_ChIP-Seq_of_CD19_Primary_Cells_trimmed_hg19_sorted.bam SRR609631_Histone_H3K27ac_ChIP-Seq_of_CD19_Primary_Cells_trimmed_hg19_sorted.bam
samtools merge H3K27_CD56_Primary_Cells_hg19_STmerged.bam SRR609634_Histone_H3K27ac_ChIP-Seq_of_CD56_Primary_Cells_trimmed_hg19_sorted.bam SRR609635_Histone_H3K27ac_ChIP-Seq_of_CD56_Primary_Cells_trimmed_hg19_sorted.bam
samtools merge Input_CD8_Primary_Cells_hg19_STmerged.bam SRR787511_ChIP-Seq_Input_of_CD8_Primary_Cells_trimmed_hg19_sorted.bam SRR787510_ChIP-Seq_Input_of_CD8_Primary_Cells_trimmed_hg19_sorted.bam SRR787512_ChIP-Seq_Input_of_CD8_Primary_Cells_trimmed_hg19_sorted.bam
samtools merge H3K27_CD8_Primary_Cells_hg19_STmerged.bam SRR787520_Histone_H3K27ac_ChIP-Seq_of_CD8_Primary_Cells_trimmed_hg19_sorted.bam SRR787521_Histone_H3K27ac_ChIP-Seq_of_CD8_Primary_Cells_trimmed_hg19_sorted.bam 

#indexing:
for a in *sorted.bam;
do 
        samtools index $a;
done

#Mark duplicates with Picard:
ml load  picard/2.1.1-Java-1.8.0_92

for a in *sorted.bam;
do 
        picard MarkDuplicates INPUT=$a \
        OUTPUT=${a%.bam}_markDup.bam \
        METRICS_FILE=${a%.bam}_markDup.txt \
        REMOVE_DUPLICATES=false \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=./ ;
done

#loop to remove duplicates and unmapped reads:
for a in *markDup.bam;
do 
        samtools view $a -b -h -F 4 -F 256 -F 1024 -F 2048 -q 15 \
        -o ${a%markDup.bam}rmdup.bam;
done

#then loop to sort:
for a in *rmdup.bam;
do 
        samtools sort $a \
        -o ${a%.bam}_resorted.bam;
done

#then loop to index:
for a in *rmdup_resorted.bam;
do 
        samtools index $a;
done

#Creating bigwig files to visualize in IGV:
ml deepTools/3.3.1-foss-2020a-Python-3.8.2

for a in *rmdup_resorted.bam;
do 
        bamCoverage -b $a \
        -o bigwig/${a%.bam}.bw -of bigwig \
        --binSize 50 --normalizeUsing CPM --extendReads 200;
done