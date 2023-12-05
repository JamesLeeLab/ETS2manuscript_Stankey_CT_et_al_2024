## Pipeline used for PAIRED-END ChIP-seq analysis with hg19 build
# Used for H3K27ac ChIP-seq data on TPP macrophages with ETS2 overexpression or KO, used for figures 4 and S2
# this pipeline works on paired-end fastq files for ChIP-seq datasets

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
folder="/Trimgalore"

a=(*.fastq.gz)
for ((i=0; i<${#a[@]}; i+=2));
do 
        trim_galore --phred33 -q 24 --illumina --length 30 --paired --fastqc \
        -o $folder/ "${a[i]}" "${a[i+1]}";
done

cd Trimgalore
multiqc .

#Align to hg19 with bwa mem (reads >70bp) in PE mode:
module load BWA/0.7.15-intel-2017a

a=(*.fq.gz)
for ((i=0; i<${#a[@]}; i+=2));
do 
        bwa mem /BWAIndex/genome.fa \
        "${a[i]}" "${a[i+1]}" \
        2> /bwa_mem_hg19/${a[i]%.fq.gz}_hg19_bwamem.log  \
        > /bwa_mem_hg19/${a[i]%.fq.gz}_hg19_bwamem.sam
done

#convert SAM into BAM files:
module load GCC/10.2.0
module load SAMtools/1.11-GCC-10.2.0
cd bwa_mem_hg19

for a in *.sam; 
do 
        samtools view -bS - > ${a%.sam}.bam; 
done

#sorting BAM:
for a in *.bam;
do 
        samtools sort $a \
        -o ${a%.bam}_sorted.bam;
done

#merging technical replicates (rerun of same libraries):
#ETS2 KO vs NTC:
samtools merge H3K27ac_TPP_NCI106NTC_hg19_STmerged.bam H3K27ac_TPP_NCI106NTC_S5_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_NCI106NTC_S5_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_NCI106KO_hg19_STmerged.bam H3K27ac_TPP_NCI106KO_S6_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_NCI106KO_S6_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_NCI107NTC_hg19_STmerged.bam H3K27ac_TPP_NCI107NTC_S7_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_NCI107NTC_S7_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_NCI107KO_hg19_STmerged.bam H3K27ac_TPP_NCI107KO_S8_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_NCI107KO_S8_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_NCI109NTC_hg19_STmerged.bam H3K27ac_TPP_NCI109NTC_S11_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_NCI109NTC_S11_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_NCI109KO_hg19_STmerged.bam H3K27ac_TPP_NCI109KO_S12_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_NCI109KO_S12_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_NCI106NTC_input_hg19_STmerged.bam H3K27ac_TPP_NCI106NTC_input_S17_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_NCI106NTC_input_S17_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_NCI106KO_input_hg19_STmerged.bam H3K27ac_TPP_NCI106KO_input_S18_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_NCI106KO_input_S18_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_NCI107NTC_input_hg19_STmerged.bam H3K27ac_TPP_NCI107NTC_input_S19_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_NCI107NTC_input_S19_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_NCI107KO_input_hg19_STmerged.bam H3K27ac_TPP_NCI107KO_input_S20_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_NCI107KO_input_S20_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_NCI109NTC_input_hg19_STmerged.bam H3K27ac_TPP_NCI109NTC_input_S23_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_NCI109NTC_input_S23_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_NCI109KO_input_hg19_STmerged.bam H3K27ac_TPP_NCI109KO_input_S24_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_NCI109KO_input_S24_L003_R1_001_val_1_hg19_bwamem_sorted.bam
#ETS2 vs REV overexpression:
samtools merge H3K27ac_TPP_D1_ETS2_500_hg19_STmerged.bam H3K27ac_TPP_D1_ETS2_500_S25_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_D1_ETS2_500_S25_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_D1_REV_500_hg19_STmerged.bam H3K27ac_TPP_D1_REV_500_S28_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_D1_REV_500_S28_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_D2_ETS2_500_hg19_STmerged.bam H3K27ac_TPP_D2_ETS2_500_S31_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_D2_ETS2_500_S31_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_D2_REV_500_hg19_STmerged.bam H3K27ac_TPP_D2_REV_500_S32_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_D2_REV_500_S32_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_D3_ETS2_500_hg19_STmerged.bam H3K27ac_TPP_D3_ETS2_500_S39_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_D3_ETS2_500_S39_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_D3_REV_500_hg19_STmerged.bam H3K27ac_TPP_D3_REV_500_S42_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_D3_REV_500_S42_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_D1_ETS2input_500_hg19_STmerged.bam H3K27ac_TPP_D1_ETS2input_500_S51_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_D1_ETS2input_500_S51_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_D1_REVinput_500_hg19_STmerged.bam H3K27ac_TPP_D1_REVinput_500_S52_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_D1_REVinput_500_S52_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_D2_ETS2input_500_hg19_STmerged.bam H3K27ac_TPP_D2_ETS2input_500_S55_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_D2_ETS2input_500_S55_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_D2_REVinput_500_hg19_STmerged.bam H3K27ac_TPP_D2_REVinput_500_S56_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_D2_REVinput_500_S56_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_D3_ETS2input_500_hg19_STmerged.bam H3K27ac_TPP_D3_ETS2input_500_S63_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_D3_ETS2input_500_S63_L003_R1_001_val_1_hg19_bwamem_sorted.bam
samtools merge H3K27ac_TPP_D3_REVinput_500_hg19_STmerged.bam H3K27ac_TPP_D3_REVinput_500_S64_L001_R1_001_val_1_hg19_bwamem_sorted.bam H3K27ac_TPP_D3_REVinput_500_S64_L003_R1_001_val_1_hg19_bwamem_sorted.bam

#indexing:
for a in *STmerged.bam;
do 
        samtools index $a;
done

#Mark duplicates with Picard:
ml load  picard/2.1.1-Java-1.8.0_92

for a in *STmerged.bam;
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
        --binSize 10 --normalizeUsing CPM --extendReads 200;
done

## Now getting number of reads from input and corresponding sample to see if downsample is necessary:
module load MACS2/2.2.7.1-foss-2021b
#example with "sample 1":
macs2 randsample -i H3K27ac_TPP_NCI106NTC_input_hg19_STmerged_rmdup_resorted.bam \
-p 100 -o H3K27ac_TPP_NCI106NTC_input_hg19.bed --outdir ./
inputTags=$(wc -l < H3K27ac_TPP_NCI106NTC_input_hg19.bed)
echo "Input size : "${inputTags}" reads"
sampleSize_1=$(samtools idxstats H3K27ac_TPP_NCI106NTC_hg19_STmerged_rmdup_resorted.bam | cut -f3 \
| awk 'BEGIN {total=0} {total += $1} END {print total}')
echo "ChIP sample size 1: "$sampleSize_1" reads"
#done for all samples but no downsampling was necessary for any sample

## Run Macs2 to call peaks: 
#Use the resorted and indexed sample bam files (post-PICARD processing).
#Change the names for samples (-t), inputs (-c) and outputs (-n) for each call:
macs2 callpeak -t sample1.bam -c input1.bam -g hs -n H3k27ac_1_hg19 -f AUTO --outdir PEAKS --qvalue 0.01 -B --nomodel --extsize=200
#example with "sample 1":
macs2 callpeak -t H3K27ac_TPP_NCI106NTC_hg19_STmerged_rmdup_resorted.bam -c H3K27ac_TPP_NCI106NTC_input_hg19_STmerged_rmdup_resorted.bam -g hs -n H3K27ac_TPP_NCI106NTC_hg19 -f AUTO --outdir PEAKS --qvalue 0.01 -B --nomodel --extsize=200 2> PEAKS/H3K27ac_TPP_NCI106NTC_hg19_macs2.log
# done for all the bam files from ETS2 KO vs NTC (no-template control), and ETS2 overexpression vs REV (reverse)