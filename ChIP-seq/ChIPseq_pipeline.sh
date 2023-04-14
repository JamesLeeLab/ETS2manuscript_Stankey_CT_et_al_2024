##Pipeline for ChIP-seq analysis with hg19:
#this pipeline works on single-end fastq files for ChIP-seq datasets. 
#It uses the genome build hg19 to be compatible with super-enhancer analysis (when looking at H3k27ac data) using ROSE (https://github.com/younglab/ROSE).

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

for a in *.fastq.gz;
do 
        trim_galore $a --phred33 -q 24 --illumina --length 30 --fastqc \
        -o $folder/ ;
done

multiqc .

#Alignment with BWA
#note: use bwa mem (instead of bwa aln) for reads >70bp
module load BWA/0.7.15-intel-2017a
cd Trimgalore/

for a in *.fastq.gz;
do 
        bwa aln /BWAIndex/genome.fa \
        $a \
        2> /bwa_aln_hg19/${a%.fq.gz}_hg19_bwaln.log \
        > /bwa_aln_hg19/${a%.fq.gz}_hg19.sam;
done

#convert SAM into BAM files:
module load GCC/10.2.0
module load SAMtools/1.11-GCC-10.2.0
cd bwa_aln_hg19/

for a in *.sam; 
do 
        samtools view -bS - > BAM_hg19/${a%.sam}.bam; 
done

#sorting BAM:
for a in *.bam;
do 
        samtools sort $a \
        -o /sortedBAM_hg19/${a%.bam}_sorted.bam;
done

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
        --binSize 10 --normalizeUsing CPM --extendReads 200;
done

##Now getting number of reads from input and corresponding sample to see if downsample is necessary:
module load MACS2/2.2.5-foss-2018b-Python-3.6.6
#example with "sample 1":
macs2 randsample -i 1_TPP_input_S6_L004_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam \
-p 100 -o 1_TPP_H3K27_hg19_input.bed --outdir ./
inputTags=$(wc -l < 1_TPP_H3K27_hg19_input.bed)
echo "Input size : "${inputTags}" reads"
sampleSize_1=$(samtools idxstats  1_TPP_H3K27ac_S2_L004_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam | cut -f3 \
| awk 'BEGIN {total=0} {total += $1} END {print total}')
echo "ChIP sample size 1: "$sampleSize_1" reads"
#if input size and sample size are different (higher read count for input) then next step is to downsample input:
macs2 randsample -i 1_TPP_input_S6_L004_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam \
-n ${sampleSize_20} -o 1_TPP_H3K27_hg19_input_downsample.bed --outdir ./
#to delete the now useless 1st bed file
rm 1_TPP_H3K27_hg19_input.bed
#do this for all samples before calling peaks

##Run Macs2 to call peaks: 
#Use the resorted and indexed sample bam files (post-PICARD processing). Use the downsampled input bed files where necessary.
#Change the names for samples (-t), inputs (-c) and outputs (-n) for each call:
macs2 callpeak -t sample1.bam -c input1.bed -g hs -n H3k27ac_1_hg19 -f AUTO --outdir PEAKS --qvalue 0.01 -B --nomodel --extsize=200

##Run ROSE to identify Super-Enhancers (cf. https://github.com/younglab/ROSE):
#first, convert the peak files (BED) to gff format using the Galaxy project website: https://usegalaxy.eu/?tool_id=CONVERTER_bed_to_gff_0&version=2.0.1
#use the same bam files as above for MACS2, and change the names for samples (-r (bam) and -i (peaks)) and inputs (-c) for each call
python ROSE_main.py -g hg19 -i /PEAKS/H3k27ac_1_hg19.gff -r /sample1.bam -o /ROSE -t 2000 -c /input1.bam 
