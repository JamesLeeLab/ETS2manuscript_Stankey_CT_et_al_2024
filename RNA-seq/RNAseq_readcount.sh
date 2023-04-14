## Pipeline used to get the read count from RNA-seq data
# This pipeline works on paired-end fastq files

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

#Trimgalore in paired-end more with extra trimmming (3 bases on each read) 
#as indicated in the library prep kit manual (SMARTer Stranded Total RNA-Seq Kit v2 - Pico Input Mammalian)
module load TrimGalore/0.6.0
folder="/Trimgalore"

a=(path/to/fastq/files/*.fastq.gz)

for ((i=0; i<${#a[@]}; i+=2));
do 
        trim_galore --phred33 -q 24 --illumina --length 20 --paired \ 
        --stringency 13 --fastqc --clip_R1 3 --clip_R2 3 --trim1 \
        -o $folder/ "${a[i]}" "${a[i+1]}";
done

#MultiQC report:
cd Trimgalore
multiqc .

#run BBSplit to get rid of ribosomal reads, using the Human ribosomal DNA complete repeating unit, GenBank: U13369.1
module load BBMap/36.20-foss-2016b-Java-1.8.0_92
mkdir bbsplit

a=(*.fq.gz)
for ((i=0; i<${#a[@]}; i+=2));
do 
        bbsplit.sh ref=humanrRNAseq.fasta \
        in1="${a[i]}" in2="${a[i+1]}" basename=bbsplit/${a[i]%.fq.gz}_output%.fq.gz \
        outu1=bbsplit/${a[i]%.fq.gz}_1_clean.fq.gz outu2=bbsplit/${a[i]%.fq.gz}_2_clean.fq.gz \
        2>bbsplit/${a[i]%.fq.gz}.summary.txt;
done

#QC of "clean" files:
cd bbsplit

for a in *.fastq.gz;
do 
        fastqc $a;
done

multiqc .

#alignment to hg38 with Hisat2:
mkdir Hisat2BAM
mkdir Hisat2BAM/fail
ml load HISAT2/2.1.0-foss-2016b
ml load SAMtools/1.11-GCC-10.2.0

a=(*clean.fq.gz)
for ((i=0; i<${#a[@]}; i+=2));
do 
        hisat2 -q --phred33 -p 8 -I 20 -X 350 --fr -x Hisat2/grch38_tran/genome_tran \
        -1 "${a[i]}" -2 "${a[i+1]}" -S Hisat2BAM/${a[i]%.fq.gz}_GRCh38.sam \
        --un-conc Hisat2BAM/fail/${a[i]%.fq.gz}_discordantpairs_%.sam \
        --un Hisat2BAM/fail/${a[i]%.fq.gz}_failalign_%.sam \
        2> Hisat2BAM/${a[i]%.fq.gz}_summarymetrics.txt ;
done

#convert SAM to BAM:
for a in *.sam; do samtools view -bS - > ${a%.sam}.bam; done
#sort BAM files:
for a in *.bam; do samtools sort -o ${a%.bam}_sorted.bam $a; done
#indexing:
for a in *sorted.bam; do samtools index $a; done
#Multiqc:
multiqc .


##Use R to get the read count:
R
library(limma)
library(Rsubread) 
save="path/to/readcountfolder"
BAM=dir(".","bam$")
r=assign(x=paste0("reads.counts",sep=""), \
value=featureCounts(BAM,annot.ext="Homo_sapiens.GRCh38.102.gtf.gz",isGTFAnnotationFile=T,\
isPairedEnd=T,requireBothEndsMapped=F,minFragLength=20,checkFragLength=F,\
strandSpecific=2,ignoreDup=F,nthreads=9),envir= .GlobalEnv)
setwd(save)
saveRDS(r,paste0("reads.counts.GRCh38.gtf102.NAME",".RDS",sep=""))
