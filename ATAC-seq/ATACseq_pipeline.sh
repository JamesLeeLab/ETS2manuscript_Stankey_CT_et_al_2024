##Pipeline for ATAC-seq analysis with hg19:
#this pipeline works on paired-end fastq files for ATAC-seq datasets. 

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

a=(/path/to/fasq/files/*.fastq.gz)
for ((i=0; i<${#a[@]}; i+=2));
do 
        trim_galore --phred33 -q 24 --nextera \
        --length 20 --paired --fastqc \
        -o $folder/ "${a[i]}" "${a[i+1]}";
done

multiqc .

#Alignment with BWA
#note: use bwa mem if reads >70bp
module load BWA/0.7.15-intel-2017a
cd Trimgalore/

a=(path/to/trimmed/files/*.fq.gz)
for ((i=0; i<${#a[@]}; i+=2));
do 
        bwa mem /BWAIndex/genome.fa \
        "${a[i]}" "${a[i+1]}" \
        2> ${a[i]%.fq.gz}_hg19_bwaln.log \
        > ${a[i]%.fq.gz}_hg19_bwaln.sam;
done

#convert SAM into BAM files:
module load GCC/10.2.0
module load SAMtools/1.11-GCC-10.2.0

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
cd sortedBAM_hg19/

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
        samtools view $a -b -h -f 2 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 \
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
        --binSize 10 --normalizeUsing CPM --extendReads;
done

##Run Macs2 to call peaks: 

for a in *rmdup_resorted.bam; 
do 
        macs2 callpeak -t $a -g hs \
        -f BAMPE --outdir PEAKS \
        -n ${a%_sorted_rmdup_resorted.bam}_macs2 \
        -B --nomodel --nolambda --keep-dup all --call-summits \
        2> PEAKS/${a%.bam}_macs2.log; 
done
