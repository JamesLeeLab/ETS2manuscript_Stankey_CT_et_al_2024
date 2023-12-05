## Analysis of TPP macrophages ChIPseq data (H3k27ac) merged by genotype at rs2836882 (cf. note below) and calling of superenhancers (SE)
# HiSeq 4000 run (50bp SE reads) aligned to hg19
# These H3K27ac data are from donors that are heterozygote, or homozygote for the minor allele, or homozygote for the major allele at rs2836882 (figures 1 and S4)
# Follow first the ChIPseq_pipeline to process the data, until the bam files have been processed using PICARD, resorted and indexed.
# The homozygotes(hom) minor or major and the heterozygotes (het) biological replicates are merged using the bam files as indicated above.
#merge BAM files from biological replicates (both samples and inputs):
module load SAMtools/1.11-GCC-10.2.0
samtools merge H3K27ac_TPP_hg19_minor_merged.bam 19_TPP_H3K27ac_S3_L004_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam 30_TPP_H3K27ac_S11_L005_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam
samtools merge H3K27ac_TPP_hg19_Inputminor_merged.bam 19_TPP_input_S7_L004_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam 30_TPP_input_S15_L005_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam
samtools merge H3K27ac_TPP_hg19_major_merged.bam 15_TPP_H3K27ac_S4_L004_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam 24_TPP_H3K27ac_S12_L005_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam
samtools merge H3K27ac_TPP_hg19_Inputmajor_merged.bam 15_TPP_input_S8_L004_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam 24_TPP_input_S16_L005_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam
samtools merge H3K27ac_TPP_hg19_HET_merged.bam 16_TPP_H3K27ac_S1_L004_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam 20_TPP_H3K27ac_S2_L004_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam 29_TPP_H3K27ac_S9_L005_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam 34_TPP_H3K27ac_S10_L005_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam
samtools merge H3K27ac_TPP_hg19_InputHET_merged.bam 16_TPP_input_S5_L004_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam 20_TPP_input_S6_L004_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam 29_TPP_input_S13_L005_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam 34_TPP_input_S14_L005_R1_001_trimmed_hg19_sorted_rmdup_resorted.bam
#index the new merged files:
for a in *or_merged.bam; do samtools index $a; done
for a in *HET_merged.bam; do samtools index $a; done
#I checked if downsample was necessary (cf. code from the ChIPseq_pipeline) but no need to downsample here.

#Run Macs2 to call peaks with merged biological replicates and inputs:
macs2 callpeak -t H3K27ac_TPP_hg19_major_merged.bam -c H3K27ac_TPP_hg19_Inputmajor.bed -g hs -n H3k27ac_TPP_hg19_major_fromMergedBAM -f AUTO --outdir PEAKSbamMerged --qvalue 0.01 -B --nomodel --extsize=200
macs2 callpeak -t H3K27ac_TPP_hg19_minor_merged.bam -c H3K27ac_TPP_hg19_Inputminor.bed -g hs -n H3k27ac_TPP_hg19_minor_fromMergedBAM -f AUTO --outdir PEAKSbamMerged --qvalue 0.01 -B --nomodel --extsize=200
macs2 callpeak -t H3K27ac_TPP_hg19_HET_merged.bam -c H3K27ac_TPP_hg19_InputHET.bed -g hs -n H3k27ac_TPP_hg19_HET_fromMergedBAM -f AUTO --outdir PEAKSbamMerged --qvalue 0.01 -B --nomodel --extsize=200

#making bigwig files from merged samples:
for a in *merged.bam; 
do 
        bamCoverage -b $a \
        -o bigwigmerged/${a%.bam}.bw -of bigwig --binSize 10 --normalizeUsing CPM --extendReads 200;
done

# run ROSE with these MACS2 peaks to call superenhancers (SE)
# ROSE described here: https://github.com/younglab/ROSE
# first, convert the peak files (BED) to gff format using the Galaxy project website: https://usegalaxy.eu/?tool_id=CONVERTER_bed_to_gff_0&version=2.0.1
#then call SE:
python ROSE_main.py -g hg19 -i /PEAKSbamMerged/H3k27ac_TPP_hg19_HET.gff -r /H3K27ac_TPP_hg19_HET_merged.bam -o /PEAKSbamMerged/ROSE -t 2000 -c /H3K27ac_TPP_hg19_InputHET_merged.bam 
python ROSE_main.py -g hg19 -i /PEAKSbamMerged/H3k27ac_TPP_hg19_major.gff -r /H3K27ac_TPP_hg19_major_merged.bam -o /PEAKSbamMerged/ROSE -t 2000 -c /H3K27ac_TPP_hg19_Inputmajor_merged.bam 
python ROSE_main.py -g hg19 -i /PEAKSbamMerged/H3k27ac_TPP_hg19_minor.gff -r /H3K27ac_TPP_hg19_minor_merged.bam -o /PEAKSbamMerged/ROSE -t 2000 -c /H3K27ac_TPP_hg19_Inputminor_merged.bam
