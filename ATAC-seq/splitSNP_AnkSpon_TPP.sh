## Code to run splitSNP to get allele-specific reads from ATAC-seq data (figure 1k and supplementary figure 4g)
# splitSNP described here: https://github.com/astatham/splitSNP
# clone from github:
git clone https://github.com/astatham/splitSNP.git
# load required modules:
ml load Pysam/0.20.0-GCC-11.3.0
ml load Python/3.5.2-foss-2016b
# input = bam files that have been processed folllowing the "ATACseq_pipeline"

# Analysis of ankylosing spondylitis patients and healthy controls CD14 ATAC-seq data from Brown et al., Cell Genomics 2023 (https://doi.org/10.1016/j.xgen.2023.100306) done following the "ATACseq_pipeline"
# running splitSNP to split reads from the 2 alleles at rs2836882 using the processed BAM files from heterozygote donors
cd /Ank_spond_Brown_Cell_Genomics/CD14/ATACseq/trim/BAM
# NOTE: the original splitSNP.py script was giving several errors. I had to add "()"" to all the "print" commands that didn't have them
# I also had to delete the part of the original script around line 87 which was checking for an index. Using now the modified python script, stored in our ATAC-seq github repository
python3 splitSNP/splitSNP.py ATAC_AS045_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam ATAC_AS045_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py ATAC_CTRL1_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam ATAC_CTRL1_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py ATAC_CTRL2_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam ATAC_CTRL2_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py ATAC_CTRL3_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam ATAC_CTRL3_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py ATAC_CTRL4_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam ATAC_CTRL4_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py FATAC_CTRL6_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam FATAC_CTRL6_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py OATAC_AS052_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam OATAC_AS052_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py OATAC_AS058_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam OATAC_AS058_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py OATAC_CTRL16_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam OATAC_CTRL16_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py OATAC_CTRL18_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam OATAC_CTRL18_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py OATAC_CTRL21_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam OATAC_CTRL21_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py OATAC_CTRL23_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam OATAC_CTRL23_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py OATAC_CTRL24_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam OATAC_CTRL24_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py OATAC_CTRL44_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam OATAC_CTRL44_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py OATAC_CTRL45_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam OATAC_CTRL45_CD14_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py OATAC_CTRL51_CD14_1_val_1_hg19_bwamem_sorted_rmdup_sorted.bam OATAC_CTRL51_CD14_rs2836882 chr21:40466570:G:A

## Now using splitSNP to split reads from the 2 alleles at rs2836882 using the processed BAM files from heterozygote donors, from TPP macrophages ATAC-seq data:
cd /TPP_ATAC/trim/BAM
python3 splitSNP/splitSNP.py d2.het_NTC_S8_L001_R1_001_val_1_hg19_bwaln_sorted_rmdup_resorted.bam ATAC_d2.het_NTC_TPP_rs2836882 chr21:40466570:G:A
python3 splitSNP/splitSNP.py d3.het_NTC_S10_L001_R1_001_val_1_hg19_bwaln_sorted_rmdup_resorted.bam ATAC_d3.het_NTC_TPP_rs2836882 chr21:40466570:G:A
