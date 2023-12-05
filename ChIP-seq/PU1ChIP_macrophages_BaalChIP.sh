## Macrophage PU1 ChIPseq data from GEO, analyzed with hg19, used for figures 1 and 4
## Fastq files from PU1 ChIP data GSM1681423 & GSM1681425, and their inputs GSM1681427 & GSM1681429, downloaded from GEO.
## Follow first the "ChIPseq_pipeline_SingleEnd" to process the data, until the bam files have been processed using PICARD, resorted and indexed.

##Call peaks with Macs2:
macs2 callpeak -t S1_GSM1681423_macro_baselinespi1_trimmed_hg19_rmdup_resorted.bam -c S1_GSM1681427_macro_input_baseline.fastqsanger.gz_trimmed_hg19_sorted_rmdup_resorted.bam -g hs -n PU1_S1_macro_Baseline_hg19 -f AUTO --outdir PEAKS --qvalue 0.05 2> PEAKS/PU1_S1_macs2.log
macs2 callpeak -t S3_GSM1681425_macro_il4spi1.fastqsanger.gz_trimmed_hg19_sorted_rmdup_resorted.bam -c S3_GSM1681429_macro_input_il4.fastqsanger.gz_trimmed_hg19_sorted_rmdup_resorted.bam -g hs -n PU1_S3_macro_il4_hg19 -f AUTO --outdir PEAKS --qvalue 0.05 2> PEAKS/PU1_S3_macs2.log

##Allele Specific ChIP analysis with BaalChIP (https://github.com/InesdeSantiago/BaalChIP):
#The "hets.txt" file contains the details of the SNPs tested in this analysis.
#First, we need to convert the peak files (.narrowpeak files) to .bed:
for a in *.narrowPeak; do cut -f 1-6 $a > ${a%.narrowPeak}.bed; done
#Then run BaalChIP in R:
R
library(BaalChIP)
path="/sortedBAM_hg19"
setwd(path)
samplesheet <- file.path("path", "PU1_ChIP_SampleSheet.txt")
#below, we could specify several groups to analyze, but here we only add "Macro" as we have only 1 group:
hets <- c("Macro"="hets.txt")
res <- BaalChIP(samplesheet=samplesheet, hets=hets)
res
BaalChIP.get(res, what="samples")
res <- alleleCounts(res, min_base_quality=10, min_mapq=15, verbose=FALSE)
#QC Filter but without "regionstokeep" since it only works for reads <50bp
data('blacklist_hg19')
data('pickrell2011cov1_hg19')
res <- QCfilter(res, RegionsToFilter=list('blacklist'=blacklist_hg19,'highcoverage'=pickrell2011cov1_hg19))
res <- mergePerGroup(res)
#retrieve mergedCounts:
counts <- BaalChIP.get(res, 'mergedCounts')
#mergedCounts are grouped by group_name:
names(counts)
sapply(counts, dim)
#check out the result for one of the groups:
head(counts[[1]])
counts
res <- filter1allele(res)
#apply RM and RAF corrections:
res <- getASB(res, Iter=5000, conf_level=0.95, cores = 4, 
              RMcorrection = TRUE, 
              RAFcorrection= TRUE)
#Getting results
result <- BaalChIP.report(res)
write.csv(result$Macro, "PU1Macro_BaalChIP_Results.csv")
write.csv(counts$Macro, "PU1Macro_BaalChIP_counts.csv")