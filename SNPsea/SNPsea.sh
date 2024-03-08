#Making input file for SNPsea using all GO BP
#James Lee, 2022-08-18
#To start, you need to download the GO biological pathways .gmt file from MSigDB. Transpose that file (so you have pathways in columns) and append any extra pathways you want to assess. Remove the second row (NA or url for source). Save as a csv.
R
library(splitstackshape)
setwd('C:/Users/leej1/Dropbox/SNPsea')
data=read.csv(file="All_GOBP_genesets.csv",header=T,na.strings=c("")) #All_GOBP_genesets.csv is a transposed gmt file downloaded from MSigDB
my_lists=lapply(data,function(col)col[!is.na(col)])
res=t(splitstackshape:::charMat(my_lists,fill=0L))
colnames(res)=colnames(data)
res2=res[res[,1]==1,]
res2=res2[,-1]
res3=as.data.frame(res2)
Description=rownames(res2)
res3$Description=Description
rows=read.csv(file="rownames.csv",header=T)#this is a csv with GeneSymbols and NCBI IDs to enable appending them
library(plyr)
gct=join(rows,res3,type='full')
## Joining by: Description
library(brio)
## 
## Attaching package: 'brio'
## The following objects are masked from 'package:base':
## 
##     readLines, writeLines
output="GO_BP_ETS2_SNPsea.gct"
writeLines(c("#1.2", paste(nrow(gct), ncol(gct) - 2, collapse="\t")), output, sep="\n")
write.table(gct, file=output, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", append = TRUE)
## Warning in write.table(gct, file = output, quote = FALSE, row.names = FALSE, :
## appending column names to file

#SNPsea - code to generate pathway file from GO Biological Pathways (downloaded as .gmt file from MSigDB) attached as RMD file
#run SNPsea v1.0.4(linux)
bin/snpsea-linux64 --args args.txt

#calculate enrichment z score from SNPsea output
library(readr)
library(ggplot2)
library(stringr)
library(ggrepel)
library(dplyr)
d <- read_tsv("condition_pvalues.txt")
colnames(d)=c("condition","pvalue","nulls_observed","nulls_tested","score","null_scores")
names=c("ETS2_DN","GOBP_ALPHA_BETA_T_CELL_ACTIVATION","GOBP_AUTOPHAGOSOME_ORGANIZATION","GOBP_CELL_ADHESION_MEDIATED_BY_INTEGRIN","GOBP_INTEGRIN_ACTIVATION","GOBP_INTERLEUKIN_10_PRODUCTION","GOBP_INTESTINAL_EPITHELIAL_STRUCTURE_MAINTENANCE","GOBP_MAINTENANCE_OF_GASTROINTESTINAL_EPITHELIUM","GOBP_NEGATIVE_REGULATION_OF_AUTOPHAGY","GOBP_POSITIVE_REGULATION_OF_AUTOPHAGY","GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_10_PRODUCTION","GOBP_POSITIVE_REGULATION_OF_T_HELPER_1_TYPE_IMMUNE_RESPONSE","GOBP_POSITIVE_REGULATION_OF_T_HELPER_17_TYPE_IMMUNE_RESPONSE","GOBP_PROCESS_UTILIZING_AUTOPHAGIC_MECHANISM","GOBP_REGULATION_OF_AUTOPHAGY","GOBP_RESPONSE_TO_MURAMYL_DIPEPTIDE","GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR","GOBP_T_HELPER_17_CELL_DIFFERENTIATION","GOBP_T_HELPER_17_TYPE_IMMUNE_RESPONSE","GOBP_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY","OE_ETS2_UP")
subset=d[d$condition %in% names,]
subset$nulls_mean <- 0
for (i in seq(nrow(subset))) {
subset$nulls_mean[i] <- mean(parse_number(str_split(subset$null_scores[i], ",")[[1]]))
}
subset$nulls_sem <- 0
for (i in seq(nrow(subset))) {
subset$nulls_sem[i] <- sd(parse_number(str_split(subset$null_scores[i], ",")[[1]])/sqrt(subset$nulls_tested[i]))
}
subset$enrich.z=(subset$score-subset$nulls_mean)/(subset$nulls_sem)
print(subset, n=21)