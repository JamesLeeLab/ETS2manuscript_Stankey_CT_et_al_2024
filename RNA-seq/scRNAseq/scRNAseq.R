# R code used for scRNAseq analysis of immune cell data from "Human CD atlas study between colon and terminal ileum" 
# available from Single Cell Portal: https://singlecell.broadinstitute.org/single_cell/study/SCP1884/human-cd-atlas-study-between-colon-and-terminal-ileum#study-download
# used for figure 5
library(dplyr)
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(magrittr)
library(dplyr)
library(limma)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
col.imm.data=readMM('CO_IMM.scp.raw.mtx')
col.imm.data2=as(col.imm.data,"CsparseMatrix")
cells=read.table(file="CO_IMM.scp.barcodes.tsv",sep='\t',header=F)
colnames(col.imm.data2)=cells$V1
features=read.table(file="CO_IMM.scp.features.tsv", sep='\t', header=F)
rownames(col.imm.data2)=features$V2
myeloid_barcodes=read.csv(file="myeloid_barcodes.csv",header=F) #list of myeloid cell barcodes based on authors' annotation
col.imm.myeloid=col.imm.data2[,colnames(col.imm.data2) %in% myeloid_barcodes$V1]
myeloid=CreateSeuratObject(counts=col.imm.myeloid)
myeloid=PercentageFeatureSet(myeloid, pattern = "^MT-", col.name = "percent.mt")
myeloid[["diagnosis"]] <- myeloid_barcodes$V2
myeloid=subset(myeloid, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 10)
myeloid <- SCTransform(myeloid, method = "glmGamPoi", vars.to.regress = "percent.mt")
myeloid <- RunPCA(myeloid)
myeloid <- RunUMAP(myeloid, dims = 1:20)
myeloid <- FindNeighbors(myeloid, dims = 1:20)
myeloid <- FindClusters(myeloid,resolution=0.1)
summary(myeloid$seurat_clusters)
myeloid.markers <- FindAllMarkers(myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pdf("ets2.vln.7clusters.res0.1.mito10.scaled.dim40.pdf") 
par(mar=c(4,4,4,4))
VlnPlot(myeloid,features=c("ETS2"),slot="scale.data", pt.size=0.02)
dev.off()
pdf("umap.7clusters.res0.1.mito10.pdf")+ xlim(-15,10)+ylim(-10,10)
par(mar=c(4,4,4,4))
DimPlot(myeloid,reduction="umap")
dev.off()


