# This code was used for analysis of CosMx SMI data which had been exported as a count matrix from AtoMx after images were segmented to obtain cell boundaries, transcripts assigned to single cells, and a transcript by cell count matrix obtained

### expMatrix = rownames(cell_id)+colnames(genes)
### metadata = rownames(cell_id) + colnames(annotations: including ifdata)
### negMatrix =  matrix of all negtive probes

#################InsituType for Cell Type Annotation##################

library(dplyr)
library(sp)
library(SeuratObject)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(InSituType)

## calculate the background
negmen <- rownames(negMatrix)
negmen.per.cell <- mean(rowMeans(negMatrix)) / mean(rowSums(expMatrix))
per.cell.bg <- rowSums(expMatrix) * negmean.per.totcount


ifdata <- as.matrix(metadata[,c("Area","AspectRatio","Mean.Membrane",
                                "Mean.PanCK","Mean.CD45","Mean.CD68","Mean.DAPI")])

## loading ioprofile and iocolors
data("ioprofiles") # Liver_HCA_profile
data("iocolors")


cohort = fastCohorting(mat = ifdata, gaussian_transform = TRUE)

# choose cluster number: (optional in semi-s
nclust = chooseClusterNumber(counts = expMatrix,
                             neg = negmean,
                             bg = per.cell.bg,
                             n_clusts = 2:15) 


sup = insitutype(  expMatrix, 
                   neg = negmean,
                   bg = per.cell.bg,
                   n_clusts = nclust$best_clust_number,
                   reference_profiles = ioprofiles, 
                   cohort = cohort,
                   n_phase1 = 1000,
                   n_phase2 = 2000,
                   n_phase3 = 3000,
                   max_iters = 20)


metadata$cell_type <- sup[["clust"]]

### plotting function are done in ggplot2

############## Enrichment Scores Analysis for Module Genes#######################
library(GSVA)
library(rstatix)

## list ETS2 module genes as a gene list
ets2_genelist <- c("ETS2","CCL5","IL1B","MMP9","TLR2",
                   "MMP14","TLR4","IL1A","NLRP3","CXCL5","PTPRC","CCL2",
                   "TLR5","SOD2","CXCL8","TLR8","PTGS2") #exclude CD163, S100A8 and S100A9 as these are used for cell typing macrophage subsets
lst_sig = list()
lst_sig[["ETS2_Module_Signatures"]] = ets2_genelist

ssgsea <- gsva( as.matrix( t(expMatrix) ), lst_sig, method="zscore")
ssgsea_result <- as.data.frame(t(ssgsea))

metadata$Enrichment_Score <- ssgsea_result$ETS2_Module_Signatures



############## Spatial metrics analysis #######################

library(phenoptr)


#metadata generated from the code above
NC_S0 <- read.csv("~/Nanostring GeoMx/CosMx/CosMx/Full Metadata_1201/S0_NC_metadata.csv")
NC_S1 <- read.csv("~/Nanostring GeoMx/CosMx/CosMx/Full Metadata_1201/S1_NC_metadata.csv")

PSC_S0 <- read.csv("~/Nanostring GeoMx/CosMx/CosMx/Full Metadata_1201/S0_PSC_metadata.csv")
PSC_S1 <- read.csv("~/Nanostring GeoMx/CosMx/CosMx/Full Metadata_1201/S1_PSC_metadata.csv")




NC_S0 <- NC_S0 %>% select(cell_id, fov, CenterX_local_px, CenterY_local_px, SemiSup, ETS2_Module_EnrichmentScore)
NC_S1 <- NC_S1 %>% select(cell_id,  fov, CenterX_local_px, CenterY_local_px, SemiSup, ETS2_Module_EnrichmentScore)

PSC_S0 <- PSC_S0 %>% select(cell_id,  fov,CenterX_local_px, CenterY_local_px, SemiSup, ETS2_Module_EnrichmentScore)
PSC_S1 <- PSC_S1 %>% select(cell_id,  fov,CenterX_local_px, CenterY_local_px, SemiSup, ETS2_Module_EnrichmentScore)







NC <- rbind(NC_S0, NC_S1)
PSC <- rbind(PSC_S0, PSC_S1)

group <- "Normal Liver Control"
NC$group <- group

group <- "PSC"
PSC$group <- group

#sum(is.na(NC$cell_id))

# Radius analysis ----

csd_NC <- NC %>%
  dplyr::rename("Cell X Position" = CenterX_local_px,
                "Cell Y Position" = CenterY_local_px,
                "Cell ID" = cell_id,
                "Phenotype" = SemiSup)



csd_NC$`Cell Y Position` <- as.numeric(csd_NC$`Cell Y Position`)
csd_NC$`Cell X Position` <- as.numeric(csd_NC$`Cell X Position`)


pairs <- list(c('Cholangiocytes', 'Inflammatory.macrophages'),
              c('Cholangiocytes', 'Non.inflammatory.macrophages'))


radii <- c(100, 200, 300, 400, 500)


Duct_NC = csd_NC %>% 
  do(count_within_many(., pairs, radius=radii))




## Stacked barchart ----
mac_col <- c("Inflammatory.macrophages" = "#AE3A7B",
             "Non.inflammatory.macrophages" = "#F194C1")


ggplot(Duct_NC, aes(x = radius, y = within_mean, fill = to)) +
  geom_bar(stat = "identity") +
  labs(x = "Radii", y = "Count", title = "Radius - Macrophage count from duct - Control") + 
  theme_minimal() + scale_fill_manual(values = mac_col) 





# PSC


csd_PSC <- PSC %>%
  dplyr::rename("Cell X Position" = CenterX_local_px,
                "Cell Y Position" = CenterY_local_px,
                "Cell ID" = cell_id,
                "Phenotype" = SemiSup)




csd_PSC$`Cell Y Position` <- as.numeric(csd_PSC$`Cell Y Position`)
csd_PSC$`Cell X Position` <- as.numeric(csd_PSC$`Cell X Position`)


# radius analysis for PSC
pairs <- list(c('Cholangiocytes', 'Inflammatory.macrophages'),
              c('Cholangiocytes', 'Non.inflammatory.macrophages'))


radii <- c(100, 200, 300, 400, 500)


Duct_PSC = csd_PSC %>% 
  do(count_within_many(., pairs, radius=radii))



## Stacked barchart ----
mac_col <- c("Inflammatory.macrophages" = "#AE3A7B",
             "Non.inflammatory.macrophages" = "#F194C1")


ggplot(Duct_PSC, aes(x = radius, y = within_mean, fill = to)) +
  geom_bar(stat = "identity") +
  labs(x = "Radii", y = "Count", title = "Radius - Macrophage count from duct - PSC") + 
  theme_minimal() + scale_fill_manual(values = mac_col) 





# Nearest neighbour -----
## keep the distance boundaries looking at as similar as possible to radius 

## Nearest neighbor analysis control -----


threshold_NC <- csd_NC %>%
  group_by(group, fov) %>% 
  do(bind_cols(., find_nearest_distance(.)))  


Average <- threshold_NC %>% group_by(Phenotype) %>% 
  select(Phenotype, starts_with('Distance to')) %>% 
  summarize_all(~round(mean(. ,na.rm = TRUE), 1))



threshold_NC_macrophages <- threshold_NC %>% select(`Cell ID`, fov, group, `Cell X Position`, `Cell Y Position`, Phenotype, `Cell ID Inflammatory.macrophages`, `Distance to Inflammatory.macrophages`,
                                                    `Cell ID Non.inflammatory.macrophages`, `Distance to Non.inflammatory.macrophages`, ETS2_Module_EnrichmentScore)




# extract duct to macrophage distances 
Average_filterd <- threshold_NC_macrophages %>% filter(Phenotype == "Cholangiocytes")




I_mac <- Average_filterd %>%  select(`Cell ID Inflammatory.macrophages`, `Distance to Inflammatory.macrophages`, `Cell X Position`,
                                        `Cell Y Position`)


NI_mac <- Average_filterd %>%  select( `Cell ID Non.inflammatory.macrophages`, `Distance to Non.inflammatory.macrophages`, `Cell X Position`,
                                     `Cell Y Position`)




Phenotype <- "Inflammatory.macrophages"
I_mac$Phenotype <- Phenotype

Phenotype <- "Non.inflammatory.macrophages"
NI_mac$Phenotype <- Phenotype


I_mac <- I_mac %>%
  dplyr::rename("Distance from duct" = `Distance to Inflammatory.macrophages`,
                "Cell ID"  = `Cell ID Inflammatory.macrophages`
                )


NI_mac <- NI_mac %>%
  dplyr::rename("Distance from duct" = `Distance to Non.inflammatory.macrophages`,
                "Cell ID"  = `Cell ID Non.inflammatory.macrophages`
  )



Duct_to_macrophage <- rbind(I_mac, NI_mac)



# extract ETS2 enrcihment score for the macrophages associated with a ductal cell 
ETS2_macrophages <- threshold_NC_macrophages %>% filter(Phenotype == 'Inflammatory.macrophages' | Phenotype == 'Non.inflammatory.macrophages')

ETS2_macrophages <- ETS2_macrophages %>% select(`Cell ID`, ETS2_Module_EnrichmentScore)


# combine datasets by cell_id

Full_control <- Duct_to_macrophage %>% right_join(ETS2_macrophages, by = "Cell ID")



group <- "Normal Liver Control"
Full_control$group <- group

Full_control <- na.omit(Full_control)


write.csv(Full_control, "S0-S1_distance_ETS2_control.csv")



# Nearest neighbor analysis PSC -----


threshold_PSC <- csd_PSC %>%
  group_by(group, fov) %>% 
  do(bind_cols(., find_nearest_distance(.)))  




Average <- threshold_PSC %>% group_by(Phenotype) %>% 
  select(Phenotype, starts_with('Distance to')) %>% 
  summarize_all(~round(mean(. ,na.rm = TRUE), 1))


# nearest neighbour barcharts 



Average_filterd <- Average %>% filter(Phenotype == "Cholangiocytes")
Average_filterd <- pivot_longer(Average_filterd, 
                                cols = 2:25,  
                                names_to = "Distance.to",
                                values_to = "Distance")

unique(Average_filterd$Distance.to)

mac_col_NN <- c("Distance to Inflammatory.macrophages" = "#AE3A7B",
                "Distance to Non.inflammatory.macrophages" = "#F194C1")

Average_filterd <- Average_filterd %>% filter(Distance.to == "Distance to Inflammatory.macrophages" 
                                              |Distance.to == "Distance to Non.inflammatory.macrophages")


ggplot(Average_filterd, aes(x = Distance.to, y = Distance, fill = Distance.to)) +
  geom_bar(stat = "identity") +
  labs(x = "Nearest Neighbour to Duct", y = "Distance (um)", title = "Nearest Neighbour per macrophage - PSC") + 
  theme_minimal() + scale_fill_manual(values = mac_col_NN) 




# extract duct to macrophage distances 
threshold_PSC_macrophages <- threshold_PSC %>% select(`Cell ID`, fov, `Cell X Position`, `Cell Y Position`, Phenotype, `Cell ID Inflammatory.macrophages`, `Distance to Inflammatory.macrophages`,
                                                    `Cell ID Non.inflammatory.macrophages`, `Distance to Non.inflammatory.macrophages`, ETS2_Module_EnrichmentScore)



Average_filterd <- threshold_PSC_macrophages %>% filter(Phenotype == "Cholangiocytes")



I_mac <- Average_filterd %>%  select(`Cell ID Inflammatory.macrophages`, `Distance to Inflammatory.macrophages`, `Cell X Position`,
                                     `Cell Y Position`)


NI_mac <- Average_filterd %>%  select( `Cell ID Non.inflammatory.macrophages`, `Distance to Non.inflammatory.macrophages`, `Cell X Position`,
                                       `Cell Y Position`)



Phenotype <- "Inflammatory.macrophages"
I_mac$Phenotype <- Phenotype

Phenotype <- "Non.inflammatory.macrophages"
NI_mac$Phenotype <- Phenotype


I_mac <- I_mac %>%
  dplyr::rename("Distance from duct" = `Distance to Inflammatory.macrophages`,
                "Cell ID"  = `Cell ID Inflammatory.macrophages`
  )


NI_mac <- NI_mac %>%
  dplyr::rename("Distance from duct" = `Distance to Non.inflammatory.macrophages`,
                "Cell ID"  = `Cell ID Non.inflammatory.macrophages`
  )



Duct_to_macrophage <- rbind(I_mac, NI_mac)



# extract ETS2 enrcihment score for the macrophages associated with a ductal cell 
ETS2_macrophages <- threshold_PSC_macrophages %>% filter(Phenotype == 'Inflammatory.macrophages' | Phenotype == 'Non.inflammatory.macrophages')

ETS2_macrophages <- ETS2_macrophages %>% select(`Cell ID`, ETS2_Module_EnrichmentScore)


# combine datasets by cell_id

Full_PSC <- Duct_to_macrophage %>% right_join(ETS2_macrophages, by = "Cell ID")



group <- "PSC"
Full_PSC$group <- group

Full_PSC <- na.omit(Full_PSC)






