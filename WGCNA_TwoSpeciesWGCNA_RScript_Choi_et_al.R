##### WGCNA =============================================================================
##### install packages ==================================================================
# BiocManager::install("impute", force = TRUE) ## for WGCNA installation
# BiocManager::install("preprocessCore", force = TRUE) ## for WGCNA installation
# install.packages("ggh4x")

##### libraries =========================================================================
library(WGCNA)
library(tidyverse)
library(dplyr) ## %>%, rename
library(stringr) ## filter
library(magrittr)
library(ggplot2)
library(readxl) ## read.xlsx
library(readr)
library(rtracklayer) ## import for gtf
library(ape) ## read.gff
library(patchwork) ## wrap_plots
library(gridExtra) ## arrangeGrob
library(flashClust)
library(pheatmap)
library(cowplot)
library(clusterProfiler)
library(ggalluvial)
library(simplifyEnrichment)
library(ggh4x)
library(ggforce)
library(reshape2)
library(igraph)
# library(xlsx)

##### Step VII. WGCNA ==========================================================
#### Step VII-a. Log2(TPM+1) Matrix with entire genes (S+M model) ==============
Data_Log2TPM <- read.csv("Data_Log2TPM.csv")

### Step VI-a. Removal of genes with 'sd=0' in calculating correlations ========
gene_vars <- apply(Data_Log2TPM[,-c(1:3)], 1, var, na.rm = TRUE)
names(gene_vars) <- Data_Log2TPM$Geneid
zero_var_genes <- names(gene_vars)[gene_vars == 0]
length(zero_var_genes) ## 3
zero_var_genes ## MtrunA17_Chr3g0097221 | MtrunA17_Chr6g0475411 | MtrunA17_Chr5g0421231

Data_Log2TPM_2 <- Data_Log2TPM %>%
  filter(!str_detect(Geneid, "MtrunA17_Chr3g0097221")) %>%
  filter(!str_detect(Geneid, "MtrunA17_Chr6g0475411")) %>%
  filter(!str_detect(Geneid, "MtrunA17_Chr5g0421231"))

WGCNA_SM <- Data_Log2TPM_2[,-c(1:3)] ## 50 libraries

#### Step VII-a-1. Running WGCNA ===============================================
### Step VII-a-2. Constructing a sample network for ourlier detection ===========
# Sample network based on squared Euclidean distances
adjmat_SM = adjacency(WGCNA_SM, type = "distance")
# this calculates the whole network connectivity
NetCon_SM = as.numeric(apply(adjmat_SM, 2, sum)) - 1
# standardized connectivity (z)
StdCon_Z_SM = scale(NetCon_SM)
# Designate samples as outlying if their Z.k value is below the threshold
thresholdZ.k = -5 ## Often use -2.5
# the color vector indicates outlyingness (red)
outlierColor_SM = ifelse(StdCon_Z_SM < thresholdZ.k, "red", "black")
# calculate the cluster tree using flahsClust or hclust
sampleTree_SM = flashClust(as.dist(1 - adjmat_SM), method = "average")
datColors_SM = data.frame(outlierC = outlierColor_SM)
# Plot the sample dendrogram with outliers =====================================
jpeg("Sample_dendrogram_SM_model.jpg", units="cm", width=12, height=6, res=600)
par(cex.main = 0.5, cex.axis = 0.5, cex.lab = 0.5,
    mar = c(0.1, 0.1, 0.1, 0.1))
plotDendroAndColors(sampleTree_SM,
                    groupLabels = names(datColors_SM),
                    colors = datColors_SM, 
                    main = "Sample dendrogram and heatmap (S+M model) (z.k = '-5')",
                    cex.colorLabels = 0.5, cex.dendroLabels = 0.5,
                    cex.main = 0.5)
dev.off()

### Step VII-a-3. Choosing the soft-threshold (beta) ===========================
# Choose a set of soft threshold powers
powers = c(c(1:20), seq(from = 22, to=30, by=2))
# Choose power based on SFT criterion
sft_SM = pickSoftThreshold(t(WGCNA_SM), powerVector = powers) 
# Plot the results; SFT index as a function of different power
jpeg("Soft-thresholing-power_SM_model.jpg", units="cm", width=10, height=5, res=600)
par(mfrow = c(1,2), cex = 0.5);
plot(sft_SM$fitIndices[,1], -sign(sft_SM$fitIndices[,3])*sft_SM$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_SM$fitIndices[,1], -sign(sft_SM$fitIndices[,3])*sft_SM$fitIndices[,2],
     labels=powers,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft_SM$fitIndices[,1], sft_SM$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft_SM$fitIndices[,1], sft_SM$fitIndices[,5], labels=powers,col="red")
dev.off() # Based on the figure, we can choose beta = 5

### Step VII-a-4. Module detection =============================================
mergingTresh = 0.25
cor <- WGCNA::cor ## to fix correlation in blockwiseModules
net_SM <- WGCNA::blockwiseModules(t(WGCNA_SM), corType = "pearson",
                                    maxBlockSize = 50000, networkType = "unsigned", 
                                    power = 5, minModuleSize = 30,
                                    mergeCutHeight = mergingTresh, numericLabels = TRUE,
                                    saveTOMs = TRUE, pamRespectsDendro = FALSE,
                                    saveTOMFileBase = "Complete_HiCor_TOM")
# cor <- stats::cor
# Convert labels to colors for plotting
moduleLabels_SM <- net_SM$colors
moduleColors_SM <- labels2colors(moduleLabels_SM)
length(moduleColors_SM) 
unique(moduleColors_SM) ## 20
# [1] "blue"         "black"        "turquoise"    "yellow"       "green"       
# [6] "grey"         "salmon"       "pink"         "lightgreen"   "magenta"     
# [11] "lightcyan"    "purple"       "red"          "brown"        "tan"         
# [16] "greenyellow"  "cyan"         "grey60"       "lightyellow"  "midnightblue"

# Select good genes for dendrogram
moduleColors_SM_trimmed <- moduleColors_SM[net_SM$goodGenes] 
unique(moduleColors_SM_trimmed) ## 20

# Select bad genes for dendrogram
badGenes <- rownames(WGCNA_SM)[!net_SM$goodGenes]

# Dendrogram
Dendro <- net_SM$dendrograms[[1]]

### Step VII-a-5. Plot the dendrogram and the module colors underneath =========
jpeg("Hierarchical_cluster_tree_SM_model.jpg", units="cm", width=12, height=6, res=600)
par(cex.main = 0.5, cex.axis = 0.5, cex.lab = 0.5)
plotDendroAndColors (Dendro,
                     colors=moduleColors_SM_trimmed,
                     groupLabels=c("Module colors"),
                     main = "Global Gene Dendrogram (All Modules)",
                     dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05,
                     cex.colorLabels = 0.4, cex.dendroLabels = 0.4)
dev.off()

### Step VII-a-6. Module membership (MM) a quantitative measure; ===============
MEs_SM <- moduleEigengenes(t(WGCNA_SM), moduleColors_SM, excludeGrey = FALSE)$eigengenes
MEs_SM <- orderMEs(MEs_SM)
names(MEs_SM)

# Plot the relationships among the eigengenes and the trait
jpeg("Eigengene_networks_SM_model.jpg", units="cm", width=8, height=12, res=600)
par(cex.main = 0.1, cex.axis = 0.5, cex = 0.6,
    mar = c(0.1, 0.1, 0.1, 0.1))
plotEigengeneNetworks (MEs_SM,"", marDendro =c(0,4,1,2),
                       marHeatmap=c(3,4,1,2), xLabelsAngle=90)
dev.off()

# Calculate the module membership values (aka. module eigengene based connectivity)
datKME_SM <- signedKME(t(WGCNA_SM), MEs_SM)
colorOfColumn = substring(names(datKME_SM), 4)

# Pairwise-corrleration of module membership between all modules
summary(datKME_SM) ## There are some NAs
kME_matrix_SM <- cor(datKME_SM, use = "pairwise.complete.obs")

# kME_plot_SM <- pheatmap(kME_matrix_SM, 
#                           cluster_rows = FALSE,  # Order module names or colors
#                           cluster_cols = FALSE, 
#                           display_numbers = TRUE,
#                           main = "Module kME correlations (S+M model)",
#                           fontsize = 7, fontsize_number = 7,
#                           cellwidth = 18, cellheight = 18,
#                           filename = "kME_correlations_SM_model.tiff",
#                           width = 8, height = 8)

kME_df <- kME_matrix_SM %>%
  as.data.frame() %>%
  rownames_to_column(var = "Module_row") %>%
  pivot_longer(
    cols = -Module_row,
    names_to = "Module_col",
    values_to = "ME"
  )

bk <- seq(-1, 1, length.out = 101)

kME_matrix_SM_df <- (kME_matrix_SM + 1) / 2

pheatmap(
  kME_matrix_SM_df,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  labels_row = gsub("^kME", "", rownames(kME_matrix_SM)),
  labels_col = gsub("^kME", "", colnames(kME_matrix_SM)),
  main = "Module Eigengene Network",
  fontsize = 10,
  fontsize_number = 8,
  cellwidth = 15,
  cellheight = 15,
  color = colorRampPalette(c("blue","white","red"))(100),
  filename = "kME_correlations_SM_model.jpg", ## tiff
  width = 10,
  height = 10
)


### Step VII-a-7. Intramodular Connectivity (kWithin, kTotal, kOut) to find hub gene ==========
adjmat_SM <- adjacency(t(WGCNA_SM), type = "distance")
IMConn_SM <- intramodularConnectivity(adjmat_SM, moduleColors_SM)

hubGenes_SM <- by(IMConn_SM, moduleColors_SM, function(moduleData) {
  moduleData[which.max(moduleData$kWithin), ] # Hub gene has the largest kWithin
})

hubGenes_SM

### Step VII-a-8. Export TOM for Cytoscape =====================================
## For each module; 
# [1] "blue"         "black"        "turquoise"    "yellow"       "green"       
# [6] "grey"         "salmon"       "pink"         "lightgreen"   "magenta"     
# [11] "lightcyan"    "purple"       "red"          "brown"        "tan"         
# [16] "greenyellow"  "cyan"         "grey60"       "lightyellow"  "midnightblue"

SM_Exp <- as.data.frame(Data_Log2TPM_2[,-c(2:3)])
rownames(SM_Exp) <- SM_Exp[,1]
SM_Exp <- SM_Exp[, -1]

datExpr_SM <- t(as.matrix(SM_Exp))
storage.mode(datExpr_SM) <- "double"

## Blue (to midnightblue)
inModule <- moduleColors_SM == "blue"
modExpr_SM_blue  <- datExpr_SM[, inModule]

adjmat_SM_blue  <- adjacency(modExpr_SM_blue, power = 5, type = "unsigned")
modTOM_SM_blue <- TOMsimilarity(adjmat_SM_blue, TOMType = "unsigned")

exportNetworkToCytoscape(modTOM_SM_blue,
                         edgeFile = "WGCNA_Edges_SM_Blue.txt",
                         nodeFile = "WGCNA_Nodes_SM_Blue.txt",
                         weighted = TRUE,
                         threshold = 0.3,
                         nodeNames = colnames(modExpr_SM_blue),
                         nodeAttr  = moduleColors_SM[inModule])

### Step VII-a-9. Adding annotations ===========================================
dim(anno_raw)
dim(datKME_SM) ## 26425 20
dim(IMConn_SM) ## 26425 4
dim(SM_Exp) ## 26425 50
ncol(SM_Exp) ## 50

# Exporting data
summary_table_SM = data.frame(
  GeneID     = Data_Log2TPM_2$Geneid,
  Module     = moduleColors_SM,
  IMConn_SM,
  datKME_SM
)

summary_table_SM <- summary_table_SM %>%
  left_join(Data_Log2TPM %>% select(-Chr, -coding_length), by = c("GeneID" = "Geneid"))

str(summary_table_SM)
# write.csv(summary_table_SM, "Summary_table_SM_model.csv")

### End of Step VII-a  =========================================================

#### Step VII-b. Log2(TPM+1) Matrix with only Sinorhizobia genes (S model) =====
Data_Log2TPM_2_Sino <- Data_Log2TPM_2 %>%
  filter(str_detect(Geneid, "CDO30"))

WGCNA_S <- Data_Log2TPM_2_Sino[,-c(1:3)]

#### Step VII-b-1. Running WGCNA ===============================================
### Step VII-b-2. Constructing a sample network for ourlier detection ==========
# Sample network based on squared Euclidean distances
adjmat_S = adjacency(WGCNA_S, type = "distance")
# this calculates the whole network connectivity
NetCon_S = as.numeric(apply(adjmat_S, 2, sum)) - 1
# standardized connectivity (z)
StdCon_Z_S = scale(NetCon_S)
# Designate samples as outlying if their Z.k value is below the threshold
thresholdZ.k = -5 ## Often use -2.5
# the color vector indicates outlyingness (red)
outlierColor_S = ifelse(StdCon_Z_S < thresholdZ.k, "red", "black")
# calculate the cluster tree using flahsClust or hclust
sampleTree_S = flashClust(as.dist(1 - adjmat_S), method = "average")
datColors_S = data.frame(outlierC = outlierColor_S)
# Plot the sample dendrogram with outliers =====================================
jpeg("Sample_dendrogram_S_model.jpg", units="cm", width=12, height=6, res=600)
par(cex.main = 0.5, cex.axis = 0.5, cex.lab = 0.5,
    mar = c(0.1, 0.1, 0.1, 0.1))
plotDendroAndColors(sampleTree_S,
                    groupLabels = names(datColors_S),
                    colors = datColors_S, 
                    main = "Sample dendrogram and heatmap (S model) (z.k = '-5')",
                    cex.colorLabels = 0.5, cex.dendroLabels = 0.5,
                    cex.main = 0.5)
dev.off()

### Step VII-b-3. Choosing the soft-threshold (beta) ===========================
# Choose a set of soft threshold powers
powers = c(c(1:20), seq(from = 22, to=30, by=2))
# Choose power based on SFT criterion
sft_S = pickSoftThreshold(t(WGCNA_S), powerVector = powers) 
# Plot the results; SFT index as a function of different power
jpeg("Soft-thresholing-power_S_model.jpg", units="cm", width=10, height=5, res=600)
par(mfrow = c(1,2), cex = 0.5);
plot(sft_S$fitIndices[,1], -sign(sft_S$fitIndices[,3])*sft_S$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_S$fitIndices[,1], -sign(sft_S$fitIndices[,3])*sft_S$fitIndices[,2],
     labels=powers,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft_S$fitIndices[,1], sft_S$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft_S$fitIndices[,1], sft_S$fitIndices[,5], labels=powers,col="red")
dev.off() # Based on the figure, we can choose beta = 5

### Step VII-b-4. Module detection =============================================
mergingTresh = 0.25
cor <- WGCNA::cor ## to fix correlation in blockwiseModules
net_S <- WGCNA::blockwiseModules(t(WGCNA_S), corType = "pearson",
                                 maxBlockSize = 50000, networkType = "unsigned", 
                                 power = 7, minModuleSize = 30,
                                 mergeCutHeight = mergingTresh, numericLabels = TRUE,
                                 saveTOMs = TRUE, pamRespectsDendro = FALSE,
                                 saveTOMFileBase = "Complete_HiCor_TOM")
# cor <- stats::cor
# Convert labels to colors for plotting
moduleLabels_S <- net_S$colors
moduleColors_S <- labels2colors(moduleLabels_S)
length(moduleColors_S) 
unique(moduleColors_S) ## 8
# "green"     "turquoise" "blue"      "brown"     "grey"      "red"       "yellow"    "black"    

# Select good genes for dendrogram
moduleColors_S_trimmed <- moduleColors_S[net_S$goodGenes] 
unique(moduleColors_S_trimmed) ## 20

# Select bad genes for dendrogram
badGenes <- rownames(WGCNA_S)[!net_S$goodGenes]

# Dendrogram
Dendro <- net_S$dendrograms[[1]]

### Step VII-b-5. Plot the dendrogram and the module colors underneath =========
jpeg("Hierarchical_cluster_tree_S_model.jpg", units="cm", width=12, height=6, res=600)
par(cex.main = 0.5, cex.axis = 0.5, cex.lab = 0.5)
plotDendroAndColors (Dendro,
                     colors=moduleColors_S_trimmed,
                     groupLabels=c("Module colors"),
                     main = "Global Gene Dendrogram (All Modules)",
                     dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05,
                     cex.colorLabels = 0.4, cex.dendroLabels = 0.4)
dev.off()

### Step VII-b-6. Module membership (MM) a quantitative measure; ===============
MEs_S <- moduleEigengenes(t(WGCNA_S), moduleColors_S, excludeGrey = FALSE)$eigengenes
MEs_S <- orderMEs(MEs_S)
names(MEs_S)

# Plot the relationships among the eigengenes and the trait
jpeg("Eigengene_networks_S_model.jpg", units="cm", width=8, height=12, res=600)
par(cex.main = 0.1, cex.axis = 0.5, cex = 0.6,
    mar = c(0.1, 0.1, 0.1, 0.1))
plotEigengeneNetworks (MEs_S,"", marDendro =c(0,4,1,2),
                       marHeatmap=c(3,4,1,2), xLabelsAngle=90)
dev.off()

# Calculate the module membership values (aka. module eigengene based connectivity)
datKME_S <- signedKME(t(WGCNA_S), MEs_S)
colorOfColumn = substring(names(datKME_S), 4)

# Pairwise-corrleration of module membership between all modules
summary(datKME_S) ## There are some NAs
kME_matrix_S <- cor(datKME_S, use = "pairwise.complete.obs")

kME_df <- kME_matrix_S %>%
  as.data.frame() %>%
  rownames_to_column(var = "Module_row") %>%
  pivot_longer(
    cols = -Module_row,
    names_to = "Module_col",
    values_to = "ME"
  )

bk <- seq(-1, 1, length.out = 101)

kME_matrix_S_df <- (kME_matrix_S + 1) / 2

pheatmap(
  kME_matrix_S_df,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  labels_row = gsub("^kME", "", rownames(kME_matrix_S)),
  labels_col = gsub("^kME", "", colnames(kME_matrix_S)),
  main = "Module Eigengene Network",
  fontsize = 10,
  fontsize_number = 8,
  cellwidth = 15,
  cellheight = 15,
  color = colorRampPalette(c("blue","white","red"))(100),
  filename = "kME_correlations_S_model.jpg", ## tiff
  width = 10,
  height = 10
)


### Step VII-b-7. Intramodular Connectivity (kWithin, kTotal, kOut) to find hub gene ==========
adjmat_S <- adjacency(t(WGCNA_S), type = "distance")
IMConn_S <- intramodularConnectivity(adjmat_S, moduleColors_S)

hubGenes_S <- by(IMConn_S, moduleColors_S, function(moduleData) {
  moduleData[which.max(moduleData$kWithin), ] # Hub gene has the largest kWithin
})

hubGenes_S

### Step VII-b-8. Export TOM for Cytoscape =====================================
## For each module; 
# [1] "blue"         "black"        "turquoise"    "yellow"       "green"       
# [6] "grey"         "salmon"       "pink"         "lightgreen"   "magenta"     
# [11] "lightcyan"    "purple"       "red"          "brown"        "tan"         
# [16] "greenyellow"  "cyan"         "grey60"       "lightyellow"  "midnightblue"

S_Exp <- as.data.frame(Data_Log2TPM_2_Sino[,-c(2:3)])
rownames(S_Exp) <- S_Exp[,1]
S_Exp <- S_Exp[, -1]

datExpr_S <- t(as.matrix(S_Exp))
storage.mode(datExpr_S) <- "double"

## Blue (to midnightblue)
inModule <- moduleColors_S == "blue"
modExpr_S_blue  <- datExpr_S[, inModule]

adjmat_S_blue  <- adjacency(modExpr_S_blue, power = 5, type = "unsigned")
modTOM_S_blue <- TOMsimilarity(adjmat_S_blue, TOMType = "unsigned")

exportNetworkToCytoscape(modTOM_S_blue,
                         edgeFile = "WGCNA_Edges_S_Blue.txt",
                         nodeFile = "WGCNA_Nodes_S_Blue.txt",
                         weighted = TRUE,
                         threshold = 0.3,
                         nodeNames = colnames(modExpr_S_blue),
                         nodeAttr  = moduleColors_S[inModule])

### Step VII-b-9. Adding annotations ===========================================
dim(anno_raw)
dim(datKME_S) ## 5194 8
dim(IMConn_S) ## 5194 4
dim(S_Exp) ## 5194 50
ncol(S_Exp) ## 50

# Exporting data
summary_table_S = data.frame(
  GeneID     = Data_Log2TPM_2_Sino$Geneid,
  Module     = moduleColors_S,
  IMConn_S,
  datKME_S
)

summary_table_S <- summary_table_S %>%
  left_join(Data_Log2TPM_2_Sino %>% select(-Chr, -coding_length), by = c("GeneID" = "Geneid"))

str(summary_table_S)
# write.csv(summary_table_S, "Summary_table_S_model.csv")

### End of Step VII-b  =========================================================

#### Step VII-c. Log2(TPM+1) Matrix with only Medicago genes (M model) =====
Data_Log2TPM_2_Medi <- Data_Log2TPM_2 %>%
  filter(str_detect(Geneid, "MtrunA"))

WGCNA_M <- Data_Log2TPM_2_Medi[,-c(1:3)]

#### Step VII-c-1. Running WGCNA ===============================================
### Step VII-c-2. Constructing a sample network for ourlier detection ==========
# Sample network based on squared Euclidean distances
adjmat_M = adjacency(WGCNA_M, type = "distance")
# this calculates the whole network connectivity
NetCon_M = as.numeric(apply(adjmat_M, 2, sum)) - 1
# standardized connectivity (z)
StdCon_Z_M = scale(NetCon_M)
# Designate samples as outlying if their Z.k value is below the threshold
thresholdZ.k = -5 ## Often use -2.5
# the color vector indicates outlyingness (red)
outlierColor_M = ifelse(StdCon_Z_M < thresholdZ.k, "red", "black")
# calculate the cluster tree using flahsClust or hclust
sampleTree_M = flashClust(as.dist(1 - adjmat_M), method = "average")
datColors_M = data.frame(outlierC = outlierColor_M)
# Plot the sample dendrogram with outliers =====================================
jpeg("Sample_dendrogram_M_model.jpg", units="cm", width=12, height=6, res=600)
par(cex.main = 0.5, cex.axis = 0.5, cex.lab = 0.5,
    mar = c(0.1, 0.1, 0.1, 0.1))
plotDendroAndColors(sampleTree_M,
                    groupLabels = names(datColors_M),
                    colors = datColors_M, 
                    main = "Sample dendrogram and heatmap (M model) (z.k = '-5')",
                    cex.colorLabels = 0.5, cex.dendroLabels = 0.5,
                    cex.main = 0.5)
dev.off()

### Step VII-c-3. Choosing the soft-threshold (beta) ===========================
# Choose a set of soft threshold powers
powers = c(c(1:20), seq(from = 22, to=30, by=2))
# Choose power based on SFT criterion
sft_M = pickSoftThreshold(t(WGCNA_M), powerVector = powers) 
# Plot the results; SFT index as a function of different power
jpeg("Soft-thresholing-power_M_model.jpg", units="cm", width=10, height=5, res=600)
par(mfrow = c(1,2), cex = 0.5);
plot(sft_M$fitIndices[,1], -sign(sft_M$fitIndices[,3])*sft_M$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_M$fitIndices[,1], -sign(sft_M$fitIndices[,3])*sft_M$fitIndices[,2],
     labels=powers,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft_M$fitIndices[,1], sft_M$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft_M$fitIndices[,1], sft_M$fitIndices[,5], labels=powers,col="red")
dev.off() # Based on the figure, we can choose beta = 5

### Step VII-c-4. Module detection =============================================
mergingTresh = 0.25
cor <- WGCNA::cor ## to fix correlation in blockwiseModules
net_M <- WGCNA::blockwiseModules(t(WGCNA_M), corType = "pearson",
                                 maxBlockSize = 50000, networkType = "unsigned", 
                                 power = 5, minModuleSize = 30,
                                 mergeCutHeight = mergingTresh, numericLabels = TRUE,
                                 saveTOMs = TRUE, pamRespectsDendro = FALSE,
                                 saveTOMFileBase = "Complete_HiCor_TOM")
# cor <- stats::cor
# Convert labels to colors for plotting
moduleLabels_M <- net_M$colors
moduleColors_M <- labels2colors(moduleLabels_M)
length(moduleColors_M) 
unique(moduleColors_M) ## 17
# "blue"         "turquoise"    "greenyellow"  "grey"         "green"        "brown"        "red"         
# "tan"          "lightcyan"    "midnightblue" "pink"         "purple"       "magenta"      "yellow"      
# "black"        "cyan"         "salmon" 

# Select good genes for dendrogram
moduleColors_M_trimmed <- moduleColors_M[net_M$goodGenes] 
unique(moduleColors_M_trimmed) ## 20

# Select bad genes for dendrogram
badGenes <- rownames(WGCNA_M)[!net_M$goodGenes]

# Dendrogram
Dendro <- net_M$dendrograms[[1]]

### Step VII-c-5. Plot the dendrogram and the module colors underneath =========
jpeg("Hierarchical_cluster_tree_M_model.jpg", units="cm", width=12, height=6, res=600)
par(cex.main = 0.5, cex.axis = 0.5, cex.lab = 0.5)
plotDendroAndColors (Dendro,
                     colors=moduleColors_M_trimmed,
                     groupLabels=c("Module colors"),
                     main = "Global Gene Dendrogram (All Modules)",
                     dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05,
                     cex.colorLabels = 0.4, cex.dendroLabels = 0.4)
dev.off()

### Step VII-c-6. Module membership (MM) a quantitative measure; ===============
MEs_M <- moduleEigengenes(t(WGCNA_M), moduleColors_M, excludeGrey = FALSE)$eigengenes
MEs_M <- orderMEs(MEs_M)
names(MEs_M)

# Plot the relationships among the eigengenes and the trait
jpeg("Eigengene_networks_M_model.jpg", units="cm", width=8, height=12, res=600)
par(cex.main = 0.1, cex.axis = 0.5, cex = 0.6,
    mar = c(0.1, 0.1, 0.1, 0.1))
plotEigengeneNetworks (MEs_M,"", marDendro =c(0,4,1,2),
                       marHeatmap=c(3,4,1,2), xLabelsAngle=90)
dev.off()

# Calculate the module membership values (aka. module eigengene based connectivity)
datKME_M <- signedKME(t(WGCNA_M), MEs_M)
colorOfColumn = substring(names(datKME_M), 4)

# Pairwise-corrleration of module membership between all modules
summary(datKME_M) ## There are some NAs
kME_matrix_M <- cor(datKME_M, use = "pairwise.complete.obs")

kME_df <- kME_matrix_M %>%
  as.data.frame() %>%
  rownames_to_column(var = "Module_row") %>%
  pivot_longer(
    cols = -Module_row,
    names_to = "Module_col",
    values_to = "ME"
  )

bk <- seq(-1, 1, length.out = 101)

kME_matrix_M_df <- (kME_matrix_M + 1) / 2

pheatmap(
  kME_matrix_M_df,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  labels_row = gsub("^kME", "", rownames(kME_matrix_M)),
  labels_col = gsub("^kME", "", colnames(kME_matrix_M)),
  main = "Module Eigengene Network",
  fontsize = 10,
  fontsize_number = 8,
  cellwidth = 15,
  cellheight = 15,
  color = colorRampPalette(c("blue","white","red"))(100),
  filename = "kME_correlations_M_model.jpg", ## tiff
  width = 10,
  height = 10
)


### Step VII-c-7. Intramodular Connectivity (kWithin, kTotal, kOut) to find hub gene ==========
adjmat_M <- adjacency(t(WGCNA_M), type = "distance")
IMConn_M <- intramodularConnectivity(adjmat_M, moduleColors_M)

hubGenes_M <- by(IMConn_M, moduleColors_M, function(moduleData) {
  moduleData[which.max(moduleData$kWithin), ] # Hub gene has the largest kWithin
})

hubGenes_M

### Step VII-c-8. Export TOM for Cytoscape =====================================
## For each module; 
# [1] "blue"         "black"        "turquoise"    "yellow"       "green"       
# [6] "grey"         "salmon"       "pink"         "lightgreen"   "magenta"     
# [11] "lightcyan"    "purple"       "red"          "brown"        "tan"         
# [16] "greenyellow"  "cyan"         "grey60"       "lightyellow"  "midnightblue"

M_Exp <- as.data.frame(Data_Log2TPM_2_Medi[,-c(2:3)])
rownames(M_Exp) <- M_Exp[,1]
M_Exp <- M_Exp[, -1]

datExpr_M <- t(as.matrix(M_Exp))
storage.mode(datExpr_M) <- "double"

## Blue (to midnightblue)
inModule <- moduleColors_M == "blue"
modExpr_M_blue  <- datExpr_M[, inModule]

adjmat_M_blue  <- adjacency(modExpr_M_blue, power = 5, type = "unsigned")
modTOM_M_blue <- TOMsimilarity(adjmat_M_blue, TOMType = "unsigned")

exportNetworkToCytoscape(modTOM_M_blue,
                         edgeFile = "WGCNA_Edges_M_Blue.txt",
                         nodeFile = "WGCNA_Nodes_M_Blue.txt",
                         weighted = TRUE,
                         threshold = 0.3,
                         nodeNames = colnames(modExpr_M_blue),
                         nodeAttr  = moduleColors_M[inModule])

### Step VII-c-9. Adding annotations ===========================================
dim(anno_raw)
dim(datKME_M) ##
dim(IMConn_M) ##
dim(M_Exp) ##
ncol(M_Exp) ## 50

# Exporting data
summary_table_M = data.frame(
  GeneID     = Data_Log2TPM_2_Medi$Geneid,
  Module     = moduleColors_M,
  IMConn_M,
  datKME_M
)

summary_table_M <- summary_table_M %>%
  left_join(Data_Log2TPM_2_Medi %>% select(-Chr, -coding_length), by = c("GeneID" = "Geneid"))

str(summary_table_M)
# write.csv(summary_table_M, "Summary_table_M_model.csv")

### End of Step VII-c  =========================================================

##### Step VIII. Co-expression patterns (Figure 3 and S4; Pie Charts) ==========
#### Step VIII-a. All genes ====================================================
### Step VIII-a-1. Importing data ==============================================
modules <- read.csv("Modules_SM_model.CSV")
df_clean <- modules %>% mutate(Total = coalesce(Total, Sino + Medi))

### Step VIII-a-2. Labelling ===================================================
# Module arrangement + facet labelling
totals <- df_clean %>% select(Module, Total)
order_vec <- totals %>% arrange(desc(Total)) %>% pull(Module)
label_map <- totals %>% mutate(strip = paste0(Module, " (n=", Total, ")")) %>%
  { setNames(.$strip, .$Module) }

df_sel <- df_clean %>%
  group_by(Module) %>%
  ungroup()

module_order <- c(
  "Turquoise", "Blue", "Grey", "Brown", "Yellow", "Green", "Red", "Black",
  "Pink", "Magenta", "Purple", "Greenyellow", "Tan", "Salmon", "Cyan",
  "Midnightblue", "Lightcyan", "Grey60", "Lightgreen", "Lightyellow"
)

df_sel$Module <- factor(df_sel$Module, levels = module_order)

# Long-form transformation with proportion calculation
df_pie <- df_sel %>%
  pivot_longer(c(Sino, Medi), names_to = "Group", values_to = "Count") %>%
  group_by(Module) %>%
  mutate(
    Total = sum(Count, na.rm = TRUE),
    Prop  = if_else(Total > 0, Count / Total, NA_real_),
    Radius = log2(Total + 1)  # diameter
  ) %>%
  ungroup() %>%
  mutate(
    Module_f = factor(Module, levels = order_vec),
    Group    = factor(Group, levels = c("Sino", "Medi"))
  ) %>%
  arrange(Module, Group) %>%
  mutate(
    angle = Prop * 2 * pi,    # angle
    start = c(0, cumsum(angle)[-n()]),
    end   = cumsum(angle),
    r0    = 0,
    r     = Radius
  ) %>%
  ungroup()

# 
prop_wide <- df_pie %>%
  select(Module, Group, Prop) %>%
  tidyr::pivot_wider(names_from = Group, values_from = Prop)

lab_tbl <- totals %>%
  left_join(prop_wide, by = "Module") %>%
  mutate(
    base_lab  = paste0(Module, " (n=", Total, ")"),
    strip_lab = paste0(base_lab))

label_map_dagger <- setNames(lab_tbl$strip_lab, lab_tbl$Module)

module_order <- c(
  "Turquoise", "Blue", "Grey", "Brown", "Yellow", "Green", "Red", "Black",
  "Pink", "Magenta", "Purple", "Greenyellow", "Tan", "Salmon", "Cyan",
  "Midnightblue", "Lightcyan", "Grey60", "Lightgreen", "Lightyellow"
)

df_pie$Module <- factor(df_pie$Module, levels = module_order)

### Step VIII-a-3. Plotting ====================================================
p <- ggplot(df_pie) +
  geom_arc_bar(aes(
    x0 = 0, y0 = 0,
    r0 = r0, r = r,
    start = start, end = end,
    fill = Group
  ), color = "white") +
  facet_wrap(~ Module, labeller = labeller(Module = label_map_dagger), nrow = 4) +
  scale_fill_manual(
    values = c(Sino = "#CC6600", Medi = "#56B4E9"),
    name = NULL
  ) +
  coord_fixed() +
  theme_void(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5,
                              margin = margin(b = 10)),
    plot.margin  = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    panel.spacing = unit(10, "pt"),
    legend.position = "bottom", 
    legend.direction = "horizontal",
    legend.background = element_rect(fill="white", color="white"),
    legend.key.size = grid::unit(10, "pt"),
    legend.text = element_text(size = 10, family = "Arial"),
    strip.text   = element_text(family = "Arial"),
    aspect.ratio = 1
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

p

# ggsave("Piechart_modules_SM.jpg", p, height = 7, width = 8, bg ="white", dpi = 480)

#### Step VIII-b. Hub genes ====================================================
### Step VIII-b-1. Importing data ==============================================
Hubmodules <- read.csv("Modules_HubGenes_SM_model.csv")
df_clean <- Hubmodules %>% mutate(Total = coalesce(Total, Sino + Medi))

### Step VIII-a-2. Labelling ===================================================
# Module arrangement + facet labelling
totals <- df_clean %>% select(Module, Total)
order_vec <- totals %>% arrange(desc(Total)) %>% pull(Module)
label_map <- totals %>% mutate(strip = paste0(Module, " (n=", Total, ")")) %>%
  { setNames(.$strip, .$Module) }

df_sel <- df_clean %>%
  group_by(Module) %>%
  ungroup()

module_order <- c(
  "Turquoise", "Blue", "Grey", "Brown", "Yellow", "Green", "Red", "Black",
  "Pink", "Magenta", "Purple", "Greenyellow", "Tan", "Salmon", "Cyan",
  "Midnightblue", "Lightcyan", "Grey60", "Lightgreen", "Lightyellow"
)

df_sel$Module <- factor(df_sel$Module, levels = module_order)

# Long-form transformation with proportion calculation
df_pie <- df_sel %>%
  pivot_longer(c(Sino, Medi), names_to = "Group", values_to = "Count") %>%
  group_by(Module) %>%
  mutate(
    Total = sum(Count, na.rm = TRUE),
    Prop  = if_else(Total > 0, Count / Total, NA_real_),
    Radius = log2(Total + 1)  # diameter
  ) %>%
  ungroup() %>%
  mutate(
    Module_f = factor(Module, levels = order_vec),
    Group    = factor(Group, levels = c("Sino", "Medi"))
  ) %>%
  arrange(Module, Group) %>%
  mutate(
    angle = Prop * 2 * pi,    # angle
    start = c(0, cumsum(angle)[-n()]),
    end   = cumsum(angle),
    r0    = 0,
    r     = Radius
  ) %>%
  ungroup()

# 
prop_wide <- df_pie %>%
  select(Module, Group, Prop) %>%
  tidyr::pivot_wider(names_from = Group, values_from = Prop)

lab_tbl <- totals %>%
  left_join(prop_wide, by = "Module") %>%
  mutate(
    base_lab  = paste0(Module, " (n=", Total, ")"),
    strip_lab = paste0(base_lab))

label_map_dagger <- setNames(lab_tbl$strip_lab, lab_tbl$Module)

module_order <- c(
  "Turquoise", "Blue", "Grey", "Brown", "Yellow", "Green", "Red", "Black",
  "Pink", "Magenta", "Purple", "Greenyellow", "Tan", "Salmon", "Cyan",
  "Midnightblue", "Lightcyan", "Grey60", "Lightgreen", "Lightyellow"
)

df_pie$Module <- factor(df_pie$Module, levels = module_order)

### Step VIII-a-3. Plotting ====================================================
p <- ggplot(df_pie) +
  geom_arc_bar(aes(
    x0 = 0, y0 = 0,
    r0 = r0, r = r,
    start = start, end = end,
    fill = Group
  ), color = "white") +
  facet_wrap(~ Module, labeller = labeller(Module = label_map_dagger), nrow = 4) +
  scale_fill_manual(
    values = c(Sino = "#CC6600", Medi = "#56B4E9"),
    name = NULL
  ) +
  coord_fixed() +
  theme_void(base_size = 12) +
  theme(
    text = element_text(family = "Arial"),
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5,
                              margin = margin(b = 10)),
    plot.margin  = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    panel.spacing = unit(10, "pt"),
    legend.position = "bottom", 
    legend.direction = "horizontal",
    legend.background = element_rect(fill="white", color="white"),
    legend.key.size = grid::unit(10, "pt"),
    legend.text = element_text(size = 10, family = "Arial"),
    strip.text   = element_text(family = "Arial"),
    aspect.ratio = 1
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

p

# ggsave("Piechart_Hubmodules_SM.jpg", p, height = 7, width = 8, bg ="white", dpi = 480)


##### End of Step VIII =========================================================
