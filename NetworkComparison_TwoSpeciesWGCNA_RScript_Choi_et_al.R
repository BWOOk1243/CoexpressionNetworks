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

##### Step X. Network comparisons (I) ==========================================
# comparison for single-species model versus two-species models in this study
summary_table_SM <- read.csv("Summary_table_SM.csv")
summary_table_S <- read.csv("Summary_table_S.csv")
summary_table_M <- read.csv("Summary_table_M.csv")
#### Step X-1. (Figure S5) S vs. S+M models ====================================
### Matching assigned modules from models
## Sinorhizobia
Compare_S_SM <- summary_table_SM %>%
  left_join(summary_table_S, by = c("GeneID")) %>%
  rename_with(~ str_replace(.x, "\\.x$", ".sm")) %>%
  rename_with(~ str_replace(.x, "\\.y$", ".s"))

Compare_S_SM_module <- Compare_S_SM %>%
  select(GeneID, Module.sm, Module.s)

## Medicago
Compare_M_SM <- summary_table_SM %>%
  left_join(summary_table_M, by = c("GeneID")) %>%
  rename_with(~ str_replace(.x, "\\.x$", ".sm")) %>%
  rename_with(~ str_replace(.x, "\\.y$", ".m"))

Compare_M_SM_module <- Compare_M_SM %>%
  select(GeneID, Module.sm, Module.m)

## Long-form for Sinorhizobia
flow_R <- Compare_S_SM_module %>%
  count(Module.sm, Module.s) %>% 
  rename(freq = n) %>%
  mutate(Species = "Sino") %>%
  rename(Module = Module.s) %>%
  filter(!is.na(Module))

## Long-form for Medicago
flow_M <- Compare_M_SM_module %>%
  count(Module.sm, Module.m) %>%
  rename(freq = n) %>%
  mutate(Species = "Medi") %>%
  rename(Module = Module.m) %>%
  filter(!is.na(Module))

## Merging two long-form data
flow_all <- bind_rows(flow_R, flow_M) %>%
  complete(Species, Module.sm, Module, fill = list(freq = 0))

## Calculating proportions
flow_all <- flow_all %>%
  group_by(Species, Module) %>%
  mutate(prop = freq / sum(freq)) %>%
  ungroup()

## counting the number
module_sm_counts <- flow_all %>%
  group_by(Species, Module.sm) %>%
  summarise(total_sm_sp = sum(freq), .groups = "drop")

## Adding the count 
flow_all <- flow_all %>%
  left_join(module_sm_counts, by = c("Species", "Module.sm")) %>%
  mutate(
    Module_sm_lab = paste0(Module.sm, " (n=", total_sm_sp, ")")
  )

## Adding the count into each module
module_counts <- bind_rows(
  Compare_S_SM_module %>% count(Species="Sino", Module = Module.s),
  Compare_M_SM_module %>% count(Species="Medi", Module = Module.m)
)

flow_all_full <- flow_all %>%
  left_join(module_counts, by = c("Species", "Module")) %>%
  filter(!is.na(n)) %>%
  mutate(Module_lab = paste0(Module, " (n=", n, ")"))

## Set the threshold to omit display each cell
thr <- 0.1

## Retain rows and columns
keep_rows <- flow_all_full %>%
  dplyr::group_by(Species, Module) %>%
  dplyr::summarise(ok = any(prop >= thr, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(ok) %>%
  dplyr::select(Species, Module)

keep_cols <- flow_all_full %>%
  dplyr::group_by(Module.sm) %>%
  dplyr::summarise(ok = any(prop >= thr, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(ok) %>%
  dplyr::pull(Module.sm)

# Set the filter
flow_filt <- flow_all_full %>%
  dplyr::filter(Module.sm %in% keep_cols) %>%
  dplyr::semi_join(keep_rows, by = c("Species", "Module"))

# Orderging rows and columns
plot_df <- flow_filt %>%
  group_by(Species) %>%
  mutate(
    Module.sm_f = reorder(Module.sm, total_sm_sp, FUN = max),
    Module_f    = reorder(Module, n, FUN = max)              
  ) %>%
  ungroup()

row_levels_Sino <- plot_df %>%
  dplyr::filter(Species == "Sino") %>%
  dplyr::distinct(Module_lab, n) %>%
  dplyr::arrange(n) %>%    
  dplyr::pull(Module_lab)

col_levels_Sino <- plot_df %>%
  filter(Species == "Sino") %>%
  distinct(Module_sm_lab, total_sm_sp) %>%
  arrange(desc(total_sm_sp)) %>%
  pull(Module_sm_lab)

plot_df_Sino <- plot_df %>%
  filter(Species == "Sino") %>%
  mutate(
    Module_sm_lab = factor(Module_sm_lab, levels = col_levels_Sino),
    Module_lab    = factor(Module_lab, levels = rev(row_levels_Sino))
  )

keep_sm_cols_Sino <- plot_df_Sino %>%
  dplyr::group_by(Module.sm) %>%
  dplyr::summarise(
    keep = any(prop >= thr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(keep) %>%
  dplyr::pull(Module.sm)

plot_df_Sino_filt <- plot_df_Sino %>%
  dplyr::filter(Module.sm %in% keep_sm_cols_Sino)

## Plotting
ggplot(plot_df_Sino_filt, aes(x = Module_sm_lab, y = Module_lab, fill = prop)) +
  geom_tile(
    color = "grey30",
    size  = 0.3,         # ← increase border thickness here
    na.rm = FALSE
  ) +
  geom_text(  data = plot_df_Sino %>% dplyr::filter(prop >= 0.1),
              aes(label = sprintf("%.2f", prop)),
              size = 2.5) +
  scale_fill_gradient(
    low      = "white",
    high     = "red",
    na.value = "white",
    name     = "Proportion"
  ) +
  labs(
    x     = "Modules (S+M model)", 
    y     = "Modules (S model)"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    # 1) remove all internal grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 2) draw a thicker border around each facet
    # panel.border     = element_rect(color = "black", fill = NA, size = .5),
    
    # keep your facet‐strip settings
    panel.spacing.y     = unit(0.5, "lines"),
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    strip.placement     = "outside",
    plot.title = element_text(size = rel(0.8)),
    strip.background.y  = element_blank(),
    strip.text.y.left   = element_text(angle = 90),
    legend.title        = element_text(size = 6),
    legend.text         = element_text(size = 6),
    legend.key.size = unit(0.5, "lines")
  )

# ggsave("Heatmap_S_SM_prop.jpg", height = 5.5, width = 6, bg ="white", dpi = 800)

#### Step X-2. (Figure S6) M vs. S+M models ====================================
# Orderging rows and columns
plot_df <- flow_filt %>%
  group_by(Species) %>%
  mutate(
    Module.sm_f = reorder(Module.sm, total_sm_sp, FUN = max),
    Module_f    = reorder(Module, n, FUN = max)              
  ) %>%
  ungroup()

row_levels_Medi <- plot_df %>%
  dplyr::filter(Species == "Medi") %>%
  dplyr::distinct(Module_lab, n) %>%
  dplyr::arrange(n) %>%
  dplyr::pull(Module_lab)

col_levels_Medi <- plot_df %>%
  filter(Species == "Medi") %>%
  distinct(Module_sm_lab, total_sm_sp) %>%
  arrange(desc(total_sm_sp)) %>%
  pull(Module_sm_lab)

plot_df_Medi <- plot_df %>%
  filter(Species == "Medi") %>%
  mutate(
    Module_sm_lab = factor(Module_sm_lab, levels = col_levels_Medi),
    Module_lab    = factor(Module_lab, levels = rev(row_levels_Medi))
  )

keep_sm_cols_Medi <- plot_df_Medi %>%
  dplyr::group_by(Module.sm) %>%
  dplyr::summarise(
    keep = any(prop >= thr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(keep) %>%
  dplyr::pull(Module.sm)

plot_df_Medi_filt <- plot_df_Medi %>%
  dplyr::filter(Module.sm %in% keep_sm_cols_Medi)

## Plotting
ggplot(plot_df_Medi_filt, aes(x = Module_sm_lab, y = Module_lab, fill = prop)) +
  geom_tile(
    color = "grey30",
    size  = 0.3,         # ← increase border thickness here
    na.rm = FALSE
  ) +
  geom_text(  data = plot_df_Medi %>% dplyr::filter(prop >= 0.1),
              aes(label = sprintf("%.2f", prop)),
              size = 2.5) +
  scale_fill_gradient(
    low      = "white",
    high     = "red",
    na.value = "white",
    name     = "Proportion"
  ) +
  labs(
    x     = "Modules (S+M model)", 
    y     = "Modules (M model)"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    # 1) remove all internal grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 2) draw a thicker border around each facet
    # panel.border     = element_rect(color = "black", fill = NA, size = .5),
    
    # keep your facet‐strip settings
    panel.spacing.y     = unit(0.5, "lines"),
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    strip.placement     = "outside",
    plot.title = element_text(size = rel(0.8)),
    strip.background.y  = element_blank(),
    strip.text.y.left   = element_text(angle = 90),
    legend.title        = element_text(size = 6),
    legend.text         = element_text(size = 6),
    legend.key.size = unit(0.5, "lines")
  )

# ggsave("Heatmap_M_SM_prop.jpg", height = 5.5, width = 6, bg ="white", dpi = 800)

#### Step X-3. (Figure S7) Eigengenes correlations b/w S vs. M models ==========
# correlations b/w module eigengenes from single-species and two-species models
MEs_SM <- read.csv("MEs_SM_model.csv")
MEs_S <- read.csv("MEs_S_model.csv")
MEs_M <- read.csv("MEs_M_model.csv")

# re-labelling
colnames(MEs_SM) <- paste0("S+M_", sub("^ME\\.", "", colnames(MEs_SM)))
colnames(MEs_S) <- paste0("M_", sub("^ME\\.", "", colnames(MEs_S)))
colnames(MEs_M) <- paste0("S_", sub("^ME\\.", "", colnames(MEs_M)))

commonSamples <- intersect(rownames(MEs_M), rownames(MEs_S))
MEs_M  <- MEs_M[commonSamples, , drop = FALSE]
MEs_S  <- MEs_S[commonSamples, , drop = FALSE]

# correlation between MEs
corMat <- cor(MEs_M[,-1], MEs_S[,-1], use = "pairwise.complete.obs")
pMat   <- corPvalueStudent(corMat, nSamples = nrow(MEs_M))

textMat <- matrix("", nrow = nrow(corMat), ncol = ncol(corMat),
                  dimnames = dimnames(corMat))
for(i in seq_len(nrow(corMat))) {
  for(j in seq_len(ncol(corMat))) {
    textMat[i,j] <- paste0("r=", sprintf("%.2f", corMat[i,j]),
                           "\np=", sprintf("%.2f", pMat[i,j]))
  }
}

# masing |r| < 0.8 or NA
mask <- abs(corMat) < 0.8
corMat_mask  <- corMat
corMat_mask[mask] <- NA

keepRows <- apply(!is.na(corMat_mask), 1, any)
keepCols <- apply(!is.na(corMat_mask), 2, any)

corMat_sub <- corMat_mask[keepRows, keepCols, drop=FALSE]
textMat_sub <- textMat[keepRows, keepCols, drop=FALSE]

# Set colors
breaksList <- seq(-1, 1, by = 0.2)
heat_colors <- colorRampPalette(c("blue","white","red"))(length(breaksList))

# Plotting heatmap
ph <- pheatmap(
  mat = corMat,
  display_numbers = textMat,   # r, p-value를 텍스트로 표시
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 10,
  fontsize_number = 8,
  number_color = "black",
  color = heat_colors,
  breaks = breaksList,
  na_col = "grey90",
  main = "Module eigengenes correlations between S and M models"
)

# png("MEheatmap_S_M.png",
#     width = 7, height = 9,
#     units = "in",
#     res = 600)
# grid::grid.draw(ph$gtable)
# dev.off()

#### Step X-4. (Figure S8) Eigengenes correlations b/w S vs. S+M models ========
commonSamples <- intersect(rownames(MEs_SM), rownames(MEs_S))
MEs_SM  <- MEs_SM[commonSamples, , drop = FALSE]
MEs_S  <- MEs_S[commonSamples, , drop = FALSE]

# correlation between MEs
corMat <- cor(MEs_SM[,-1], MEs_S[,-1], use = "pairwise.complete.obs")
pMat   <- corPvalueStudent(corMat, nSamples = nrow(MEs_SM))

textMat <- matrix("", nrow = nrow(corMat), ncol = ncol(corMat),
                  dimnames = dimnames(corMat))
for(i in seq_len(nrow(corMat))) {
  for(j in seq_len(ncol(corMat))) {
    textMat[i,j] <- paste0("r=", sprintf("%.2f", corMat[i,j]),
                           "\np=", sprintf("%.2f", pMat[i,j]))
  }
}

# masing |r| < 0.8 or NA
mask <- abs(corMat) < 0.8
corMat_mask  <- corMat
corMat_mask[mask] <- NA

keepRows <- apply(!is.na(corMat_mask), 1, any)
keepCols <- apply(!is.na(corMat_mask), 2, any)

corMat_sub <- corMat_mask[keepRows, keepCols, drop=FALSE]
textMat_sub <- textMat[keepRows, keepCols, drop=FALSE]

# Set colors
breaksList <- seq(-1, 1, by = 0.2)
heat_colors <- colorRampPalette(c("blue","white","red"))(length(breaksList))

# Plotting heatmap
ph <- pheatmap(
  mat = corMat,
  display_numbers = textMat,   # r, p-value를 텍스트로 표시
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 10,
  fontsize_number = 8,
  number_color = "black",
  color = heat_colors,
  breaks = breaksList,
  na_col = "grey90",
  main = "Module eigengenes correlations between S and S+M models"
)

# png("MEheatmap_S_SM.png",
#     width = 7, height = 9,
#     units = "in",
#     res = 600)
# grid::grid.draw(ph$gtable)
# dev.off()

#### Step X-5. (Figure S9) Eigengenes correlations b/w M vs. S+M models ========
commonSamples <- intersect(rownames(MEs_SM), rownames(MEs_M))
MEs_SM  <- MEs_SM[commonSamples, , drop = FALSE]
MEs_M  <- MEs_M[commonSamples, , drop = FALSE]

# correlation between MEs
corMat <- cor(MEs_SM[,-1], MEs_M[,-1], use = "pairwise.complete.obs")
pMat   <- corPvalueStudent(corMat, nSamples = nrow(MEs_SM))

textMat <- matrix("", nrow = nrow(corMat), ncol = ncol(corMat),
                  dimnames = dimnames(corMat))
for(i in seq_len(nrow(corMat))) {
  for(j in seq_len(ncol(corMat))) {
    textMat[i,j] <- paste0("r=", sprintf("%.2f", corMat[i,j]),
                           "\np=", sprintf("%.2f", pMat[i,j]))
  }
}

# masing |r| < 0.8 or NA
mask <- abs(corMat) < 0.8
corMat_mask  <- corMat
corMat_mask[mask] <- NA

keepRows <- apply(!is.na(corMat_mask), 1, any)
keepCols <- apply(!is.na(corMat_mask), 2, any)

corMat_sub <- corMat_mask[keepRows, keepCols, drop=FALSE]
textMat_sub <- textMat[keepRows, keepCols, drop=FALSE]

# Set colors
breaksList <- seq(-1, 1, by = 0.2)
heat_colors <- colorRampPalette(c("blue","white","red"))(length(breaksList))

# Plotting heatmap
ph <- pheatmap(
  mat = corMat,
  display_numbers = textMat,   # r, p-value를 텍스트로 표시
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 10,
  fontsize_number = 8,
  number_color = "black",
  color = heat_colors,
  breaks = breaksList,
  na_col = "grey90",
  main = "Module eigengenes correlations between M and S+M models"
)

# png("MEheatmap_M_SM.png",
#     width = 10, height = 9,
#     units = "in",
#     res = 600)
# grid::grid.draw(ph$gtable)
# dev.off()

##### End of Step X ============================================================

##### Step XI. Network comparisons (II) ========================================
# comparison for different networks b/w Riaz et al (2025) and this study
# S' and M' are from Riaz study
# S, M and S+M are from this study
#### Step XI-1. Building new dataset with different networks ===================
### Step XI-1-a. Raw data ======================================================
## 1) Raw gene/gene cluster counts from Riaz ===================================
R <- readxl::read_xlsx("Riaz_raw.xlsx")
str(R)
dim(R) ## 62035 88

R2 <- R %>%
  select(-length, -chr, -type, -def, -Strain_MABNR56) %>%
  rename(Panid = gene_id, Species = species, Geneid = Strain_1021)

dim(R2) ## 62035 83

## Medicago
R_Medi <- R2 %>%
  filter(str_detect(Species, "Medi")) %>%
  select(-Species, -Geneid) ## 49,928

## Sinorhizobia
R_Sino <- R2 %>%
  filter(str_detect(Species, "Sino")) %>%
  select(-Species) ## 12,107

## 2) LogCPM of Riaz, Sino (geneID (SM_) + proteinID (WP)) =====================
A <- read.csv("GSE212235_Salmon_counts_all_samples.csv") 
dim(A) ## 56121 161

## 3) Read counts from Choi ====================================================
B <- read.csv("WGCNA_Readcount.csv")
dim(B) ## 26428 68

# Sinorhizobia
B_Sino <- B %>%
  select(-X, -gene_biotype, -name, -Description, -PMID, -GO.annotation, -Sense_Anti,
         -attributes.Sino_only., -Chr, -chr_Start, -chr_End, -Strand,
         -coding_length, -max, -average, -total) %>%
  filter(str_detect(Geneid, "CDO30")) %>%
  rename(Proteinid = EC._proteinID)

dim(B_Sino) ## 5194

# Medicago
B_Medi <- B %>%
  select(-X, -gene_biotype, -name, -Description, -PMID, -GO.annotation, -Sense_Anti,
         -attributes.Sino_only., -Chr, -chr_Start, -chr_End, -Strand, -EC._proteinID,
         -coding_length, -max, -average, -total) %>%
  filter(str_detect(Geneid, "MtrunA"))

dim(B_Medi) ## 21234

## 4) Medicago Module information of Riaz, including Grey module ===============
C <- read.csv("suppl_file_S8 (plant module genes).csv") 
dim(C) ## 30962 109

# Medicago with Grey
C_Medi <- C %>%
  select(gene, AssignedModule) %>%
  rename(Panid = gene, Module_kh = AssignedModule)

## Counting module 
C_Medi %>%
  group_by(Module_kh) %>%
  count() %>%
  print(n=100) ## 53 modules (35 to 7167, M0 (grey) of 6701)

## 5) Pan-genome Sinorhizobia Module information of Riaz, including Grey module =====================
D <- read.csv("suppl_file_S6 (pangenome-rhizobia module genes).csv") 
dim(D) ## 5939 66

# Sinorhizobia with Grey
D_Sino <- D %>%
  select(gene, AssignedModule) %>%
  rename(Panid = gene, Module_kh = AssignedModule)

# Counting modules
D_Sino %>%
  group_by(Module_kh) %>%
  count() %>%
  print(n=100) ## 32 modules (10 to 2145, M0 (grey) of 1051)

## 6) Sino + Medi Module information of Choi, including Grey module ============
E <- read.csv("Summary_table_SM_model.csv") 
dim(E) ## 26428 31

# Sinorhizobia with Grey
E_Sino <- E %>%
  select(GeneID, Module) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  rename(Geneid = GeneID, Module_pt = Module)

# Counting modules
E_Sino %>%
  group_by(Module_pt) %>%
  count() %>%
  print(n=100) ## 19 modules (3 to 1720, M0 (grey) of 798)

# Medicago with Grey
E_Medi <- E %>%
  select(GeneID, Module) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  rename(Geneid = GeneID, Module_pt = Module)

# Counting modules
E_Medi %>%
  group_by(Module_pt) %>%
  count() %>%
  print(n=100) ## 20 modules (21 to 7482, M0 (grey) of 1622)

## 7) Sinorhizobia reference in Choi ===========================================
G <- import("GCF_002197065.1_ASM219706v1_genomic.gtf.gz") 
G <- as.data.frame(G)
head(G)

G_pt_WP <- G %>%
  select(gene_id, protein_id) %>%
  mutate(row_id = row_number()) %>%
  filter(row_id %% 4 == 2) %>%
  select(-row_id) %>%
  rename(Geneid = gene_id, Proteinid = protein_id) %>%
  filter(!is.na(Proteinid)) ## GeneID plus proteinID for Choi

## 8) Sinorhizobia reference in Riaz ===========================================
H <- import("GCF_000006965.1_ASM696v1_genomic.gtf.gz") 
H <- as.data.frame(H)
head(H)

H_kh_WP <- H %>%
  select(gene_id, protein_id) %>%
  mutate(row_id = row_number()) %>%
  filter(row_id %% 4 == 2) %>%
  select(-row_id) %>%
  rename(Geneid = gene_id, Proteinid = protein_id) %>%
  filter(!is.na(Proteinid)) ## GeneID plus proteinID for Riaz

gtf_R_WP <- G_pt_WP %>%
  full_join(H_kh_WP, by = "Proteinid") %>%
  rename(Geneid_pt = Geneid.x, Geneid_kh = Geneid.y) %>%
  select(Proteinid, Geneid_pt, Geneid_kh) %>%
  ungroup()

dim(gtf_R_WP) ## 6094 3 ## Choi-made reference for proteinID comparison

## 9) Module information of Choi, including Grey module From S Model network =====
# Sinorhizobia with Grey
I <- read.csv("Summary_table_S_model.csv") 

I_Sino <- I %>%
  select(GeneID, Module) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  rename(Geneid = GeneID, Module_pt = Module)

# Counting modules
I_Sino %>%
  group_by(Module_pt) %>%
  count() %>%
  print(n=100) ## 8 modules (91 to 1998, grey of 1127)

## 10) Module information of Choi, including Grey module From M Model network =====
# Medicago with Grey
J <- read.csv("Summary_table_M_model.csv") 

# Medicago with Grey
J_Medi <- J %>%
  select(GeneID, Module) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  rename(Geneid = GeneID, Module_pt = Module)

# Counting modules
J_Medi %>%
  group_by(Module_pt) %>%
  count() %>%
  print(n=100) ## 17 modules (64 to 8037, M0 (grey) of 1846)

### Step XI-1-b. Modifyng data for Riaz dataset ================================
## Step XI-1-b-1. Sinorhizobia ================================================
# 1) D+R
dim(D_Sino) # 5939
dim(R_Sino) # 12107

D_Sino <- D_Sino %>% mutate(Panid = as.character(Panid))
R_Sino <- R_Sino %>% mutate(Panid = as.character(Panid))

DR_Sino <- D_Sino %>%
  left_join(R_Sino, by = c("Panid")) ## 5939

dim(DR_Sino)

# 2) DR+A
# 2-1). Sorting A id to make A_Sino
# Genes labelled with Geneid (like "RS_")
A_Sino1 <- A %>%
  filter(str_detect(X, "RS")) %>%
  rename(Geneid = X)
dim(A_Sino1) ## 285 161

# Check genes in A_Sino1 didn't have protein ID?
## 271 genes don't have protein ID, 14 genes 
A_Sino1p <- A_Sino1 %>%
  left_join(H_kh_WP, by = c("Geneid")) ## 285

# Genes labelled with proteinID (like "WP_")
A_Sino2 <- A %>%
  filter(str_detect(X, "WP")) %>%
  rename(Proteinid = X)
dim(A_Sino2) ## 5908 161

A_Sino2p <- A_Sino2 %>%
  left_join(H_kh_WP, by = c("Proteinid")) ## 5939

A_Sino <- bind_rows(A_Sino1p, A_Sino2p) %>%
  distinct()
dim(A_Sino) ## 6224 162

A_Sinop <- A_Sino %>%
  select(Geneid, Proteinid)
dim(A_Sinop) ## 6225 2

# 3) Merging DR+A
# Unmatched genes having multiple geneid and pseudo genes (= couldn't match proteinIDs)
DR_Sino %>%
  filter(str_detect(Geneid, ";") | str_detect(Geneid, "pseudo") | str_detect(Panid, "none")) ## 71 genes

DR_Sino2 <- DR_Sino %>%
  select(Panid, Module_kh, Geneid) %>%
  mutate(Matching = case_when(
    is.na(Panid) | Panid == "NA"              ~ "NoPanID",
    tolower(Panid) == "none"                  ~ "NoName",
    is.na(Geneid) | Geneid == "NA"            ~ "Undefined",
    str_detect(Geneid, "refound")                   ~ "Refound",
    str_detect(Geneid, ";")                   ~ "Duplicated",
    str_detect(tolower(Geneid), "pseudo")     ~ "Pseudo",
    str_detect(Panid, "^group") & str_detect(Geneid, "^SM_") ~ "Defined"
  )) 

dim(DR_Sino2) ## 5939 4

# Check the number of Geneid
DR_Sino2 %>%
  group_by(Module_kh, Matching)%>%
  count() %>%
  print(n=100)

DRA_Sino <- DR_Sino2 %>%
  drop_na(Geneid) %>%
  left_join(A_Sino, by = c("Geneid"))

summary(unique(DRA_Sino$Geneid)) ## 5465
summary(unique(DRA_Sino$Panid)) ## 5510
summary(unique(DRA_Sino$Proteinid)) ## 3575

DRA_Sino2 <- DRA_Sino %>%
  mutate(Matching2 = case_when(
    is.na(Proteinid) | Proteinid == "NA" ~ "NoProteinID",
    str_detect(Proteinid, "^WP_") ~ "ProteinID"))

DRA_Sino2 %>%
  select(Panid, Module_kh, Geneid, Proteinid, Matching, Matching2) %>%
  group_by(Module_kh, Matching, Matching2) %>%
  count() %>%
  filter(str_detect(Matching, "Defined")) %>%
  print(n=100)

DRA_Sino2 ## Riaz' Sinorhizobia dataset for network comparison


## Step XI-1-b-2. Medicago ================================================
# 1) C+A
dim(A) ## 56121 161
dim(C) ## 30962 109

A_Medi <- A %>%
  rename(Geneid = X)

C_Medi <- C %>%
  rename(Geneid = gene, Module_kh = AssignedModule) %>%
  select(Geneid, Module_kh)

CA_Medi <- C_Medi %>%
  left_join(A_Medi, by = c("Geneid")) ## 30962

CA_Medi %>%
  group_by(Module_kh) %>%
  count() %>%
  print(n=100)

CA_Medi ## Riaz's Medicago dataset for network comparison

### Step XI-1-c. Modifyng data for Choi dataset ================================
## Step XI-1-c-1. Sinorhizobia =================================================
# 1) B + E
dim(B_Sino) ## 5194 52
dim(E_Sino) ## 5194 2

BE_Sino <- E_Sino %>%
  left_join(B_Sino, by = c("Geneid"))

# Check the number of Geneid
BE_Sino2 <- BE_Sino %>%
  mutate(Matching3 = case_when(
    is.na(Proteinid) | Proteinid == "NA" ~ "NoProteinID",
    str_detect(Proteinid, "^WP_") ~ "ProteinID"))

## Step XI-1-c-2. Medicago =====================================================
# 1) B + E
dim(B_Medi) 
dim(E_Medi)

BE_Medi <- E_Medi %>%
  left_join(B_Medi, by = c("Geneid"))

### Step XI-1-d. Combining Riaz and Choi ========================================
## Step XI-1-d-1. Sinorhizobia ==================================================
# 1) DRA_Sino + BE_Sino
dim(DRA_Sino2) ## 5510 165
dim(BE_Sino) ## 5194 53

DRABE_Sino <- DRA_Sino2 %>%
  full_join(BE_Sino2, by = c("Proteinid")) %>%
  drop_na(Proteinid) %>% 
  rename(Geneid_kh = Geneid.x, Geneid_pt = Geneid.y)

dim(DRABE_Sino) ## 5441 219

# Labelling and Transforming
DRABE_Sino_module <- DRABE_Sino %>%
  select(Panid, Geneid_kh, Module_kh, Geneid_pt, Module_pt, Proteinid) %>%
  mutate(Module_kh = as.character(Module_kh))

## Step XI-1-d-2. Medicago ======================================================
# 1) CA_Medi + BE_Medi
dim(CA_Medi) ## 30962 162
dim(BE_Medi) ## 21231 52

CABE_Medi <- CA_Medi %>%
  full_join(BE_Medi, by = c("Geneid"))

dim(CABE_Medi) ## 31833 213

# Labelling and Transforming
CABE_Medi_module <- CABE_Medi %>%
  select(Geneid, Module_kh, Module_pt) %>%
  mutate(Module_kh = as.character(Module_kh))

#### Step XI-2. (Figure S10) S' vs. S+M models =================================
## With DRABE_Sino, DRABE_Sino_module
# Some are Not expressed in only multiple genotypes (S+M)
# Some are Not expressed in only single genotype (S')
DRABE_list_Sino <- DRABE_Sino_module %>%
  replace_na(list(Module_pt = "Not Expressed")) %>%
  mutate(Module_kh = as.character(Module_kh)) %>%
  mutate(Module_kh = if_else(Module_kh %in% c("NA", ""), NA_character_, Module_kh)) %>%
  replace_na(list(Module_kh = "Not Expressed"))

DRABE_list_Sino_NEC <- DRABE_list_Sino %>%
  filter(str_detect(Module_pt, c("Not Expressed")))
# %>%
#   count(is.na(Geneid_pt)) ## 353 These are not expressed in Choi

DRABE_list_Sino_NEC %>%
  group_by(Module_kh) %>%
  count() %>%
  print(n=100)

##
DRABE_list_Sino_NER <- DRABE_list_Sino %>%
  filter(str_detect(Module_kh, c("Not Expressed")))
# %>%
#   count(is.na(Geneid_kh)) ## 1855 These are not expressed in Riaz

DRABE_list_Sino_NER %>%
  group_by(Module_pt) %>%
  count() %>%
  print(n=100)

# Counting 
# Choi
DRABE_countlist_Sino_Choi <- DRABE_list_Sino %>%
  group_by(Module_pt) %>%
  group_modify(~ {
    df_sub <- .x  
    group_key <- .y
    
    if(group_key$Module_pt == "Not Expressed") {
      tibble(total_pt = nrow(df_sub))
    } else {
      tibble(total_pt = n_distinct(df_sub$Geneid_pt))
    }
  }) %>%
  ungroup() %>%
  mutate(Module_pt_lab = paste0(Module_pt, " (n=", total_pt, ")"))

# Riaz
DRABE_countlist_Sino_Riaz <- DRABE_list_Sino %>%
  group_by(Module_kh) %>%
  group_modify(~ {
    df_sub <- .x  
    group_key <- .y
    
    if(group_key$Module_kh == "Not Expressed") {
      tibble(total_kh = nrow(df_sub))
    } else {
      tibble(total_kh = n_distinct(df_sub$Geneid_kh))
    }
  }) %>%
  ungroup() %>%
  mutate(Module_kh_lab = paste0(Module_kh, " (n=", total_kh, ")"))

DRABE_Sino_labeled <- DRABE_list_Sino %>%
  left_join(DRABE_countlist_Sino_Choi %>% select(Module_pt, Module_pt_lab, total_pt), by = "Module_pt") %>%
  left_join(DRABE_countlist_Sino_Riaz %>% select(Module_kh, Module_kh_lab, total_kh), by = "Module_kh")

DRABE_Sino_Counted <- DRABE_Sino_labeled %>%
  count(Module_pt, Module_kh, Module_pt_lab, Module_kh_lab, name = "freq") %>%
  tidyr::complete(
    Module_pt_lab,
    Module_kh_lab,
    fill = list(freq = 0)
  )

DRABE_Sino_All <- DRABE_Sino_Counted %>%
  group_by(Module_kh_lab) %>% 
  mutate(prop = if (sum(freq) > 0) freq / sum(freq) else 0
  ) %>%
  ungroup()

## Ordering name tags;  
kh_levels <- unique(DRABE_Sino_All$Module_kh_lab)
kh_levels_sorted <- kh_levels[order(as.integer(str_extract(kh_levels, "\\d+")))]
DRABE_Sino_All$Module_kh_lab <- factor(DRABE_Sino_All$Module_kh_lab, levels = kh_levels_sorted)

# removal letters on each cell with proportion < 0.1
thr <- 0.1

keep_pt <- DRABE_Sino_All %>%
  dplyr::group_by(Module_pt_lab) %>%
  dplyr::summarise(ok = any(prop >= thr, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(ok) %>%
  dplyr::pull(Module_pt_lab)

## 
keep_kh <- DRABE_Sino_All %>%
  dplyr::group_by(Module_kh_lab) %>%
  dplyr::summarise(ok = any(prop >= thr, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(ok) %>%
  dplyr::pull(Module_kh_lab)

flow_filt_R <- DRABE_Sino_All %>%
  dplyr::filter(Module_pt_lab %in% keep_pt,
                Module_kh_lab %in% keep_kh)

## Ordering by Module gene size
# x axis: total_pt
pt_levels <- DRABE_Sino_labeled %>%
  dplyr::distinct(Module_pt_lab, total_pt) %>%
  dplyr::filter(Module_pt_lab %in% keep_pt) %>%
  dplyr::arrange(dplyr::desc(total_pt)) %>% # bigger on left
  dplyr::pull(Module_pt_lab)

#
kh_order_tbl <- DRABE_Sino_labeled %>%
  dplyr::distinct(Module_kh_lab, total_kh) %>%
  dplyr::filter(Module_kh_lab %in% keep_kh)

kh_unmatched <- kh_order_tbl %>%
  dplyr::filter(stringr::str_detect(Module_kh_lab, "^Not Expressed")) %>%
  dplyr::pull(Module_kh_lab) %>%
  unique()

kh_non_unmatched <- kh_order_tbl %>%
  dplyr::filter(!stringr::str_detect(Module_kh_lab, "^Not Expressed")) %>%
  dplyr::pull(Module_kh_lab) %>%
  unique()

pt_order_tbl <- DRABE_Sino_labeled %>%
  dplyr::distinct(Module_pt_lab, total_pt) %>%
  dplyr::filter(Module_pt_lab %in% keep_pt)

# Unmatched genes
pt_unmatched <- pt_order_tbl %>%
  dplyr::filter(stringr::str_detect(Module_pt_lab, "^Not Expressed")) %>%
  dplyr::pull(Module_pt_lab) %>%
  unique()

# Ordering by size
pt_non_unmatched <- pt_order_tbl %>%
  dplyr::filter(!stringr::str_detect(Module_pt_lab, "^Not Expressed")) %>%
  dplyr::arrange(dplyr::desc(total_pt)) %>% 
  dplyr::pull(Module_pt_lab) %>%
  unique()

# positioning unmatched genes
pt_levels <- c(pt_non_unmatched, pt_unmatched)
kh_levels <- c(kh_non_unmatched, kh_unmatched)

DRABE_Sino_All <- DRABE_Sino_All %>%
  dplyr::mutate(
    Module_pt_lab = factor(Module_pt_lab, levels = pt_levels),
    Module_kh_lab = factor(Module_kh_lab, levels = kh_levels)
  )

# Removal of NA module? of Choi data
DRABE_Sino_All <- DRABE_Sino_All %>%
  drop_na(Module_pt_lab)

## Plotting
ggplot(DRABE_Sino_All, aes(x = Module_pt_lab, y = Module_kh_lab, fill = prop)) +
  geom_tile(
    color = "grey30",
    size  = 0.5,         # ← increase border thickness here
    na.rm = FALSE
  ) +
  geom_text(
    data = flow_filt_R %>% dplyr::filter(prop >= 0.1),
    aes(label = sprintf("%.2f", prop)),
    size = 2.3
  ) +
  scale_fill_gradient(
    low      = "white",
    high     = "red",
    na.value = "white",
    name     = "Proportion"
  ) +
  labs(
    x     = "Modules (Multiple genotypes) (S+M Model)", 
    y     = "Modules (Single genotype)",
    title = expression(italic("S. meliloti"))
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y     = unit(0.5, "lines"),
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 7,  vjust = 10),
    axis.title.y = element_text(size = 7, vjust = -10),
    strip.placement     = "outside",
    plot.title = element_text(size = 8, vjust = -1),
    strip.background.y  = element_blank(),
    strip.text.y.left   = element_text(angle = 90),
    legend.title        = element_text(size = 6),
    legend.text         = element_text(size = 6),
    legend.key.size = unit(0.5, "lines"),
    legend.position = "right",
    legend.box.margin = margin(l = -10)
  )

# ggsave("Heatmap_Riza(S)_Choi(SM)_Sino.jpg", height = 10, width = 6, bg ="white", dpi = 800)

#### Step XI-3. (Figure S11) S' vs. S models ===================================
### Step XI-3-a. DRABI_Sino dataset ============================================
## 1) I+E
dim(I_Sino)
dim(B_Sino)

BI_Sino <- I_Sino %>%
  left_join(B_Sino, by = c("Geneid"))

# Check the number of Geneid
BI_Sino2 <- BI_Sino %>%
  mutate(Matching3 = case_when(
    is.na(Proteinid) | Proteinid == "NA" ~ "NoProteinID",
    str_detect(Proteinid, "^WP_") ~ "ProteinID"))

## 2) DRA_Sino + BE_Sino
dim(DRA_Sino2)
dim(BI_Sino2) 

DRABI_Sino <- DRA_Sino2 %>%
  full_join(BI_Sino2, by = c("Proteinid")) %>%
  drop_na(Proteinid) %>% 
  rename(Geneid_kh = Geneid.x, Geneid_pt = Geneid.y)

dim(DRABI_Sino)

# Labelling and Transforming
DRABI_Sino_module <- DRABI_Sino %>%
  select(Panid, Geneid_kh, Module_kh, Geneid_pt, Module_pt, Proteinid) %>%
  mutate(Module_kh = as.character(Module_kh))

## With DRABI_Sino and DRABI_Sino_module
DRABI_list_Sino <- DRABI_Sino_module %>%
  replace_na(list(Module_pt = "Not Expressed")) %>%
  mutate(Module_kh = as.character(Module_kh)) %>%
  mutate(Module_kh = if_else(Module_kh %in% c("NA", ""), NA_character_, Module_kh)) %>%
  replace_na(list(Module_kh = "Not Expressed"))

DRABI_list_Sino_NEC <- DRABI_list_Sino %>%
  filter(str_detect(Module_pt, c("Not Expressed")))
# %>%
#   count(is.na(Geneid_pt)) ## 353 These are not expressed in Choi

DRABI_list_Sino_NEC %>%
  group_by(Module_kh) %>%
  count() %>%
  print(n=100)

DRABI_list_Sino_NER <- DRABI_list_Sino %>%
  filter(str_detect(Module_kh, c("Not Expressed")))
# %>%
#   count(is.na(Geneid_kh)) ## 1855 These are not expressed in Riaz

DRABI_list_Sino_NER %>%
  group_by(Module_pt) %>%
  count() %>%
  print(n=100)

# Counting 
# Choi
DRABI_countlist_Sino_Choi <- DRABI_list_Sino %>%
  group_by(Module_pt) %>%
  group_modify(~ {
    df_sub <- .x  
    group_key <- .y
    
    if(group_key$Module_pt == "Not Expressed") {
      tibble(total_pt = nrow(df_sub))
    } else {
      tibble(total_pt = n_distinct(df_sub$Geneid_pt))
    }
  }) %>%
  ungroup() %>%
  mutate(Module_pt_lab = paste0(Module_pt, " (n=", total_pt, ")"))

# Riaz
DRABI_countlist_Sino_Riaz <- DRABI_list_Sino %>%
  group_by(Module_kh) %>%
  group_modify(~ {
    df_sub <- .x 
    group_key <- .y 
    
    if(group_key$Module_kh == "Not Expressed") {
      tibble(total_kh = nrow(df_sub))
    } else {
      tibble(total_kh = n_distinct(df_sub$Geneid_kh))
    }
  }) %>%
  ungroup() %>%
  mutate(Module_kh_lab = paste0(Module_kh, " (n=", total_kh, ")"))

DRABI_Sino_labeled <- DRABI_list_Sino %>%
  left_join(DRABI_countlist_Sino_Choi %>% select(Module_pt, Module_pt_lab, total_pt), by = "Module_pt") %>%
  left_join(DRABI_countlist_Sino_Riaz %>% select(Module_kh, Module_kh_lab, total_kh), by = "Module_kh")

DRABI_Sino_Counted <- DRABI_Sino_labeled %>%
  count(Module_pt, Module_kh, Module_pt_lab, Module_kh_lab, name = "freq") %>%
  tidyr::complete(
    Module_pt_lab,
    Module_kh_lab,
    fill = list(freq = 0)
  )

DRABI_Sino_All <- DRABI_Sino_Counted %>%
  group_by(Module_kh_lab) %>% 
  mutate(prop = if (sum(freq) > 0) freq / sum(freq) else 0
  ) %>%
  ungroup()

## Ordering name tags;  
kh_levels <- unique(DRABI_Sino_All$Module_kh_lab)
kh_levels_sorted <- kh_levels[order(as.integer(str_extract(kh_levels, "\\d+")))]
DRABI_Sino_All$Module_kh_lab <- factor(DRABI_Sino_All$Module_kh_lab, levels = kh_levels_sorted)

#
thr <- 0.1

keep_pt <- DRABI_Sino_All %>%
  dplyr::group_by(Module_pt_lab) %>%
  dplyr::summarise(ok = any(prop >= thr, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(ok) %>%
  dplyr::pull(Module_pt_lab)

keep_kh <- DRABI_Sino_All %>%
  dplyr::group_by(Module_kh_lab) %>%
  dplyr::summarise(ok = any(prop >= thr, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(ok) %>%
  dplyr::pull(Module_kh_lab)

flow_filt_R <- DRABI_Sino_All %>%
  dplyr::filter(Module_pt_lab %in% keep_pt,
                Module_kh_lab %in% keep_kh)

pt_levels <- DRABI_Sino_labeled %>%
  dplyr::distinct(Module_pt_lab, total_pt) %>%
  dplyr::filter(Module_pt_lab %in% keep_pt) %>%
  dplyr::arrange(dplyr::desc(total_pt)) %>%    
  dplyr::pull(Module_pt_lab)

#
kh_order_tbl <- DRABI_Sino_labeled %>%
  dplyr::distinct(Module_kh_lab, total_kh) %>%
  dplyr::filter(Module_kh_lab %in% keep_kh)

kh_unmatched <- kh_order_tbl %>%
  dplyr::filter(stringr::str_detect(Module_kh_lab, "^Not Expressed")) %>%
  dplyr::pull(Module_kh_lab) %>%
  unique()

kh_non_unmatched <- kh_order_tbl %>%
  dplyr::filter(!stringr::str_detect(Module_kh_lab, "^Not Expressed")) %>%
  dplyr::pull(Module_kh_lab) %>%
  unique()

pt_order_tbl <- DRABI_Sino_labeled %>%
  dplyr::distinct(Module_pt_lab, total_pt) %>%
  dplyr::filter(Module_pt_lab %in% keep_pt)

pt_unmatched <- pt_order_tbl %>%
  dplyr::filter(stringr::str_detect(Module_pt_lab, "^Not Expressed")) %>%
  dplyr::pull(Module_pt_lab) %>%
  unique()

pt_non_unmatched <- pt_order_tbl %>%
  dplyr::filter(!stringr::str_detect(Module_pt_lab, "^Not Expressed")) %>%
  dplyr::arrange(dplyr::desc(total_pt)) %>% 
  dplyr::pull(Module_pt_lab) %>%
  unique()

pt_levels <- c(pt_non_unmatched, pt_unmatched)
kh_levels <- c(kh_non_unmatched, kh_unmatched)

DRABI_Sino_All <- DRABI_Sino_All %>%
  dplyr::mutate(
    Module_pt_lab = factor(Module_pt_lab, levels = pt_levels),
    Module_kh_lab = factor(Module_kh_lab, levels = kh_levels)
  )

# Removal of NA module? of Choi data
DRABI_Sino_All <- DRABI_Sino_All %>%
  drop_na(Module_pt_lab)

##
ggplot(DRABI_Sino_All, aes(x = Module_pt_lab, y = Module_kh_lab, fill = prop)) +
  geom_tile(
    color = "grey30",
    size  = 0.5,       
    na.rm = FALSE
  ) +
  geom_text(
    data = flow_filt_R %>% dplyr::filter(prop >= 0.1),
    aes(label = sprintf("%.2f", prop)),
    size = 2.3
  ) +
  scale_fill_gradient(
    low      = "white",
    high     = "red",
    na.value = "white",
    name     = "Proportion"
  ) +
  labs(
    x     = "Modules (Multiple genotypes) (S Model)", 
    y     = "Modules (Single genotype)",
    title = expression(italic("S. meliloti"))
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y     = unit(0.5, "lines"),
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 7,  vjust = 10),
    axis.title.y = element_text(size = 7, vjust = -10),
    strip.placement     = "outside",
    plot.title = element_text(size = 8, vjust = -1),
    strip.background.y  = element_blank(),
    strip.text.y.left   = element_text(angle = 90),
    legend.title        = element_text(size = 6),
    legend.text         = element_text(size = 6),
    legend.key.size = unit(0.5, "lines"),
    legend.position = "right",
    legend.box.margin = margin(l = -10)
  )

# ggsave("Heatmap_Riza(S)_Choi(S)_Sino.jpg", height = 10, width = 5, bg ="white", dpi = 800)

#### Step XI-4. (Figure S12) M' vs. S+M models =================================
## With CABE_Medi, CABE_Medi_module
# Some are Not expressed in only multiple genotypes (S+M, Choi)
# Some are Not expressed in only single genotype (M', Riaz)
CABE_list_Medi <- CABE_Medi_module %>%
  replace_na(list(Module_pt = "Not Expressed")) %>%
  mutate(Module_kh = as.character(Module_kh)) %>%
  mutate(Module_kh = if_else(Module_kh %in% c("NA", ""), NA_character_, Module_kh)) %>%
  replace_na(list(Module_kh = "Not Expressed"))

CABE_list_Medi_NEC <- CABE_list_Medi %>%
  filter(str_detect(Module_pt, c("Not Expressed")))
# %>%
#   count(is.na(Geneid)) ## 10602 These are not expressed in Choi

CABE_list_Medi_NEC %>%
  group_by(Module_kh) %>%
  count() %>%
  print(n=100)

# List of genes not expressed in Riaz, but expressed in Choi
CABE_list_Medi_NER <- CABE_list_Medi %>%
  filter(str_detect(Module_kh, c("Not Expressed"))) ## 871 genes not expressed in Riaz

##
CABE_countlist_Medi_Choi <- CABE_list_Medi %>%
  group_by(Module_pt) %>%
  group_modify(~ {
    df_sub <- .x  
    group_key <- .y
    
    if(group_key$Module_pt == "Not Expressed") {
      tibble(total_pt = nrow(df_sub))
    } else {
      tibble(total_pt = n_distinct(df_sub$Geneid))
    }
  }) %>%
  ungroup() %>%
  mutate(Module_pt_lab = paste0(Module_pt, " (n=", total_pt, ")"))

# Riaz
CABE_countlist_Medi_Riaz <- CABE_list_Medi %>%
  group_by(Module_kh) %>%
  group_modify(~ {
    df_sub <- .x  
    group_key <- .y
    
    if(group_key$Module_kh == "Not Expressed") {
      tibble(total_kh = nrow(df_sub))
    } else {
      tibble(total_kh = n_distinct(df_sub$Geneid))
    }
  }) %>%
  ungroup() %>%
  mutate(Module_kh_lab = paste0(Module_kh, " (n=", total_kh, ")"))

CABE_Medi_labeled <- CABE_list_Medi %>%
  left_join(CABE_countlist_Medi_Choi %>% select(Module_pt, Module_pt_lab, total_pt), by = "Module_pt") %>%
  left_join(CABE_countlist_Medi_Riaz %>% select(Module_kh, Module_kh_lab, total_kh), by = "Module_kh")

CABE_Medi_Counted <- CABE_Medi_labeled %>%
  count(Module_pt, Module_kh, Module_pt_lab, Module_kh_lab, name = "freq") %>%
  tidyr::complete(
    Module_pt_lab,
    Module_kh_lab,
    fill = list(freq = 0)
  )

CABE_Medi_All <- CABE_Medi_Counted %>%
  group_by(Module_kh_lab) %>% 
  mutate(prop = if (sum(freq) > 0) freq / sum(freq) else 0
  ) %>%
  ungroup()

## Ordering name tags;  
kh_levels <- unique(CABE_Medi_All$Module_kh_lab)
kh_levels_sorted <- kh_levels[order(as.integer(str_extract(kh_levels, "\\d+")))]
CABE_Medi_All$Module_kh_lab <- factor(CABE_Medi_All$Module_kh_lab, levels = kh_levels_sorted)

# removal letters on each cell with proportion < 0.1
thr <- 0.1

keep_pt <- CABE_Medi_All %>%
  dplyr::group_by(Module_pt_lab) %>%
  dplyr::summarise(ok = any(prop >= thr, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(ok) %>%
  dplyr::pull(Module_pt_lab)

keep_kh <- CABE_Medi_All %>%
  dplyr::group_by(Module_kh_lab) %>%
  dplyr::summarise(ok = any(prop >= thr, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(ok) %>%
  dplyr::pull(Module_kh_lab)

flow_filt_R <- CABE_Medi_All %>%
  dplyr::filter(Module_pt_lab %in% keep_pt,
                Module_kh_lab %in% keep_kh)

## Ordering by Module gene size
# x axis: total_pt
pt_levels <- CABE_Medi_labeled %>%
  dplyr::distinct(Module_pt_lab, total_pt) %>%
  dplyr::filter(Module_pt_lab %in% keep_pt) %>%
  dplyr::arrange(dplyr::desc(total_pt)) %>%   
  dplyr::pull(Module_pt_lab)

#
kh_order_tbl <- CABE_Medi_labeled %>%
  dplyr::distinct(Module_kh_lab, total_kh) %>%
  dplyr::filter(Module_kh_lab %in% keep_kh)

kh_unmatched <- kh_order_tbl %>%
  dplyr::filter(stringr::str_detect(Module_kh_lab, "^Not Expressed")) %>%
  dplyr::pull(Module_kh_lab) %>%
  unique()

kh_non_unmatched <- kh_order_tbl %>%
  dplyr::filter(!stringr::str_detect(Module_kh_lab, "^Not Expressed")) %>%
  dplyr::pull(Module_kh_lab) %>%
  unique()

pt_order_tbl <- CABE_Medi_labeled %>%
  dplyr::distinct(Module_pt_lab, total_pt) %>%
  dplyr::filter(Module_pt_lab %in% keep_pt)

pt_unmatched <- pt_order_tbl %>%
  dplyr::filter(stringr::str_detect(Module_pt_lab, "^Not Expressed")) %>%
  dplyr::pull(Module_pt_lab) %>%
  unique()

pt_non_unmatched <- pt_order_tbl %>%
  dplyr::filter(!stringr::str_detect(Module_pt_lab, "^Not Expressed")) %>%
  dplyr::arrange(dplyr::desc(total_pt)) %>% 
  dplyr::pull(Module_pt_lab) %>%
  unique()

pt_levels <- c(pt_non_unmatched, pt_unmatched)
kh_levels <- c(kh_non_unmatched, kh_unmatched)

CABE_Medi_All <- CABE_Medi_All %>%
  dplyr::mutate(
    Module_pt_lab = factor(Module_pt_lab, levels = pt_levels),
    Module_kh_lab = factor(Module_kh_lab, levels = kh_levels)
  )

# Removal of NA module? of Choi data
CABE_Medi_All <- CABE_Medi_All %>%
  drop_na(Module_pt_lab)

##
ggplot(CABE_Medi_All, aes(x = Module_pt_lab, y = Module_kh_lab, fill = prop)) +
  geom_tile(
    color = "grey30",
    size  = 0.5,         # ← increase border thickness here
    na.rm = FALSE
  ) +
  geom_text(
    data = flow_filt_R %>% dplyr::filter(prop >= 0.1),
    aes(label = sprintf("%.2f", prop)),
    size = 2.3
  ) +
  scale_fill_gradient(
    low      = "white",
    high     = "red",
    na.value = "white",
    name     = "Proportion"
  ) +
  labs(
    x     = "Modules (Multiple genotypes) (S+M Model)", 
    y     = "Modules (Single genotype)",
    title = expression(italic("M. truncatula"))
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y     = unit(0.5, "lines"),
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 7,  vjust = 10),
    axis.title.y = element_text(size = 7, vjust = -10),
    strip.placement     = "outside",
    plot.title = element_text(size = 8, vjust = -1),
    strip.background.y  = element_blank(),
    strip.text.y.left   = element_text(angle = 90),
    legend.title        = element_text(size = 6),
    legend.text         = element_text(size = 6),
    legend.key.size = unit(0.5, "lines"),
    legend.position = "right",
    legend.box.margin = margin(l = -10)
  )

# ggsave("Heatmap_Riza(M)_Choi(SM)_Medi.jpg", height = 10, width = 6, bg ="white", dpi = 800)

#### Step XI-5. (Figure S13) M' vs. M models ===================================
### Step XI-5-a. CABJ_Medi dataset =============================================
## 1) J+E
dim(J_Medi)
dim(B_Medi)

BJ_Medi <- J_Medi %>%
  left_join(B_Medi, by = c("Geneid"))

## 2) CA_Medi + BJ_Medi
dim(CA_Medi)
dim(BJ_Medi)

CABJ_Medi <- CA_Medi %>%
  full_join(BJ_Medi, by = c("Geneid"))

# Labelling and Transforming
CABJ_Medi_module <- CABJ_Medi %>%
  select(Geneid, Module_kh, Module_pt) %>%
  mutate(Module_kh = as.character(Module_kh))

CABJ_list_Medi <- CABJ_Medi_module %>%
  replace_na(list(Module_pt = "Not Expressed")) %>%
  mutate(Module_kh = as.character(Module_kh)) %>%
  mutate(Module_kh = if_else(Module_kh %in% c("NA", ""), NA_character_, Module_kh)) %>%
  replace_na(list(Module_kh = "Not Expressed"))

CABJ_list_Medi_NEC <- CABJ_list_Medi %>%
  filter(str_detect(Module_pt, c("Not Expressed")))
# %>%
#   count(is.na(Geneid)) ## 10602 These are not expressed in Choi

##
CABJ_countlist_Medi_Choi <- CABJ_list_Medi %>%
  group_by(Module_pt) %>%
  group_modify(~ {
    df_sub <- .x 
    group_key <- .y 
    
    if(group_key$Module_pt == "Not Expressed") {
      tibble(total_pt = nrow(df_sub))
    } else {
      tibble(total_pt = n_distinct(df_sub$Geneid))
    }
  }) %>%
  ungroup() %>%
  mutate(Module_pt_lab = paste0(Module_pt, " (n=", total_pt, ")"))

# Riaz
CABJ_countlist_Medi_Riaz <- CABJ_list_Medi %>%
  group_by(Module_kh) %>%
  group_modify(~ {
    df_sub <- .x  
    group_key <- .y
    
    if(group_key$Module_kh == "Not Expressed") {
      tibble(total_kh = nrow(df_sub))
    } else {
      tibble(total_kh = n_distinct(df_sub$Geneid))
    }
  }) %>%
  ungroup() %>%
  mutate(Module_kh_lab = paste0(Module_kh, " (n=", total_kh, ")"))

CABJ_Medi_labeled <- CABJ_list_Medi %>%
  left_join(CABJ_countlist_Medi_Choi %>% select(Module_pt, Module_pt_lab, total_pt), by = "Module_pt") %>%
  left_join(CABJ_countlist_Medi_Riaz %>% select(Module_kh, Module_kh_lab, total_kh), by = "Module_kh")

CABJ_Medi_Counted <- CABJ_Medi_labeled %>%
  count(Module_pt, Module_kh, Module_pt_lab, Module_kh_lab, name = "freq") %>%
  tidyr::complete(
    Module_pt_lab,
    Module_kh_lab,
    fill = list(freq = 0)
  )

CABJ_Medi_All <- CABJ_Medi_Counted %>%
  group_by(Module_kh_lab) %>% 
  mutate(prop = if (sum(freq) > 0) freq / sum(freq) else 0
  ) %>%
  ungroup()

## Ordering name tags;  
kh_levels <- unique(CABJ_Medi_All$Module_kh_lab)
kh_levels_sorted <- kh_levels[order(as.integer(str_extract(kh_levels, "\\d+")))]
CABJ_Medi_All$Module_kh_lab <- factor(CABJ_Medi_All$Module_kh_lab, levels = kh_levels_sorted)

#
thr <- 0.1

keep_pt <- CABJ_Medi_All %>%
  dplyr::group_by(Module_pt_lab) %>%
  dplyr::summarise(ok = any(prop >= thr, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(ok) %>%
  dplyr::pull(Module_pt_lab)

keep_kh <- CABJ_Medi_All %>%
  dplyr::group_by(Module_kh_lab) %>%
  dplyr::summarise(ok = any(prop >= thr, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(ok) %>%
  dplyr::pull(Module_kh_lab)

flow_filt_R <- CABJ_Medi_All %>%
  dplyr::filter(Module_pt_lab %in% keep_pt,
                Module_kh_lab %in% keep_kh)

pt_levels <- CABJ_Medi_labeled %>%
  dplyr::distinct(Module_pt_lab, total_pt) %>%
  dplyr::filter(Module_pt_lab %in% keep_pt) %>%
  dplyr::arrange(dplyr::desc(total_pt)) %>%   
  dplyr::pull(Module_pt_lab)

#
kh_order_tbl <- CABJ_Medi_labeled %>%
  dplyr::distinct(Module_kh_lab, total_kh) %>%
  dplyr::filter(Module_kh_lab %in% keep_kh)

kh_unmatched <- kh_order_tbl %>%
  dplyr::filter(stringr::str_detect(Module_kh_lab, "^Not Expressed")) %>%
  dplyr::pull(Module_kh_lab) %>%
  unique()

kh_non_unmatched <- kh_order_tbl %>%
  dplyr::filter(!stringr::str_detect(Module_kh_lab, "^Not Expressed")) %>%
  dplyr::pull(Module_kh_lab) %>%
  unique()

pt_order_tbl <- CABJ_Medi_labeled %>%
  dplyr::distinct(Module_pt_lab, total_pt) %>%
  dplyr::filter(Module_pt_lab %in% keep_pt)

pt_unmatched <- pt_order_tbl %>%
  dplyr::filter(stringr::str_detect(Module_pt_lab, "^Not Expressed")) %>%
  dplyr::pull(Module_pt_lab) %>%
  unique()

pt_non_unmatched <- pt_order_tbl %>%
  dplyr::filter(!stringr::str_detect(Module_pt_lab, "^Not Expressed")) %>%
  dplyr::arrange(dplyr::desc(total_pt)) %>% 
  dplyr::pull(Module_pt_lab) %>%
  unique()

pt_levels <- c(pt_non_unmatched, pt_unmatched)
kh_levels <- c(kh_non_unmatched, kh_unmatched)

CABJ_Medi_All <- CABJ_Medi_All %>%
  dplyr::mutate(
    Module_pt_lab = factor(Module_pt_lab, levels = pt_levels),
    Module_kh_lab = factor(Module_kh_lab, levels = kh_levels)
  )

# Removal of NA module? of Choi data
CABJ_Medi_All <- CABJ_Medi_All %>%
  drop_na(Module_pt_lab)

##
ggplot(CABJ_Medi_All, aes(x = Module_pt_lab, y = Module_kh_lab, fill = prop)) +
  geom_tile(
    color = "grey30",
    size  = 0.5,   
    na.rm = FALSE
  ) +
  geom_text(
    data = flow_filt_R %>% dplyr::filter(prop >= 0.1),
    aes(label = sprintf("%.2f", prop)),
    size = 2.3
  ) +
  scale_fill_gradient(
    low      = "white",
    high     = "red",
    na.value = "white",
    name     = "Proportion"
  ) +
  labs(
    x     = "Modules (Multiple genotypes (M Model))", 
    y     = "Modules (Single genotype)",
    title = expression(italic("M. truncatula"))
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y     = unit(0.5, "lines"),
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 7,  vjust = 10),
    axis.title.y = element_text(size = 7, vjust = -10),
    strip.placement     = "outside",
    plot.title = element_text(size = 8, vjust = -1),
    strip.background.y  = element_blank(),
    strip.text.y.left   = element_text(angle = 90),
    legend.title        = element_text(size = 6),
    legend.text         = element_text(size = 6),
    legend.key.size = unit(0.5, "lines"),
    legend.position = "right",
    legend.box.margin = margin(l = -10)
  )

# ggsave("Heatmap_Riza(M)_Choi(M)_Medi.jpg", height = 10, width = 6, bg ="white", dpi = 800)
##### End of Step XI ===========================================================