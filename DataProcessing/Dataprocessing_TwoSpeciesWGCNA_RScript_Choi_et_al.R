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

##### Pre-Steps. Importing dataset =============================================
### Kono_reads_set2_march20.xlsx
# This file is newly created on March 20, 2025 (by Peter Tiffin)
# Most genes were also removed based on their annotations;
# such as tRNA (974), rRNA (62), pre-miRNA (799), and repeat-region (245,073) for Medicago;
# RNA (12), 6s rRNA (2), tRNA (65) for Sinorhizobia

### Filtering was performed based on the mean Reads at least 10 across 54 libraries
### Gene filtered based on their read counts (the average number of mapped reads < 10)
# Deletion of genes; ncRNA (4,064), mRNA (42,040), and Sino (3,287) for Antisense
# Deletion of genes; ncRNA (4,839), mRNA (23,143), and Sino (1,029) for Sense
### Some of the remaining genes were nested in the other gene models

##### Step I. Data Information  ================================================
#### Step I-a. Importing dataset ==============================================
Data_raw <- read.csv("WGCNA_Readcount.csv")
Data_library_50 <- read.csv("Data_library_50.csv")

#### Step I-b. (Dataset S1) Read counts for each library =======================
Data_raw_R <- Data_raw %>%
  filter(str_detect(Geneid, "CDO30_"))
dim(Data_raw_R) ## 5194 68
Data_raw_M <- Data_raw %>%
  filter(str_detect(Geneid, "MtrunA"))
dim(Data_raw_M) ## 21234 68

### Step I-b-1. For Sinorhizobia  ==============================================
ColSum_Sn_R <- colSums(Data_raw_R[,-c(1:17)], na.rm = TRUE)
ColSum_Sn_R <- data.frame(ColSum_Sn_R)
ColSum_Sn_R <- rownames_to_column(ColSum_Sn_R, var = "library")


### Step I-b-2. For Medicago  ==================================================
ColSum_Sn_M <- colSums(Data_raw_M[,-c(1:17)], na.rm = TRUE)
ColSum_Sn_M <- data.frame(ColSum_Sn_M)
ColSum_Sn_M <- rownames_to_column(ColSum_Sn_M, var = "library")

## Merging columns from Sinorhizobia and Medicago
ColSum_Sn <- cbind(ColSum_Sn_R, ColSum_Sn_M)
ColSum_Sn <- ColSum_Sn[, !duplicated(colnames(ColSum_Sn))]

ColSum_Sn$Total_Sn <- rowSums(ColSum_Sn[2:3],)

ColSum_Sn$ColSum_Sn_R <- paste0(
  ColSum_Sn$ColSum_Sn_R,
  " (", formatC((ColSum_Sn$ColSum_Sn_R / ColSum_Sn$Total_Sn) * 100, format = "f", digits = 1), "%)"
)

ColSum_Sn$ColSum_Sn_M <- paste0(
  ColSum_Sn$ColSum_Sn_M,
  " (", formatC((ColSum_Sn$ColSum_Sn_M / ColSum_Sn$Total_Sn) * 100, format = "f", digits = 1), "%)"
)

## Making tables
ColSum_All <- ColSum_Sn %>%
  left_join(Data_library_50, by = "library") %>%
  data.frame()

### End of Step I  =============================================================

##### Step II. TPM calculation =================================================
Data_raw_R ## 5,194
Data_raw_M ## 21,234

#### Step II-a. Sinorhizobia ====================================================
## converting into long-form
Data_raw_R_long <- Data_raw_R %>%
  select(-X, -name, -Description, -PMID, -gene_biotype, -Sense_Anti, 
         -`GO.annotation`, -`EC._proteinID`, -`attributes.Sino_only.`,
         -chr_Start, -chr_End, -Strand, -max, -average, -total) %>%
  pivot_longer(cols = SRR18707692:SRR18707778,
               names_to = ("library"), values_to = "Reads") 

# calculating TPM, log2(TPM+1), and TPM/coding length
TPMExp_R <- Data_raw_R_long %>%
  mutate(RPK = Reads / (coding_length)) %>%
  group_by(library) %>%
  mutate(
    TPM = (RPK / sum(RPK, na.rm =  TRUE)) * 1e6,
    Log2TPM = log2(TPM+1),
    TPM_cdlg = (TPM / (coding_length)), 
  ) %>%
  ungroup() 

mean_TPMExp_R <- TPMExp_R %>%
  group_by(Geneid) %>%
  summarise(
    meanTPM = mean(TPM, na.rm = TRUE),
    meanLog2TPM = mean(Log2TPM, na.rm = TRUE),
    meanTPM_cdlg = mean(TPM_cdlg, na.rm = TRUE),
    .groups = "drop"
  )

## converting into wide-form
TPM_R <- TPMExp_R %>%
  select(Geneid, library, Chr, coding_length, TPM) %>%
  pivot_wider(names_from = library, values_from = TPM)

Log2TPM_R <- TPMExp_R %>%
  select(Geneid, library, Chr, coding_length, Log2TPM) %>%
  pivot_wider(names_from = library, values_from = Log2TPM)

TPM_cdlg_R <- TPMExp_R %>%
  select(Geneid, library, Chr, coding_length, TPM_cdlg) %>%
  pivot_wider(names_from = library, values_from = TPM_cdlg)

#### Step II-b. Medicago =======================================================
Data_raw_M_long <- Data_raw_M %>%
  select(-X, -name, -Description, -PMID, -gene_biotype, -Sense_Anti, 
         -`GO.annotation`, -`EC._proteinID`, -`attributes.Sino_only.`,
         -chr_Start, -chr_End, -Strand, -max, -average, -total) %>%
  pivot_longer(cols = SRR18707692:SRR18707778,
               names_to = ("library"), values_to = "Reads") 

# calculating TPM, log2(TPM+1), and TPM/coding length
TPMExp_M <- Data_raw_M_long %>%
  mutate(RPK = Reads / (coding_length)) %>%
  group_by(library) %>%
  mutate(
    TPM = (RPK / sum(RPK, na.rm =  TRUE)) * 1e6,
    Log2TPM = log2(TPM+1),
    TPM_cdlg = (TPM / (coding_length)), 
  ) %>%
  ungroup() 

mean_TPMExp_M <- TPMExp_M %>%
  group_by(Geneid) %>%
  summarise(
    meanTPM = mean(TPM, na.rm = TRUE),
    meanLog2TPM = mean(Log2TPM, na.rm = TRUE),
    meanTPM_cdlg = mean(TPM_cdlg, na.rm = TRUE),
    .groups = "drop"
  )

## converting into wide-form
TPM_M <- TPMExp_M %>%
  select(Geneid, library, Chr, coding_length, TPM) %>%
  pivot_wider(names_from = library, values_from = TPM)

Log2TPM_M <- TPMExp_M %>%
  select(Geneid, library, Chr, coding_length, Log2TPM) %>%
  pivot_wider(names_from = library, values_from = Log2TPM)

TPM_cdlg_M <- TPMExp_M %>%
  select(Geneid, library, Chr, coding_length, TPM_cdlg) %>%
  pivot_wider(names_from = library, values_from = TPM_cdlg)

#### Step II-c. Exporting TPM, Log2(TPM+1), and TPM/coding length profiles ====
### Read counts 
Data_raw_R
Data_raw_M

Data_readcounts <- rbind(Data_raw_R, Data_raw_R) ##

# write.csv(Data_readcounts, "Data_readcounts.csv")

TPM_R 
TPM_M 

Data_TPM <- rbind(TPM_R, TPM_M) ##

# write.csv(Data_TPM, "Data_TPM.csv")

Log2TPM_R #
Log2TPM_M #

Data_Log2TPM <- rbind(Log2TPM_R, Log2TPM_M)

# write.csv(Data_Log2TPM, "Data_Log2TPM.csv")

TPM_cdlg_R #
TPM_cdlg_M #

Data_TPM_cdlg <- rbind(TPM_cdlg_R, TPM_cdlg_M)

# write.csv(Data_TPM_cdlg, "Data_TPM_cdlg.csv")

### End of Step II  ============================================================

##### Step III. Calculation of mean and CV of TPM, Log2(TPM+1), and TPM/coding length =====
##### For Suppl. tables S4, S5, S6, S7, S8 =====================================
#### Step III-a. Importing datasets ============================================
Data_TPM_cdlg <- read.csv("Data_TPM_cdlg.csv")
Data_readcounts <- read.csv("Data_readcounts.csv")
Data_TPM <- read.csv("Data_TPM.csv")
Data_Log2TPM <- read.csv("Data_Log2TPM.csv")

#### Step III-b. TPM/coding length =============================================
### Calculating mean, sd, and CV
Stat.TPM_cdlg <- Data_TPM_cdlg %>%
  pivot_longer(cols = SRR18707692:SRR18707778,
               names_to = ("library"), values_to = "TPM_cdlg") %>%
  group_by(Geneid) %>%
  summarise(
    mean = mean(TPM_cdlg, na.rm = TRUE),
    sd = sd(TPM_cdlg, na.rm = TRUE),
    CV = (sd/mean) * 100,
    .groups = "drop"
  )

## Select annotation information
anno_raw <- Data_raw %>%
  select(Geneid, Description, Chr, gene_biotype, EC._proteinID) %>%
  rename(Proteinid = EC._proteinID)

## Extract GeneID and add annotations
Sum.TPM_cdlg <- Stat.TPM_cdlg %>%
  left_join(anno_raw, by = c("Geneid"))

### Step III-b-1. Selecting Top 20 highly expressed genes ======================
## For Sinorhizobia ============================================================
Top20_TPM_cdlg_R <- Sum.TPM_cdlg %>%
  filter(str_detect(Geneid, "CDO30")) %>%
  select(Geneid, Proteinid, Chr, gene_biotype, Description, mean) %>%
  slice_max(mean, n=20) %>%
  data.frame()

Top20_TPM_cdlg_R

## For Medicago ================================================================
Top20_TPM_cdlg_M <- Sum.TPM_cdlg %>%
  filter(str_detect(Geneid, "CDO30")) %>%
  select(Geneid, Proteinid, Chr, gene_biotype, Description, mean) %>%
  slice_max(mean, n=25) %>%
  data.frame()

Top20_TPM_cdlg_M <- Sum.TPM_cdlg %>%
  filter(str_detect(Geneid, "MtrunA")) %>%
  filter(!str_detect(gene_biotype, "ncRNA")) %>%
  filter(!grepl("MtrunA17_Chr5g0422301", Geneid)) %>% ## from BLAST
  select(Geneid, Proteinid, Chr, gene_biotype, Description, mean) %>%
  slice_max(mean, n=20) %>%
  data.frame()

#### Step III-c. Read counts ===================================================
### Calculating mean, sd, and CV
Stat.Reads <- Data_readcounts %>%
  pivot_longer(cols = SRR18707692:SRR18707778,
               names_to = ("library"), values_to = "Reads") %>%
  group_by(Geneid) %>%
  summarise(
    mean = mean(Reads, na.rm = TRUE),
    sd = sd(Reads, na.rm = TRUE),
    CV = (sd/mean) * 100,
    .groups = "drop"
  )

#### Step III-d. TPM ===========================================================
### Calculating mean, sd, and CV
Stat.TPM <- Data_TPM %>%
  pivot_longer(cols = SRR18707692:SRR18707778,
               names_to = ("library"), values_to = "TPM") %>%
  group_by(Geneid) %>%
  summarise(
    mean = mean(TPM, na.rm = TRUE),
    sd = sd(TPM, na.rm = TRUE),
    CV = (sd/mean) * 100,
    .groups = "drop"
  )

#### Step III-d. Log2(TPM+1) ===================================================
### Calculating mean, sd, and CV
Stat.Log2TPM <- Data_Log2TPM %>%
  pivot_longer(cols = SRR18707692:SRR18707778,
               names_to = ("library"), values_to = "Log2TPM") %>%
  group_by(Geneid) %>%
  summarise(
    mean = mean(Log2TPM, na.rm = TRUE),
    sd = sd(Log2TPM, na.rm = TRUE),
    CV = (sd/mean) * 100,
    .groups = "drop"
  )

#### Step III-e. Combining TPM/coding length + TPM =============================
## (Table S4) For Sinorhizobia  ================================================
Sum_Top20_R <- Top20_TPM_cdlg_R %>%
  left_join(Stat.TPM %>% select (Geneid, mean) %>% rename(meanTPM = mean),
            by = c("Geneid"))

## (Table S6) For Medicago================================================
Sum_Top20_M <- Top20_TPM_cdlg_M %>%
  left_join(Stat.TPM %>% select (Geneid, mean) %>% rename(meanTPM = mean),
            by = c("Geneid"))

### End of Step III  ===========================================================

##### Step IV. Highly variable genes (CV) based on Log2(TPM+1) =================
Sum.Log2TPM <- Stat.Log2TPM %>%
  left_join(anno_raw, by = c("Geneid"))
  
## (Table S7) For Sinorhizobia =================================================
CV_Top20_R <- Sum.Log2TPM %>%
  filter(str_detect(Geneid, "CDO30")) %>%
  select(Geneid, Chr, gene_biotype, Description, mean, CV) %>%
  slice_max(mean, n=1000) %>%
  slice_max(CV, n=20) %>%
  left_join(Stat.TPM_cdlg %>% select (Geneid, mean) %>% rename(meanTPM_cdlg = mean),
            by = "Geneid") %>%
  data.frame()

CV_Top20_R

## (Table S8) For Medicago =====================================================
CV_Top20_M <- Sum.Log2TPM %>%
  filter(str_detect(Geneid, "MtrunA")) %>%
  filter(!str_detect(gene_biotype, "ncRNA")) %>%
  filter(!grepl("MtrunA17_Chr6g0488421", Geneid)) %>% ## from BLAST
  select(Geneid, Chr, gene_biotype, Description, mean, CV) %>%
  slice_max(mean, n=1000) %>%
  slice_max(CV, n=20) %>%
  left_join(Stat.TPM_cdlg %>% select (Geneid, mean) %>% rename(meanTPM_cdlg = mean),
            by = "Geneid") %>%
  data.frame()

CV_Top20_M

### End of Step IV  ============================================================

##### Step V. Relationships between Log2(TPM+1) and CV =========================
Sum.Log2TPM_2 <- Sum.Log2TPM %>%
  left_join(Stat.TPM %>% select(Geneid, mean) %>%
              rename(meanTPM = mean), by = c("Geneid"))

#### (Figure S1-a) Distribution of expression for all expressed genes ==========
Sum.Log2TPM_2 %>% 
  mutate(
    Species = case_when(
      str_starts(Geneid, "^CDO30")  ~ "Sinorhizobium",
      str_starts(Geneid, "^MtrunA") ~ "Medicago",
      TRUE                       ~ NA_character_
    )
  ) %>%
  ggplot(aes(x = mean, fill = Species)) +
  geom_density(aes(y = ..count..), alpha = 0.5, linewidth = 0.2) +
  scale_fill_manual(
    values = c("Sinorhizobium" = "#CC6600", "Medicago" = "#56B4E9"),
    na.value = "grey80"      # optional: color for any NA/missing Geneid patterns
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    y = "Counts",
    x = "Log2(TPM+1)",
    title = "a) All genes",
    subtitle = "",
    fill = "Species"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9),
    # legend.position = c(0.85, 0.5),
    legend.justification = "right",
    legend.position = c(0.9, 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.6, "lines")
  )

# ggsave("Log2TPM_Distribution.jpg", height = 4, width = 4, bg ="white", dpi = 650)

#### (Figure S1-b) Distribution of expression for highly expressed 1k genes ====
Sum.Log2TPM_2 %>% 
  slice_max(order_by = meanTPM, n=1000) %>%
  mutate(
    Species = case_when(
      str_starts(Geneid, "^CDO30")  ~ "Sinorhizobium",
      str_starts(Geneid, "^MtrunA") ~ "Medicago",
      TRUE                       ~ NA_character_
    )
  ) %>%
  ggplot(aes(x = mean, fill = Species)) +
  geom_density(aes(y = ..count..), alpha = 0.5, linewidth = 0.2) +
  scale_fill_manual(
    values = c("Sinorhizobium" = "#CC6600", "Medicago" = "#56B4E9"),
    na.value = "grey80"      # optional: color for any NA/missing Geneid patterns
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    y = "Counts",
    x = "Log2(TPM+1)",
    title = "b) Highly expressed 1k genes",
    subtitle = "",
    fill = "Species"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9),
    # legend.position = c(0.85, 0.5),
    legend.justification = "right",
    legend.position = c(0.9, 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.6, "lines")
  )

# ggsave("Log2TPM_Distribution_1k.jpg", height = 4, width = 4, bg ="white", dpi = 650)

#### (Figure S1-c) Distribution of CV for all expressed genes ==========
Sum.Log2TPM_2 %>% 
  mutate(Species = case_when(
    str_detect(Geneid, "^CDO30") ~ "Sinorhizobium",
    str_detect(Geneid, "^MtrunA") ~ "Medicago"
  )) %>%
  ggplot(aes(x = CV, fill = Species, color = Species)) +
  geom_density(aes(y = ..count..), alpha = 0.5, size = 0.3)+
  scale_fill_manual(
    values = c("Sinorhizobium" = "#CC6600", "Medicago" = "#56B4E9"),
    na.value = "grey80"      # optional: color for any NA/missing Geneid patterns
  ) +
  scale_color_manual(
    values = c("sense" = "#CC6600", "antisense" = "#56B4E9")) +
  labs(y = "Counts", x = "CV",
       title = "c) All genes",
       subtitle = " ") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9),
    # legend.position = c(0.85, 0.5),
    legend.justification = "right",
    legend.position = c(0.9, 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.6, "lines")
  )

# ggsave("CV_Distribution.jpg", height = 4, width = 4, bg ="white", dpi = 650)

#### (Figure S1-d) Distribution of CV for for highly expressed 1k genes =========
Sum.Log2TPM_2 %>% 
  slice_max(order_by = meanTPM, n=1000) %>%
  mutate(Species = case_when(
    str_detect(Geneid, "^CDO30") ~ "Sinorhizobium",
    str_detect(Geneid, "^MtrunA") ~ "Medicago"
  )) %>%
  ggplot(aes(x = CV, fill = Species, color = Species)) +
  geom_density(aes(y = ..count..), alpha = 0.5, size = 0.3)+
  scale_fill_manual(
    values = c("Sinorhizobium" = "#CC6600", "Medicago" = "#56B4E9"),
    na.value = "grey80"      # optional: color for any NA/missing Geneid patterns
  ) +
  scale_color_manual(
    values = c("sense" = "#CC6600", "antisense" = "#56B4E9")) +
  labs(y = "Counts", x = "Log2(TPM+1)",
       title = "d) Highly expressed 1k genes",
       subtitle = " ") +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9),
    # legend.position = c(0.85, 0.5),
    legend.justification = "right",
    legend.position = c(0.9, 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.6, "lines")
  )

# ggsave("CV_Distribution_1k.jpg", height = 4, width = 4, bg ="white", dpi = 650)

#### (Figure S1-e) Scatter plots between CV and Log(TPM+1) for all genes =======
df <- Sum.Log2TPM_2 %>%
  mutate(
    Species = case_when(
      str_detect(Geneid, "^CDO30")  ~ "Sinorhizobium",
      str_detect(Geneid, "^MtrunA")  ~ "Medicago",
      TRUE                           ~ NA_character_
    )) %>%
  filter(!is.na(Species))

ggplot() +
  geom_point(
    data      = filter(df, Species == "Medicago"),
    aes(x = CV, y = mean, fill = Species),
    shape     = 21,
    color     = "grey80",   
    size      = 1.0,
    alpha     = 0.4,
  ) +
  geom_point(
    data      = filter(df, Species == "Sinorhizobium"),
    aes(x = CV, y = mean, fill = Species), 
    shape     = 21,
    color     = "grey80",
    size      = 1.0,
    alpha     = 0.4
  ) +
  scale_fill_manual(
    values   = c(
      "Sinorhizobium" = "#CC6600",
      "Medicago"      = "#56B4E9"
    ),
    na.value = "grey80"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    y     = "Log2(TPM+1)",
    x     = "CV",
    title = "e) All genes",
    subtitle = " ",
    fill  = "Species" 
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9),
    # legend.position = c(0.85, 0.5),
    legend.justification = "right",
    legend.position = c(0.9, 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.6, "lines")
  )

# ggsave("mean_CV_Distribution.jpg", height = 4, width = 4, bg ="white", dpi = 650)


#### (Figure S1-f) Scatter plots between CV and Log(TPM+1) for HE 1k genes ====
df_1k <- Sum.Log2TPM_2 %>%
  slice_max(order_by = mean, n=1000) %>%
  mutate(
    Species = case_when(
      str_detect(Geneid, "^CDO30")  ~ "Sinorhizobium",
      str_detect(Geneid, "^MtrunA")  ~ "Medicago",
      TRUE                           ~ NA_character_
    )) %>%
  filter(!is.na(Species))

ggplot() +
  geom_point(
    data      = filter(df_1k, Species == "Medicago"),
    aes(x = CV, y = mean, fill = Species),
    shape     = 21,
    color     = "grey80", 
    size      = 1.0,
    alpha     = 0.4
  ) +
  geom_point(
    data      = filter(df_1k, Species == "Sinorhizobium"),
    aes(x = CV, y = mean, fill = Species),
    shape     = 21,
    color     = "grey80",
    size      = 1.0,
    alpha     = 0.4
  ) +
  scale_fill_manual(
    values   = c(
      "Sinorhizobium" = "#CC6600",
      "Medicago"      = "#56B4E9"
    ),
    na.value = "grey80"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    y     = "Log2(TPM+1)",
    x     = "CV",
    title = "f) Highly expressed 1k genes",
    subtitle = "",
    fill  = "Species"  
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9),
    # legend.position = c(0.85, 0.5),
    legend.justification = "right",
    legend.position = c(0.9, 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.6, "lines")
  )

ggsave("mean_CV_Distribution_1k.jpg", height = 4, width = 4, bg ="white", dpi = 650)

### End of Step V  =============================================================

##### Step VI. Correlations between genes ======================================
dim(Data_raw) ## 26,428 68
dim(Data_Log2TPM) ## 26,428 53

cor <- stats::cor
CorrMat <- cor(t(Data_Log2TPM[apply(Data_Log2TPM[,-c(1:3)], 1, sd) != 0,
                              -c(1:3)]), use = "pairwise.complete.obs")

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

### Labelling gene names on the correlation matrix
gene_names <- Data_Log2TPM_2$Geneid
length(gene_names)

CorrMat_tri <- CorrMat
CorrMat_tri[lower.tri(CorrMat_tri)] <- NA
rownames(CorrMat_tri) <- gene_names
colnames(CorrMat_tri) <- gene_names

### Step VI-b. Calculating t, pvalue, fdr for every gene correlation pairs =====
# number of library
number_of_library <- ncol(Data_Log2TPM_2[,-c(1:3)])
number_of_library ### 50

# Calculation t, pvalue, fdr
CorrMat_table <- CorrMat_tri %>% 
  as.data.frame() %>% 
  rownames_to_column("from") %>%
  pivot_longer(cols = -from,
               names_to = "to",
               values_to = "r") %>% 
  filter(is.na(r) == F, from != to) %>% 
  mutate(t = r*sqrt((number_of_library-2)/(1-r^2))) %>% 
  mutate(p.value = case_when(
    t > 0 ~ pt(t, df = number_of_library-2, lower.tail = F),
    TRUE ~ pt(t, df = number_of_library-2, lower.tail = T)
  ),
  FDR = p.adjust(p.value, method = "fdr")
  ) %>%
  data.frame()

head(CorrMat_table)

### Step VI-c. Searching the highly correlated interspecies genes ==============
## Step VI-c-1. Intraspecies for Sinorhizobia (Within species) =================
CorrMat_table_wSn <- CorrMat_table %>%
  filter(str_detect(from, "^CDO30")) %>%
  filter(str_detect(to, "^CDO30"))
CorrMat_table_wSn

## Step VI-c-2. Intraspecies for Medicago (Within species) =====================
CorrMat_table_wMd <- CorrMat_table %>%
  filter(str_detect(from, "^MtrunA")) %>%
  filter(str_detect(to, "^MtrunA"))
CorrMat_table_wMd

## Step VI-c-3. Interspecies gene pairs (between species) ======================
CorrMat_table_bSM <- CorrMat_table %>%
  filter(str_detect(from, "^CDO30")) %>%
  filter(str_detect(to, "^MtrunA"))
CorrMat_table_bSM

## Step VI-c-4. (Dataset S2) Highly correlated interspecieis genes =============
## Positive
CorrMat_bSM_Pos <- CorrMat_table_bSM %>%
  left_join(anno_raw %>%
              filter(str_detect(Geneid, "CDO30_")) %>%
                       select(Geneid, Proteinid, Description, Chr) %>%
              rename(Description_from = Description,
                     Replicon = Chr),
            by = c("from" = "Geneid")) %>%
  left_join(anno_raw %>%
              filter(str_detect(Geneid, "MtrunA")) %>%
                       select (Geneid, Description) %>%
                       rename(Description_to = Description),
            by = c("to" = "Geneid")) %>%
  select(-t, -p.value, -FDR) %>%
  slice_max(order_by = r, n=100) %>%
  arrange(-r) %>%
  mutate(r = sprintf("%.2f", r)) %>%
  data.frame()

# Counting the number of interspecies correlation pairs
CorrMat_table_bSM %>%
  select(-t, -p.value, -FDR) %>%
  slice_max(order_by = r, n=100) %>%
  group_by(from) %>%
  count() %>%
  arrange(-n) %>%
  data.frame()

## Negative
CorrMat_bSM_Neg <- CorrMat_table_bSM %>%
  left_join(anno_raw %>%
              filter(str_detect(Geneid, "CDO30_")) %>%
              select(Geneid, Proteinid, Description, Chr) %>%
              rename(Description_from = Description,
                     Replicon = Chr),
            by = c("from" = "Geneid")) %>%
  left_join(anno_raw %>%
              filter(str_detect(Geneid, "MtrunA")) %>%
              select (Geneid, Description) %>%
              rename(Description_to = Description),
            by = c("to" = "Geneid")) %>%
  select(-t, -p.value, -FDR) %>%
  slice_max(order_by = -r, n=100) %>%
  arrange(r) %>%
  mutate(r = sprintf("%.2f", r)) %>%
  data.frame()

# Counting the number of interspecies correlation pairs
CorrMat_table_bSM %>%
  select(-t, -p.value, -FDR) %>%
  slice_max(order_by = -r, n=100) %>%
  group_by(from) %>%
  count() %>%
  arrange(-n) %>%
  data.frame()


## Step VI-c-5. Intraspecies correlations in the interspecies correlations =====
# Listing genes from each correlation and species
# Gene list 
# For Sinorhizobia
Genes_Ps_bSM_Sino <- CorrMat_bSM_Pos %>%
  group_by(from) %>%
  count() %>%
  arrange(-n)

Genes_Ng_bSM_Sino <- CorrMat_bSM_Neg %>%
  group_by(from) %>%
  count() %>%
  arrange(-n)

# For Medicago
Genes_Ps_bSM_Medi <- CorrMat_bSM_Pos %>%
  group_by(to) %>%
  count() %>%
  arrange(-n)

Genes_Ng_bSM_Medi <- CorrMat_bSM_Neg %>%
  group_by(to) %>%
  count() %>%
  arrange(-n)

## Step VI-c-6. (Table 1 and Dataset S2) Within Sinorhizobia ===================
Corr_bSM_Sino <- rbind(Genes_Ps_bSM_Sino, Genes_Ng_bSM_Sino)
genes_of_interest <- Corr_bSM_Sino[[1]] %>% unique()

## Making tables for Intraspecies correlations
Intra_Corr_Sino <- CorrMat_table_wSn %>%
  select(from, to, r) %>%
  dplyr::filter(from %in% genes_of_interest & to %in% genes_of_interest)

# write.csv(Intra_Corr_Sino, "Table_correlations_among_SinoGenes.csv")

## Step VI-c-7. (Table 2 and Dataset S2) Within Sinorhizobia ===================
Corr_bSM_Medi <- rbind(Genes_Ps_bSM_Medi, Genes_Ng_bSM_Medi)
genes_of_interest <- Corr_bSM_Medi[[1]] %>% unique()

## Making tables for Intraspecies correlations
Intra_Corr_Medi <- CorrMat_table_wMd %>%
  select(from, to, r) %>%
  dplyr::filter(from %in% genes_of_interest & to %in% genes_of_interest)

# write.csv(Intra_Corr_Medi, "Table_correlations_among_MediGenes.csv")

### End of Step VI  ============================================================