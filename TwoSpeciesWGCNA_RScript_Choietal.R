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

##### Step VII. WGCNA ==========================================================
#### Step VII-a. Log2(TPM+1) Matrix with entire genes (S+M model) ==============
Data_Log2TPM_2
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

##### Step IX. GO enrichment analysis ==========================================
#### Step IX-a. GO References for Sinorhizobia =================================
gtf_R <- import("GCF_002197065.1_ASM219706v1_genomic.gtf.gz")
gtf_R <- as.data.frame(gtf_R)
head(gtf_R, n=6)

gtf_R <- gtf_R %>%
  select("gene_id", "gene", "gene_biotype", "locus_tag", "Ontology_term", "go_component", "go_function", "go_process", "product")
head(gtf_R, n=8)

### Step IX.a-1. Component =====================================================
gtf_R_Comp1 <- gtf_R %>%
  select("gene_id","go_component") %>%
  mutate(go_component = str_extract(go_component, "^[^|]+"))
colnames(gtf_R_Comp1) <- c("GeneID", "GO_CC_Description")
gtf_R_Comp1$GO_component <- gtf_R$go_component

gtf_R_Comp1 <- gtf_R_Comp1 %>%
  mutate(GO_component = sub(".*\\|(\\d+)\\|.*", "\\1", GO_component)) %>%
  mutate(GO_component = ifelse(is.na(GO_component), NA, paste0("GO:", GO_component))) %>%
  select("GeneID", "GO_component", "GO_CC_Description")
head(gtf_R_Comp1)

gtf_R_Comp <- gtf_R_Comp1 %>%
  dplyr::slice(seq(2, nrow(gtf_R_Comp1), 5)) 

str(gtf_R_Comp) ## 5,053

### Step IX.a-2. Functions =====================================================
gtf_R_Func1 <- gtf_R %>%
  select("gene_id","go_function") %>%
  mutate(go_function = str_extract(go_function, "^[^|]+"))
colnames(gtf_R_Func1) <- c("GeneID", "GO_MF_Description")
gtf_R_Func1$GO_function <- gtf_R$go_function

gtf_R_Func1 <- gtf_R_Func1 %>%
  mutate(GO_function = sub(".*\\|(\\d+)\\|.*", "\\1", GO_function)) %>%
  mutate(GO_function = ifelse(is.na(GO_function), NA, paste0("GO:", GO_function))) %>%
  select("GeneID", "GO_function", "GO_MF_Description")
head(gtf_R_Func1)

gtf_R_Func <- gtf_R_Func1 %>%
  dplyr::slice(seq(2, nrow(gtf_R_Func1), 5)) 

str(gtf_R_Func) ## 5,053

### Step IX.a-3. Processes =====================================================
gtf_R_Proc1 <- gtf_R %>%
  select("gene_id","go_process") %>%
  mutate(go_process = str_extract(go_process, "^[^|]+"))
colnames(gtf_R_Proc1) <- c("GeneID", "GO_BP_Description")
gtf_R_Proc1$GO_process <- gtf_R$go_process

gtf_R_Proc1 <- gtf_R_Proc1 %>%
  mutate(GO_process = sub(".*\\|(\\d+)\\|.*", "\\1", GO_process)) %>%
  mutate(GO_process = ifelse(is.na(GO_process), NA, paste0("GO:", GO_process))) %>%
  select("GeneID", "GO_process", "GO_BP_Description")
head(gtf_R_Proc1)

gtf_R_Proc <- gtf_R_Proc1 %>%
  dplyr::slice(seq(2, nrow(gtf_R_Proc1), 5)) 

str(gtf_R_Proc) ## 5,053


#### Step IX-b. GO References for Medicago =====================================
gtf <- import("MtrunA17r5.0-ANR-EGN-r1.9.gtf")
ann <- read_tsv("MtrunA17r5.0-ANR-EGN-r1.9.b2g.tsv")
xls <- readxl::read_excel("MtrunA17r5.0-ANR-EGN-r1.9.xlsx")

gtf <- as.data.frame(gtf)
ann <- as.data.frame(ann)
xls <- as.data.frame(xls)

head(gtf, n=6) ## gene_id, locus_tag
head(ann, n=6) ## Sequence Name, Annotation GO ID, Sequence Description
head(xls, n=6) ## seqname, description, go_ids, go_names

gtf_sort <- gtf %>%
  dplyr::select("gene_id", "locus_tag")

ann_sort <- ann %>%
  dplyr::select("Sequence Name", "Annotation GO ID", "Sequence Description")

xls_sort <- xls %>%
  dplyr::select("seqname", "description", "go_ids", "go_names")

head(gtf_sort, n=4)
head(ann_sort, n=4)
head(xls_sort, n=4)

### Step IX-b-1. Total =========================================================
gtf_M <- xls_sort %>%
  drop_na(go_ids) %>%
  separate_rows(go_ids, go_names, sep = ";") %>%
  separate(go_ids, into = c("category", "GO"), sep = ":GO:", remove = FALSE) %>%
  separate(go_names, into = c("names", "extra"), sep = ":", remove = FALSE) %>%
  mutate(names = coalesce(names, extra)) %>%
  dplyr::select(seqname, description, category, GO, names, extra) %>%
  mutate(GO=paste0("GO:", GO)) %>%
  filter(!is.na(GO) & !is.na(names)) %>%
  distinct(seqname, description, category, GO, names, .keep_all = TRUE) %>%
  dplyr::select(seqname, category, GO, extra)

colnames(gtf_M) <- c("seqname", "category", "GO", "description")
gtf_M <- as.data.frame(gtf_M)
head(gtf_M, n=6)

### Step IX-b-2. Component =====================================================
gtf_M_Comp <- gtf_M %>%
  filter(str_detect(category, "C"))
colnames(gtf_M_Comp) <- c("GeneID", "category", "GO_component", "GO_CC_Description")
str(gtf_M_Comp) ## 42676

### Step IX-b-3. Functions =====================================================
gtf_M_Func <- gtf_M %>%
  filter(str_detect(category, "F"))
colnames(gtf_M_Func) <- c("GeneID", "category", "GO_function", "GO_MF_Description")
str(gtf_M_Func) ## 64336

### Step IX-b-4. Processes =====================================================
gtf_M_Proc <- gtf_M %>%
  filter(str_detect(category, "P"))
colnames(gtf_M_Proc) <- c("GeneID", "category", "GO_process", "GO_BP_Description")
str(gtf_M_Proc) ## 44981

#### Step IX-c. Making TERM2GENE and TERM2NAME =================================
### TERM2GENE for mapping between GO terms and genes ===========================
### TERM2NAME for mapping between GO terms and their description ===============
### Step IX-c-1. Sino_Component ================================================
gtf_R_GO <- gtf_R_Comp[!is.na(gtf_R_Comp$GO_component), ]
TERM2GENE_R_CC <- gtf_R_GO[, c("GO_component", "GeneID")]
colnames(TERM2GENE_R_CC) <- c("ID", "gene_id")
TERM2NAME_R_CC <- unique(gtf_R_GO[, c("GO_component", "GO_CC_Description")])
colnames(TERM2NAME_R_CC) <- c("ID", "gene_name")

### Step IX-c-2. Sino_Functions ================================================
gtf_R_GO <- gtf_R_Func[!is.na(gtf_R_Func$GO_function), ]
TERM2GENE_R_MF <- gtf_R_GO[, c("GO_function", "GeneID")]
colnames(TERM2GENE_R_MF) <- c("ID", "gene_id")
TERM2NAME_R_MF <- unique(gtf_R_GO[, c("GO_function", "GO_MF_Description")])
colnames(TERM2NAME_R_MF) <- c("ID", "gene_name")

### Step IX-c-3. Sino_Processes ================================================
gtf_R_GO <- gtf_R_Proc[!is.na(gtf_R_Proc$GO_process), ]
TERM2GENE_R_BP <- gtf_R_GO[, c("GO_process", "GeneID")]
colnames(TERM2GENE_R_BP) <- c("ID", "gene_id")
TERM2NAME_R_BP <- unique(gtf_R_GO[, c("GO_process", "GO_BP_Description")])
colnames(TERM2NAME_R_BP) <- c("ID", "gene_name")

### Step IX-c-4. Medi_Component ================================================
gtf_M_GO <- gtf_M_Comp[!is.na(gtf_M_Comp$GO_component), ]
TERM2GENE_M_CC <- gtf_M_GO[, c("GO_component", "GeneID")]
colnames(TERM2GENE_M_CC) <- c("ID", "gene_id")
TERM2NAME_M_CC <- unique(gtf_M_GO[, c("GO_component", "GO_CC_Description")])
colnames(TERM2NAME_M_CC) <- c("ID", "gene_name")

### Step IX-c-5. Medi_Functions ================================================
gtf_M_GO <- gtf_M_Func[!is.na(gtf_M_Func$GO_function), ]
TERM2GENE_M_MF <- gtf_M_GO[, c("GO_function", "GeneID")]
colnames(TERM2GENE_M_MF) <- c("ID", "gene_id")
TERM2NAME_M_MF <- unique(gtf_M_GO[, c("GO_function", "GO_MF_Description")])
colnames(TERM2NAME_M_MF) <- c("ID", "gene_name")

### Step IX-c-6. Medi_Processes ================================================
gtf_M_GO <- gtf_M_Proc[!is.na(gtf_M_Proc$GO_process), ]
TERM2GENE_M_BP <- gtf_M_GO[, c("GO_process", "GeneID")]
colnames(TERM2GENE_M_BP) <- c("ID", "gene_id")
TERM2NAME_M_BP <- unique(gtf_M_GO[, c("GO_process", "GO_BP_Description")])
colnames(TERM2NAME_M_BP) <- c("ID", "gene_name")

#### Step IX-d. Universal genes ================================================
### Step IX-d-1. Sinorhizobia genes ============================================
universe_list_R <- summary_table_Set7 %>%
  filter(str_detect(GeneID, "CDO30"))
universe_list_R <- as.vector(universe_list_R$GeneID)

### Step IX-d-2. Medicago genes ================================================
universe_list_M <- summary_table_Set7 %>%
  filter(str_detect(GeneID, "MtrunA"))
universe_list_M <- as.vector(universe_list_M$GeneID)

#### Step IX-e. Highly correlated interspecies genes ===========================
### Step IX-e-1. Sinorhizobia genes ============================================
## Step IX-e-1-a. Positive 200 genes ===========================================
Genes_R_top200_pos <- CorrMat_table_bSM %>%
  select(-t, -p.value, -FDR) %>%
  slice_max(order_by = r, n = 2000) %>%
  arrange(-r) %>%
  mutate(r = sprintf("%.2f", r)) %>%
  filter(str_detect(from, "CDO30")) %>%
  select(from) %>%
  distinct() %>%
  slice_head(n = 200)

Genes_R_top200_pos <- Genes_R_top200_pos[,1]

## Step IX-e-1-b. Negative 200 genes ===========================================
Genes_R_top200_neg <- CorrMat_table_bSM %>%
  select(-t, -p.value, -FDR) %>%
  slice_max(order_by = -r, n=2000) %>%
  arrange(r) %>%
  mutate(r = sprintf("%.2f", r)) %>%
  data.frame()  %>%
  filter(str_detect(from, "CDO30")) %>%
  select(from) %>%
  distinct() %>%
  slice_head(n = 200)

Genes_R_top200_neg <- Genes_R_top200_neg[,1]

### Step IX-e-2. Mediago genes =================================================
## Step IX-e-2-a. Positive 200 genes ===========================================
Genes_M_top200_pos <- CorrMat_table_bSM %>%
  select(-t, -p.value, -FDR) %>%
  slice_max(order_by = r, n=2000) %>%
  arrange(-r) %>%
  mutate(r = sprintf("%.2f", r)) %>%
  data.frame() %>%
  filter(str_detect(to, "Mtrun"))  %>%
  select(to) %>%
  distinct() %>%
  slice_head(n = 200)

Genes_M_top200_pos <- Genes_M_top200_pos[,1]

## Step IX-e-2-b. Negative 200 genes ===========================================
Genes_M_top200_neg <- CorrMat_table_bSM %>%
  select(-t, -p.value, -FDR) %>%
  slice_max(order_by = -r, n=2000) %>%
  arrange(r) %>%
  mutate(r = sprintf("%.2f", r)) %>%
  data.frame() %>%
  filter(str_detect(to, "MtrunA")) %>%
  select(to) %>%
  distinct() %>%
  slice_head(n = 200)

Genes_M_top200_neg <- Genes_M_top200_neg[,1]

### Step IX-e-3. GO enrichment for highly correlated interspecies genes ========
## Step IX-e-3-a. Sinorhizobia =================================================
## Genes_R_top200_pos and Genes_R_top200_neg
# Step IX-e-3-a-1. Component ===================================================
Enriched_R_CC <- enricher(
  gene = Genes_R_top200_pos,
  TERM2GENE = TERM2GENE_R_CC,
  TERM2NAME = TERM2NAME_R_CC,
  universe = universe_list_R,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
Enriched_R_CC@ontology <- "CC"
simplified_R_CC <- simplify(Enriched_R_CC, cutoff = 0.7, by = "p.adjust", select_fun = min)
simplified_R_CC %>% data.frame()

# write.csv(simplified_R_CC, "simplified_R_CC_top200_pos.csv")

# Step IX-e-3-a-2. Functions ===================================================
Enriched_R_MF <- enricher(
  gene = Genes_R_top200_pos,
  TERM2GENE = TERM2GENE_R_MF,
  TERM2NAME = TERM2NAME_R_MF,
  universe = universe_list_R,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
Enriched_R_MF@ontology <- "MF"
simplified_R_MF <- simplify(Enriched_R_MF, cutoff = 0.7, by = "p.adjust", select_fun = min)
simplified_R_MF %>% data.frame()

# write.csv(simplified_R_MF, "simplified_R_MF_top200_pos.csv")

# Step IX-e-3-a-3. Processes ===================================================
Enriched_R_BP <- enricher(
  gene = Genes_R_top200_pos,
  TERM2GENE = TERM2GENE_R_BP,
  TERM2NAME = TERM2NAME_R_BP,
  universe = universe_list_R,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
Enriched_R_BP@ontology <- "BP"
simplified_R_BP <- simplify(Enriched_R_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
simplified_R_BP %>% data.frame()

# write.csv(simplified_R_BP, "simplified_R_BP_top200_pos.csv")

## Step IX-e-3-b. Medicago =====================================================
## Genes_M_top200_pos and Genes_M_top200_neg
# Step IX-e-3-b-1. Component ===================================================
Enriched_M_CC <- enricher(
  gene = Genes_M_top200_pos,
  TERM2GENE = TERM2GENE_M_CC,
  TERM2NAME = TERM2NAME_M_CC,
  universe = universe_list_M,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
Enriched_M_CC@ontology <- "CC"
simplified_M_CC <- simplify(Enriched_M_CC, cutoff = 0.7, by = "p.adjust", select_fun = min)
simplified_M_CC %>% data.frame()

# write.csv(simplified_M_CC, "simplified_M_CC_top200_pos.csv")

# Step IX-e-3-b-2. Functions ===================================================
Enriched_M_MF <- enricher(
  gene = Genes_M_top200_pos,
  TERM2GENE = TERM2GENE_M_MF,
  TERM2NAME = TERM2NAME_M_MF,
  universe = universe_list_M,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
Enriched_M_MF@ontology <- "MF"
simplified_M_MF <- simplify(Enriched_M_MF, cutoff = 0.7, by = "p.adjust", select_fun = min)
simplified_M_MF %>% data.frame()

# write.csv(simplified_M_MF, "simplified_M_MF_top200_pos.csv")

# Step IX-e-3-b-3. Processes ===================================================
Enriched_M_BP <- enricher(
  gene = Genes_M_top200_pos,
  TERM2GENE = TERM2GENE_M_BP,
  TERM2NAME = TERM2NAME_M_BP,
  universe = universe_list_M,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
Enriched_M_BP@ontology <- "BP"
simplified_M_BP <- simplify(Enriched_M_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
simplified_M_BP %>% data.frame()

# write.csv(simplified_M_BP, "simplified_M_BP_top200_pos.csv")

#### Step IX-f. Module-wise genes ==============================================
### Step IX-f-1. Sinorhizobia ==================================================
## Step IX-f-1-a. Module selections ============================================
# turquoise
Genes_R_tur <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "turquoise"))
Genes_R_tur <- as.vector(Genes_R_tur$GeneID)

# blue
Genes_R_blu <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "blue"))
Genes_R_blu <- as.vector(Genes_R_blu$GeneID)

# grey
Genes_R_gre <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "grey"))
Genes_R_gre <- as.vector(Genes_R_gre$GeneID)

# brown
Genes_R_bro <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "brown"))
Genes_R_bro <- as.vector(Genes_R_bro$GeneID)

# yellow
Genes_R_yel <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "yellow"))
Genes_R_yel <- as.vector(Genes_R_yel$GeneID)

# green 
Genes_R_grn <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "green"))
Genes_R_grn <- as.vector(Genes_R_grn$GeneID)

# red
Genes_R_red <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "red"))
Genes_R_red <- as.vector(Genes_R_red$GeneID)

# black
Genes_R_bla <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "black"))
Genes_R_bla <- as.vector(Genes_R_bla$GeneID)

# pink
Genes_R_pin <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "pink"))
Genes_R_pin <- as.vector(Genes_R_pin$GeneID)

# magenta
Genes_R_mag <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "magenta"))
Genes_R_mag <- as.vector(Genes_R_mag$GeneID)

# purple 
Genes_R_pur <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "purple"))
Genes_R_pur <- as.vector(Genes_R_pur$GeneID)

# greenyellow
Genes_R_gyl <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "greenyellow"))
Genes_R_gyl <- as.vector(Genes_R_gyl$GeneID)

# tan 
Genes_R_tan <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "tan"))
Genes_R_tan <- as.vector(Genes_R_tan$GeneID)

# salmon 
Genes_R_sal <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "salmon"))
Genes_R_sal <- as.vector(Genes_R_sal$GeneID)

# cyan
Genes_R_cya <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "cyan"))
Genes_R_cya <- as.vector(Genes_R_cya$GeneID)

# midnightblue 
Genes_R_mid <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "midnightblue"))
Genes_R_mid <- as.vector(Genes_R_mid$GeneID)

# lightcyan 
Genes_R_lcy <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "lightcyan"))
Genes_R_lcy <- as.vector(Genes_R_lcy$GeneID)

# grey60 
Genes_R_g60 <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "grey60"))
Genes_R_g60 <- as.vector(Genes_R_g60$GeneID)

# lightgreen
Genes_R_lgr <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "lightgreen"))
Genes_R_lgr <- as.vector(Genes_R_lgr$GeneID)

# lightyellow 
Genes_R_lyr <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "lightyellow"))
Genes_R_lyr <- as.vector(Genes_R_lyr$GeneID)


### Step IX-f-2. Medicago ======================================================
## Step IX-f-2-a. Module selections ============================================
Genes_M_tur <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "turquoise"))
Genes_M_tur <- as.vector(Genes_M_tur$GeneID)

# blue
Genes_M_blu <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "blue"))
Genes_M_blu <- as.vector(Genes_M_blu$GeneID)

# grey
Genes_M_gre <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "grey"))
Genes_M_gre <- as.vector(Genes_M_gre$GeneID)

# brown
Genes_M_bro <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "brown"))
Genes_M_bro <- as.vector(Genes_M_bro$GeneID)

# yellow
Genes_M_yel <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "yellow"))
Genes_M_yel <- as.vector(Genes_M_yel$GeneID)

# green 
Genes_M_grn <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "green"))
Genes_M_grn <- as.vector(Genes_M_grn$GeneID)

# red
Genes_M_red <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "red"))
Genes_M_red <- as.vector(Genes_M_red$GeneID)

# black
Genes_M_bla <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "black"))
Genes_M_bla <- as.vector(Genes_M_bla$GeneID)

# pink
Genes_M_pin <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "pink"))
Genes_M_pin <- as.vector(Genes_M_pin$GeneID)

# magenta
Genes_M_mag <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "magenta"))
Genes_M_mag <- as.vector(Genes_M_mag$GeneID)

# purple 
Genes_M_pur <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "purple"))
Genes_M_pur <- as.vector(Genes_M_pur$GeneID)

# greenyellow
Genes_M_gyl <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "greenyellow"))
Genes_M_gyl <- as.vector(Genes_M_gyl$GeneID)

# tan 
Genes_M_tan <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "tan"))
Genes_M_tan <- as.vector(Genes_M_tan$GeneID)

# salmon 
Genes_M_sal <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "salmon"))
Genes_M_sal <- as.vector(Genes_M_sal$GeneID)

# cyan
Genes_M_cya <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "cyan"))
Genes_M_cya <- as.vector(Genes_M_cya$GeneID)

# midnightblue 
Genes_M_mid <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "midnightblue"))
Genes_M_mid <- as.vector(Genes_M_mid$GeneID)

# lightcyan 
Genes_M_lcy <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "lightcyan"))
Genes_M_lcy <- as.vector(Genes_M_lcy$GeneID)

# grey60 
Genes_M_g60 <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "grey60"))
Genes_M_g60 <- as.vector(Genes_M_g60$GeneID)

# lightgreen
Genes_M_lgr <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "lightgreen"))
Genes_M_lgr <- as.vector(Genes_M_lgr$GeneID)

# lightyellow 
Genes_M_lyr <- summary_table_Set7 %>%
  mutate(GeneID = str_remove(GeneID, "_s$")) %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "lightyellow"))
Genes_M_lyr <- as.vector(Genes_M_lyr$GeneID)






















### Step IX-f-3. GO enrichment for Module-wise genes ===========================
## Step IX-f-3-a. Sinorhizobia =================================================
## Genes_R_tur... to Genes_R_(other modules)
# Step IX-f-3-a-1. Component ===================================================
Enriched_R_CC <- enricher(
  gene = Genes_R_tur,
  TERM2GENE = TERM2GENE_R_CC,
  TERM2NAME = TERM2NAME_R_CC,
  universe = universe_list_R,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
Enriched_R_CC@ontology <- "CC"
simplified_R_CC <- simplify(Enriched_R_CC, cutoff = 0.7, by = "p.adjust", select_fun = min)
simplified_R_CC %>% data.frame()

# write.csv(simplified_R_CC, "simplified_R_CC_tur.csv")

# Step IX-f-3-a-2. Functions ===================================================
Enriched_R_MF <- enricher(
  gene = Genes_R_tur,
  TERM2GENE = TERM2GENE_R_MF,
  TERM2NAME = TERM2NAME_R_MF,
  universe = universe_list_R,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
Enriched_R_MF@ontology <- "MF"
simplified_R_MF <- simplify(Enriched_R_MF, cutoff = 0.7, by = "p.adjust", select_fun = min)
simplified_R_MF %>% data.frame()

# write.csv(simplified_R_MF, "simplified_R_MF_tur.csv")

# Step IX-f-3-a-3. Processes ===================================================
Enriched_R_BP <- enricher(
  gene = Genes_R_tur,
  TERM2GENE = TERM2GENE_R_BP,
  TERM2NAME = TERM2NAME_R_BP,
  universe = universe_list_R,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
Enriched_R_BP@ontology <- "BP"
simplified_R_BP <- simplify(Enriched_R_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
simplified_R_BP %>% data.frame()

# write.csv(simplified_R_BP, "simplified_R_BP_tur.csv")

## Step IX-f-3-b. Medicago =====================================================
## Genes_M_tur... to Genes_M_(other modules)
# Step IX-f-3-b-1. Component ===================================================
Enriched_M_CC <- enricher(
  gene = Genes_M_tur,
  TERM2GENE = TERM2GENE_M_CC,
  TERM2NAME = TERM2NAME_M_CC,
  universe = universe_list_M,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
Enriched_M_CC@ontology <- "CC"
simplified_M_CC <- simplify(Enriched_M_CC, cutoff = 0.7, by = "p.adjust", select_fun = min)
simplified_M_CC %>% data.frame()

# write.csv(simplified_M_CC, "simplified_M_CC_tur.csv")

# Step IX-f-3-b-2. Functions ===================================================
Enriched_M_MF <- enricher(
  gene = Genes_M_tur,
  TERM2GENE = TERM2GENE_M_MF,
  TERM2NAME = TERM2NAME_M_MF,
  universe = universe_list_M,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
Enriched_M_MF@ontology <- "MF"
simplified_M_MF <- simplify(Enriched_M_MF, cutoff = 0.7, by = "p.adjust", select_fun = min)
simplified_M_MF %>% data.frame()

# write.csv(simplified_M_MF, "simplified_M_MF_tur.csv")

# Step IX-f-3-b-3. Processes ===================================================
Enriched_M_BP <- enricher(
  gene = Genes_M_tur,
  TERM2GENE = TERM2GENE_M_BP,
  TERM2NAME = TERM2NAME_M_BP,
  universe = universe_list_M,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 5000,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
Enriched_M_BP@ontology <- "BP"
simplified_M_BP <- simplify(Enriched_M_BP, cutoff = 0.7, by = "p.adjust", select_fun = min)
simplified_M_BP %>% data.frame()

# write.csv(simplified_M_BP, "simplified_M_BP_tur.csv")


### Step IX-f-4. Plotting ======================================================
df_cor <- read.csv("EnrichedGO_highcorr.CSV") ## summary of above GO results
df_go <- read.csv("EnrichedGO_modules_SM_model.CSV") ## summary of above GO results

## Step IX-f-4-a. (Figure 1) Sinorhizobia ======================================
df_cor_R <- df_cor %>%
  drop_na(Species, Correlations, Category, ID, Description, Count) %>%
  filter(str_detect(Species, "Sino"))

df_clean <- df_cor_R %>%
  mutate(
    Category = factor(Category, levels = c("CC","MF","BP")),
    Description_wrapped = str_wrap(Description, width = 25))

df_clean$Correlations <- factor(
  df_clean$Correlations,
  levels = c("Positive", "Negative"))

p1 <- ggplot(
  df_clean, aes(x = Count, y = Description_wrapped, fill = Category)) +
  geom_col(width = 0.7) +
  facet_grid(rows = vars(Correlations), cols = vars(Category), 
             scales = "free", space = "free_y") +
  labs(x = "Gene Counts", y = NULL,
       subtitle = expression(bold(italic("S. meliloti")))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = scales::breaks_extended(n = 3),
                     labels = scales::label_number(big.mark = "", accuracy = 1)) + 
  theme_bw(base_size = 10) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    legend.position = "none",
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    strip.text.y = element_text(size = 5.5, face = "bold"),
    strip.text.x = element_text(size = 7, face = "bold"),
    panel.spacing.y = unit(1.5, "pt"),
    panel.spacing.x = unit(1.5, "pt")
  )

p1

# ggsave("barplot_GOs_highcorr_Sino.jpg", p1, height = 3.8, width = 6, bg ="white", dpi = 800)

## Step IX-f-4-b. (Figure 2) Medicago ==========================================
df_cor_M <- df_cor %>%
  drop_na(Species, Correlations, Category, ID, Description, Count) %>%
  filter(str_detect(Species, "Medi"))

df_clean <- df_cor_M %>%
  mutate(
    Category = factor(Category, levels = c("CC","MF","BP")),
    Description_wrapped = str_wrap(Description, width = 30)
  )

df_clean$Correlations <- factor(
  df_clean$Correlations,
  levels = c("Positive", "Negative"))

p2 <- ggplot(
  df_clean, aes(x = Count, y = Description_wrapped, fill = Category)) +
  geom_col(width = 0.8) +
  facet_grid(rows = vars(Correlations), cols = vars(Category), 
             scales = "free", space = "free_y") +
  labs(x = "Gene Counts", y = NULL,
       subtitle = expression(bold(italic("M. truncatula")))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = scales::breaks_extended(n = 3),
                     labels = scales::label_number(big.mark = "", accuracy = 1)) + 
  scale_fill_manual(
    values = c(
      "CC" = "#F8766D",
      "MF" = "#00BA38",
      "BP" = "#619CFF" 
    )
  ) +
  theme_bw(base_size = 9) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    legend.position = "none",
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    strip.text.y = element_text(size = 7, face = "bold"),
    strip.text.x = element_text(size = 7, face = "bold"),
    panel.spacing.y = unit(1.5, "pt"),
    panel.spacing.x = unit(1.5, "pt")
  )

p2

# ggsave("barplot_GOs_highcorr_Medi.jpg", p2, height = 6, width = 6, bg ="white", dpi = 800)

## Step IX-f-4-c. (Figure 4) Sinorhizobia ======================================
df_go_R <- df_go %>%
  drop_na(Species, Category, ID, Description, Count) %>%
  filter(str_detect(Species, "Sino"))

df_clean <- df_go_R %>%
  mutate(
    GO = str_wrap(str_squish(Description), width = 10),
    Category = factor(Category, levels = c("CC","MF","BP"))
  ) %>%
  drop_na(Module, Category, GO, Count)

df_clean <- df_clean %>%
  mutate(
    Module = str_to_title(Module),
    GO = str_wrap(str_squish(Description), width = 10),
    Description_wrapped = str_wrap(Description, width = 20),
    Category = factor(Category, levels = c("CC","MF","BP"))
  ) %>%
  drop_na(Module, Category, GO, Count, p.adjust)

df_top <- df_clean %>%
  group_by(Module) %>%
  slice_max(order_by = Count, n = 5, with_ties = FALSE) %>%
  ungroup()

module_order <- c(
  "Turquoise", "Blue", "Grey", "Brown", "Yellow", "Green", "Black",
  "Salmon", "Cyan")

df_top$Module <- factor(df_top$Module, levels = module_order)

p1 <- ggplot(
  df_top, aes(x = Count, y = Description_wrapped, fill = Category)) +
  geom_col(width = 0.7) +
  facet_grid(rows = vars(Category), cols = vars(Module), 
             scales = "free", space = "free_y") +
  labs(x = "Gene Counts", y = NULL,
       subtitle = expression(bold(italic("S. meliloti")))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = scales::breaks_extended(n = 4),
                     labels = scales::label_number(big.mark = "", accuracy = 1)) + 
  theme_bw(base_size = 9) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    legend.position = "none",
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 5.5),
    axis.text.y = element_text(size = 5.5),
    strip.text.y = element_text(size = 7, face = "bold"),
    strip.text.x = element_text(size = 6.5, face = "bold"),
    panel.spacing.y = unit(1.5, "pt"),
    panel.spacing.x = unit(1.5, "pt")
  )
p1

# ggsave("barplot_GOs_modules_Sino.jpg", p1, height = 3.5, width = 5.5, bg ="white", dpi = 650)

## Step IX-f-4-d. (Figure 5) Medicago ==========================================
df_go_M <- df_go %>%
  drop_na(Species, Category, ID, Description, Count) %>%
  filter(str_detect(Species, "Medi"))

df_clean <- df_go_M %>%
  mutate(
    GO = str_wrap(str_squish(Description), width = 40),
    Category = factor(Category, levels = c("CC","MF","BP"))
  ) %>%
  drop_na(Module, Category, GO, Count)

df_clean <- df_clean %>%
  mutate(
    Module = str_to_title(Module),
    GO = str_wrap(str_squish(Description), width = 40),
    Description_wrapped = str_wrap(Description, width = 35),
    Category = factor(Category, levels = c("CC","MF","BP"))
  ) %>%
  drop_na(Module, Category, GO, Count, p.adjust)

df_sel <- df_clean %>%
  group_by(Module) %>%  
  ungroup()

unique(df_sel$GO)
module_order <- c(
  "Turquoise", "Blue", "Grey", "Brown", "Yellow", "Green", "Red", "Black",
  "Pink", "Magenta", "Purple", "Greenyellow", "Tan", "Salmon", "Cyan",
  "Midnightblue", "Lightcyan", "Grey60", "Lightgreen", "Lightyellow"
)

df_sel$Module <- factor(df_sel$Module, levels = module_order)

df_sel2 <- df_sel %>%
  filter(Module %in% c("Blue", "Grey", "Yellow", "Black", "Pink", "Magenta", "Salmon", "Lightgreen"))

p3 <- ggplot(
  df_sel2, aes(x = Count, y = Description_wrapped, fill = Category)) +
  geom_col(width = 0.7) +
  facet_grid(rows = vars(Category), cols = vars(Module), 
             scales = "free", space = "free_y") +
  labs(x = "Gene Counts", y = NULL,
       subtitle = expression(bold(italic("M. truncatula")))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = scales::breaks_extended(n = 3),
                     labels = scales::label_number(big.mark = "", accuracy = 1)) + 
  theme_bw(base_size = 8) +
  theme(
    text = element_text(face = "bold", colour = "black"),
    legend.position = "none",
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 5),
    axis.text.y = element_text(size = 6),
    strip.text.y = element_text(size = 7, face = "bold"),
    strip.text.x = element_text(size = 7, face = "bold"),
    panel.spacing.y = unit(3, "pt"),
    panel.spacing.x = unit(3, "pt")
  )

p3

ggsave("barplot_GOs_modules_Medi.jpg", p3, height = 10, width = 12, bg ="white", dpi = 480)



##### End of Step XI ===========================================================

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