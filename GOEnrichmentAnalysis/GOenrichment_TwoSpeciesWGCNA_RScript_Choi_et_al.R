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

##### Step IX. GO enrichment analysis ==========================================
Summary_SM <- read.csv("Summary_Module_SM_Model.csv")

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
universe_list_R <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30"))
universe_list_R <- as.vector(universe_list_R$GeneID)

### Step IX-d-2. Medicago genes ================================================
universe_list_M <- Summary_SM %>%
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
Genes_R_tur <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "turquoise"))
Genes_R_tur <- as.vector(Genes_R_tur$GeneID)

# blue
Genes_R_blu <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "blue"))
Genes_R_blu <- as.vector(Genes_R_blu$GeneID)

# grey
Genes_R_gre <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "grey"))
Genes_R_gre <- as.vector(Genes_R_gre$GeneID)

# brown
Genes_R_bro <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "brown"))
Genes_R_bro <- as.vector(Genes_R_bro$GeneID)

# yellow
Genes_R_yel <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "yellow"))
Genes_R_yel <- as.vector(Genes_R_yel$GeneID)

# green 
Genes_R_grn <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "green"))
Genes_R_grn <- as.vector(Genes_R_grn$GeneID)

# red
Genes_R_red <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "red"))
Genes_R_red <- as.vector(Genes_R_red$GeneID)

# black
Genes_R_bla <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "black"))
Genes_R_bla <- as.vector(Genes_R_bla$GeneID)

# pink
Genes_R_pin <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "pink"))
Genes_R_pin <- as.vector(Genes_R_pin$GeneID)

# magenta
Genes_R_mag <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "magenta"))
Genes_R_mag <- as.vector(Genes_R_mag$GeneID)

# purple 
Genes_R_pur <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "purple"))
Genes_R_pur <- as.vector(Genes_R_pur$GeneID)

# greenyellow
Genes_R_gyl <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "greenyellow"))
Genes_R_gyl <- as.vector(Genes_R_gyl$GeneID)

# tan 
Genes_R_tan <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "tan"))
Genes_R_tan <- as.vector(Genes_R_tan$GeneID)

# salmon 
Genes_R_sal <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "salmon"))
Genes_R_sal <- as.vector(Genes_R_sal$GeneID)

# cyan
Genes_R_cya <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "cyan"))
Genes_R_cya <- as.vector(Genes_R_cya$GeneID)

# midnightblue 
Genes_R_mid <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "midnightblue"))
Genes_R_mid <- as.vector(Genes_R_mid$GeneID)

# lightcyan 
Genes_R_lcy <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "lightcyan"))
Genes_R_lcy <- as.vector(Genes_R_lcy$GeneID)

# grey60 
Genes_R_g60 <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "grey60"))
Genes_R_g60 <- as.vector(Genes_R_g60$GeneID)

# lightgreen
Genes_R_lgr <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "lightgreen"))
Genes_R_lgr <- as.vector(Genes_R_lgr$GeneID)

# lightyellow 
Genes_R_lyr <- Summary_SM %>%
  filter(str_detect(GeneID, "CDO30")) %>%
  filter(str_detect(Module, "lightyellow"))
Genes_R_lyr <- as.vector(Genes_R_lyr$GeneID)


### Step IX-f-2. Medicago ======================================================
## Step IX-f-2-a. Module selections ============================================
Genes_M_tur <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "turquoise"))
Genes_M_tur <- as.vector(Genes_M_tur$GeneID)

# blue
Genes_M_blu <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "blue"))
Genes_M_blu <- as.vector(Genes_M_blu$GeneID)

# grey
Genes_M_gre <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "grey"))
Genes_M_gre <- as.vector(Genes_M_gre$GeneID)

# brown
Genes_M_bro <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "brown"))
Genes_M_bro <- as.vector(Genes_M_bro$GeneID)

# yellow
Genes_M_yel <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "yellow"))
Genes_M_yel <- as.vector(Genes_M_yel$GeneID)

# green 
Genes_M_grn <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "green"))
Genes_M_grn <- as.vector(Genes_M_grn$GeneID)

# red
Genes_M_red <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "red"))
Genes_M_red <- as.vector(Genes_M_red$GeneID)

# black
Genes_M_bla <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "black"))
Genes_M_bla <- as.vector(Genes_M_bla$GeneID)

# pink
Genes_M_pin <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "pink"))
Genes_M_pin <- as.vector(Genes_M_pin$GeneID)

# magenta
Genes_M_mag <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "magenta"))
Genes_M_mag <- as.vector(Genes_M_mag$GeneID)

# purple 
Genes_M_pur <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "purple"))
Genes_M_pur <- as.vector(Genes_M_pur$GeneID)

# greenyellow
Genes_M_gyl <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "greenyellow"))
Genes_M_gyl <- as.vector(Genes_M_gyl$GeneID)

# tan 
Genes_M_tan <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "tan"))
Genes_M_tan <- as.vector(Genes_M_tan$GeneID)

# salmon 
Genes_M_sal <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "salmon"))
Genes_M_sal <- as.vector(Genes_M_sal$GeneID)

# cyan
Genes_M_cya <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "cyan"))
Genes_M_cya <- as.vector(Genes_M_cya$GeneID)

# midnightblue 
Genes_M_mid <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "midnightblue"))
Genes_M_mid <- as.vector(Genes_M_mid$GeneID)

# lightcyan 
Genes_M_lcy <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "lightcyan"))
Genes_M_lcy <- as.vector(Genes_M_lcy$GeneID)

# grey60 
Genes_M_g60 <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "grey60"))
Genes_M_g60 <- as.vector(Genes_M_g60$GeneID)

# lightgreen
Genes_M_lgr <- Summary_SM %>%
  filter(str_detect(GeneID, "MtrunA")) %>%
  filter(str_detect(Module, "lightgreen"))
Genes_M_lgr <- as.vector(Genes_M_lgr$GeneID)

# lightyellow 
Genes_M_lyr <- Summary_SM %>%
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
