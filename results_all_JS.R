library("tidyverse")
library("DESeq2")
library("pheatmap")

local <- 0
prefix <- ifelse(local, "~", "")

results_dir <- file.path(prefix, "data/bioinf/projects/data/2022_MLseq_KP/results")

dir_cervicalD <- file.path(prefix, "data/bioinf/projects/data/2022_MLseq_KP/cervical_DESeq2/results/")
dir_cervicalM <- file.path(prefix, "data/bioinf/projects/data/2022_MLseq_KP/cervical_MLSeq/results/")
dir_lungD <- file.path(prefix, "data/bioinf/projects/data/2022_MLseq_KP/lung_DESeq2/results/")
dir_lungM <- file.path(prefix, "data/bioinf/projects/data/2022_MLseq_KP/lung_MLSeq/results/")
dir_MAITD <- file.path(prefix, "data/bioinf/projects/data/2022_MLseq_KP/MAIT_DESeq2/results/")
dir_MAITM <- file.path(prefix, "data/bioinf/projects/data/2022_MLseq_KP/MAIT_MLSeq/results/")

results <- read_csv(file.path(dir_cervicalD, "07_results_cervical_all.csv")) %>%
  mutate("gene_name" = `...1`) %>%
  select(-`...1`)
results %>% head()

load(file.path(dir_cervicalM, "07_pam.rda"))
load(file.path(dir_cervicalM, "08_plda.rda"))
load(file.path(dir_cervicalM, "09_plda2.rda"))
load(file.path(dir_cervicalM, "10_vnsc.rda"))

pam <- tibble("gene_name" = genes.pam) %>%
  mutate("pam" = 1)
plda <- tibble("gene_name" = genes.plda) %>%
  mutate("plda" = 1)
plda2 <- tibble("gene_name" = genes.plda2) %>%
  mutate("plda2" = 1)
vnsc <- tibble("gene_name" = genes.vnsc) %>%
  mutate("vnsc" = 1)

mlseq <- full_join(pam, full_join(plda, full_join(plda2, vnsc)))

deseq2_mlseq <- left_join(results, mlseq, by = "gene_name") %>%
  replace_na(list(pam = 0, plda = 0, plda2 = 0, vnsc = 0))


##########################

load(file = file.path(dir_cervicalD, "02_rld.Rdata"))
rld
# load( file= file.path(dir_cervicalD, "04_dds.Rdata"))
load(file = file.path(dir_cervicalD, "03_dds_pre_processed.Rdata"))
dds <- estimateSizeFactors(dds)
select <- order(rowMeans(counts(dds, normalized = TRUE)),
                decreasing = TRUE)[1:20]


df <- as.data.frame(colData(dds)[, c("condition")])
row <- results
rownames(row) <- row[, 1]
#row <- (row)[select,]
temp <- as.data.frame(assay(rld[select]))
temp$MLSeq <- 0

a <- 0
list <- row.names(temp)
a <- 0
for (gene in list) {
  a <- a + 1
  pos <- which(results == gene)
  print(gene)
  print(results[pos, 8])
  temp[a, 59] <- results[pos, 8]
  #use variable that increases
  #copy value from results to temp
}
#temp <- temp[,59]

assay(rld) %>% head()

pam_data <- deseq2_mlseq %>%
  filter(gene_name %in% (assay(rld)[select,] %>% rownames())) %>%
  select(gene_name, pam)
row_anno_df <- data.frame(pam_data$pam, row.names = pam_data$gene_name)

#How do I get the MLSeq data into the dds?
pheatmap(assay(rld)[select,],
         cluster_rows = TRUE,
         show_rownames = TRUE,
         cluster_cols = TRUE,
         # annotation_col = df,
         annotation_row = row_anno_df
)
