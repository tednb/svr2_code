load("~/data/Liver_scRNA/E-MTAB-10553/hep+chol.Rdata")

library(DESeq2)
library(SingleCellExperiment)
library(magrittr)
library(tidyverse)
library(scuttle)
library(Matrix.utils)
library(pheatmap)
# quality control
# hc.o <- perCellQCMetrics(hc.o)
# 1. meta_data
counts <- hc.o@assays$RNA@counts
meta <- hc.o@meta.data
idx <- na.omit(match(hc_gene,rownames(counts)))
mat <- counts[idx,]
mat <- mat[rowSums(mat) > 0,]
sce <- SingleCellExperiment(assays = list(counts= mat),
          colData = meta)
dim(sce)
sce$Cell.type <- as.factor(as.character(sce$Cell.type))
sce$age <- as.factor(sce$age)
sce$Sample <- sub("_",".",sce$Sample)
sce$Sample <- as.factor(sce$Sample)

kids <- purrr::set_names(levels(sce$Cell.type))
nk <- length(kids)
sids <- purrr::set_names(levels(sce$Sample))
ns <- length(sids)

table(sce$Sample)
n_cells <- as.numeric(table(sce$Sample))
m <- match(sids, sce$Sample)
ei <- data.frame(colData(sce)[m,],
                 n_cells, row.names = NULL) %>%
                select(-"Cell.type")
# 2. QC #
# Computing per-cell QC metrics
#is.mito <- grep("MT\\.", rownames(sce))
#per.cell <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
# sce$percent.mt < 10% for human and 5% for mouse
summary(sce$percent.mt)
#outliers: > median absolute deviations (MADs) 
low.total <- isOutlier(per.cell$sum, nmads = 3, type="lower", log=TRUE)
summary(low.total) # no outlier
# per-gene QC
# remove lowly expressed genes which has less than 10 cells with any counts
#sce <- sce[rowSums(counts(sce)>1)>=10,] 

# 3. aggregation
groups <- colData(sce)[,c("Cell.type","Sample")]
pb <- aggregate.Matrix(t(counts(sce)),
                           groupings = groups, fun ="sum")
dim(pb)

splitf <- sapply(str_split(rownames(pb), pattern = "_",  n = 2), `[`, 1)
pb <- split.data.frame(pb,
                       factor(splitf)) %>% 
  lapply(function(u) 
    set_colnames(t(u), 
                 stringr::str_extract(rownames(u), "(?<=_)[a-zA-Z][\\w.]+")))
get_sample_ids <- function(x){
  pb[[x]] %>% colnames()
  } 
de_samples <- map(1:length(kids), get_sample_ids) %>% unlist()
samples_list <- map(1:length(kids), get_sample_ids) 
get_cluster_ids <- function(x){ 
  rep(names(pb)[x], each = length(samples_list[[x]]))
  } 
de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>% unlist()
gg_df <- data.frame(Cell.type = de_cluster_ids, Sample = de_samples) 
gg_df <- left_join(gg_df, ei[, c("Sample", "age")]) 
metadata <- gg_df %>% dplyr::select(Cell.type, Sample, age) 
metadata

clusters <- levels(as.factor(metadata$Cell.type))
# run deseq2 on chol
cluster_metadata <- metadata[which(metadata$Cell.type == clusters[2]), ]
head(cluster_metadata)
# Assign the rownames of the metadata to be the sample IDs
rownames(cluster_metadata) <- cluster_metadata$Sample
head(cluster_metadata)
# Subset the counts to only the Chol cells
count <- pb[[clusters[2]]]
cluster_counts <- data.frame(count[, which(colnames(count) %in% rownames(cluster_metadata))])
# Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(rownames(cluster_metadata) == colnames(cluster_counts))

#deseq2
out <- c(1,4,5,6)
matrix <- cluster_counts[,-out]
cluster_metadata$age <- factor(as.numeric(as.character(cluster_metadata$age)))
cluster_metadata$group <- factor(c("old","mid","mid","old","old","old","mid1","mid1"))
cluster_metadata$Sample <- as.factor(cluster_metadata$Sample)
me <- cluster_metadata[-out,]
dds <- DESeqDataSetFromMatrix(matrix, 
                              colData = me, 
                              design = ~ group)
rld <- rlog(dds, blind=TRUE)
# Plot PCA
DESeq2::plotPCA(rld, intgroup = "group")
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
# Plot heatmap
pheatmap(rld_cor, annotation = cluster_metadata[, c("group"), drop=F])
dds <- DESeq(dds)
res <- results(dds,alpha = 0.05)
plotMA(res)
resultsNames(dds)
res_s <- lfcShrink(dds,res = res,type = "apeglm",coef =  "group_old_vs_mid")
plotMA(res_s)

res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

chol_stat <- res_tbl$stat
hep_stat <- res_tbl$stat

stat <- data.frame(chol = chol_stat,hep = hep_stat)
rownames(stat) <- rownames(cluster_counts)
stat <- stat[complete.cases(stat), ]

















