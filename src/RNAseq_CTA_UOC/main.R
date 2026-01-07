# Análisis del RNAseq genes CTA
# Este script usa renv para asegurar la reproducibilidad


if (!require('renv')) {
  install.packages('renv')
}else{
    "renv ya instalado"
}

#renv::init(bioconductor="3.21")
#  # restaura el rev usando el fichero renv.lock

# Cargamos las librerías 
library(here)
library(rtracklayer)
library(rhdf5)
library(tximport)
library(stringr)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(apeglm)
library(svglite)
library(patchwork)
library(org.Dr.eg.db)
library(clusterProfiler)
library(pheatmap)
library(enrichplot)
# library(readr)
# library(ComplexHeatmap)
# library(RColorBrewer)
# library(circlize)


# library(EnhancedVolcano)


here::i_am("main.R")

t2g_filename <- "results/t2g.csv" # change path if yours is different!
# kallisto_dir <- "kallisto_aug21"
#kallisto_dir <- "kallisto_dec21"
kallisto_dirs <- c("kallisto_aug21", "kallisto_dec21")

# generate output dirs for figures and results
dir.create("results")
dir.create("figures")


# generate t2g
if (!file.exists(t2g_filename)) {
  message("File not detected")
  source(here::here("R", "01-t2g.R"))
  t2g <- generate_t2g(file_name = t2g_filename, write_t2g = T)
} else{
  t2g <- read.csv(file = t2g_filename)
}

# define files (kallisto's output)
#kallisto_abundance_files <- Sys.glob(here::here("data", kallisto_dir, "*", "*.tsv"))

kallisto_abundance_files <- unlist(
  lapply(kallisto_dirs, function(d) {
    Sys.glob(here::here("data", d, "*", "abundance.tsv"))
  })
)

length(kallisto_abundance_files)

samples <- data.frame(
  file   = kallisto_abundance_files,
  batch  = gsub("kallisto_", "", basename(dirname(dirname(kallisto_abundance_files)))),   # kallisto_aug21 / kallisto_dec21
  sample = basename(dirname(kallisto_abundance_files)),            # Control_1, dazl_1, etc.
  stringsAsFactors = FALSE
)

# condición: Control, dazl, ddx4, elavl2, zar1l
samples$condition <- gsub("_[0-9]+$", "", samples$sample)

samples_filt <- subset(samples,!(batch == "aug21" & condition %in% c("dazl", "zar1l")))

table(samples_filt$batch, samples_filt$condition)

samples_filt$sample_id <- paste(samples_filt$batch, samples_filt$sample, sep = "_")
rownames(samples_filt) <- samples_filt$sample_id

# read abundance.tsv (kallisto quant)
source(here::here("R", "02-tximport.R"))
# txi <- read_quant_files(kallisto_files = kallisto_abundance_files, tx2gene = t2g)
txi <- tximport(files = samples_filt$file, type="kallisto", tx2gene = t2g, ignoreTxVersion = T)
colnames(txi$counts) <- samples_filt$sample_id
colnames(txi$abundance) <- samples_filt$sample_id
colnames(txi$length) <- samples_filt$sample_id

rm(t2g);gc()

# define metadata
# condition <- factor((gsub("_[0-9]", "",x=colnames(txi$counts))))
# coldata <- data.frame(row.names = colnames(txi$counts), condition)
# coldata$samples <- factor(rownames(coldata))
# coldata
coldata <- samples_filt[, c("batch", "condition")]
coldata$batch     <- factor(coldata$batch)
coldata$condition <- factor(coldata$condition)
coldata$condition <- relevel(coldata$condition, ref = "Control")
coldata$batch     <- relevel(coldata$batch, ref = "aug21")

# tpm matrix
tpm_mat <- txi$abundance
colnames(tpm_mat) <- colnames(txi$counts)
# head(tpm_mat)
write.table(tpm_mat, here::here("results", "tpm_rnaseq.tsv"), row.names = T, quote = F)

# genera DESeq object
source(here::here("R", "03-deseq2.R"))
# dds <- build_deseq(txi = txi, coldata = coldata, condition = condition, ref="Control")
dds <- build_deseq(txi = txi, coldata = coldata, condition = condition, ref="Control", batch = T)

# filtra objeto deseq
dds <- filter_dds(dds, min_counts = 10)

# modelado estadístico y normalización
dds <- DESeq(dds, test = "Wald")

# VST calculado una sola vez para reutilizar
vsd <- DESeq2::vst(dds, blind = FALSE)


# plot dispersion
plot_SparDisp(dds)
# boxplot normalized counts
boxplot_counts(dds)

# VST (variance stabilizing transformation) + PCA
# vst_pca(dds)
vst_pca(dds, check_batch = T, vsd = vsd)

pca_batch_condition(dds, vsd = vsd)

# Heatmap VST (figura suplementaria: top genes variables)
source(here::here("R", "07-heatmap.R"))
heatmap_top_n <- 50
heatmap_remove_batch <- TRUE
plot_vst_heatmap(vsd = vsd, top_n = heatmap_top_n,
                 out_file = here::here("figures", "VST_topvar_heatmap.svg"),
                 remove_batch = heatmap_remove_batch,
                 show_rownames = FALSE)


genes_of_interest <- c(
  dazl   = "ENSDARG00000036214",
  ddx4   = "ENSDARG00000014373",
  elavl2 = "ENSDARG00000040732")
# zar1l discarded

# Visualize gene of interest
plot_genes_boxplots(dds, genes = genes_of_interest)  # normalized counts
plot_genes_boxplots_tpm(tpm_mat, genes = genes_of_interest, sample_info = colData(dds))  # tpm



# dazl vs Control
res_dazl <- results(dds, contrast = c("condition", "dazl", "Control"), alpha = 0.05)
summary(res_dazl)
# ddx4 vs Control
res_ddx4 <- results(dds, contrast = c("condition", "ddx4", "Control"), alpha = 0.05)
summary(res_ddx4)
# elavl2 vs Control
res_elavl2 <- results(dds, contrast = c("condition", "elavl2", "Control"), alpha = 0.05)
summary(res_elavl2)
# zar1l vs Control
# res_zar1l <- results(dds, contrast = c("condition", "zar1l", "Control"))

# log2FC shrinkage
res_dazl_shrink   <- lfcShrink(dds, coef = "condition_dazl_vs_Control",   type = "apeglm")
summary(res_dazl_shrink)
res_ddx4_shrink   <- lfcShrink(dds, coef = "condition_ddx4_vs_Control",   type = "apeglm")
summary(res_ddx4_shrink)
res_elavl2_shrink <- lfcShrink(dds, coef = "condition_elavl2_vs_Control", type = "apeglm")
summary(res_elavl2_shrink)
# res_zar1l_shrink  <- lfcShrink(dds, coef = "condition_zar1l_vs_Control",  type = "apeglm")


# añadir la columna de LFC shrunk a las tablas originales
res_dazl$log2FC_shrunk   <- res_dazl_shrink$log2FoldChange
res_ddx4$log2FC_shrunk   <- res_ddx4_shrink$log2FoldChange
res_elavl2$log2FC_shrunk <- res_elavl2_shrink$log2FoldChange

# Heatmaps de DEGs
deg_heatmap_top_n <- 50
deg_heatmap_padj <- 0.05
deg_heatmap_remove_batch <- TRUE
deg_heatmap_control <- "Control"
plot_deg_heatmap(vsd = vsd, deseq_res = res_dazl, padj_th = deg_heatmap_padj,
                 top_n = deg_heatmap_top_n,
                 out_file = here::here("figures", "DEG_dazl_vs_control_heatmap.svg"),
                 remove_batch = deg_heatmap_remove_batch,
                 show_rownames = TRUE,
                 condition_value = "dazl", control_value = deg_heatmap_control)
plot_deg_heatmap(vsd = vsd, deseq_res = res_ddx4, padj_th = deg_heatmap_padj,
                 top_n = deg_heatmap_top_n,
                 out_file = here::here("figures", "DEG_ddx4_vs_control_heatmap.svg"),
                 remove_batch = deg_heatmap_remove_batch,
                 show_rownames = TRUE,
                 condition_value = "ddx4", control_value = deg_heatmap_control)
plot_deg_heatmap(vsd = vsd, deseq_res = res_elavl2, padj_th = deg_heatmap_padj,
                 top_n = deg_heatmap_top_n,
                 out_file = here::here("figures", "DEG_elavl2_vs_control_heatmap.svg"),
                 remove_batch = deg_heatmap_remove_batch,
                 show_rownames = TRUE,
                 condition_value = "elavl2", control_value = deg_heatmap_control)



# MAplots
svg(file=here::here("figures", "MAplots.svg"), width = 12, height = 9)
par(mfrow=c(3,2))
plotMA(res_dazl, main="dazl vs Control", colSig="steelblue", ylim=c(-3,3))
plotMA(res_dazl_shrink, main="dazl vs Control (lfcShrink)", colSig="blue3", ylim=c(-3,3))
plotMA(res_ddx4, main="ddx4 vs Control", colSig="steelblue", ylim=c(-3,3))
plotMA(res_ddx4_shrink, main="ddx4 vs Control (lfcShrink)", colSig="blue3", ylim=c(-3,3))
plotMA(res_elavl2, main="elavl2 vs Control", colSig="steelblue", ylim=c(-3,3))
plotMA(res_elavl2_shrink, main="elavl2 vs Control (lfcShrink)", colSig="blue3", ylim=c(-3,3))
# plotMA(res_zar1l, main="zar1l vs Control", colSig="steelblue", ylim=c(-3,3))
# plotMA(res_zar1l_shrink, main="zar1l vs Control (lfcShrink)", colSig="blue3", ylim=c(-3,3))
par(mfrow=c(1,1))
dev.off()

# histogram p-values
svg(file=here::here("figures", "pval_histograms.svg"), width = 12, height = 7)
par(mfrow=c(3,1))
hist(res_dazl$pvalue, main="dazl vs Control", xlab = "p-value")
hist(res_ddx4$pvalue, main="ddx4 vs Control",  xlab = "p-value")
hist(res_elavl2$pvalue, main="elavl2 vs Control", xlab = "p-value")
# hist(res_zar1l$pvalue, main="zar1l vs Control",  xlab = "p-value")
par(mfrow=c(1,1))
dev.off()

rm(res_dazl_shrink);gc()
rm(res_ddx4_shrink);gc()
rm(res_elavl2_shrink);gc()

#---------------------------------------
## Annotation ##
#---------------------------------------

columns(org.Dr.eg.db)

source(here::here("R", "04-annotation.R"))
res_dazl_annot <- annot(res_dazl)
res_ddx4_annot <- annot(res_ddx4)
res_elavl2_annot <- annot(res_elavl2)

write_deseq_table <- function(res_obj, out_file) {
  res_df <- as.data.frame(res_obj)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]
  write.table(res_df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

write_deseq_table(res_dazl_annot, here::here("results", "deseq2_dazl_vs_control.tsv"))
write_deseq_table(res_ddx4_annot, here::here("results", "deseq2_ddx4_vs_control.tsv"))
write_deseq_table(res_elavl2_annot, here::here("results", "deseq2_elavl2_vs_control.tsv"))

#---------------------------------------
## Volcano Plots ##
#---------------------------------------
source(here::here("R", "10-volcano.R"))
plot_volcano(res_dazl_annot,   "dazl_vs_control")
plot_volcano(res_ddx4_annot,   "ddx4_vs_control")
plot_volcano(res_elavl2_annot, "elavl2_vs_control")

#########################################
# Enriquecimiento funcional
#########################################

# Analisis de sobre-representacion (ORA)
source(here::here("R", "05-ORA.R"))

dazl_ora <- run_GO_ORA(deseq_res = res_dazl_annot, save_tables = T, 
                       universe = T, contrast_name = "dazl_vs_control")
ddx4_ora <- run_GO_ORA(deseq_res = res_ddx4_annot, save_tables = T, 
                       universe = T, contrast_name = "ddx4_vs_control")
elavl2_ora <- run_GO_ORA(deseq_res = res_elavl2_annot, save_tables = T, 
                       universe = T, contrast_name = "elavl2_vs_control")
# Analisis de rutas con GSEA
source(here::here("R", "06-GSEA.R"))
# Ajustes para GSEA (fgsea)
gsea_rank_col <- "stat"
gsea_eps <- 0
gsea_nPermSimple <- 10000
gsea_tie_break <- TRUE
gsea_simplify <- FALSE


gsea_dazl <- run_GO_GSEA(deseq_res = res_dazl, contrast_name = "dazl_vs_control",
                         eps = gsea_eps, nPermSimple = gsea_nPermSimple,
                         tie_break = gsea_tie_break, rank_col = gsea_rank_col,
                         simplify_res = gsea_simplify)
gsea_ddx4 <- run_GO_GSEA(deseq_res = res_ddx4, contrast_name = "ddx4_vs_control",
                         eps = gsea_eps, nPermSimple = gsea_nPermSimple,
                         tie_break = gsea_tie_break, rank_col = gsea_rank_col,
                         simplify_res = gsea_simplify)
gsea_elavl2 <- run_GO_GSEA(deseq_res = res_elavl2, contrast_name = "elavl2_vs_control",
                          eps = gsea_eps, nPermSimple = gsea_nPermSimple,
                          tie_break = gsea_tie_break, rank_col = gsea_rank_col,
                          simplify_res = gsea_simplify)

# Analisis de enriquecimiento en rutas (MSigDB)
source(here::here("R", "08-PEA.R"))
pea_lfc_col <- "log2FC_shrunk"

pea_dazl <- run_pea(deseq_res = res_dazl, contrast_name = "dazl_vs_control",
                    lfc_col = pea_lfc_col)
pea_ddx4 <- run_pea(deseq_res = res_ddx4, contrast_name = "ddx4_vs_control",
                    lfc_col = pea_lfc_col)
pea_elavl2 <- run_pea(deseq_res = res_elavl2, contrast_name = "elavl2_vs_control",
                     lfc_col = pea_lfc_col)

#---------------------------------------
## Visualizaciones avanzadas (enrichplot) ##
#---------------------------------------
source(here::here("R", "09-plots_enrich.R"))

# 1. ORA plots
lapply(names(dazl_ora),  function(ont) plot_ora_advanced(dazl_ora[[ont]],  "dazl_vs_control",  ont))
lapply(names(ddx4_ora),  function(ont) plot_ora_advanced(ddx4_ora[[ont]],  "ddx4_vs_control",  ont))
lapply(names(elavl2_ora), function(ont) plot_ora_advanced(elavl2_ora[[ont]], "elavl2_vs_control", ont))

# 2. GSEA plots
lapply(names(gsea_dazl),  function(ont) plot_gsea_advanced(gsea_dazl[[ont]],  "dazl_vs_control",  ont))
lapply(names(gsea_ddx4),  function(ont) plot_gsea_advanced(gsea_ddx4[[ont]],  "ddx4_vs_control",  ont))
lapply(names(gsea_elavl2), function(ont) plot_gsea_advanced(gsea_elavl2[[ont]], "elavl2_vs_control", ont))

# 3. PEA plots (custom MSigDB ORA)
lapply(names(pea_dazl),  function(id) plot_ora_advanced(pea_dazl[[id]],  "dazl_vs_control",  id))
lapply(names(pea_ddx4),  function(id) plot_ora_advanced(pea_ddx4[[id]],  "ddx4_vs_control",  id))
lapply(names(pea_elavl2), function(id) plot_ora_advanced(pea_elavl2[[id]], "elavl2_vs_control", id))

#---------------------------------------
## Análisis Comparativo (Venn / UpSet) ##
#---------------------------------------
source(here::here("R", "11-comparative.R"))

# Lista nombrada con los resultados anotados
res_list_all <- list(
  dazl   = res_dazl_annot,
  ddx4   = res_ddx4_annot,
  elavl2 = res_elavl2_annot
)

run_comparative_analysis(res_list_all)

#---------------------------------------
## Core UP genes annotated table ##
#---------------------------------------
core_file <- here::here("figures", "comparative", "Core_Germinal_Genes_UP.txt")
core_symbols <- character(0)
if (file.exists(core_file)) {
  core_symbols <- readLines(core_file)
}
core_symbols <- core_symbols[core_symbols != ""]
if (length(core_symbols) == 0) {
  up_lists <- get_gene_lists(res_list_all, direction = "UP", lfc_th = 1, padj_th = 0.05)
  if (length(up_lists) > 0) {
    core_symbols <- Reduce(intersect, up_lists)
  }
}

if (length(core_symbols) > 0) {
  build_symbol_map <- function(res_obj) {
    df <- as.data.frame(res_obj)
    df$gene_id <- rownames(df)
    if (!"symbol" %in% colnames(df)) {
      df$symbol <- df$gene_id
    }
    if (!"genename" %in% colnames(df)) {
      df$genename <- NA_character_
    }
    df <- df[, c("symbol", "gene_id", "genename"), drop = FALSE]
    df <- df[!duplicated(df$symbol), , drop = FALSE]
    df
  }

  build_stat_table <- function(res_obj, suffix) {
    df <- as.data.frame(res_obj)
    df$gene_id <- rownames(df)
    if (!"symbol" %in% colnames(df)) {
      df$symbol <- df$gene_id
    }
    df <- df[df$symbol %in% core_symbols, c("symbol", "log2FC_shrunk", "padj"), drop = FALSE]
    if (nrow(df) > 0) {
      df <- df[order(df$symbol, df$padj), , drop = FALSE]
      df <- df[!duplicated(df$symbol), , drop = FALSE]
    }
    colnames(df)[colnames(df) == "log2FC_shrunk"] <- paste0("log2FC_shrunk_", suffix)
    colnames(df)[colnames(df) == "padj"] <- paste0("padj_", suffix)
    df
  }

  symbol_map <- build_symbol_map(res_dazl_annot)
  symbol_map <- rbind(symbol_map, build_symbol_map(res_ddx4_annot))
  symbol_map <- rbind(symbol_map, build_symbol_map(res_elavl2_annot))
  symbol_map <- symbol_map[!duplicated(symbol_map$symbol), , drop = FALSE]

  core_tbl <- data.frame(symbol = core_symbols, stringsAsFactors = FALSE)
  core_tbl <- merge(core_tbl, symbol_map, by = "symbol", all.x = TRUE)

  core_tbl <- merge(core_tbl, build_stat_table(res_dazl_annot, "dazl"), by = "symbol", all.x = TRUE)
  core_tbl <- merge(core_tbl, build_stat_table(res_ddx4_annot, "ddx4"), by = "symbol", all.x = TRUE)
  core_tbl <- merge(core_tbl, build_stat_table(res_elavl2_annot, "elavl2"), by = "symbol", all.x = TRUE)

  core_tbl <- core_tbl[match(core_symbols, core_tbl$symbol), , drop = FALSE]
  core_tbl$mean_log2FC_shrunk <- rowMeans(core_tbl[, c("log2FC_shrunk_dazl",
                                                      "log2FC_shrunk_ddx4",
                                                      "log2FC_shrunk_elavl2")],
                                          na.rm = TRUE)

  out_dir <- here::here("results", "comparative")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.table(core_tbl, here::here(out_dir, "core_up_genes_annotated.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}
