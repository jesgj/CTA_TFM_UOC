# 03-deseq2.R

build_deseq <- function(txi, coldata, condition, ref, batch=F) {
  # generamos el objeto DeSeq desde txi
  if (batch == T) {
    dds <- DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~ batch + condition)
  } else {
    dds <- DESeqDataSetFromTximport(txi = txi, colData = coldata, design = ~ condition)
  }
  dds$condition <- relevel(dds$condition, ref=ref)  # set the reference
  print(dds$condition)
  return(dds)
}

filter_dds <- function(deseq_obj, min_counts) {
  # Filtramos el objeto DESeq
  keep <- rowSums(counts(deseq_obj)) >= min_counts
  deseq_obj <- deseq_obj[keep,]
  print(table(keep))
  print(dim(deseq_obj))
  return(deseq_obj)
}

plot_SparDisp <- function(deseq_obj, pdf_w=12, pdf_h=7) {
  #  Sparsity and DispEsts 
  svg(file = here::here("figures", "dispersion_plot.svg"), width = pdf_w, height = pdf_h, )
  par(mfrow=c(1,2))
  DESeq2::plotSparsity(deseq_obj)  # raw dispersion
  DESeq2::plotDispEsts(deseq_obj)  # after correction
  par(mfrow=c(1,1))
  dev.off()
}

boxplot_counts <- function(deseq_obj, pdf_w=12, pdf_h=7) {
  #  boxplots raw and normalized counts
  log_raw_counts <- log2(counts(deseq_obj, normalized=F)+1)
  log_norm_counts <- log2(counts(deseq_obj, normalized=T)+1)
  svg(file = here::here("figures", "counts_boxplots.svg"), width = pdf_w, height = pdf_h, )
  par(mfrow=c(1,2))
  boxplot(log_raw_counts, notch=T, outline=F, las=2, col="grey50", main="log2(raw counts + 1)")
  boxplot(log_norm_counts, notch=T, outline=F, las=2, col="lightblue", main="log2(normalized counts + 1)")
  par(mfrow=c(1,1))
  dev.off()
}

vst_pca <- function(deseq_obj, pdf_w=12, pdf_h=7, intg="condition", check_batch=F, vsd=NULL) {
  # genera pca plot de la matriz vst
  if (is.null(vsd)) {
    vsd <- vst(deseq_obj, blind = FALSE)
  }
  if (check_batch == F) {
    vst_pca_pl <- DESeq2::plotPCA(vsd, intgroup=intg) + theme_bw()
  } else {
    pca_cond <- DESeq2::plotPCA(vsd, intgroup=intg) + theme_bw() + scale_color_brewer(palette = "Dark2")
    pca_batch <- DESeq2::plotPCA(vsd, intgroup='batch') + theme_bw() + scale_color_brewer(palette = "Set1")
    vst_pca_pl <- pca_cond + pca_batch + plot_annotation(tag_levels = "A")
  }
  ggsave(here::here("figures", "PCA.svg"), vst_pca_pl, dpi="print", width = pdf_w, height = pdf_h)
}
pca_batch_condition <- function(dds, width = 12, height = 7, vsd=NULL) {
  # VST con informacion de diseño
  if (is.null(vsd)) {
    vsd <- vst(dds, blind = FALSE)
  }
  # matriz VST
  mat <- assay(vsd)
  batch <- vsd$batch
  design <- model.matrix(~ vsd$condition)  # protege la señal de la condicion
  # correccion de efecto de batch en la matriz VST
  mat_corr <- limma::removeBatchEffect(mat, batch = batch, design = design)
  vsd_corr <- vsd
  assay(vsd_corr) <- mat_corr
  p_cond <- DESeq2::plotPCA(vsd_corr, intgroup = "condition") +
    theme_bw() +  scale_color_brewer(palette = "Dark2")
  
  # PCA coloreado por batch
  p_batch <- DESeq2::plotPCA(vsd_corr, intgroup = "batch") +
    theme_bw() + scale_color_brewer(palette = "Set1")
  
  # combinar en una sola figura A|B
  p_comb <- p_cond + p_batch + plot_annotation(tag_levels = "A")
  ggsave(here::here("figures", "PCA_batch_condition.svg"), p_comb, width = width, height = height, dpi = 300)
}
plot_genes_boxplots <- function(deseq_obj, genes, group_var="condition") {
  all_counts <- lapply(names(genes),
                       function(x) {
                         ens <- genes[[x]]
                         d <- plotCounts(deseq_obj, gene = ens, intgroup = group_var, returnData = T)
                         d$gene <- x
                         d}) %>% bind_rows()
  all_counts$gene <- factor(all_counts$gene, levels = names(genes))
  pl <- ggplot(all_counts, aes_string(x = group_var, y = "count")) +
    geom_boxplot(outlier.shape = NA, fill = "#F5EFFF") + geom_jitter(width = 0.1, size = 1.8) +
    facet_wrap(~ gene, scales = "free_y") +
    labs(x = "", y = "Normalized counts") + theme_bw()
  ggsave(filename = here::here("figures", "selected_genes_counts_boxplots.svg"), pl, width = 12, height = 7, dpi="print")
}

plot_genes_boxplots_tpm <- function(tpm_mat,genes,sample_info,group_var = "condition") {
  # plot log2(TPM + 1) de los genes seleccionados
  sample_info <- sample_info %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample")
  
  all_tpm <- lapply(names(genes), function(x) {
    ens <- genes[[x]]
    if (!ens %in% rownames(tpm_mat)) {
      stop(paste("El gen", ens, "no está en rownames(tpm_mat)"))
    }
    vals <- tpm_mat[ens, ]
    d <- data.frame(
      sample = colnames(tpm_mat),
      TPM    = log2(as.numeric(vals)+1),
      gene   = x,
      stringsAsFactors = FALSE
    )
    d
  }) %>% bind_rows()
  
  # añadimos la información de la condición u otro factor
  all_tpm <- all_tpm %>%
    left_join(sample_info, by = "sample")
  
  all_tpm$gene <- factor(all_tpm$gene, levels = names(genes))
  
  pl <- ggplot(all_tpm, aes(x = .data[[group_var]], y = TPM)) +
    geom_boxplot(outlier.shape = NA, fill = "lavender") +
    geom_jitter(width = 0.1, size = 1.8) +
    facet_wrap(~ gene, scales = "free_y") +
    labs(x = "", y = "log2(TPM + 1)") +
    theme_bw()
  
  ggsave(
    filename = here::here("figures", "selected_genes_tpm_boxplots.svg"),
    plot = pl,
    width = 12,
    height = 7,
    dpi = "print"
  )
}

