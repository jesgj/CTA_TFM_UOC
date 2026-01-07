# 07-heatmap.R

save_pheatmap_svg <- function(ph, out_file, width = 11.7, height = 8.3) {
  svglite::svglite(out_file, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  dev.off()
}

plot_vst_heatmap <- function(dds = NULL, vsd = NULL, top_n = 50, out_file = "VST_topvar_heatmap.svg",
                             blind = FALSE, scale_rows = TRUE, show_rownames = FALSE,
                             cluster_rows = TRUE, cluster_cols = TRUE,
                             remove_batch = FALSE, batch_col = "batch",
                             condition_col = "condition",
                             width = 11.7, height = 8.3) {
  if (is.null(vsd)) {
    if (is.null(dds)) {
      stop("dds o vsd son requeridos")
    }
    # Transformacion VST
    vsd <- DESeq2::vst(dds, blind = blind)
  }
  mat <- SummarizedExperiment::assay(vsd)

  # Correccion de batch opcional para visualizacion
  if (remove_batch && batch_col %in% colnames(SummarizedExperiment::colData(vsd))) {
    batch <- SummarizedExperiment::colData(vsd)[[batch_col]]
    if (condition_col %in% colnames(SummarizedExperiment::colData(vsd))) {
      design <- model.matrix(reformulate(condition_col),
                              data = as.data.frame(SummarizedExperiment::colData(vsd)))
      mat <- limma::removeBatchEffect(mat, batch = batch, design = design)
    } else {
      mat <- limma::removeBatchEffect(mat, batch = batch)
    }
  }

  # Seleccion de genes mas variables
  vars <- apply(mat, 1, var, na.rm = TRUE)
  vars <- vars[!is.na(vars)]
  if (length(vars) == 0) {
    stop("No hay varianzas validas para seleccionar genes")
  }
  top_n <- min(top_n, length(vars))
  top_genes <- names(sort(vars, decreasing = TRUE))[seq_len(top_n)]
  mat_top <- mat[top_genes, , drop = FALSE]

  # Escalado por filas (z-score)
  if (scale_rows) {
    row_sd <- apply(mat_top, 1, sd, na.rm = TRUE)
    mat_top <- mat_top[row_sd > 0, , drop = FALSE]
    if (nrow(mat_top) == 0) {
      stop("No hay genes con varianza > 0 para el heatmap")
    }
    mat_top <- t(scale(t(mat_top)))
  }

  # Anotaciones de columnas (muestras)
  annot <- as.data.frame(SummarizedExperiment::colData(vsd))
  if (!is.null(rownames(annot))) {
    annot <- annot[colnames(mat_top), , drop = FALSE]
  }

  # Paleta divergente
  heat_cols <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

  ph <- pheatmap::pheatmap(mat_top,
                           color = heat_cols,
                           annotation_col = annot,
                           show_rownames = show_rownames,
                           cluster_rows = cluster_rows,
                           cluster_cols = cluster_cols,
                           fontsize_col = 9,
                           silent = TRUE)
  save_pheatmap_svg(ph, out_file, width = width, height = height)
}

plot_deg_heatmap <- function(dds = NULL, deseq_res, vsd = NULL, padj_th = 0.05, top_n = 50,
                             lfc_col = "log2FC_shrunk", use_symbol = TRUE, symbol_col = "symbol", OrgDb = org.Dr.eg.db,
                             out_file = "DEG_heatmap.svg",
                             blind = FALSE, scale_rows = TRUE, show_rownames = FALSE,
                             cluster_rows = TRUE, cluster_cols = TRUE,
                             remove_batch = FALSE, batch_col = "batch",
                             condition_col = "condition",
                             condition_value = NULL, control_value = "Control",
                             width = 11.7, height = 8.3) {
  if (is.null(deseq_res)) {
    stop("deseq_res no puede ser NULL")
  }
  if (is.null(condition_value)) {
    stop("condition_value es obligatorio para comparar vs control")
  }
  if (!"padj" %in% colnames(deseq_res)) {
    stop("deseq_res no contiene columna padj")
  }

  # Filtrar DEGs por padj
  degs <- deseq_res[!is.na(deseq_res$padj) & deseq_res$padj < padj_th, , drop = FALSE]
  if (nrow(degs) == 0) {
    stop("No hay DEGs con el umbral de padj")
  }

  # Columna de ranking para ordenar DEGs
  if (!lfc_col %in% colnames(degs)) {
    if ("log2FoldChange" %in% colnames(degs)) {
      lfc_col <- "log2FoldChange"
    } else {
      stop("No se encontro una columna de LFC en deseq_res")
    }
  }
  degs <- degs[order(abs(degs[[lfc_col]]), decreasing = TRUE), , drop = FALSE]
  top_n <- min(top_n, nrow(degs))
  top_genes <- rownames(degs)[seq_len(top_n)]

  # Transformacion VST
  if (is.null(vsd)) {
    if (is.null(dds)) {
      stop("dds o vsd son requeridos")
    }
    vsd <- DESeq2::vst(dds, blind = blind)
  }

  # Subset de muestras: condicion vs control
  if (!condition_col %in% colnames(SummarizedExperiment::colData(vsd))) {
    stop("condition_col no existe en colData")
  }
  cond <- SummarizedExperiment::colData(vsd)[[condition_col]]
  if (!condition_value %in% cond) {
    stop("condition_value no existe en las muestras")
  }
  if (!control_value %in% cond) {
    stop("control_value no existe en las muestras")
  }
  keep <- cond %in% c(condition_value, control_value)
  vsd <- vsd[, keep]
  mat <- SummarizedExperiment::assay(vsd)

  # Correccion de batch opcional para visualizacion (solo en el subset)
  if (remove_batch && batch_col %in% colnames(SummarizedExperiment::colData(vsd))) {
    batch <- SummarizedExperiment::colData(vsd)[[batch_col]]
    if (condition_col %in% colnames(SummarizedExperiment::colData(vsd))) {
      design <- model.matrix(reformulate(condition_col),
                              data = as.data.frame(SummarizedExperiment::colData(vsd)))
      mat <- limma::removeBatchEffect(mat, batch = batch, design = design)
    } else {
      mat <- limma::removeBatchEffect(mat, batch = batch)
    }
  }

  # Subset de genes en matriz
  top_genes <- intersect(top_genes, rownames(mat))
  if (length(top_genes) == 0) {
    stop("Ningun gen DEG esta presente en la matriz VST")
  }
  mat_top <- mat[top_genes, , drop = FALSE]

  # Etiquetas de filas: usar SYMBOL si existe o mapear desde Ensembl
  if (use_symbol) {
    symbols <- NULL
    if (!is.null(symbol_col) && symbol_col %in% colnames(deseq_res)) {
      symbols <- deseq_res[rownames(mat_top), symbol_col]
    } else {
      suppressMessages({
        symbols <- AnnotationDbi::mapIds(OrgDb,
                                         keys = rownames(mat_top),
                                         column = "SYMBOL",
                                         keytype = "ENSEMBL",
                                         multiVals = "first")
      })
    }
    symbols <- as.character(symbols)
    missing <- is.na(symbols) | symbols == ""
    symbols[missing] <- rownames(mat_top)[missing]
    rownames(mat_top) <- make.unique(symbols)
  }

  # Escalado por filas (z-score)
  if (scale_rows) {
    row_sd <- apply(mat_top, 1, sd, na.rm = TRUE)
    mat_top <- mat_top[row_sd > 0, , drop = FALSE]
    if (nrow(mat_top) == 0) {
      stop("No hay genes con varianza > 0 para el heatmap")
    }
    mat_top <- t(scale(t(mat_top)))
  }

  # Anotaciones de columnas (muestras)
  annot <- as.data.frame(SummarizedExperiment::colData(vsd))
  if (!is.null(rownames(annot))) {
    annot <- annot[colnames(mat_top), , drop = FALSE]
  }

  # Paleta divergente
  heat_cols <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

  ph <- pheatmap::pheatmap(mat_top,
                           color = heat_cols,
                           annotation_col = annot,
                           show_rownames = show_rownames,
                           cluster_rows = cluster_rows,
                           cluster_cols = cluster_cols,
                           fontsize_col = 9,
                           silent = TRUE)
  save_pheatmap_svg(ph, out_file, width = width, height = height)
}

