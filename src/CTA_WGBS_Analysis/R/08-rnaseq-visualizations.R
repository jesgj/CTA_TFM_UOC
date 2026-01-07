# Visualizaciones para la integracion RNA-seq/DMR

plot_rnaseq_integration <- function(config) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Paquete faltante para graficas: ggplot2")
    return(NULL)
  }

  suppressPackageStartupMessages({
    library(ggplot2)
  })

  out_dir <- file.path(config$dirs$figures, "rnaseq_integration")
  ensure_dir(out_dir)

  summary_path <- file.path(config$rnaseq$output_dir, "summary_rnaseq_dmr_overlap.tsv")
  if (!file.exists(summary_path)) {
    message("No se encontro summary_rnaseq_dmr_overlap.tsv; omitiendo graficas.")
    return(NULL)
  }

  summary_df <- read.table(
    summary_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (nrow(summary_df) == 0) {
    message("Resumen de integracion vacio; omitiendo graficas.")
    return(NULL)
  }

  comp_col <- if ("Comparison_DMR" %in% colnames(summary_df)) {
    "Comparison_DMR"
  } else if ("Comparison_RNAseq" %in% colnames(summary_df)) {
    "Comparison_RNAseq"
  } else {
    colnames(summary_df)[1]
  }

  # Opcion 1: barras de conteos (DEGs, genes con DMR, solape)
  required_counts <- c("DEGs", "DMR_genes", "Overlap_genes")
  if (all(required_counts %in% colnames(summary_df))) {
    counts_df <- data.frame(
      comparison = rep(summary_df[[comp_col]], times = 3),
      metric = rep(required_counts, each = nrow(summary_df)),
      value = c(summary_df$DEGs, summary_df$DMR_genes, summary_df$Overlap_genes),
      stringsAsFactors = FALSE
    )

    metric_labels <- c(
      DEGs = "DEGs (RNA-seq)",
      DMR_genes = "Genes con DMR",
      Overlap_genes = "Solape DEGs+DMR"
    )
    counts_df$metric <- factor(
      counts_df$metric,
      levels = names(metric_labels),
      labels = unname(metric_labels)
    )

    p_counts <- ggplot(counts_df, aes(x = comparison, y = value, fill = metric)) +
      geom_col(position = "dodge", color = "black", linewidth = 0.3) +
      labs(
        x = "Comparacion",
        y = "Numero de genes",
        fill = "Categoria"
      ) +
      theme_minimal()

    ggsave(
      file.path(out_dir, "rnaseq_dmr_overlap_counts.pdf"),
      p_counts,
      width = 8,
      height = 5
    )
    ggsave(
      file.path(out_dir, "rnaseq_dmr_overlap_counts.svg"),
      p_counts,
      width = 6.3,
      height = 4.5,
      device = grDevices::svg
    )
  } else {
    warning("Faltan columnas en summary_rnaseq_dmr_overlap.tsv para graficas de conteo.")
  }

  overlap_files <- sort(list.files(
    config$rnaseq$output_dir,
    pattern = "^rnaseq_dmr_overlap_.*\\.tsv$",
    full.names = TRUE
  ))
  overlap_files <- overlap_files[!grepl("_all\\.tsv$", overlap_files)]

  invisible(TRUE)
}
