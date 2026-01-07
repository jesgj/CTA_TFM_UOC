# Analisis de DMRs con DSS por pares

run_pca_bsseq <- function(BSobj, config) {
  message("Iniciando PCA sobre BSseq object")
  suppressPackageStartupMessages({
    library(bsseq)
    library(ggplot2)
    library(RColorBrewer)
  })

  output_fig <- config$dirs$qc_figures
  ensure_dir(output_fig)

  # 1. Obtener matriz de metilación (beta values)
  # Usamos type="raw" para obtener los ratios de metilación directamente
  beta_matrix <- getMeth(BSobj, type = "raw")
  
  # 2. Filtrar sitios completos (sin NAs en ninguna muestra)
  complete_sites <- complete.cases(beta_matrix)
  beta_complete <- beta_matrix[complete_sites, ]
  
  message(sprintf("   Sitios totales: %d | Sitios completos para PCA: %d", 
                  nrow(beta_matrix), nrow(beta_complete)))

  if (nrow(beta_complete) < 100) {
    warning("Muy pocos sitios completos para un PCA fiable.")
    return(NULL)
  }

  # 3. Ejecutar PCA
  # t(beta_complete) porque prcomp espera muestras en filas
  pca_res <- prcomp(t(beta_complete), center = TRUE, scale. = FALSE)
  
  # 4. Calcular varianza explicada
  var_explained <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
  
  # 5. Preparar datos para ggplot
  pca_data <- as.data.frame(pca_res$x)
  pca_data$sample_id <- rownames(pca_data)
  
  # Añadir metadata desde pData(BSobj)
  metadata <- as.data.frame(pData(BSobj))
  # Asumimos que la columna de grupo se llama 'group' o 'tissue'
  # Intentamos 'tissue' primero (usado en run_dss_pairwise) o 'group'
  if ("tissue" %in% colnames(metadata)) {
    pca_data$group <- metadata$tissue
  } else if ("group" %in% colnames(metadata)) {
    pca_data$group <- metadata$group
  } else {
    pca_data$group <- "Unknown"
  }

  # 6. Generar Gráficos
  svg_width <- 6.3
  svg_height <- 5
  
  # PC1 vs PC2
  p12 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = sample_id)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text(vjust = -1, size = 3, show.legend=F) +
    labs(
      title = NULL,
      subtitle = sprintf("Basado en %s sitios CpG completos", format(nrow(beta_complete), big.mark=".")),
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2])
    ) +
    theme_minimal() +
    scale_color_brewer(palette = "Set1")

  ggsave(file.path(output_fig, "pca_bsseq_pc1_pc2.pdf"), p12, width = 10, height = 8)
  ggsave(
    file.path(output_fig, "pca_bsseq_pc1_pc2.svg"),
    p12,
    width = svg_width,
    height = svg_height,
    device = grDevices::svg
  )

  # PC1 vs PC3 (útil si PC2 no separa bien)
  if (ncol(pca_res$x) >= 3) {
    p13 <- ggplot(pca_data, aes(x = PC1, y = PC3, color = group, label = sample_id)) +
      geom_point(size = 4, alpha = 0.8) +
      geom_text(vjust = -1, size = 3) +
      labs(
        title = NULL,
        subtitle = sprintf("Basado en %s sitios CpG completos", format(nrow(beta_complete), big.mark=".")),
        x = sprintf("PC1 (%.1f%%)", var_explained[1]),
        y = sprintf("PC3 (%.1f%%)", var_explained[3])
      ) +
      theme_minimal() +
      scale_color_brewer(palette = "Set1")
    ggsave(file.path(output_fig, "pca_bsseq_pc1_pc3.pdf"), p13, width = 10, height = 8)
    ggsave(
      file.path(output_fig, "pca_bsseq_pc1_pc3.svg"),
      p13,
      width = svg_width,
      height = svg_height,
      device = grDevices::svg
    )
  }

  # Scree Plot
  pdf(file.path(output_fig, "pca_bsseq_screeplot.pdf"), width = 8, height = 6)
  barplot(
    var_explained[1:min(10, length(var_explained))],
    names.arg = paste0("PC", 1:min(10, length(var_explained))),
    main = "Varianza explicada por Componentes Principales",
    ylab = "% Varianza",
    col = "steelblue",
    las = 2
  )
  dev.off()
  
  svg(file.path(output_fig, "pca_bsseq_screeplot.svg"), width = svg_width, height = 4.5)
  barplot(
    var_explained[1:min(10, length(var_explained))],
    names.arg = paste0("PC", 1:min(10, length(var_explained))),
    main = "Varianza explicada por Componentes Principales",
    ylab = "% Varianza",
    col = "steelblue",
    las = 2
  )
  dev.off()

  message("PCA completado exitosamente")
  invisible(pca_res)
}

run_dss_pairwise <- function(BSobj, sample_info, config) {
  suppressPackageStartupMessages({
    library(DSS)
  })

  output_dir <- config$dirs$dss_output
  ensure_dir(output_dir)

  if (!"tissue" %in% colnames(pData(BSobj))) {
    idx <- match(sampleNames(BSobj), sample_info$sample_id)
    if (any(is.na(idx))) {
      stop("No se pudo mapear la metadata a BSobj.")
    }
    pData(BSobj)$tissue <- factor(as.character(sample_info$group[idx]))
  }
  if (!"individual" %in% colnames(pData(BSobj))) {
    idx <- match(sampleNames(BSobj), sample_info$sample_id)
    if (any(is.na(idx))) {
      stop("No se pudo mapear la metadata a BSobj.")
    }
    pData(BSobj)$individual <- factor(as.character(sample_info$replicate[idx]))
  }

  tissues <- sort(unique(as.character(pData(BSobj)$tissue)))
  if (length(tissues) < 2) {
    stop("Se necesitan al menos dos grupos para DSS.")
  }

  resolve_group_name <- function(name) {
    if (name %in% tissues) {
      return(name)
    }
    idx <- which(tolower(tissues) == tolower(name))
    if (length(idx) == 1) {
      return(tissues[idx])
    }
    stop("Grupo no encontrado: ", name)
  }

  comparison_mode <- if (!is.null(config$dss$comparison_mode)) {
    config$dss$comparison_mode
  } else {
    "all_pairs"
  }

  if (identical(comparison_mode, "control_vs_groups")) {
    if (is.null(config$dss$reference_group)) {
      stop("Falta configurar dss$reference_group.")
    }
    ref_group <- resolve_group_name(config$dss$reference_group)

    target_groups <- config$dss$target_groups
    if (is.null(target_groups) || length(target_groups) == 0) {
      target_groups <- setdiff(tissues, ref_group)
    }
    target_groups <- unique(vapply(target_groups, resolve_group_name, character(1)))
    if (length(target_groups) == 0) {
      stop("No hay grupos objetivos para comparar contra ", ref_group)
    }

    # Comparaciones control vs cada grupo objetivo (Tratamiento vs Control)
    # Ponemos ref_group segundo para que diff = Treat - Control (si DSS hace G1 - G2)
    pair_list <- lapply(target_groups, function(g) c(g, ref_group))
  } else {
    pair_list <- combn(tissues, 2, simplify = FALSE)
  }
  summary_table <- data.frame(
    Comparison = character(),
    DMRs_total = integer(),
    stringsAsFactors = FALSE
  )

  for (pair in pair_list) {
    t1 <- pair[1]
    t2 <- pair[2]
    message("Comparacion: ", t1, " vs ", t2)

    keep_samples <- rownames(pData(BSobj))[pData(BSobj)$tissue %in% c(t1, t2)]
    BS.sub <- BSobj[, keep_samples]
    pData(BS.sub)$tissue <- factor(as.character(pData(BS.sub)$tissue))

    samples1 <- rownames(pData(BS.sub))[pData(BS.sub)$tissue == t1]
    samples2 <- rownames(pData(BS.sub))[pData(BS.sub)$tissue == t2]

    if (length(samples1) < 2 || length(samples2) < 2) {
      warning("Hay grupos con menos de 2 muestras: ", t1, " (", length(samples1), "), ",
              t2, " (", length(samples2), ")")
    }

    dml_test <- DMLtest(
      BS.sub,
      group1 = samples1,
      group2 = samples2,
      smoothing = TRUE,
      smoothing.span = config$dss$smoothing_span,
      ncores = config$dss$ncores
    )

    dmrs <- callDMR(
      dml_test,
      delta = config$dss$delta,
      p.threshold = config$dss$p_threshold,
      minCG = config$dss$minCG,
      minlen = config$dss$minlen,
      dis.merge = config$dss$dis_merge
    )

    n_dmrs <- if (is.null(dmrs)) 0L else nrow(dmrs)
    summary_table <- rbind(
      summary_table,
      data.frame(
        Comparison = paste(t1, "vs", t2),
        DMRs_total = n_dmrs,
        stringsAsFactors = FALSE
      )
    )

    if (!is.null(dmrs) && nrow(dmrs) > 0) {
      out_dmrs <- file.path(output_dir, paste0("DMR_results_", t1, "_vs_", t2, ".tsv"))
      write.table(dmrs, file = out_dmrs, sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
      message("No se detectaron DMRs para ", t1, " vs ", t2)
    }
  }

  summary_path <- file.path(output_dir, "summary_dmrs.tsv")
  write.table(summary_table, summary_path, sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(summary_table)
}
