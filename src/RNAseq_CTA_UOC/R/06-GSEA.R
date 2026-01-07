# 06-GSEA.R

run_GO_GSEA <- function(deseq_res,
                        rank_col = "stat",
                        OrgDb = org.Dr.eg.db,
                        keyType = "ENTREZID",
                        ont = c("BP", "MF"),
                        minGSSize = 10,
                        maxGSSize = 500,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        eps = 1e-10,
                        nPermSimple = NULL,
                        contrast_name = "gsea",
                        showCategory = 10,
                        save_tables = TRUE,
                        simplify_res = FALSE,
                        tie_break = FALSE,
                        seed = 1234) {
  # Opción A: Usar la columna 'stat' directamente (Muy recomendado)
  # rank_col = "stat"

  # Fijar semilla para reproducibilidad (fgsea usa permutaciones aleatorias)
  set.seed(seed)

  # Opción B: Calcular signed p-value manualmente con LFC Shrunken
  if (rank_col == "signed_pvalue") {
    sign_col <- "log2FC_shrunk"
    if (!sign_col %in% colnames(deseq_res)) {
      if ("log2FoldChange" %in% colnames(deseq_res)) {
        sign_col <- "log2FoldChange"
      } else {
        stop("signed_pvalue requiere log2FC_shrunk o log2FoldChange")
      }
    }
    # Aseguramos que no haya p-valores 0 para el log10
    deseq_res$pvalue[deseq_res$pvalue == 0] <- .Machine$double.xmin

    # Usamos LFC para determinar la dirección (signo)
    deseq_res$signed_pvalue <- -log10(deseq_res$pvalue) * sign(deseq_res[[sign_col]])
  }
  if (!rank_col %in% colnames(deseq_res)) {
    stop(paste("Columna no encontrada:", rank_col))
  }
  if (!"EntrezID" %in% colnames(deseq_res)) {
    suppressMessages({
      deseq_res$EntrezID <- mapIds(OrgDb,
                                  keys = rownames(deseq_res),
                                  column = "ENTREZID",
                                  keytype = "ENSEMBL",
                                  multiVals = "first")
    })
  }

  gene_df <- deseq_res[, c("EntrezID", rank_col)]
  gene_df <- gene_df[!is.na(gene_df$EntrezID) & !is.na(gene_df[[rank_col]]), , drop = FALSE]

  if (nrow(gene_df) == 0) {
    stop("No hay valores validos de EntrezID/ranking para GSEA")
  }

  # Colapsar IDs duplicados conservando el mayor cambio absoluto
  gene_list <- tapply(gene_df[[rank_col]], gene_df$EntrezID, function(x) x[which.max(abs(x))])
  gene_names <- names(gene_list)
  if (is.null(gene_names)) {
    gene_names <- dimnames(gene_list)[[1]]
  }
  gene_list <- as.numeric(gene_list)
  names(gene_list) <- gene_names
  # Si hay empates, opcionalmente se añade un jitter minimo para ordenar
  if (tie_break) {
    set.seed(1)
    gene_list <- gene_list + stats::rnorm(length(gene_list), mean = 0, sd = 1e-6)
  }
  gene_list <- sort(gene_list, decreasing = TRUE)

  # Validación del ranking
  # 1. Chequeo de infinitos (pueden romper fgsea)
  if (any(is.infinite(gene_list))) {
    message("Advertencia: Se han detectado valores Infinitos en el ranking. Se reemplazarán por los valores finitos extremos.")
    max_finite <- max(gene_list[is.finite(gene_list)], na.rm = TRUE)
    min_finite <- min(gene_list[is.finite(gene_list)], na.rm = TRUE)
    gene_list[gene_list == Inf]  <- max_finite * 1.01
    gene_list[gene_list == -Inf] <- min_finite * 1.01
  }

  # 2. Chequeo de cantidad minima
  if (length(gene_list) < 10) {
    stop("Error: El ranking contiene menos de 10 genes válidos. Insuficiente para GSEA.")
  }

  # 3. Chequeo de varianza cero
  if (stats::var(gene_list) == 0) {
    stop("Error: El ranking tiene varianza 0 (todos los valores son iguales). Imposible calcular enriquecimiento.")
  }

  # 4. Verificar que hay genes con valores positivos y negativos
  if (all(gene_list > 0) || all(gene_list < 0)) {
    warning("Atención: El ranking es unidireccional (todos los valores tienen el mismo signo). El algoritmo GSEA puede no ser óptimo en este escenario.")
  }

  res_list <- list()
  for (o in ont) {
    # Re-fijar semilla antes de cada ontología para consistencia interna si se desea
    set.seed(seed)
    
    # Construir argumentos para gseGO, pasando nPermSimple solo si es soportado
    gse_args <- list(
      geneList = gene_list,
      OrgDb = OrgDb,
      keyType = keyType,
      ont = o,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      eps = eps,
      verbose = FALSE,
      by = "fgsea"
    )
    if (!is.null(nPermSimple)) {
      gse_formals <- names(formals(clusterProfiler::gseGO))
      if (("..." %in% gse_formals) || ("nPermSimple" %in% gse_formals)) {
        gse_args$nPermSimple <- nPermSimple
      }
    }
    gsea_res <- do.call(gseGO, gse_args)

    # Convertir IDs a Simbolos para visualizacion (cnetplot)
    if (!is.null(gsea_res)) {
      try({
        gsea_res <- setReadable(gsea_res, OrgDb = OrgDb, keyType = keyType)
      }, silent = TRUE)
    }

    gsea_df <- as.data.frame(gsea_res)
    if (simplify_res && nrow(gsea_df) > 0) {
      gsea_res <- simplify(gsea_res)
      gsea_df <- as.data.frame(gsea_res)
    }
    if (nrow(gsea_df) > 0 && "p.adjust" %in% colnames(gsea_df)) {
      gsea_df$neg_log10_padj <- -log10(pmax(gsea_df$p.adjust, .Machine$double.xmin))
      
      # NOTA TÉCNICA SOBRE GENERATIO EN GSEA:
      # A diferencia de ORA (donde GeneRatio = genes_significativos_en_termino / total_genes_significativos),
      # en GSEA el "GeneRatio" se define habitualmente como la proporción de genes del "Core Enrichment"
      # (leading edge) respecto al tamaño total del conjunto (setSize).
      # Esto indica qué parte de la vía está realmente contribuyendo al enriquecimiento.
      
      core_counts <- sapply(strsplit(as.character(gsea_df$core_enrichment), "/"), length)
      gsea_df$GeneRatio <- core_counts / gsea_df$setSize
      
      gsea_res@result <- gsea_df
    }

    res_list[[o]] <- gsea_res

    if (nrow(gsea_df) > 0) {
      title <- paste0("GSEA ", o, " (", contrast_name, ")")
      gsea_dot <- enrichplot::dotplot(
        gsea_res,
        showCategory = showCategory,
        title = title,
        x = "GeneRatio",
        color = "NES",
        size = "neg_log10_padj"
      ) +
        ggplot2::scale_color_gradient2(
          low = "royalblue3",
          mid = "white",
          high = "firebrick3",
          midpoint = 0,
          name = "NES"
        ) +
        ggplot2::scale_size_continuous(
          name = "-log10(padj)",
          range = c(2, 6)
        ) +
        theme(plot.title = element_text(face = "bold",
                                        hjust = 0.5,
                                        size = 16),
              axis.text = element_text(face = "bold", size = 11),
              axis.title.x = element_text(face = "bold"),
              legend.text = element_text(face = "bold", size = 11),
              legend.title = element_text(face = "bold", size = 13))

      ggsave(filename = here::here("figures", paste0(contrast_name, "_GSEA_", o, "_dotplot.svg")),
             plot = gsea_dot,
             dpi = "print",
             width = 12,
             height = 9)

      top_id <- gsea_df$ID[1]
      top_desc <- gsea_df$Description[1]
      gsea_top <- enrichplot::gseaplot2(gsea_res, geneSetID = top_id, title = top_desc)
      ggsave(filename = here::here("figures", paste0(contrast_name, "_GSEA_", o, "_top_term.svg")),
             plot = gsea_top,
             dpi = "print",
             width = 9,
             height = 6)
             
      # -----------------------------------------------------------------------
      # Gráfico Multi-GSEA (Top 3 UP + Top 3 DOWN)
      # -----------------------------------------------------------------------
      # Ordenamos por NES
      gsea_sorted <- gsea_df[order(gsea_df$NES, decreasing = TRUE), ]
      
      # Top 3 Positivos (Mayor NES > 0)
      top_pos_ids <- head(gsea_sorted[gsea_sorted$NES > 0, "ID"], 3)
      
      # Top 3 Negativos (Menor NES < 0, es decir, más negativo)
      # Al estar ordenado decreciente, los negativos más fuertes están al final
      top_neg_ids <- tail(gsea_sorted[gsea_sorted$NES < 0, "ID"], 3)
      
      multi_ids <- c(top_pos_ids, top_neg_ids)
      
      if (length(multi_ids) > 0) {
        # Título descriptivo
        multi_title <- paste0("Top GSEA Pathways (", contrast_name, " - ", o, ")")
        
        # Generar plot con múltiples curvas
        # pvalue_table = FALSE para evitar saturar el gráfico con tablas grandes
        gsea_multi <- enrichplot::gseaplot2(gsea_res, 
                                            geneSetID = multi_ids, 
                                            title = multi_title,
                                            pvalue_table = FALSE,
                                            ES_geom = "line") # ES_geom="line" es más limpio para múltiples curvas
        
        ggsave(filename = here::here("figures", paste0(contrast_name, "_GSEA_", o, "_multi_top.svg")),
               plot = gsea_multi,
               dpi = "print",
               width = 12,
               height = 8)
      }
      # -----------------------------------------------------------------------
    }

    if (save_tables) {
      gsea_df_out <- gsea_df
      if ("neg_log10_padj" %in% colnames(gsea_df_out)) {
        gsea_df_out$neg_log10_padj <- NULL
      }
      write.table(x = gsea_df_out,
                  file = here::here("results", paste0(contrast_name, "_GSEA_", o, ".txt")),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE)
    }
  }

  return(res_list)
}
