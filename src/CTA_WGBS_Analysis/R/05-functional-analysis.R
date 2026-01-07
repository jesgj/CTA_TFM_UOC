# Analisis Funcional y Visualizacion Avanzada (Basado en NGS101 Parte 3)

load_annotation_packages <- function(config) {
  pkgs <- c(
    "ChIPseeker",
    "clusterProfiler",
    "GenomeInfoDb",
    "GenomicRanges",
    "GenomicFeatures", # Necesario para makeTxDbFromGFF
    "txdbmaker",       # Nueva ubicacion de makeTxDbFromGFF
    "AnnotationDbi",   # Necesario para loadDb/saveDb
    "pheatmap",
    "ggplot2",
    "svglite",
    config$annotation$org_package
  )
  
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    warning("Paquetes faltantes para analisis funcional: ", paste(missing, collapse = ", "))
    return(FALSE)
  }
  
  suppressPackageStartupMessages({
    for (pkg in pkgs) {
      library(pkg, character.only = TRUE)
    }
  })
  return(TRUE)
}

get_txdb_from_gtf <- function(config) {
  gtf_path <- config$annotation$gtf_file
  sqlite_path <- config$annotation$txdb_sqlite
  
  if (file.exists(sqlite_path)) {
    message("Cargando TxDb desde cache: ", sqlite_path)
    return(AnnotationDbi::loadDb(sqlite_path))
  }
  
  if (!file.exists(gtf_path)) {
    stop("No se encuentra el archivo GTF: ", gtf_path)
  }
  
  message("Generando TxDb desde GTF (esto puede tardar unos minutos): ", gtf_path)
  txdb <- txdbmaker::makeTxDbFromGFF(
    file = gtf_path,
    format = "gtf",
    organism = "Danio rerio",
    taxonomyId = 7955 # TaxID para Zebrafish
  )
  
  message("Guardando TxDb compilado en: ", sqlite_path)
  AnnotationDbi::saveDb(txdb, file = sqlite_path)
  return(txdb)
}

normalizar_seqlevels <- function(gr, txdb) {
  # Igualar estilo de cromosomas entre DMRs y anotacion
  txdb_style <- GenomeInfoDb::seqlevelsStyle(txdb)
  if (length(txdb_style) > 0) {
    tryCatch({
      GenomeInfoDb::seqlevelsStyle(gr) <- txdb_style[1]
    }, error = function(e) {
      warning("No se pudo convertir el estilo de cromosomas: ", conditionMessage(e))
    })
  }

  comunes <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(txdb))
  if (length(comunes) == 0) {
    warning("No hay cromosomas en comun entre DMRs y anotacion.")
    # Intento fallback manual comun: quitar "chr" o poner "chr"
    style_txdb <- seqlevelsStyle(txdb)[1]
    if (style_txdb == "Ensembl") { # 1, 2, ...
       new_levels <- sub("^chr", "", seqlevels(gr))
       # Solo renombramos si los niveles coinciden
       if (any(new_levels %in% seqlevels(txdb))) {
          message("Aplicando correccion manual de 'chr' -> ''")
          seqlevels(gr) <- sub("^chr", "", seqlevels(gr))
       }
    } else if (style_txdb == "UCSC") { # chr1, chr2 ...
       if (!any(grepl("^chr", seqlevels(gr)))) {
          message("Aplicando correccion manual de '' -> 'chr'")
           seqlevels(gr) <- paste0("chr", seqlevels(gr))
       }
    }
    
    # Re-chequear interseccion
    comunes <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(txdb))
    if (length(comunes) == 0) {
        warning("CRITICO: Siguen sin coincidir cromosomas incluso tras intento de arreglo.")
        return(gr)
    }
  }

  GenomeInfoDb::keepSeqlevels(gr, comunes, pruning.mode = "coarse")
}

annotate_dmrs_workflow <- function(BSobj, config) {
  if (!load_annotation_packages(config)) {
    return(NULL)
  }

  has_bsobj <- !is.null(BSobj)
  if (!has_bsobj) {
    message("BSobj no disponible: se omite el heatmap de DMRs")
  }
  
  dss_dir <- config$dirs$dss_output
  dmr_files <- list.files(dss_dir, pattern = "^DMR_results_.*\\.tsv$", full.names = TRUE)
  
  if (length(dmr_files) == 0) {
    message("No se encontraron ficheros de DMRs para anotar.")
    return(NULL)
  }
  
  # Cargar TxDb desde GTF o Cache
  txdb <- get_txdb_from_gtf(config)
  org_db <- get(config$annotation$org_package)
  
  # Definir el universo de genes (background) para el análisis de enriquecimiento
  # Usamos todos los genes presentes en la anotación (TxDb) como fondo para WGBS
  message("Definiendo universo de genes para enriquecimiento (WGBS background)")
  universe_genes <- tryCatch({
    AnnotationDbi::keys(txdb, keytype = "GENEID")
  }, error = function(e) {
    warning("No se pudo extraer keys del TxDb, se usará el genoma completo por defecto.")
    NULL
  })
  
  if (!is.null(universe_genes)) {
    message("   Total genes en el universo (TxDb): ", length(universe_genes))
  }
  
  # Directorio para figuras de anotacion
  anno_fig_dir <- file.path(config$dirs$figures, "annotation")
  ensure_dir(anno_fig_dir)
  anno_res_dir <- file.path(config$dirs$results, "annotation")
  ensure_dir(anno_res_dir)
  
  for (fpath in dmr_files) {
    # Limpiar nombre de la comparacion de forma segura
    # Primero quitamos el prefijo, luego la extension
    filename <- basename(fpath)
    comp_name <- sub("^DMR_results_", "", filename)
    comp_name <- sub("\\.tsv$", "", comp_name)
    
    message("Anotando DMRs para: ", comp_name)
    
    dmrs <- read.table(fpath, header = TRUE, stringsAsFactors = FALSE)
    if (nrow(dmrs) == 0) next
    
    if (!all(c("chr", "start", "end") %in% colnames(dmrs))) {
      warning("La tabla no tiene columnas chr/start/end: ", basename(fpath))
      next
    }

    # Crear GRanges
    gr <- GRanges(
      seqnames = dmrs$chr,
      ranges = IRanges(start = dmrs$start, end = dmrs$end),
      diff_meth = dmrs$diff.Methy
    )

    gr <- normalizar_seqlevels(gr, txdb)
    
    # Anotar con ChIPseeker
    anno <- annotatePeak(
      gr,
      tssRegion = c(-config$annotation$promoter_upstream, config$annotation$promoter_downstream),
      TxDb = txdb,
      annoDb = config$annotation$org_package
    )
    
    # Guardar tabla anotada
    anno_df <- as.data.frame(anno)
    out_table <- file.path(anno_res_dir, paste0("Annotated_DMRs_", comp_name, ".tsv"))
    write.table(anno_df, out_table, sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Plots de anotacion
    pdf(file.path(anno_fig_dir, paste0("Annotation_Pie_", comp_name, ".pdf")))
    plotAnnoPie(anno)
    dev.off()
    
    # Ancho 6.3 in para A4 con margenes
    svg(file.path(anno_fig_dir, paste0("Annotation_Pie_", comp_name, ".svg")), width = 6.3, height = 5)
    plotAnnoPie(anno)
    dev.off()
    
    # Analisis de Enriquecimiento (GO y KEGG)
    # Pasamos el universo calculado para mejorar la precisión estadística
    run_enrichment(anno_df, comp_name, anno_res_dir, anno_fig_dir, config, org_db, universe = universe_genes)
    
    # Heatmap de top DMRs
    if (has_bsobj) {
      plot_dmr_heatmap(BSobj, dmrs, comp_name, anno_fig_dir, config)
    }
    
    # Volcano Plot
    plot_dmr_volcano(dmrs, comp_name, anno_fig_dir, config)
  }
}

run_enrichment <- function(anno_df, comp_name, res_dir, fig_dir, config, org_db, universe = NULL) {
  gene_id_type <- config$annotation$gene_id_type
  if (is.null(gene_id_type) || !nzchar(gene_id_type)) {
    gene_id_type <- "ENSEMBL"
  }

  genes_raw <- anno_df$geneId[!is.na(anno_df$geneId)]
  # A veces vienen multiples genes separados por ;
  genes_raw <- unlist(strsplit(as.character(genes_raw), "[/;]"))
  
  # Limpieza de IDs (quitar versiones tipo .1, .2)
  # Como TxDb viene de Ensembl, asumimos que genes_raw son ENSEMBL
  if (toupper(gene_id_type) == "ENSEMBL") {
    genes_raw <- sub("\\.\\d+$", "", genes_raw)
    if (!is.null(universe)) {
       universe <- sub("\\.\\d+$", "", universe)
    }
  }
  genes <- unique(trimws(genes_raw))

  suppress_bitr_warnings <- isTRUE(config$annotation$bitr_suppress_warnings)

  map_ids_with_summary <- function(ids, from_type, to_type, label) {
    if (is.null(ids)) {
      return(NULL)
    }

    ids <- ids[!is.na(ids)]
    ids <- ids[nzchar(ids)]
    ids <- unique(ids)
    if (length(ids) == 0) {
      return(NULL)
    }

    map_call <- function() {
      bitr(ids, fromType = from_type, toType = to_type, OrgDb = org_db)
    }

    label_suffix <- if (!is.null(label) && nzchar(label)) {
      paste0(" (", label, ")")
    } else {
      ""
    }

    mapped <- tryCatch(
      if (suppress_bitr_warnings) {
        suppressWarnings(map_call())
      } else {
        map_call()
      },
      error = function(e) {
        warning("Fallo mapeando IDs", label_suffix, ": ", e$message)
        NULL
      }
    )

    mapped_input <- 0L
    if (!is.null(mapped) && nrow(mapped) > 0 && from_type %in% colnames(mapped)) {
      mapped_input <- length(unique(mapped[[from_type]]))
    }

    pct <- (mapped_input / length(ids)) * 100
    msg_label <- if (!is.null(label) && nzchar(label)) label else paste(from_type, "->", to_type)
    message(sprintf(
      "  Mapeo %s: %s/%s (%.1f%%)",
      msg_label,
      format(mapped_input, big.mark = ","),
      format(length(ids), big.mark = ","),
      pct
    ))

    mapped
  }

  if (length(genes) < 10) {
    message("  Pocos genes (<10) para enriquecimiento en ", comp_name)
    return()
  }
  
  message("  Ejecutando enriquecimiento GO/KEGG con ", length(genes), " genes...")
  
  # --- Armonización de IDs para GO ---
  # El TxDb devuelve IDs ENSEMBL. Si config$gene_id_type es distinto (ej. ENTREZID),
  # debemos mapear tanto los genes como el universo al nuevo espacio para que enrichGO funcione.
  
  effective_genes_go <- genes
  effective_universe_go <- universe
  
  if (toupper(gene_id_type) != "ENSEMBL") {
      message("  [INFO] Mapeando genes y universo de ENSEMBL a ", gene_id_type, " para consistencia en GO...")
      
      # Mapear genes
      mapped_genes <- map_ids_with_summary(genes, "ENSEMBL", gene_id_type, "GO genes")
      
      if (!is.null(mapped_genes) && nrow(mapped_genes) > 0) {
          effective_genes_go <- unique(mapped_genes[[gene_id_type]])
      }
      
      # Mapear universo
      if (!is.null(universe)) {
          mapped_univ <- map_ids_with_summary(universe, "ENSEMBL", gene_id_type, "GO universo")
          
          if (!is.null(mapped_univ) && nrow(mapped_univ) > 0) {
              effective_universe_go <- unique(mapped_univ[[gene_id_type]])
          }
      }
  }

  # GO (Gene Ontology) - Biological Process (BP)
  go <- enrichGO(
    gene = effective_genes_go,
    universe = effective_universe_go,
    OrgDb = org_db,
    keyType = gene_id_type,
    ont = "BP", # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff = config$annotation$enrichment_pvalue,
    readable = TRUE
  )
  
  go_df <- as.data.frame(go)
  if (!is.null(go_df) && nrow(go_df) > 0) {
    write.table(go_df, file.path(res_dir, paste0("GO_BP_", comp_name, ".tsv")), sep = "\t", quote=F, row.names=F)
    
    p_go <- dotplot(go, showCategory=20, title = paste("GO:BP -", comp_name))
    ggsave(file.path(fig_dir, paste0("GO_Dotplot_", comp_name, ".pdf")), p_go, width=10, height=8)
    ggsave(file.path(fig_dir, paste0("GO_Dotplot_", comp_name, ".svg")), p_go, width=6.3, height=6.5)
  }
}

plot_dmr_heatmap <- function(BSobj, dmr_table, comp_name, fig_dir, config) {
  suppressPackageStartupMessages({
    library(bsseq)
    library(pheatmap)
  })
  
  top_n <- config$annotation$heatmap_top_n
  if (nrow(dmr_table) < 2) return()
  
  # Ordenar todas las DMRs por diferencia absoluta de metilacion
  all_dmrs <- dmr_table[order(abs(dmr_table$diff.Methy), decreasing = TRUE), ]
  
  # --- Filtrado de muestras segun comparacion ---
  # ... (resto del filtrado de muestras) ...
  pd <- pData(BSobj)
  group_col <- if ("tissue" %in% names(pd)) "tissue" else "group"
  
  parts <- strsplit(comp_name, "_vs_")[[1]]
  if (length(parts) == 2) {
    g1_raw <- parts[1]
    g2_raw <- parts[2]
    
    # Helper para buscar grupo case-insensitive
    find_group <- function(g, available) {
      if (g %in% available) return(g)
      idx <- which(tolower(available) == tolower(g))
      if (length(idx) > 0) return(available[idx[1]])
      return(NULL)
    }
    
    available_groups <- as.character(unique(pd[[group_col]]))
    g1 <- find_group(g1_raw, available_groups)
    g2 <- find_group(g2_raw, available_groups)
    
    if (!is.null(g1) && !is.null(g2)) {
      keep_mask <- pd[[group_col]] %in% c(g1, g2)
      samples_to_keep <- rownames(pd)[keep_mask]
      message(sprintf("  Heatmap: filtrando muestras para %s vs %s (%d muestras)", g1, g2, length(samples_to_keep)))
      
      if (length(samples_to_keep) > 0) {
         BSobj <- BSobj[, samples_to_keep]
      }
    }
  }
  
  # Optimizacion: Iterar hasta encontrar top_n DMRs validas (sin NAs)
  mat_list <- list()
  row_names <- c()
  
  message(sprintf("  Buscando top %d DMRs sin valores faltantes...", top_n))
  
  for (i in seq_len(nrow(all_dmrs))) {
    if (length(mat_list) >= top_n) break
    
    dmr <- all_dmrs[i,]
    gr <- GRanges(seqnames = dmr$chr, ranges = IRanges(dmr$start, dmr$end))
    
    bs_sub <- subsetByOverlaps(BSobj, gr)
    
    if (nrow(bs_sub) > 0) {
      betas <- getMeth(bs_sub, type = "raw")
      mean_vals <- colMeans(betas, na.rm = TRUE)
      
      # Solo aceptamos la region si no hay NAs en ninguna de las muestras seleccionadas
      if (all(!is.na(mean_vals))) {
        mat_list[[length(mat_list)+1]] <- mean_vals
        row_names <- c(row_names, paste0(dmr$chr, ":", dmr$start))
      }
    }
  }
  
  if (length(mat_list) < 2) {
    warning("No se encontraron suficientes DMRs sin NAs para el heatmap.")
    return()
  }
  
  mat_res <- do.call(rbind, mat_list)
  rownames(mat_res) <- row_names
  colnames(mat_res) <- colnames(BSobj)
  
  # Ya no es necesario el filtrado de NAs posterior porque lo hemos hecho en el bucle
  
  # Anotaciones para el heatmap (Grupos)
  pd_final <- pData(BSobj)
  annot_col <- data.frame(
    Group = pd_final[[group_col]],
    row.names = colnames(mat_res)
  )
  
  # Plot
  pdf(file.path(fig_dir, paste0("Heatmap_DMRs_", comp_name, ".pdf")), width=8, height=10)
  pheatmap(
    mat_res,
    annotation_col = annot_col,
    show_rownames = TRUE,
    show_colnames = TRUE,
    scale = "row", # Escalar por fila resalta diferencias relativas
    main = paste("Top", nrow(mat_res), "DMRs (by delta) -", comp_name),
    fontsize_row = 8
  )
  dev.off()
  
  svg(file.path(fig_dir, paste0("Heatmap_DMRs_", comp_name, ".svg")), width=6.3, height=8.5)
  pheatmap(
    mat_res,
    annotation_col = annot_col,
    show_rownames = TRUE,
    show_colnames = TRUE,
    scale = "row", # Escalar por fila resalta diferencias relativas
    main = paste("Top", nrow(mat_res), "DMRs (by delta) -", comp_name),
    fontsize_row = 8
  )
  dev.off()
}

plot_dmr_volcano <- function(dmr_table, comp_name, fig_dir, config) {
  # Asumimos que dmr_table viene de DSS: diff.Methy, areaStat, etc.
  
  df <- dmr_table
  if (!"areaStat" %in% names(df)) {
    # Si por alguna razon no esta areaStat, usamos diff como proxy visual
    df$areaStat <- abs(df$diff.Methy)
  }
  
  df$Direction <- ifelse(df$diff.Methy > 0, "Hyper", "Hypo")
  
  # Conteos para subtitulo
  n_hyper <- sum(df$Direction == "Hyper")
  n_hypo <- sum(df$Direction == "Hypo")
  
  # Umbral delta para lineas de referencia
  delta <- if (!is.null(config$dss$delta)) config$dss$delta else 0.1
  
  p <- ggplot(df, aes(x = diff.Methy, y = abs(areaStat), color = Direction)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(
      values = c("Hyper" = "#E74C3C", "Hypo" = "#3498DB"),
      name = "Methylation Change"
    ) +
    labs(
      subtitle = sprintf("Hyper: %s | Hypo: %s (Delta threshold: %s)",
                         format(n_hyper, big.mark = ","),
                         format(n_hypo, big.mark = ","),
                         delta),
      x = "Methylation Difference",
      y = "Abs(AreaStat)"
    ) +
    geom_vline(xintercept = c(-delta, delta), 
               linetype = "dashed", alpha = 0.7, color = "gray50") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )
  
  ggsave(file.path(fig_dir, paste0("Volcano_DMR_", comp_name, ".pdf")), p, width = 8, height = 6)
  ggsave(file.path(fig_dir, paste0("Volcano_DMR_", comp_name, ".svg")), p, width = 6.3, height = 5)
}
