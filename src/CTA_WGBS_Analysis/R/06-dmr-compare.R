# Comparacion de DMRs entre condiciones

read_dmrs_as_granges <- function(path) {
  dmrs <- read.table(path, header = TRUE, stringsAsFactors = FALSE)
  if (nrow(dmrs) == 0) {
    return(GenomicRanges::GRanges())
  }

  required <- c("chr", "start", "end")
  if (!all(required %in% colnames(dmrs))) {
    stop("La tabla de DMRs no tiene columnas chr/start/end: ", basename(path))
  }

  dmrs <- dmrs[complete.cases(dmrs[, required]), , drop = FALSE]

  GenomicRanges::GRanges(
    seqnames = as.character(dmrs$chr),
    ranges = IRanges::IRanges(start = dmrs$start, end = dmrs$end)
  )
}

compare_dmrs_across_conditions <- function(config) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE)) {
    warning("Paquetes faltantes para comparar DMRs: GenomicRanges, IRanges")
    return(NULL)
  }

  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(IRanges)
  })

  dss_dir <- config$dirs$dss_output
  dmr_files <- sort(list.files(dss_dir, pattern = "^DMR_results_.*\\.tsv$", full.names = TRUE))

  if (length(dmr_files) == 0) {
    message("No se encontraron ficheros de DMRs para comparar en ", dss_dir)
    return(NULL)
  }

  comp_names <- vapply(dmr_files, function(path) {
    name <- sub("^DMR_results_", "", basename(path))
    sub("\\.tsv$", "", name)
  }, character(1))
  names(dmr_files) <- comp_names

  gr_list <- lapply(dmr_files, read_dmrs_as_granges)
  names(gr_list) <- comp_names

  if (length(gr_list) < 2) {
    message("Se necesitan al menos 2 comparaciones para calcular la interseccion.")
    return(invisible(NULL))
  }

  if (length(gr_list) != 4) {
    message("Se encontraron ", length(gr_list), " comparaciones. Se usaran todas.")
  }

  dmr_counts <- vapply(gr_list, length, integer(1))
  empty_sets <- names(gr_list)[dmr_counts == 0]
  if (length(empty_sets) > 0) {
    warning("Comparaciones sin DMRs: ", paste(empty_sets, collapse = ", "),
            ". La interseccion sera vacia.")
  }

  common_gr <- gr_list[[1]]
  if (length(gr_list) > 1) {
    for (i in 2:length(gr_list)) {
      common_gr <- GenomicRanges::intersect(common_gr, gr_list[[i]])
      if (length(common_gr) == 0) {
        break
      }
    }
  }

  if (length(common_gr) > 0) {
    common_gr <- unique(common_gr)
    common_gr <- sort(common_gr)
  }

  shared_counts <- vapply(gr_list, function(gr) {
    if (length(common_gr) == 0) {
      0L
    } else {
      sum(GenomicRanges::countOverlaps(gr, common_gr) > 0L)
    }
  }, integer(1))

  out_dir <- file.path(config$dirs$results, "dmr_compare")
  ensure_dir(out_dir)

  summary_table <- data.frame(
    Comparison = names(gr_list),
    DMRs_total = dmr_counts,
    DMRs_shared_all = shared_counts,
    stringsAsFactors = FALSE
  )
  write.table(
    summary_table,
    file.path(out_dir, "summary_dmrs_common.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  common_df <- if (length(common_gr) == 0) {
    data.frame(
      chr = character(),
      start = integer(),
      end = integer(),
      width = integer(),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      chr = as.character(GenomicRanges::seqnames(common_gr)),
      start = GenomicRanges::start(common_gr),
      end = GenomicRanges::end(common_gr),
      width = GenomicRanges::width(common_gr),
      stringsAsFactors = FALSE
    )
  }
  write.table(
    common_df,
    file.path(out_dir, "dmrs_common_all.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  stats_table <- data.frame(
    n_comparisons = length(gr_list),
    n_common_regions = nrow(common_df),
    stringsAsFactors = FALSE
  )
  write.table(
    stats_table,
    file.path(out_dir, "summary_common_stats.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message("Comparacion de DMRs")
  message("Comparaciones: ", paste(names(gr_list), collapse = ", "))
  message("DMRs comunes en todas las condiciones: ", nrow(common_df))

  # --- Generacion de Graficos de Interseccion ---
  
  # Directorio para figuras de comparacion
  fig_dir <- file.path(config$dirs$figures, "dmr_compare")
  ensure_dir(fig_dir)
  
  # Crear matriz binaria de presencia/ausencia basada en regiones 'reduce' (consenso)
  # Unimos todos los GRanges y reducimos para obtener regiones unicas contiguas
  all_gr <- unlist(GRangesList(gr_list))
  consensus_gr <- GenomicRanges::reduce(all_gr)
  
  if (length(consensus_gr) > 0) {
    message("Generando matriz de solapamiento para ", length(consensus_gr), " regiones de consenso...")
    
    # Matriz logica
    mat <- matrix(0L, nrow = length(consensus_gr), ncol = length(gr_list))
    colnames(mat) <- names(gr_list)
    
    for (i in seq_along(gr_list)) {
      # findOverlaps devuelve indices de query (consensus) que solapan con subject (lista original)
      hits <- GenomicRanges::findOverlaps(consensus_gr, gr_list[[i]])
      mat[unique(S4Vectors::queryHits(hits)), i] <- 1L
    }
    
    mat_df <- as.data.frame(mat)
    
    # --- UpSet Plot ---
    if (requireNamespace("UpSetR", quietly = TRUE)) {
      message("Generando UpSet plot...")
      
      pdf(file.path(fig_dir, "UpSet_DMRs.pdf"), width = 10, height = 7, onefile = FALSE) 
      # Usamos print() porque UpSetR a veces no autoplotea en dispositivos no interactivos correctamente sin ello
      print(UpSetR::upset(
        mat_df, 
        sets = names(gr_list), 
        order.by = "freq", 
        main.bar.color = "steelblue",
        sets.bar.color = "tomato",
        text.scale = c(1.3, 1.3, 1, 1, 1.3, 1)
      ))
      dev.off()
      
      # Version SVG
      svg(file.path(fig_dir, "UpSet_DMRs.svg"), width = 10, height = 7)
      print(UpSetR::upset(
        mat_df, 
        sets = names(gr_list), 
        order.by = "freq",
        main.bar.color = "steelblue",
        sets.bar.color = "tomato"
      ))
      dev.off()
    } else {
      warning("Paquete 'UpSetR' no instalado. Se omite el grafico UpSet.")
    }
    
    # --- Venn Diagram ---
    # Venn es util para <= 4 conjuntos. Mas se vuelve ilegible.
    if (length(gr_list) <= 4 && length(gr_list) >= 2) {
      if (requireNamespace("ggVennDiagram", quietly = TRUE)) {
        message("Generando diagrama de Venn (ggVennDiagram)...")
        library(ggplot2)
        
        # ggVennDiagram espera una lista de vectores (IDs). 
        # Convertimos la matriz binaria a lista de indices de regiones de consenso.
        x_list <- apply(mat, 2, function(col) which(col == 1))
        
        p_venn <- ggVennDiagram::ggVennDiagram(x_list, label_alpha = 0) + 
          scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
          theme_void()
        
        ggsave(file.path(fig_dir, "Venn_DMRs.pdf"), p_venn, width = 8, height = 8)
        ggsave(file.path(fig_dir, "Venn_DMRs.svg"), p_venn, width = 8, height = 8)
        
      } else if (requireNamespace("VennDiagram", quietly = TRUE)) {
        message("Generando diagrama de Venn (VennDiagram)...")
        
        # VennDiagram genera logs por defecto, intentamos silenciarlo
        futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
        
        venn_plot <- VennDiagram::venn.diagram(
          x = apply(mat, 2, function(col) which(col == 1)),
          filename = NULL, # No guardar directo, queremos objeto grid
          category.names = names(gr_list),
          fill = RColorBrewer::brewer.pal(length(gr_list), "Pastel1"),
          output = TRUE
        )
        
        pdf(file.path(fig_dir, "Venn_DMRs.pdf"), width = 8, height = 8)
        grid::grid.draw(venn_plot)
        dev.off()
        
        svg(file.path(fig_dir, "Venn_DMRs.svg"), width = 8, height = 8)
        grid::grid.draw(venn_plot)
        dev.off()
      } else {
        warning("Ni 'ggVennDiagram' ni 'VennDiagram' instalados. Se omite el grafico Venn.")
      }
    } else if (length(gr_list) > 4) {
      message("Mas de 4 grupos, se omite diagrama de Venn (usar UpSet).")
    }
  }

  invisible(list(summary = summary_table, common = common_gr))
}
