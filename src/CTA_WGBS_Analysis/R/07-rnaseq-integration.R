# Integracion simple de RNA-seq (DESeq2) con DMRs anotadas

strip_version_ids <- function(ids) {
  sub("\\.\\d+$", "", ids)
}

format_numeric_columns <- function(df, digits = 4) {
  numeric_cols <- vapply(df, is.numeric, logical(1))
  if (!any(numeric_cols)) {
    return(df)
  }

  df[numeric_cols] <- lapply(df[numeric_cols], function(col) {
    formatted <- sprintf(paste0("%.", digits, "f"), round(col, digits))
    formatted[is.na(col)] <- NA_character_
    formatted
  })

  df
}

resolve_dmr_comparison <- function(comp_rna, dmr_map, config) {
  comp_map <- config$rnaseq$comparison_map
  if (!is.null(comp_map) && comp_rna %in% names(comp_map)) {
    candidate <- comp_map[[comp_rna]]
  } else {
    candidate <- sub("_vs_control$", "_vs_Ctr", comp_rna, ignore.case = TRUE)
  }

  if (tolower(candidate) %in% names(dmr_map)) {
    return(candidate)
  }
  if (tolower(comp_rna) %in% names(dmr_map)) {
    return(comp_rna)
  }

  NULL
}

read_deseq2_results <- function(path, config) {
  df <- read.table(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  gene_col <- config$rnaseq$gene_id_column
  if (!is.null(gene_col) && nzchar(gene_col) && !gene_col %in% colnames(df)) {
    warning("gene_id_column no encontrado (", gene_col, "). Se intentara inferir.")
    gene_col <- NULL
  }

  if (is.null(gene_col) || !nzchar(gene_col)) {
    candidates <- c("gene_id", "row.names", "X", "gene", "Gene", "GENE", "")
    hit <- candidates[candidates %in% colnames(df)]
    if (length(hit) > 0) {
      gene_col <- hit[1]
    }
  }

  if (is.null(gene_col) || !nzchar(gene_col)) {
    first_col <- colnames(df)[1]
    first_vals <- as.character(df[[1]])
    ensembl_ratio <- mean(grepl("^ENS", first_vals))
    if (!is.na(ensembl_ratio) && ensembl_ratio > 0.5) {
      gene_col <- first_col
    }
  }

  if (is.null(gene_col) || !nzchar(gene_col)) {
    stop("No se pudo identificar columna de genes en ", basename(path),
         ". Ajusta rnaseq$gene_id_column.")
  }

  if (!"padj" %in% colnames(df)) {
    stop("Falta columna padj en ", basename(path))
  }

  lfc_col <- config$rnaseq$lfc_column
  if (is.null(lfc_col) || !lfc_col %in% colnames(df)) {
    if ("log2FC_shrunk" %in% colnames(df)) {
      lfc_col <- "log2FC_shrunk"
    } else if ("log2FoldChange" %in% colnames(df)) {
      lfc_col <- "log2FoldChange"
    } else {
      stop("No se encontro columna de log2FC en ", basename(path))
    }
  }

  df$gene_id <- strip_version_ids(df[[gene_col]])
  df$log2fc <- df[[lfc_col]]

  padj_threshold <- config$rnaseq$padj_threshold
  if (is.null(padj_threshold)) {
    padj_threshold <- 0.05
  }

  keep <- !is.na(df$padj) & df$padj <= padj_threshold
  lfc_threshold <- config$rnaseq$lfc_threshold
  if (!is.null(lfc_threshold) && lfc_threshold > 0) {
    keep <- keep & !is.na(df$log2fc) & abs(df$log2fc) >= lfc_threshold
  }

  df <- df[keep, , drop = FALSE]
  df <- df[!is.na(df$gene_id) & nzchar(df$gene_id), , drop = FALSE]

  symbol_col <- if ("symbol" %in% colnames(df)) {
    "symbol"
  } else if ("SYMBOL" %in% colnames(df)) {
    "SYMBOL"
  } else {
    NULL
  }

  out <- data.frame(
    gene_id = df$gene_id,
    log2fc = df$log2fc,
    padj = df$padj,
    stringsAsFactors = FALSE
  )

  if (!is.null(symbol_col)) {
    out$symbol <- df[[symbol_col]]
    out <- out[, c("gene_id", "symbol", "log2fc", "padj")]
  }

  out
}

read_annotated_dmrs <- function(path) {
  dmrs <- read.table(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE,
    quote = "",
    comment.char = ""
  )

  if (nrow(dmrs) == 0) {
    return(list(by_gene = data.frame(), dmr_gene_count = 0L))
  }

  if (!"geneId" %in% colnames(dmrs)) {
    stop("No se encontro columna geneId en ", basename(path))
  }

  diff_col <- if ("diff_meth" %in% colnames(dmrs)) {
    "diff_meth"
  } else if ("diff.Methy" %in% colnames(dmrs)) {
    "diff.Methy"
  } else {
    NULL
  }

  dmrs$geneId <- strip_version_ids(dmrs$geneId)
  dmrs <- dmrs[!is.na(dmrs$geneId) & nzchar(dmrs$geneId), , drop = FALSE]

  if (nrow(dmrs) == 0) {
    return(list(by_gene = data.frame(), dmr_gene_count = 0L))
  }

  diff_values <- if (!is.null(diff_col)) dmrs[[diff_col]] else rep(NA_real_, nrow(dmrs))
  split_diff <- split(diff_values, dmrs$geneId)

  dmr_summary <- data.frame(
    gene_id = names(split_diff),
    dmr_count = vapply(split_diff, length, integer(1)),
    dmr_mean_diff_meth = vapply(split_diff, function(x) mean(x, na.rm = TRUE), numeric(1)),
    stringsAsFactors = FALSE
  )

  dmr_summary$dmr_direction <- ifelse(
    is.na(dmr_summary$dmr_mean_diff_meth),
    NA_character_,
    ifelse(
      dmr_summary$dmr_mean_diff_meth > 0,
      "hyper",
      ifelse(dmr_summary$dmr_mean_diff_meth < 0, "hypo", "neutral")
    )
  )

  list(by_gene = dmr_summary, dmr_gene_count = nrow(dmr_summary))
}

integrate_rnaseq_dmrs <- function(config) {
  rnaseq_dir <- config$rnaseq$input_dir
  dmr_dir <- config$rnaseq$annotation_dir
  out_dir <- config$rnaseq$output_dir

  if (!dir.exists(rnaseq_dir)) {
    warning("No existe el directorio RNA-seq: ", rnaseq_dir)
    return(NULL)
  }
  if (!dir.exists(dmr_dir)) {
    warning("No existe el directorio de anotaciones DMR: ", dmr_dir)
    return(NULL)
  }

  rnaseq_files <- sort(list.files(rnaseq_dir, pattern = "\\.tsv$", full.names = TRUE))
  dmr_files <- sort(list.files(dmr_dir, pattern = "^Annotated_DMRs_.*\\.tsv$", full.names = TRUE))

  if (length(rnaseq_files) == 0) {
    message("No se encontraron ficheros DESeq2 en ", rnaseq_dir)
    return(NULL)
  }
  if (length(dmr_files) == 0) {
    message("No se encontraron DMRs anotadas en ", dmr_dir)
    return(NULL)
  }

  dmr_names <- vapply(dmr_files, function(path) {
    name <- sub("^Annotated_DMRs_", "", basename(path))
    sub("\\.tsv$", "", name)
  }, character(1))
  dmr_map <- setNames(dmr_files, tolower(dmr_names))

  ensure_dir(out_dir)

  summary_list <- list()
  overlap_list <- list()

  for (path in rnaseq_files) {
    comp_rna <- sub("^deseq2_", "", basename(path))
    comp_rna <- sub("\\.tsv$", "", comp_rna)

    comp_dmr <- resolve_dmr_comparison(comp_rna, dmr_map, config)
    if (is.null(comp_dmr)) {
      warning(
        "No se encontro anotacion DMR para ",
        comp_rna,
        ". Disponibles: ",
        paste(dmr_names, collapse = ", ")
      )
      next
    }

    dmr_path <- dmr_map[[tolower(comp_dmr)]]

    deg_df <- read_deseq2_results(path, config)
    dmr_info <- read_annotated_dmrs(dmr_path)
    dmr_summary <- dmr_info$by_gene

    overlap <- merge(deg_df, dmr_summary, by = "gene_id")
    overlap$comparison <- comp_dmr

    desired <- c(
      "comparison",
      "gene_id",
      "symbol",
      "log2fc",
      "padj",
      "dmr_count",
      "dmr_mean_diff_meth",
      "dmr_direction"
    )
    overlap <- overlap[, intersect(desired, colnames(overlap)), drop = FALSE]

    if (nrow(overlap) > 0 && "padj" %in% colnames(overlap) && "log2fc" %in% colnames(overlap)) {
      overlap <- overlap[order(overlap$padj, -abs(overlap$log2fc)), , drop = FALSE]
    }

    out_file <- file.path(out_dir, paste0("rnaseq_dmr_overlap_", comp_dmr, ".tsv"))
    write.table(overlap, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

    summary_list[[comp_rna]] <- data.frame(
      Comparison_RNAseq = comp_rna,
      Comparison_DMR = comp_dmr,
      DEGs = nrow(deg_df),
      DMR_genes = nrow(dmr_summary),
      Overlap_genes = nrow(overlap),
      stringsAsFactors = FALSE
    )

    if (nrow(overlap) > 0) {
      overlap_list[[comp_rna]] <- overlap
    }
  }

  if (length(summary_list) == 0) {
    message("No se generaron comparaciones RNA-seq/DMR.")
    return(NULL)
  }

  summary_table <- do.call(rbind, summary_list)
  write.table(
    summary_table,
    file.path(out_dir, "summary_rnaseq_dmr_overlap.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  if (length(overlap_list) > 0) {
    combined <- do.call(rbind, overlap_list)
    write.table(
      combined,
      file.path(out_dir, "rnaseq_dmr_overlap_all.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    combined_rounded <- format_numeric_columns(combined, digits = 4)
    write.table(
      combined_rounded,
      file.path(out_dir, "rnaseq_dmr_overlap_all_rounded4.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  message("Integracion RNA-seq/DMR completada")
  invisible(summary_table)
}
