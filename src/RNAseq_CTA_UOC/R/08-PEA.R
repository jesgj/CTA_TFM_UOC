# 08-PEA.R

msigdb_collections <- function() {
  data.frame(
    collection = c("H", "REACTOME", "KEGG"),
    category = c("H", "C2", "C2"),
    subcategory = c(NA, "CP:REACTOME", "CP:KEGG_LEGACY"),
    stringsAsFactors = FALSE
  )
}

fetch_msigdb <- function(species, category, subcategory = NA) {
  if (is.na(subcategory) || is.null(subcategory) || subcategory == "") {
    msigdbr::msigdbr(species = species, collection = category)
  } else {
    msigdbr::msigdbr(species = species, collection = category, subcollection = subcategory)
  }
}

resolve_msigdb_ensembl_col <- function(msig_df) {
  candidates <- c("ensembl_gene", "ensembl_gene_id", "ensembl", "ensembl_id")
  available <- intersect(candidates, colnames(msig_df))
  if (length(available) == 0) {
    stop("No ENSEMBL column found in msigdbr output.")
  }
  available[[1]]
}

build_term2gene <- function(msig_df, genes_in_data) {
  id_col <- resolve_msigdb_ensembl_col(msig_df)
  term2gene <- data.frame(
    term = msig_df$gs_name,
    gene = as.character(msig_df[[id_col]]),
    stringsAsFactors = FALSE
  )
  term2gene <- term2gene[!is.na(term2gene$term) & !is.na(term2gene$gene), , drop = FALSE]
  term2gene <- term2gene[term2gene$gene %in% genes_in_data, , drop = FALSE]
  term2gene
}

split_deg_by_direction <- function(deseq_res, lfc_col, padj_col, id_col,
                                   log2fc_th, padj_th) {
  df <- as.data.frame(deseq_res)
  req_cols <- c(lfc_col, padj_col, id_col)
  if (!all(req_cols %in% colnames(df))) {
    stop("Missing required columns: ", paste(setdiff(req_cols, colnames(df)), collapse = ", "))
  }
  df <- df[, req_cols, drop = FALSE]
  df <- df[!is.na(df[[id_col]]) & !is.na(df[[lfc_col]]) & !is.na(df[[padj_col]]), , drop = FALSE]
  df$direction <- "NO"
  df$direction[df[[lfc_col]] > log2fc_th & df[[padj_col]] < padj_th] <- "UP"
  df$direction[df[[lfc_col]] < -log2fc_th & df[[padj_col]] < padj_th] <- "DOWN"
  df <- df[df$direction != "NO", , drop = FALSE]
  if (nrow(df) == 0) {
    return(list())
  }
  split(df, df$direction)
}

run_pea <- function(deseq_res,
                    contrast_name = "contrast",
                    species = "Danio rerio",
                    OrgDb = org.Dr.eg.db,
                    collections = msigdb_collections(),
                    log2fc_th = 0.5,
                    padj_th = 0.05,
                    lfc_col = "log2FoldChange",
                    padj_col = "padj",
                    ensembl_col = NULL,
                    padj_cutoff = 0.05,
                    gene_count_cutoff = 5,
                    save_background = TRUE,
                    out_dir = here::here("results", "pea"),
                    background_dir = here::here("results", "pea", "background")) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (save_background) {
    dir.create(background_dir, showWarnings = FALSE, recursive = TRUE)
  }

  deseq_df <- as.data.frame(deseq_res)
  if (!is.null(ensembl_col)) {
    if (!ensembl_col %in% colnames(deseq_df)) {
      stop("ensembl_col not found in deseq_res: ", ensembl_col)
    }
  } else {
    if (is.null(rownames(deseq_df))) {
      stop("No rownames found to derive ENSEMBL IDs.")
    }
    deseq_df$ensembl_id <- rownames(deseq_df)
    ensembl_col <- "ensembl_id"
  }
  deseq_df[[ensembl_col]] <- as.character(deseq_df[[ensembl_col]])
  deseq_df[[ensembl_col]][deseq_df[[ensembl_col]] == ""] <- NA
  genes_in_data <- unique(na.omit(deseq_df[[ensembl_col]]))
  if (length(genes_in_data) == 0) {
    stop("No valid ENSEMBL IDs found in deseq_res.")
  }

  term2gene_list <- list()
  for (i in seq_len(nrow(collections))) {
    collection <- collections$collection[i]
    category <- collections$category[i]
    subcategory <- collections$subcategory[i]
    msig_df <- fetch_msigdb(species = species, category = category, subcategory = subcategory)
    term2gene <- build_term2gene(msig_df, genes_in_data)
    term2gene_list[[collection]] <- term2gene
    if (save_background) {
      saveRDS(msig_df, file = here::here(background_dir, paste0(collection, ".RDS")))
    }
  }

  deg_list <- split_deg_by_direction(deseq_df, lfc_col, padj_col, ensembl_col, log2fc_th, padj_th)
  if (length(deg_list) == 0) {
    stop("No UP/DOWN genes after filtering.")
  }

  res_list <- list()
  for (collection in names(term2gene_list)) {
    term2gene <- term2gene_list[[collection]]
    if (nrow(term2gene) == 0) {
      next
    }
    for (direction in names(deg_list)) {
      gene_ids <- unique(deg_list[[direction]][[ensembl_col]])
      if (length(gene_ids) == 0) {
        next
      }
      enrich_res <- clusterProfiler::enricher(
        gene = gene_ids,
        TERM2GENE = term2gene,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        pAdjustMethod = "BH"
      )

      # Convertir IDs a Simbolos para visualizacion (cnetplot)
      try({
        enrich_res <- setReadable(enrich_res, OrgDb = OrgDb, keyType = "ENSEMBL")
      }, silent = TRUE)

      enrich_df <- as.data.frame(enrich_res)
      if (nrow(enrich_df) == 0) {
        next
      }
      enrich_df <- enrich_df[enrich_df$p.adjust < padj_cutoff &
                               enrich_df$Count > gene_count_cutoff, , drop = FALSE]
      if (nrow(enrich_df) == 0) {
        next
      }
      enrich_df$contrast <- contrast_name
      enrich_df$direction <- direction
      enrich_df$collection <- collection
      out_file <- here::here(out_dir, paste0(contrast_name, "_", collection, "_", direction, "_enrichment.csv"))
      write.csv(enrich_df, out_file, row.names = FALSE)
      
      # Retenemos el objeto de enriquecimiento filtrado por los mismos criterios
      res_list[[paste(collection, direction, sep = "_")]] <- enrich_res
    }
  }

  return(res_list)
}
