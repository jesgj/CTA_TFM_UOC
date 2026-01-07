# Construccion de BSseq a partir de bedGraph

read_bedgraph_as_dss <- function(path, remove_regex = NULL) {
  suppressPackageStartupMessages({
    library(data.table)
  })

  first_line <- readLines(path, n = 1, warn = FALSE)
  skip_lines <- if (length(first_line) > 0 && grepl("^track", first_line)) 1 else 0

  dt <- fread(path, skip = skip_lines, header = FALSE, data.table = TRUE)
  setnames(dt, c("chr", "start", "end", "percent", "meth", "unmeth"))
  dt[, chr := as.character(chr)]

  if (!is.null(remove_regex) && nzchar(remove_regex)) {
    dt <- dt[!grepl(remove_regex, chr)]
  }

  dss_df <- dt[, .(
    chr = chr,
    pos = start + 1,
    N = meth + unmeth,
    X = meth
  )]
  dss_df <- dss_df[N > 0]

  dss_df
}

filter_to_common_sites <- function(df_list) {
  get_site_id <- function(df) paste(df$chr, df$pos, sep = ":")
  common_sites <- Reduce(intersect, lapply(df_list, get_site_id))
  lapply(df_list, function(df) {
    site_id <- get_site_id(df)
    df[site_id %in% common_sites, ]
  })
}

build_bsseq_object <- function(files, sample_info, config) {
  suppressPackageStartupMessages({
    library(DSS)	  
    library(bsseq)
  })

  bsobj_path <- config$dss$bsobj_path
  if (file.exists(bsobj_path) && !isTRUE(config$run$rebuild_bsobj)) {
    message("Cargando BSobj existente")
    load(bsobj_path)
    if (!exists("BSobj")) {
      stop("No se encontro el objeto BSobj en ", bsobj_path)
    }
    existing_samples <- colnames(BSobj)
    current_samples <- sample_info$sample_id
    if (!setequal(existing_samples, current_samples)) {
      missing_current <- setdiff(current_samples, existing_samples)
      missing_existing <- setdiff(existing_samples, current_samples)
      msg <- c("El BSobj existente no coincide con las muestras actuales.",
               "Activa run$rebuild_bsobj = TRUE para regenerarlo.",
               if (length(missing_current) > 0) paste("Faltan en BSobj:", paste(missing_current, collapse = ", ")),
               if (length(missing_existing) > 0) paste("Sobrantes en BSobj:", paste(missing_existing, collapse = ", ")))
      stop(paste(msg, collapse = " "))
    }
    if (!identical(existing_samples, current_samples)) {
      order_idx <- match(current_samples, existing_samples)
      if (any(is.na(order_idx))) {
        stop("No se pudo reordenar BSobj para coincidir con las muestras actuales.")
      }
      BSobj <- BSobj[, order_idx]
      message("BSobj reordenado para coincidir con el orden de muestras actual.")
    }
    return(BSobj)
  }

  if (length(files) == 0) {
    stop("No hay archivos bedGraph para construir BSobj.")
  }

  message("Leyendo bedGraph y creando BSobj")
  df_list <- lapply(
    sample_info$file,
    read_bedgraph_as_dss,
    remove_regex = config$dss$remove_chr_regex
  )

  standard_chr <- config$dss$standard_chromosomes
  if (is.null(standard_chr)) {
    standard_chr <- infer_standard_chromosomes(
      df_list[[1]]$chr,
      max_chr = config$dss$standard_chr_max
    )
  }
  df_list <- lapply(df_list, function(df) df[df$chr %in% standard_chr, ])

  if (isTRUE(config$dss$require_common_cpgs)) {
    df_list <- filter_to_common_sites(df_list)
  }

  BSobj <- makeBSseqData(df_list, sample_info$sample_id)

  meta <- data.frame(
    tissue = factor(sample_info$group),
    individual = factor(sample_info$replicate),
    row.names = sample_info$sample_id,
    stringsAsFactors = FALSE
  )
  pData(BSobj) <- meta

  save(BSobj, file = bsobj_path)
  message("BSobj guardado en ", bsobj_path)

  BSobj
}
