# Utilidades comunes para el pipeline

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

list_input_files <- function(input_dir, pattern) {
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  files <- sort(files)
  if (length(files) == 0) {
    stop("No se encontraron archivos en ", input_dir, " con patron ", pattern)
  }
  files
}

parse_sample_info <- function(files) {
  base <- basename(files)
  base_no_ext <- sub("\\..*$", "", base)

  match <- regexec("^([^_]+)_([0-9]+)_CpG$", base_no_ext)
  parts <- regmatches(base_no_ext, match)
  ok <- lengths(parts) == 3
  if (!all(ok)) {
    stop("Nombres de muestra inesperados: ", paste(base_no_ext[!ok], collapse = ", "))
  }

  group <- vapply(parts, function(x) x[2], character(1))
  replicate <- vapply(parts, function(x) x[3], character(1))
  sample_id <- paste(group, replicate, sep = "_")

  df <- data.frame(
    sample_id = sample_id,
    group = group,
    replicate = replicate,
    file = files,
    stringsAsFactors = FALSE
  )

  if (anyDuplicated(df$sample_id)) {
    dupes <- unique(df$sample_id[duplicated(df$sample_id)])
    stop("Muestras duplicadas detectadas: ", paste(dupes, collapse = ", "))
  }

  df[order(df$sample_id), , drop = FALSE]
}

check_sample_consistency <- function(info_a, info_b) {
  if (is.null(info_a) || is.null(info_b)) {
    return(invisible(FALSE))
  }

  if (!setequal(info_a$sample_id, info_b$sample_id)) {
    missing_a <- setdiff(info_b$sample_id, info_a$sample_id)
    missing_b <- setdiff(info_a$sample_id, info_b$sample_id)
    msg <- c("Las muestras de methylKit y bedGraph no coinciden.",
             if (length(missing_a) > 0) paste("Faltan en methylKit:", paste(missing_a, collapse = ", ")),
             if (length(missing_b) > 0) paste("Faltan en bedGraph:", paste(missing_b, collapse = ", ")))
    stop(paste(msg, collapse = " "))
  }

  merged <- merge(info_a, info_b, by = "sample_id", suffixes = c(".a", ".b"))
  if (any(merged$group.a != merged$group.b)) {
    stop("Las etiquetas de grupo no coinciden entre fuentes.")
  }
  if (any(merged$replicate.a != merged$replicate.b)) {
    stop("Los replicados no coinciden entre fuentes.")
  }

  invisible(TRUE)
}

infer_standard_chromosomes <- function(chr_values, max_chr = 25) {
  chr_values <- unique(as.character(chr_values))
  has_prefix <- any(grepl("^chr", chr_values))
  if (has_prefix) {
    paste0("chr", seq_len(max_chr))
  } else {
    as.character(seq_len(max_chr))
  }
}
