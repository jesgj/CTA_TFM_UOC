# QC con methylKit

filter_methylkit_chromosomes <- function(myobj, config) {
  standard_chr <- config$methylkit$standard_chromosomes
  if (is.null(standard_chr)) {
    standard_chr <- infer_standard_chromosomes(
      myobj[[1]]$chr,
      max_chr = config$methylkit$standard_chr_max
    )
  }

  remove_regex <- config$methylkit$remove_chr_regex

  for (i in seq_along(myobj)) {
    chr <- as.character(myobj[[i]]$chr)
    keep <- chr %in% standard_chr
    if (!is.null(remove_regex) && nzchar(remove_regex)) {
      keep <- keep & !grepl(remove_regex, chr)
    }
    myobj[[i]] <- myobj[[i]][keep, ]
  }

  myobj
}

run_methylkit_qc <- function(files, sample_info, config) {
  suppressPackageStartupMessages({
    library(methylKit)
  })

  if (length(files) == 0) {
    stop("No hay archivos methylKit para QC.")
  }
  if (is.null(sample_info) || nrow(sample_info) == 0) {
    stop("No hay metadata de muestras para QC.")
  }

  output_fig <- config$dirs$qc_figures
  output_res <- config$dirs$qc_results
  ensure_dir(output_fig)
  ensure_dir(output_res)

  sample_ids <- sample_info$sample_id
  treatment <- rep(0L, length(sample_ids))

  message("Cargando methylKit para QC")
  myobj <- methRead(
    as.list(files),
    sample.id = as.list(sample_ids),
    assembly = config$methylkit$assembly,
    treatment = treatment,
    context = config$methylkit$context,
    mincov = config$methylkit$mincov
  )
  names(myobj) <- sample_ids

  myobj <- filter_methylkit_chromosomes(myobj, config)

  message("Estadisticas de metilacion y cobertura")
  pdf(file = file.path(output_fig, "MethStats.pdf"), width = 12, height = 9)
  par(mfrow = c(2, 2))
  for (i in seq_along(myobj)) {
    getMethylationStats(myobj[[i]], plot = TRUE, both.strands = FALSE)
  }
  par(mfrow = c(1, 1))
  dev.off()

  pdf(file = file.path(output_fig, "MethCoverage.pdf"), width = 12, height = 9)
  par(mfrow = c(2, 2))
  for (i in seq_along(myobj)) {
    getCoverageStats(myobj[[i]], plot = TRUE, both.strands = FALSE)
  }
  par(mfrow = c(1, 1))
  dev.off()

  filtered <- myobj
  if (!is.null(config$methylkit$filter_by_coverage)) {
    filtered <- filterByCoverage(
      myobj,
      lo.count = config$methylkit$filter_by_coverage$lo.count,
      lo.perc = NULL,
      hi.count = NULL,
      hi.perc = config$methylkit$filter_by_coverage$hi.perc
    )
  }

  if (isTRUE(config$methylkit$normalize)) {
    filtered <- normalizeCoverage(filtered, method = "median")
  }

  coverage_median <- sapply(filtered, function(x) median(x$coverage, na.rm = TRUE))
  coverage_df <- data.frame(
    sample_id = names(coverage_median),
    median_coverage = as.numeric(coverage_median),
    stringsAsFactors = FALSE
  )
  write.table(
    coverage_df,
    file = file.path(output_res, "median_coverage.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message("Uniendo CpGs para analisis exploratorio")
  meth <- unite(
    filtered,
    destrand = TRUE,
    min.per.group = config$methylkit$min_per_group,
    mc.cores = config$methylkit$cores
  )

  message("Correlacion y clustering")
  correlation_matrix <- getCorrelation(meth, plot = FALSE)
  write.table(
    correlation_matrix,
    file = file.path(output_res, "correlation_methylkit.tsv"),
    sep = "\t",
    quote = FALSE
  )

  write.table(
    sample_info,
    file = file.path(output_res, "sample_info.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  invisible(list(methyl_raw = myobj, methyl_base = meth))
}
