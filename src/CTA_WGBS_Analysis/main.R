# Punto de entrada principal

options(stringsAsFactors = FALSE)

source("R/00-config.R")
source("R/01-utils.R")
source("R/02-methylkit-qc.R")
source("R/03-bsseq-build.R")
source("R/04-dss-analysis.R")
source("R/05-functional-analysis.R")
source("R/06-dmr-compare.R")
source("R/07-rnaseq-integration.R")
source("R/08-rnaseq-visualizations.R")
source("R/09-pca-visualization.R")

config <- get_config()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  args_lower <- tolower(args)
  valid_args <- c("methylkit", "dss", "annotation", "annot", "compare", "rnaseq", "integration", "pca", "all")
  unknown_args <- setdiff(args_lower, valid_args)
  run_all <- "all" %in% args_lower
  run_methylkit <- "methylkit" %in% args_lower
  run_dss <- "dss" %in% args_lower
  run_annotation <- any(args_lower %in% c("annotation", "annot"))
  run_compare <- "compare" %in% args_lower
  run_rnaseq <- any(args_lower %in% c("rnaseq", "integration"))
  run_pca <- "pca" %in% args_lower
  
  if (length(unknown_args) > 0) {
    stop(
      "Argumentos no reconocidos: ", paste(unknown_args, collapse = ", "),
      "\nUso: Rscript main.R [methylkit|dss|annotation|compare|rnaseq|pca|all]"
    )
  }

  if (run_all) {
    config$run$methylkit_qc <- TRUE
    config$run$dss <- TRUE
    config$run$annotation <- TRUE
    config$run$compare <- TRUE
    config$run$rnaseq <- TRUE
    # PCA se ejecuta normalmente dentro de dss, pero lo activamos explícitamente si se pide
    config$run$pca <- TRUE
    message("Modo de ejecucion: ALL (methylkit + dss + annotation + compare + rnaseq + pca)")
  } else {
    config$run$methylkit_qc <- run_methylkit
    config$run$dss <- run_dss
    config$run$annotation <- run_annotation
    config$run$compare <- run_compare
    config$run$rnaseq <- run_rnaseq
    config$run$pca <- run_pca
    
    selected <- c()
    if (run_methylkit) selected <- c(selected, "methylkit")
    if (run_dss) selected <- c(selected, "dss")
    if (run_annotation) selected <- c(selected, "annotation")
    if (run_compare) selected <- c(selected, "compare")
    if (run_rnaseq) selected <- c(selected, "integration")
    if (run_pca) selected <- c(selected, "pca")
    message("Modo de ejecucion: ", paste(selected, collapse = " + "), " (argumentos)")
  }
}

ensure_dir(config$dirs$results)
ensure_dir(config$dirs$figures)
ensure_dir(config$dirs$qc_results)
ensure_dir(config$dirs$qc_figures)
ensure_dir(config$dirs$dss_output)

methylkit_files <- NULL
bedgraph_files <- NULL
samples_methylkit <- NULL
samples_bedgraph <- NULL
BSobj <- NULL

if (isTRUE(config$run$methylkit_qc)) {
  methylkit_files <- list_input_files(
    config$dirs$methylkit_input,
    pattern = "\\.methylKit$"
  )
  samples_methylkit <- parse_sample_info(methylkit_files)
}

if (isTRUE(config$run$dss)) {
  bedgraph_files <- list_input_files(
    config$dirs$dss_input,
    pattern = "\\.bedGraph$"
  )
  samples_bedgraph <- parse_sample_info(bedgraph_files)
}

check_sample_consistency(samples_methylkit, samples_bedgraph)

if (isTRUE(config$run$methylkit_qc)) {
  run_methylkit_qc(methylkit_files, samples_methylkit, config)
  rm(methylkit_files, samples_methylkit)
  gc()
}

if (isTRUE(config$run$dss)) {
  BSobj <- build_bsseq_object(bedgraph_files, samples_bedgraph, config)

  if (!"tissue" %in% colnames(pData(BSobj))) {
    idx <- match(sampleNames(BSobj), samples_bedgraph$sample_id)
    if (any(is.na(idx))) {
      stop("No se pudo mapear la metadata a BSobj.")
    }
    pData(BSobj)$tissue <- factor(as.character(samples_bedgraph$group[idx]))
  }
  if (!"individual" %in% colnames(pData(BSobj))) {
    idx <- match(sampleNames(BSobj), samples_bedgraph$sample_id)
    if (any(is.na(idx))) {
      stop("No se pudo mapear la metadata a BSobj.")
    }
    pData(BSobj)$individual <- factor(as.character(samples_bedgraph$replicate[idx]))
  }
  
  # Ejecutar PCA sobre el objeto BSseq (siguiendo metodología NGS101)
  run_pca_bsseq(BSobj, config)

  run_dss_pairwise(BSobj, samples_bedgraph, config)

  if (!isTRUE(config$run$annotation)) {
    rm(BSobj, bedgraph_files, samples_bedgraph)
    BSobj <- NULL
    gc()
  }
}

if (isTRUE(config$run$annotation)) {
  if (is.null(BSobj)) {
    bsobj_path <- config$dss$bsobj_path
    if (file.exists(bsobj_path)) {
      message("Cargando BSobj para anotacion")
      loaded <- load(bsobj_path)
      if (!"BSobj" %in% loaded) {
        warning("No se encontro el objeto BSobj en ", bsobj_path, ". Se omitira el heatmap.")
        BSobj <- NULL
      }
    } else {
      warning("No se encontro ", bsobj_path, ". Se omitira el heatmap.")
    }
  }

  # Analisis funcional post-DSS
  message("Iniciando analisis funcional y anotacion")
  annotate_dmrs_workflow(BSobj, config)
  if (!is.null(BSobj)) {
    rm(BSobj)
  }
  rm(bedgraph_files, samples_bedgraph)
  gc()
}

if (isTRUE(config$run$compare)) {
  message("Comparando DMRs entre condiciones")
  compare_dmrs_across_conditions(config)
  gc()
}

if (isTRUE(config$run$rnaseq)) {
  message("Integrando RNA-seq con DMRs")
  integrate_rnaseq_dmrs(config)
  plot_rnaseq_integration(config)
  gc()
}

if (isTRUE(config$run$pca)) {
  message("Ejecutando PCA")
  run_pca_standalone(config)
  gc()
}
