# Configuracion principal del pipeline
get_config <- function() {
  list(
    dirs = list(
      methylkit_input = "../WGBS_methyldackel_results/methyldackel",
      dss_input = "../WGBS_methyldackel_results/methyldackel_mergecontext",
      results = "results",
      figures = "figures",
      qc_results = file.path("results", "qc"),
      qc_figures = file.path("figures", "qc"),
      dss_output = "dmrs_span500_delta0.2"
    ),
    run = list(
      methylkit_qc = TRUE,
      dss = TRUE,
      annotation = TRUE,
      compare = FALSE,
      rnaseq = FALSE,
      rebuild_bsobj = FALSE
    ),
    methylkit = list(
      # Ensamblado de referencia Ensembl para pez cebra
      assembly = "GRCz11",
      context = "CpG",
      mincov = 10,
      filter_by_coverage = list(lo.count = 10, hi.perc = 99.9),
      normalize = TRUE,
      min_per_group = 1L,
      cores = 1,
      standard_chr_max = 25,
      standard_chromosomes = NULL,
      remove_chr_regex = "^(Lambda|chrLambda|lambda)$"
    ),
    dss = list(
      bsobj_path = "BSobj_CTA.RData",
      standard_chr_max = 25,
      standard_chromosomes = NULL,
      remove_chr_regex = "^(Lambda|chrLambda|lambda)$",
      comparison_mode = "control_vs_groups",
      reference_group = "Ctr",
      target_groups = NULL,
      require_common_cpgs = FALSE,
      smoothing_span = 500,
      delta = 0.2,
      p_threshold = 0.05,
      minCG = 10,
      minlen = 100,
      dis_merge = 100,
      ncores = 8
    ),
    annotation = list(
      # Anotacion desde fichero GTF local (Ensembl 115)
      gtf_file = "ref/Danio_rerio.GRCz11.115.gtf.gz",
      txdb_sqlite = "ref/Danio_rerio.GRCz11.115.sqlite", # Cache compilado
      org_package = "org.Dr.eg.db",
      gene_id_type = "ENTREZID",
      bitr_suppress_warnings = TRUE,
      kegg_organism = "dre", # Codigo KEGG para Danio rerio
      promoter_upstream = 2000,
      promoter_downstream = 500,
      enrichment_pvalue = 0.05,
      heatmap_top_n = 100
    ),
    rnaseq = list(
      input_dir = file.path("results", "deseq2"),
      annotation_dir = file.path("results", "annotation"),
      output_dir = file.path("results", "rnaseq_integration"),
      gene_id_column = NULL,
      padj_threshold = 0.05,
      lfc_threshold = 0,
      lfc_column = NULL,
      comparison_map = c(elavl2_vs_control = "elav2_vs_Ctr")
    )
  )
}
