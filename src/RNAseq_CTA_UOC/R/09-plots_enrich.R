# 09-plots_enrich.R

# Visualizaciones avanzadas para resultados de enriquecimiento (ORA y GSEA)
# usando el paquete enrichplot.

library(enrichplot)
library(ggplot2)

strip_term_prefixes <- function(labels) {
  if (is.null(labels)) {
    return(labels)
  }
  labels <- as.character(labels)
  patterns <- c("^REACTOME_", "^KEGG_", "^HALLMARK_")
  for (pat in patterns) {
    labels <- gsub(pat, "", labels)
  }
  labels
}

sanitize_term_descriptions <- function(enrich_obj) {
  if (is.null(enrich_obj)) {
    return(enrich_obj)
  }
  clean_df <- function(df) {
    if ("Description" %in% colnames(df)) {
      df$Description <- strip_term_prefixes(df$Description)
    }
    df
  }
  if (methods::is(enrich_obj, "compareClusterResult")) {
    enrich_obj@compareClusterResult <- clean_df(enrich_obj@compareClusterResult)
  } else if (methods::is(enrich_obj, "enrichResult") || methods::is(enrich_obj, "gseaResult")) {
    enrich_obj@result <- clean_df(enrich_obj@result)
  } else if ("result" %in% methods::slotNames(enrich_obj)) {
    enrich_obj@result <- clean_df(enrich_obj@result)
  }
  enrich_obj
}

plot_ora_advanced <- function(enrich_obj, contrast_name, ont_name, out_dir = here::here("figures")) {
  # Verificar si hay resultados
  if (is.null(enrich_obj)) return(NULL)
  res_df <- as.data.frame(enrich_obj)
  if (nrow(res_df) == 0) {
    message("No hay resultados significativos para ", contrast_name, " ", ont_name)
    return(NULL)
  }

  message("Generando visualizaciones avanzadas ORA para: ", contrast_name, " ", ont_name)
  enrich_obj <- sanitize_term_descriptions(enrich_obj)

  # 1. Cnetplot: Relación entre genes y términos (Top 5)
  try({
    p_cnet <- cnetplot(enrich_obj, showCategory = 5, cex_label_gene = 0.8) + 
      ggtitle(paste("Cnetplot:", contrast_name, ont_name))
    ggsave(file.path(out_dir, paste0(contrast_name, "_", ont_name, "_cnetplot.svg")), 
           p_cnet, width = 12, height = 9)
  }, silent = TRUE)

      # 2. Treeplot: Clustering jerárquico de términos enriquecidos

      try({
        enrich_obj_sim <- pairwise_termsim(enrich_obj)
        p_tree <- treeplot(enrich_obj_sim, label_format = 30, fontsize_tiplab=2, fontsize_cladelab=4) + 
          ggtitle(paste("Treeplot:", contrast_name, ont_name))
        ggsave(file.path(out_dir, paste0(contrast_name, "_", ont_name, "_treeplot.svg")), 
               p_tree, width = 30, height = 12)
      }, silent = TRUE)
    }

    # Generar visualizaciones avanzadas para GSEA
    plot_gsea_advanced <- function(gsea_obj, contrast_name, ont_name, out_dir = here::here("figures")) {
      if (is.null(gsea_obj)) return(NULL)
      res_df <- as.data.frame(gsea_obj)
      if (nrow(res_df) == 0) {
        message("No hay resultados significativos GSEA para ", contrast_name, " ", ont_name)
        return(NULL)
      }
      message("Generando visualizaciones avanzadas GSEA para: ", contrast_name, " ", ont_name)
      gsea_obj <- sanitize_term_descriptions(gsea_obj)
      # 1. Cnetplot: Relación entre genes y términos (Top 5)
      # Desactivado para GSEA por redundancia con Ridgeplot
      # try({
      #   p_cnet <- cnetplot(gsea_obj, showCategory = 5, cex_label_gene = 0.8) + 
      #   ggtitle(paste("Cnetplot GSEA:", contrast_name, ont_name))
      #   ggsave(file.path(out_dir, paste0(contrast_name, "_", ont_name, "_GSEA_cnetplot.svg")), 
      #          p_cnet, width = 12, height = 9)
      # }, silent = TRUE)

      # 2. Treeplot: Clustering jerárquico de términos enriquecidos

      try({
        gsea_obj_sim <- pairwise_termsim(gsea_obj)
        p_tree <- treeplot(gsea_obj_sim, label_format = 30, fontsize_tiplab=2, fontsize_cladelab=4) + 
          ggtitle(paste("Treeplot GSEA:", contrast_name, ont_name))
        ggsave(file.path(out_dir, paste0(contrast_name, "_", ont_name, "_GSEA_treeplot.svg")), 
               p_tree, width = 30, height = 12)
      }, silent = TRUE)

      # 3. Ridgeplot: Distribución de Fold Changes por término
      try({
        p_ridge <- ridgeplot(gsea_obj, showCategory = 15) + 
          labs(title = paste("Ridgeplot GSEA:", contrast_name, ont_name),
               x = "stat (Distribution)") +
          theme(plot.title = element_text(face = "bold", hjust = 0.5),
                axis.text.y = element_text(size = 9))
        ggsave(file.path(out_dir, paste0(contrast_name, "_", ont_name, "_GSEA_ridgeplot.svg")), 
               p_ridge, width = 15, height = 9)
      }, silent = TRUE)

    }

    

  
