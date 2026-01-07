# 11-comparative.R

# Este módulo realiza un análisis comparativo entre los diferentes contrastes
# para identificar genes compartidos (core response) usando diagramas de Venn y UpSet plots.

library(dplyr)
library(ggplot2)
library(ggvenn)   # Para diagramas de Venn clásicos (hasta 4 conjuntos)
library(ggupset)  # Para intersecciones complejas

#' Extraer listas de genes significativos
#'
#' @param res_list Lista nombrada de dataframes de resultados (annotados)
#' @param padj_th Umbral de p-valor ajustado
#' @param lfc_th Umbral de log2FoldChange (valor absoluto)
#' @param direction Dirección del cambio: "UP", "DOWN" o "BOTH"
#'
#' @return Lista nombrada de vectores de símbolos de genes
get_gene_lists <- function(res_list, padj_th = 0.05, lfc_th = 1, direction = "BOTH") {
  gene_lists <- list()
  
  for (name in names(res_list)) {
    df <- as.data.frame(res_list[[name]])
    
    # Filtrar por significancia
    df_sig <- df %>% filter(!is.na(padj) & padj < padj_th)
    
    # Filtrar por dirección y magnitud
    if (direction == "UP") {
      df_sig <- df_sig %>% filter(log2FoldChange > lfc_th)
    } else if (direction == "DOWN") {
      df_sig <- df_sig %>% filter(log2FoldChange < -lfc_th)
    } else {
      df_sig <- df_sig %>% filter(abs(log2FoldChange) > lfc_th)
    }
    
    # Extraer símbolos (o IDs si símbolo no existe/es NA)
    if ("symbol" %in% colnames(df_sig)) {
      genes <- na.omit(df_sig$symbol)
    } else {
      genes <- rownames(df_sig)
    }
    
    gene_lists[[name]] <- unique(as.character(genes))
  }
  return(gene_lists)
}

#' Plot Venn Diagram
#'
#' @param gene_lists Lista nombrada de vectores de genes
#' @param title Título del gráfico
#' @param out_file Ruta de salida (SVG)
plot_venn_custom <- function(gene_lists, title, out_file) {
  # ggvenn soporta hasta 4 conjuntos, perfecto para 3 contrastes
  p <- ggvenn(
    gene_lists, 
    fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
    stroke_size = 0.5, set_name_size = 4
  ) + ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  ggsave(out_file, p, width = 8, height = 8)
}

#' Plot UpSet Diagram (Alternativa moderna al Venn)
#'
#' @param gene_lists Lista nombrada de vectores de genes
#' @param title Título del gráfico
#' @param out_file Ruta de salida (SVG)
plot_upset_custom <- function(gene_lists, title, out_file) {
  # Transformar lista a formato compatible con ggupset (tibble)
  # Creamos un dataframe con todos los genes únicos y una columna lista con los sets a los que pertenecen
  all_genes <- unique(unlist(gene_lists))
  
  dat <- tibble(gene = all_genes) %>%
    mutate(sets = lapply(gene, function(g) {
      names(gene_lists)[sapply(gene_lists, function(s) g %in% s)]
    }))
  
  p <- ggplot(dat, aes(x = sets)) +
    geom_bar() +
    scale_x_upset() +
    labs(title = title, x = "Intersección de Conjuntos", y = "Número de Genes") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face="bold"))
  
  ggsave(out_file, p, width = 10, height = 7)
}

#' Ejecutar análisis comparativo completo
#'
#' @param res_list Lista con los resultados de DESeq2 (e.g., list(dazl=res_dazl, ddx4=res_ddx4...))
#' @param out_dir Directorio de salida
run_comparative_analysis <- function(res_list, out_dir = here::here("figures", "comparative")) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. Genes UPREGULATED (Activación ectópica)
  up_genes <- get_gene_lists(res_list, direction = "UP", lfc_th = 1, padj_th = 0.05)
  plot_venn_custom(up_genes, "Genes UP-regulados Compartidos", file.path(out_dir, "Venn_UP.svg"))
  plot_upset_custom(up_genes, "Intersección Genes UP", file.path(out_dir, "UpSet_UP.svg"))
  
  # 2. Genes DOWNREGULATED (Represión)
  down_genes <- get_gene_lists(res_list, direction = "DOWN", lfc_th = 1, padj_th = 0.05)
  plot_venn_custom(down_genes, "Genes DOWN-regulados Compartidos", file.path(out_dir, "Venn_DOWN.svg"))
  
  # 3. Guardar las listas de intersección (Core Germinal Program)
  # Intersección de los tres (si existen 3 elementos)
  if (length(up_genes) >= 2) {
    common_up <- Reduce(intersect, up_genes)
    write.table(common_up, file.path(out_dir, "Core_Germinal_Genes_UP.txt"), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    message("Núcleo de genes UP compartido (", length(common_up), " genes) guardado.")
  }
}
