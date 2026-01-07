# 10-volcano.R

# Este módulo se encarga de generar Volcano Plots
# Utiliza ggplot2 para la base gráfica y ggrepel para evitar el solapamiento de nombres de genes.

library(ggplot2)
library(ggrepel)
library(dplyr)

plot_volcano <- function(res_df, contrast_name, log2FC_th = 1, padj_th = 0.05, top_n = 15, out_dir = here::here("figures")) {
  
  # Convertimos a data frame por si viene como objeto de DESeq2
  df <- as.data.frame(res_df)
  
  # Verificación de columnas necesarias
  required_cols <- c("log2FoldChange", "padj")
  if (!all(required_cols %in% colnames(df))) {
    stop("Faltan columnas obligatorias en el dataframe (log2FoldChange, padj)")
  }
  
  # Si no existe la columna 'symbol', intentamos usar los nombres de las filas
  if (!"symbol" %in% colnames(df)) {
    df$symbol <- rownames(df)
  }
  
  # Clasificamos los genes según su expresión diferencial
  df <- df %>%
    mutate(status = case_when(
      padj < padj_th & log2FoldChange > log2FC_th ~ "Sobreexpresado (UP)",
      padj < padj_th & log2FoldChange < -log2FC_th ~ "Subexpresado (DOWN)",
      TRUE ~ "No significativo (NS)"
    ))
  
  # Filtramos genes con padj NA (suelen ser genes con pocos conteos)
  df <- df %>% filter(!is.na(padj))
  df <- df %>% mutate(padj_plot = pmax(padj, .Machine$double.xmin))
  
  # Seleccionamos los genes top para etiquetar (ordenados por padj)
  top_genes <- df %>%
    filter(status != "No significativo (NS)") %>%
    arrange(padj) %>%
    head(top_n)
  
  # Construcción del gráfico con ggplot2
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj_plot), color = status)) +
    # Capa de puntos (los NS con mayor transparencia)
    geom_point(alpha = 0.5, size = 1.2) +
    # Colores personalizados (Estilo tesis)
    scale_color_manual(values = c("Sobreexpresado (UP)" = "#E41A1C", 
                                  "Subexpresado (DOWN)" = "#377EB8", 
                                  "No significativo (NS)" = "grey80")) +
    # Líneas de umbral
    geom_vline(xintercept = c(-log2FC_th, log2FC_th), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_th), linetype = "dashed", color = "grey40") +
    # Etiquetas inteligentes con ggrepel
    geom_text_repel(data = top_genes, aes(label = symbol),
                    size = 3.5, 
                    color = "black", 
                    fontface = "italic",
                    box.padding = 0.6,
                    point.padding = 0.3,
                    segment.color = "grey50",
                    max.overlaps = 20) +
    # Estética y etiquetas
    theme_minimal() +
    labs(title = contrast_name,
         subtitle = paste0("Umbrales: |log2FC| > ", log2FC_th, " y padj < ", padj_th),
         x = "Log2FC",
         y = "-log10adjPval",
         color = "Estado") +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          panel.grid.minor = element_blank())

  # Guardado en formato vectorial SVG para máxima resolución en la tesis
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_dir, paste0("Volcano_", contrast_name, ".svg"))
  ggsave(out_file, p, width = 8, height = 7, dpi = 300)
  
  message("Volcano Plot guardado en: ", out_file)
  
  return(p)
}
