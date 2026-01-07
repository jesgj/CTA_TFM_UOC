# remotes::install_github("giocomai/ganttrify")

# Cargar las librerías necesarias
library(dplyr)
library(lubridate)
library(ganttrify)
library(ggplot2)
library(svglite)

proyecto_df <- data.frame(
  etapa = factor(c(
    rep("ETAPA 1: Preparación del análisis y los datos", 3),
    rep("ETAPA 2: Procesamiento de los datos en entorno HPC", 3),
    rep("ETAPA 3: Análisis e integración", 4),
    rep("ETAPA 4: Finalización y Entrega", 3)
  ), levels = c(
    "ETAPA 1: Preparación del análisis y los datos",
    "ETAPA 2: Procesamiento de los datos en entorno HPC",
    "ETAPA 3: Análisis e integración",
    "ETAPA 4: Finalización y Entrega"
  )),
  tarea = c(
    # ETAPA 1
    "Tarea 1.1 - Descarga de datos crudos alojados en bases de datos públicas",
    "Tarea 1.2 - Configuración de entornos",
    "Tarea 1.3 - Preparación de pipeline de Snakemake",

    # ETAPA 2
    "Tarea 2.1 - Procesamiento de datos de RNA-seq",
    "Tarea 2.2 - Procesamiento de datos de WGBS",
    "Tarea 2.3 - Visualización de los datos públicos frente a los CTA",

    # ETAPA 3
    "Tarea 3.1 - Análisis diferencial y enriquecimiento funcional (RNAseq)",
    "Tarea 3.2 - Análisis de metilación: DMRs (WGBS)",
    "Tarea 3.3 - Integración de los resultados de RNAseq y WGBS",
    "Tarea 3.4 - Preparación de figuras finales",

    # ETAPA 4
    "Tarea 4.1 - Terminar redacción de la memoria",
    "Tarea 4.2 - Preparación de la presentación",
    "Tarea 4.3 - Defensa"
  ),
  fecha_inicio = as.Date(c(
    "2025-09-25", "2025-09-25", "2025-10-02",
    "2025-10-09", "2025-10-16", "2025-10-30",
    "2025-11-13", "2025-11-20", "2025-12-04", "2025-12-18",
    "2025-12-18", "2026-01-07", "2026-01-19"
  )),
  fecha_fin = as.Date(c(
    "2025-10-31",
    "2025-10-01", "2025-10-08",
    "2025-10-15", "2025-10-29", "2025-11-19",
    "2025-12-09", "2025-12-09", "2025-12-17", "2026-01-07",
    "2026-01-07", "2026-01-13", "2026-01-20"
  ))
)

hitos_df <- data.frame(
  fecha = as.Date(c("2025-10-08", "2025-11-05", "2025-12-03",
                    "2026-01-07", "2026-01-13", "2026-01-20")),
  etiqueta = c("Entrega PEC1", "Entrega PEC2", "Entrega PEC3",
               "Entrega PEC4", "Entrega PEC5", "Defensa Final")
)

proyecto_gantt <- proyecto_df %>%
  rename(wp = etapa, activity = tarea, start_date = fecha_inicio, end_date = fecha_fin)

gantt_plot <- ganttrify(
  project = proyecto_gantt,
  project_start_date = "2025-09-25", by_date = TRUE, exact_date = TRUE,
  font_family = "sans", mark_years = TRUE, month_breaks = 1,
  wp_label_bold = FALSE, show_vertical_lines = FALSE, colour_stripe = "grey90"
)

for (i in 1:nrow(hitos_df)) {
  gantt_plot <- gantt_plot +
    annotate(
      geom = "segment",
      x = hitos_df$fecha[i], xend = hitos_df$fecha[i],
      y = 0, yend = Inf,
      linetype = "dashed", color = "purple"
    ) +
    annotate(
      geom = "text",
      x = hitos_df$fecha[i], y = Inf,
      label = hitos_df$etiqueta[i],
      color = "purple", angle = 90,
      vjust = 1.5, hjust = 1.1,
      size = 3.5, fontface = "bold"
    )
}

gantt_plot <- gantt_plot +
  scale_x_date(
    limits = as.Date(c("2025-09-25", "2026-01-25")),
    date_labels = "%d %b",
    date_breaks = "1 week"
  ) +
  labs(
    title = "",
    subtitle = "(25/09/2025 - 20/01/2026)"
  ) +
  theme_minimal(base_size = 12) +
  xlab(NULL) + ylab(NULL) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12, margin = margin(b = 20)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

print(gantt_plot)
ggsave("gantt_diagram_TFM_JGJ.svg", plot = gantt_plot, width = 15, height = 9)
