# Visualizacion PCA independiente

run_pca_standalone <- function(config) {
  suppressPackageStartupMessages({
    library(bsseq)
    library(ggplot2)
    library(RColorBrewer)
  })

  # Cargar BSobj
  bsobj_path <- config$dss$bsobj_path
  if (!file.exists(bsobj_path)) {
    stop("No se encuentra el objeto BSobj en: ", bsobj_path, 
         "\nAsegurate de haber ejecutado el paso 'dss' al menos una vez para crearlo.")
  }
  
  message("Cargando BSobj desde: ", bsobj_path)
  load(bsobj_path)
  
  if (!exists("BSobj")) {
    stop("El archivo RData no contiene un objeto llamado 'BSobj'.")
  }

  # Reutilizamos la logica de PCA que ya definimos en dss-analysis.R
  # pero nos aseguramos de que este disponible o la redefinimos aqui si es necesario.
  # Como run_pca_bsseq esta en 04-dss-analysis.R y ese script se carga en main.R,
  # podemos llamarla directamente.
  
  message("Generando PCA actualizado...")
  run_pca_bsseq(BSobj, config)
}

