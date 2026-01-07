# 01-t2g.R

# Generamos el fichero t2g.csv si no se detecta en el directorio de trabajo
# a partir del gtf del ensamblaje 115 de Danio rerio (Ensembl)

generate_t2g <- function(file_name, write_t2g=F) {
  # Función que genera el fichero t2g.csv si no se detecta
    message("Generating t2g...")
    gtf_url <- "https://ftp.ensembl.org/pub/release-115/gtf/danio_rerio/Danio_rerio.GRCz11.115.gtf.gz"
    gtf <- import(gtf_url, format = "gtf")
    tx <- gtf[mcols(gtf)$type == "transcript"]
    m <- mcols(tx)
    t2g <- data.frame(
      transcript_id = as.character(m$transcript_id),
      gene_id       = as.character(m$gene_id),
      stringsAsFactors = FALSE)
    if (write_t2g == T) {
      out_file <- here::here("results", paste0(strsplit(basename(file_name), split = "\\.")[[1]][1],".csv"))
      write.csv(t2g, out_file, row.names = FALSE)
      message("t2g table written to: ", out_file)
    }
    return(t2g)
  } 
