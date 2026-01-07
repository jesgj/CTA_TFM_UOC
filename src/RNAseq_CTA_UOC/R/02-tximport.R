# 02-tximport.R

# Cargamos los ficheros resultantes de la ejecución de kallisto quant

read_quant_files <- function(kallisto_files, tx2gene) {
  txi <- tximport(kallisto_files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)
  sample_names <-  basename(dirname(kallisto_files))
  colnames(txi$counts) <- sample_names
  return(txi)
}
