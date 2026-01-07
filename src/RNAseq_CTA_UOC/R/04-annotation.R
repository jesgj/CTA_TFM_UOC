# 04-annotation.R

annot <- function(deseq_res) {
  # add column "ENTREZID"
  deseq_res$EntrezID <- mapIds(org.Dr.eg.db, keys=rownames(deseq_res), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
  # add column "SYMBOL"
  deseq_res$symbol <- mapIds(org.Dr.eg.db, keys=rownames(deseq_res), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  # add column "GENENAME"
  deseq_res$genename <- mapIds(org.Dr.eg.db, keys=rownames(deseq_res), column="GENENAME", keytype="ENSEMBL", multiVals="first")

  #  si el símbolo es NA, usamos el Ensembl ID (rownames)
  missing_sym <- is.na(deseq_res$symbol)
  deseq_res$symbol[missing_sym] <- rownames(deseq_res)[missing_sym]

  return(deseq_res)
}