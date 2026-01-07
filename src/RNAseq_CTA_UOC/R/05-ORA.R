# 05-ORA.R

run_GO_ORA <- function(deseq_res,
                       log2fc_th = 1,
                       padj_th = 0.05,
                       OrgDb         = org.Dr.eg.db,
                       contrast_name = "ora",
                       showCategory  = 10,
                       save_tables   = TRUE,
                       universe      = FALSE) {
  # Ordenados la tabla de resultados por el p-value ajustado por FDR
  deseq_res <- deseq_res[order(deseq_res$padj),]
  # Fijamos previamente dos umbrales para determinar los genes up- o down- regulados
  deg <- list(up = deseq_res$EntrezID[which(deseq_res$log2FC_shrunk > log2fc_th & deseq_res$padj < padj_th)],
                   down = deseq_res$EntrezID[which(deseq_res$log2FC_shrunk < -log2fc_th & deseq_res$padj < padj_th)])
  # Eliminamos de la lista los términos "NA"
  deg <- lapply(deg, na.omit)
  
  if (universe == FALSE) {
    universe <- NULL
  } else {
    universe <- na.omit(deseq_res$EntrezID[!is.na(deseq_res$padj)])
  }
  print(head(deseq_res))
  # 1) GO: BP
  fea_GO_BP <- simplify(compareCluster(deg,
                                ont   = "BP",
                                fun   = "enrichGO",
                                OrgDb = OrgDb,
                                readable = TRUE,
                                universe = universe))
  BP_title <- paste0("BP: Biological process (", contrast_name, ")")
  BPp <- dotplot(fea_GO_BP,
                 showCategory = showCategory,
                 title = BP_title) +
    theme(plot.title   = element_text(face = "bold",
                                      hjust = 0.5,
                                      size  = 16),
          axis.text    = element_text(face = "bold", size = 11),
          axis.title.x = element_text(face = "bold"),
          legend.text  = element_text(face = "bold", size = 11),
          legend.title = element_text(face = "bold", size = 13))
  
  ggsave(filename = here::here("figures", paste0(contrast_name, "_BP_cluster.svg")),
         plot     = BPp,
         dpi      = "print",
         width    = 12,
         height   = 9)
  
  if (save_tables) {
    write.table(x         = as.data.frame(fea_GO_BP),
                file      = here::here("results", paste0(contrast_name, "_GO_BP_enriquecimiento.txt")),
                row.names = FALSE,
                sep       = "\t",
                quote     = FALSE)
  }
  
  # 2) GO: MF
  fea_GO_MF <- simplify(compareCluster(deg,
                                ont   = "MF",
                                fun   = "enrichGO",
                                OrgDb = OrgDb,
                                readable = TRUE,
                                universe = universe))
  MF_title <- paste0("MF: Molecular function (", contrast_name, ")")
  MFp <- dotplot(fea_GO_MF,
                 showCategory = showCategory,
                 title = MF_title) +
    theme(plot.title   = element_text(face = "bold",
                                      hjust = 0.5,
                                      size  = 16),
          axis.text    = element_text(face = "bold", size = 11),
          axis.title.x = element_text(face = "bold"),
          legend.text  = element_text(face = "bold", size = 11),
          legend.title = element_text(face = "bold", size = 13))
  ggsave(filename = here::here("figures", paste0(contrast_name, "_MF_cluster.svg")),
         plot     = MFp,
         dpi      = "print",
         width    = 12,
         height   = 9)
  if (save_tables) {
    write.table(x         = as.data.frame(fea_GO_MF),
                file      = here::here("results", paste0(contrast_name, "_GO_MF_enriquecimiento.txt")),
                row.names = FALSE,
                sep       = "\t",
                quote     = FALSE)
  }
  return(list(BP = fea_GO_BP, MF = fea_GO_MF))
}
