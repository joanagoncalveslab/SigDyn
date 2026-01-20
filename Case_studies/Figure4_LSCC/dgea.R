

#Sara Costa 21/04/2025

#install the libraries required 

install.packages("reticulate")
if (!requireNamespace("zellkonverter", quietly = TRUE)) {
  BiocManager::install("zellkonverter")
}
if (!requireNamespace("MAST", quietly = TRUE)) {
  BiocManager::install("MAST")
}

#import libraries
library(zellkonverter)
library(MAST)


# Load adata from .h5ad file
adata_neck <- readH5AD("adata_sbs32_all_26_sep.h5ad")

#define the function that will do the DGEA 
find_de_MAST_RE_simple <- function(adata_){
  sca <- SceToSingleCellAssay(adata_, class = "SingleCellAssay")
  print("Dimensions before subsetting:")
  print(dim(sca)); cat("\n")
  
  sca <- sca[freq(sca) > 0.1, ]
  print("Dimensions after subsetting:")
  print(dim(sca)); cat("\n")
  
  cdr2 <- colSums(assay(sca) > 0)
  colData(sca)$ngeneson <- scale(cdr2)
  
  colData(sca)$label <- factor(colData(sca)$dge_group)
  
  colData(sca)$label <- relevel(factor(colData(sca)$label), ref = "group1")
  colData(sca)$group <- factor(colData(sca)$label)
  
  zlmCond <- zlm(formula = ~ngeneson + group, 
                 sca = sca, 
                 method = 'glm', 
                 ebayes = FALSE, 
                 strictConvergence = FALSE)
  
  lrt_term <- levels(colData(sca)$group)[2]
  summaryCond <- summary(zlmCond, doLRT = paste0("group", lrt_term))
  summaryDt <- summaryCond$datatable
  
  result <- merge(
    summaryDt[contrast == paste0("group", lrt_term) & component == "H", .(primerid, `Pr(>Chisq)`)],
    summaryDt[contrast == paste0("group", lrt_term) & component == "logFC", .(primerid, coef)],
    by = "primerid"
  )
  
  result[, coef := coef / log(2)]
  result[, FDR := p.adjust(`Pr(>Chisq)`, method = "fdr")]
  #result <- result[FDR < 0.05, , drop = FALSE]
  
  result <- stats::na.omit(as.data.frame(result))
  return(result)
}


res <- find_de_MAST_RE_simple(adata_neck)

write.csv(res, "mast_results_26_sep.csv", row.names = FALSE)