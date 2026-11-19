# Fig 4c/d over the missingness levels, step 2: MAST DGEA for one level.
# find_de_MAST_RE_simple is copied unchanged from Case_studies/Figure4_LSCC/dgea.R (Sara Costa 21/04/2025).
#
#   Rscript dgea_mast.R --level 0                       # neck, default --out levels
#   Rscript dgea_mast.R --level 30 --out test --cores 4
#
# Reads  <out>/<name>/dgea/expression.h5ad and <out>/<name>/dgea/pXX/groups.csv (dgea_prepare.py)
# Writes <out>/<name>/dgea/pXX/mast_results.csv  (primerid, Pr(>Chisq), coef = log2FC H32 vs L32, FDR)
# --cores sets options(mc.cores) for MAST's zlm (same results, more memory).

library(zellkonverter)
library(MAST)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  args[i + 1]
}
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
name  <- get_arg("--name", "neck")
level <- get_arg("--level")
out   <- get_arg("--out", file.path(script_dir, "levels"))
cores <- as.integer(get_arg("--cores", "1"))
if (is.null(level)) stop("--level is required (e.g. 0, 30, p60)")
level <- if (grepl("^p", level)) level else sprintf("p%02d", as.integer(level))
options(mc.cores = cores)

dgea_dir <- file.path(out, name, "dgea")
groups_file <- file.path(dgea_dir, level, "groups.csv")
result_file <- file.path(dgea_dir, level, "mast_results.csv")
if (!file.exists(groups_file)) stop("missing ", groups_file, " (run dgea_prepare.py first)")

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

cat(format(Sys.time()), level, ": reading", file.path(dgea_dir, "expression.h5ad"), "\n")
adata_neck <- readH5AD(file.path(dgea_dir, "expression.h5ad"))
groups <- read.csv(groups_file, stringsAsFactors = FALSE)
stopifnot(setequal(groups$cell, colnames(adata_neck)))
colData(adata_neck)$dge_group <- groups$dge_group[match(colnames(adata_neck), groups$cell)]
print(table(colData(adata_neck)$dge_group))

cat(format(Sys.time()), level, ": MAST\n")
res <- find_de_MAST_RE_simple(adata_neck)

write.csv(res, paste0(result_file, ".tmp"), row.names = FALSE)
file.rename(paste0(result_file, ".tmp"), result_file)
cat(format(Sys.time()), level, ": wrote", result_file, "\n")
