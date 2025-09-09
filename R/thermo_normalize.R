#' Within-cancer Modified Z-score normalization
#'
#' Performs per-cancer-type Modified Z-score normalization on a 
#' gene-by-sample expression matrix.
#'
#' Z = 0.6745 * (x - median) / MAD
#'
#' @param expr Numeric matrix, rows = genes, cols = samples. Row/col names required.
#' @param groups Character/factor vector of cancer types per sample.
#'   Names of \code{groups} must match \code{colnames(expr)}.
#' @param eps Small positive constant to avoid division by zero when MAD = 0.
#' @return A numeric matrix of the same dimension as \code{expr}, normalized within cancer types.
#' @export
thermo_normalize <- function(expr, groups, eps = 1e-8) {
  stopifnot(is.matrix(expr), !is.null(rownames(expr)), !is.null(colnames(expr)))
  stopifnot(all(colnames(expr) %in% names(groups)))
  groups <- factor(groups[colnames(expr)])
  
  out <- expr
  for (lev in levels(groups)) {
    idx <- which(groups == lev)
    block <- expr[, idx, drop = FALSE]
    med <- apply(block, 1, median, na.rm = TRUE)
    centered <- sweep(block, 1, med, "-")
    madv <- apply(centered, 1, function(x) median(abs(x), na.rm = TRUE))
    madv[is.na(madv) | madv < eps] <- eps
    z <- 0.6745 * sweep(centered, 1, madv, "/")
    out[, idx] <- z
  }
  out
}
