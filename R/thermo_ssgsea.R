#' Immune or stroma signature gene-set score for one set
#'
#' @param ranks numeric matrix, rows = genes, cols = samples; entries are ranks per gene across samples
#' @param genes character vector of gene symbols in the gene set
#' @return numeric named vector of length ncol(ranks), names = sample IDs
#' @keywords internal
thermo_ssgsea_score <- function(ranks, genes) {
  common_genes <- intersect(genes, rownames(ranks))
  if (length(common_genes) == 0L) {
    return(setNames(rep(0, ncol(ranks)), colnames(ranks)))
  }
  sranks <- ranks[common_genes, , drop = FALSE]
  num <- colSums(sranks ^ 1.25)
  den <- colSums(sranks ^ 0.25)
  base <- (nrow(ranks) - length(common_genes) + 1) / 2
  out <- (num / den) - base
  setNames(out, colnames(ranks))
}

#' Immune and stroma signatures ssGSEAscores
#'
#' @description
#' ranks = data.T.rank(method=rank_method, na_option='bottom')
#' score = (sranks**1.25).sum()/(sranks**0.25).sum() - (N - |S| + 1)/2
#'
#' @param expr numeric matrix with rows = genes, cols = samples
#' @param gene_sets named list: each element is a character vector of genes
#' @param rank_method one of c("max","min") for tie handling
#' @param na_to_bottom if TRUE, NAs are assigned the smallest rank (0)
#' @return data.frame with rows = samples and columns = gene set names
#' @examples
#' set.seed(1)
#' expr <- matrix(runif(100), nrow = 10,
#'                dimnames = list(paste0("G",1:10), paste0("S",1:10)))
#' gsets <- list(SetA = paste0("G", c(1,3,5,7)), SetB = paste0("G", c(2,4,6)))
#' thermo_ssgsea(expr, gsets)
#' @export
thermo_ssgsea <- function(expr,
                          gene_sets,
                          rank_method = c("max", "min"),
                          na_to_bottom = TRUE) {
  stopifnot(is.matrix(expr), !is.null(rownames(expr)), !is.null(colnames(expr)))
  rank_method <- match.arg(rank_method)
  
  # ranks per gene across samples (rows = genes, cols = samples)
  ranks <- t(apply(expr, 1, function(v) {
    r <- rank(v, ties.method = rank_method, na.last = "keep")
    if (na_to_bottom) r[is.na(r)] <- 0
    r
  }))
  rownames(ranks) <- rownames(expr)
  colnames(ranks) <- colnames(expr)
  
  scores_list <- lapply(names(gene_sets), function(gs) {
    thermo_ssgsea_score(ranks, gene_sets[[gs]])
  })
  names(scores_list) <- names(gene_sets)
  
  scores_df <- as.data.frame(scores_list, check.names = FALSE)
  rownames(scores_df) <- colnames(expr)  # rows = samples
  scores_df
}
