#' Demo expression matrix (genes x samples, 500 x 24)
#' @format Numeric matrix with genes as rows and samples as columns.
#' @examples
#' data(expr_demo); dim(expr_demo)
"expr_demo"

#' Demo sample metadata (24 samples)
#' @format Data frame with columns: sample, group.
#' @examples
#' data(sample_info_demo); head(sample_info_demo)
"sample_info_demo"

#' Demo gene sets (5 sets, sizes 30–60)
#' @format Named list of character vectors (gene symbols).
#' @examples
#' data(gene_sets_demo); names(gene_sets_demo)
"gene_sets_demo"
