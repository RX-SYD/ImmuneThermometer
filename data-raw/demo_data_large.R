set.seed(2025)

# Expression matrix: 500 genes x 24 samples (rows = genes, cols = samples)
genes <- paste0("GENE", sprintf("%04d", 1:500))
samples <- paste0("S", sprintf("%02d", 1:24))

group <- rep(c("GroupA","GroupB"), each = 12)

expr <- matrix(rnorm(500 * 24, mean = 0, sd = 1),
               nrow = 500, dimnames = list(genes, samples))

shift_genes <- sample(genes, 100)
expr[shift_genes, group == "GroupA"] <- expr[shift_genes, group == "GroupA"] + 0.6

sample_info_demo <- data.frame(
  sample = samples,
  group  = group,
  stringsAsFactors = FALSE
)

# Five demo gene sets with realistic sizes (30–60 genes each)
set.seed(2025)
gene_sets_demo <- list(
  Immune_Sig_1  = sample(genes, 40),
  Immune_Sig_2  = sample(genes, 50),
  Stromal_Sig_1 = sample(genes, 30),
  Stromal_Sig_2 = sample(genes, 60),
  Mixed_Sig     = sample(genes, 45)
)
expr_demo <- expr
usethis::use_data(expr_demo = expr, sample_info_demo, gene_sets_demo, overwrite = TRUE)
