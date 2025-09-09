cancer <- data.frame(
  TAX  = c(1,2,3,4,5,6,7,8),
  OPTN = c(6,7,8,9,10,11,1,0.5),
  P62  = c(16,5,4,8,3,5,99,0.3),
  row.names = paste0("TCGA", 1:8)
)
cancer_expr <- t(cancer)

cancerinfo <- data.frame(
  sample = paste0("TCGA", 1:8),
  cancer_type = c("breast","breast","breast","breast",
                  "liver","liver","liver","liver"),
  stringsAsFactors = FALSE
)

usethis::use_data(cancer_expr, cancerinfo, overwrite = TRUE)
