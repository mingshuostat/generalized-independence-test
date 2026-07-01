# Real data example for the generalized independence test

source("R/GIT_test.R")

breast_file <- "data/gtex/Breast_Mammary_Tissue_normalized_expression_final.txt"
ebv_file <- "data/gtex/Cells_EBV-transformed_lymphocytes_normalized_expression_final.txt"

Breast_Mammary_Tissue_normalized_expression_final <- read.delim(breast_file)
Cells_EBV.transformed_lymphocytes_normalized_expression_final <- read.delim(ebv_file)

inter.ebv.breast <- intersect(
  colnames(Breast_Mammary_Tissue_normalized_expression_final),
  colnames(Cells_EBV.transformed_lymphocytes_normalized_expression_final)
)

breast <- Breast_Mammary_Tissue_normalized_expression_final[, inter.ebv.breast][, 5:length(inter.ebv.breast)]
ebv <- Cells_EBV.transformed_lymphocytes_normalized_expression_final[, inter.ebv.breast][, 5:length(inter.ebv.breast)]

X <- t(breast)
Y <- t(ebv)

set.seed(1)
result <- GIT_test(X, Y)

print(result)
