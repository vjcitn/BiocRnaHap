library(BiocRnaHap)

context("test look1kg")
test_that("retrieve VCF data working", {
 l1 = look1kg()
 dd = dim(l1)
 expect_true(all(dd == c(3, 2504)))
})
