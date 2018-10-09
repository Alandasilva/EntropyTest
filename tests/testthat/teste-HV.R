library(testthat)
library(EntropyTest)
context("Tamanho da saída")
test_that("Testa se HVmn(x) retorna um argumento unico", {
  expect_equal(length(HVmn(rnorm(10))), 1)
})

