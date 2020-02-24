library(tidyverse)
library(testthat)
source("~/bin/comparativeSRA/R/buildReferenceLookup.R")


ref_test <- data.frame(start.a = c(1, 10, 30, 39,  139, 160, 170, 181),
                       end.a = c(9, 29, 38, 138, 157, 169, 180, 190),
                       start.b = c(1, 10, 31, 0, 40, 60, -70, -81),
                       end.b = c(9, 29, 39, 0, 58, 69, -80, -89))

ref_test_correct <- data.frame(start.a = c(1, 30,  139, 160, 181),
                           end.a = c(29, 38, 157, 180, 190),
                           strand.a = c("+", "+", "+", "+", "+"),
                           start.b = c(1, 31, 40, 60, 81),
                           end.b = c( 29, 39, 58, 80, 89),
                           strand.b = c("+", "+", "+", "+", "-"),
                           diff = c(0,-1, 99, 100, 100))



test_that("Correct Diff", {

  ref_test_out <- buildReferenceLookup(reference = ref_test, seqA = 1, seqB = 2, remove.missing = T, collapse.alignment = T, quiet = T)


  expect_equal(ref_test_out$diff, ref_test_correct$diff)

})
