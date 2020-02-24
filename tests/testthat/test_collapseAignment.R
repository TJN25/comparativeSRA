library(tidyverse)
library(testthat)
source("~/bin/comparativeSRA/R/collapseAlignment.R")


ref_test <- data.frame(start.a = c(1, 10, 30, 39, 60, 70, 81),
                  end.a = c(9, 29, 38, 57, 69, 80, 90),
                  strand.a = c("+", "+", "+", "+", "+", "+", "+"),
                  start.b = c(1, 10, 31, 40, 60, 70, 81),
                  end.b = c(9, 29, 39, 58, 69, 80, 89),
                  strand.b = c("+", "+", "+", "+", "+", "-", "+"))

ref_test_2 <- data.frame(start.a = c(1, 10, 30, 39, 60, 70, 81),
                       end.a = c(9, 29, 38, 57, 69, 80, 90),
                       strand.a = c("+", "+", "+", "+", "+", "+", "+"),
                       start.b = c(1, 10, 31, 40, 60, 70, 81),
                       end.b = c(9, 29, 39, 58, 69, 80, 90),
                       strand.b = c("+", "+", "+", "+", "+", "-", "+"))

ref_test_correct <- data.frame(start.a = c(1, 30, 60, 81),
                                 end.a = c(29, 57, 80, 90),
                                 strand.a = c("+", "+", "+", "+"),
                                 start.b = c(1, 31, 60, 81),
                                 end.b = c(29, 58, 80, 89),
                                 strand.b = c("+", "+", "+", "+"))

ref_test_2_correct <- data.frame(start.a = c(1, 30, 60),
                       end.a = c(29, 57, 90),
                       strand.a = c("+", "+", "+"),
                       start.b = c(1, 31, 60),
                       end.b = c(29, 58, 90),
                       strand.b = c("+", "+", "+"))



test_that("Correct Merging", {

  ref_test_out <- collapseAlignment(ref_test)
  ref_test_2_out <- collapseAlignment(ref_test_2)

  ##ignoring strand as this is not used
  expect_equal(ref_test_out[,c(1,2,4,5)], ref_test_correct[,c(1,2,4,5)])
  expect_equal(ref_test_2_out[,c(1,2,4,5)], ref_test_2_correct[,c(1,2,4,5)])

})
