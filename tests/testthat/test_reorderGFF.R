library(testthat)
library(tidyverse)
source("~/bin/comparativeSRA/R/reorderGFF.R")


load("~/bin/r_git/R/reorderGFFData.Rda")
load("~/bin/r_git/R/reorderGFFDataCorrect.Rda")

reference <- reorderGFFData[['reference']]
reference <- buildReferenceLookup(reference = reference, seqA = 5, seqB = 6)

gff1 <- reorderGFFData[['gff1']]
gff2 <- reorderGFFData[['gff2']][1:50,]

gff1_correct <- reorderGFFDataCorrect[['gff1_correct']]
gff2_correct <- reorderGFFDataCorrect[['gff2_correct']]





test_that("Correct reordering", {

  gff1_out <- reorderGFF(ref = reference, gff = gff1, reference.genome = T, quiet = T)
  gff2_out <- reorderGFF(ref = reference, gff = gff2, reference.genome = F, quiet = T)

  ##ignoring strand as this is not used
  expect_equal(gff1_out$start, gff1_correct$start)
  expect_equal(gff2_out$start, gff2_correct$start)
  expect_equal(gff1_out$end, gff1_correct$end)
  expect_equal(gff2_out$end, gff2_correct$end)
})
