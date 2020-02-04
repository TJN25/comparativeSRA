#' alignAndCombine Function
#'
#' @param reference The file reference positons
#' @param gff1 The gff file which is the reference gff
#' @param gff2 the gff file that needs rearranging
#' @param filenum1 number or name for the first gff file
#' @param filenum2 number or name for the second gff file
#' @keywords combines multiple gff files from different genomes
#' @export
#' @examples
#' alignAndCombine()


alignAndCombine <- function(reference, gff1, gff2, time.it = T, quiet = F, filenum1 = "1", filenum2 = "2", seqA = 1, seqB = 2){

  referenceEsch1Serr1 <- read.table(reference, header = T, as.is = T)
  referenceEsch1Serr1Built <- buildReferenceLookup(reference = referenceEsch1Serr1,
                                                   as.numeric(seqA), seqB = as.numeric(seqB),
                                                   collapse.alignment = T,
                                                   quiet = quiet)

  esch1 <- gff1
  serr1 <- gff2

  esch1 <- esch1%>%mutate(changed = F)
  serr1 <- serr1%>%mutate(changed = F)



  serr1b <- reorderGFF(ref = referenceEsch1Serr1Built, gff = serr1, time.it = time.it, quiet = quiet)
  esch1b <- reorderGFF(ref = referenceEsch1Serr1Built, gff = esch1, reference.genome = T, time.it = time.it, quiet = quiet)

  serr1b <- serr1b%>%mutate(file_id = filenum2)
  esch1b <- esch1b%>%mutate(file_id = filenum1)
  ncRNAgff <- esch1b%>%bind_rows(serr1b)

  serr1 <- serr1%>%mutate(file_id = filenum2)
  esch1 <- esch1%>%mutate(file_id = filenum1)
  ncRNAgffUnchanged <- esch1%>%bind_rows(serr1)
  tmp1 <- ncRNAgff %>% select(id) %>% mutate(found = T)
  tmp2 <- ncRNAgffUnchanged %>% select(id)
  tmp1 <- tmp1 %>% full_join(tmp2, by = "id") %>% filter(is.na(found)) %>% mutate(found = F)

  otherncRNA <- ncRNAgffUnchanged %>% full_join(tmp1, by = "id") %>% filter(found == F) %>% select(-found) %>%
    mutate(start = -start, end = -end)

  ncRNAgff <- ncRNAgff %>% bind_rows(otherncRNA)

  return(ncRNAgff)

}
