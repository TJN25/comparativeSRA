#' buildReferenceLookup Function
#'
#' @param reference The file reference positons
#' @param seqA which columns to use for the first sequence
#' @param seqB which columns to use for the second sequence
#' @param remove.missing remove rows that are not conserved
#' @param collapse.alignment combines segments which can be combines once other alignments are removed
#' @keywords Takes a full .backbone input and produces a file with two sequences, sorted and all with +ve values and strand info
#' @export
#' @examples
#' buildReferenceLookup()

buildReferenceLookup <- function(reference, seqA = 1, seqB = 2, remove.missing = T, collapse.alignment = T, quiet = F){
  refA <- reference[,(seqA*2 - 1):(seqA*2)]
  refB <- reference[,(seqB*2 - 1):(seqB*2)]
  colnames(refA) <- c("start.a", "end.a")
  colnames(refB) <- c("start.b", "end.b")


  refA <- refA%>%
    mutate(strand.a = ifelse(start.a < 0, "-", "+"))%>%
    mutate(start.a = abs(start.a))%>%
    mutate(end.a = abs(end.a))

  refB <- refB%>%
    mutate(strand.b = ifelse(start.b < 0, "-", "+"))%>%
    mutate(start.b = abs(start.b))%>%
    mutate(end.b = abs(end.b))


  ref <- refA%>%bind_cols(refB)

  if(remove.missing == T){
    ref <- ref%>%filter((start.a - end.a) != 0, (start.b - end.b) != 0)
  }
  ref <- ref%>%arrange(start.a)



  if(collapse.alignment == T){
    ref <- collapseAlignment(ref = ref, quiet = quiet)
  }

  ref <- ref%>%mutate(diff = start.a - start.b)

  return(ref)

}
