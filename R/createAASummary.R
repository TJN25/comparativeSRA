#' createAASummary Function
#'
#' @param aa The amino acid frequency in the bacterial genomes
#' @param cl Codons mapping to amino acids
#' @keywords reformat the phage data
#' @export
#' @examples
#' createAASummary()

createAASummary <- function(aa = aminoAcidFreq, cl = codonLookup){
  colnames(cl) <- c("codon", "amino.acid")
  aa <- aa%>%
    filter(Genome != "Genome")%>%
    filter(Genome != "A")%>%select(-X.1)

  aminoSummary <- data.frame(amino.acid = colnames(aa)[2:ncol(aa)],
                             freq = rep(NA, ncol(aa) - 1))


  for(i in 2:ncol(aa)){
    aminoSummary[(i - 1), 2] <- mean(as.numeric(aa[,i]))
  }

  aminoSummary <- aminoSummary%>%arrange(-freq)

  aminoSummary <- aminoSummary%>%full_join(cl, by = "amino.acid")
  return(aminoSummary)
}
