#' codonSummary Function
#'
#' @param bac The bacterial dataset formated by bacterialCodonsFormat
#' @param ph The phage dataset formated by phageCodonsFormat
#' @param codon_order Determines whether to reorder the codons by amino acid frequency or not
#' @keywords reformat the phage data
#' @export
#' @examples
#' codonSummary()


codonSummary <- function(ph = phageCodons, bac = bacCodons, codon_order = F){

  test_list <- c(match("codon", colnames(ph)),match("codon", colnames(bac)))

  if(anyNA(test_list) == FALSE){

  codons <- ph%>%left_join(bac, by = "codon")
  codons <- codons%>%mutate(log_odds = log2(prop/bac_prop))


  logOddsSummary <- codons%>%select(codon, prop, bac_prop, log_odds)%>%unique()%>%
    group_by(codon)%>%summarise(score = sum(abs(log_odds)))


  if(codon_order == T){
  codons <- codons%>%left_join(aminoSummary, by = "codon")

  codons <- codons%>%arrange(-freq)

  }
  codons <- codons
  codons <- codons%>%left_join(logOddsSummary, by = "codon")

  if(codon_order ==T){
  codons$codon <- factor(codons$codon, levels = unique(codons$codon))
  }
  return(codons)
  }else{
    print("Error: Files not formated correctly. The 'codon' column does not exist.")
  }
}




