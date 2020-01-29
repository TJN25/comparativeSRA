#' dinucleotideSummary Function
#'
#' @param bac The bacterial dataset formated by bacterialDinucleotideFormat
#' @param ph The phage dataset formated by phageDinucleotideFormat
#' @keywords reformat the phage data
#' @export
#' @examples
#' dinucleotideSummary()


dinucleotideSummary <- function(ph = phageDinucleotides, bac = bacDinucleotides){

  test_list <- c(match("dinucleotide", colnames(ph)),match("dinucleotide", colnames(bac)))

  if(anyNA(test_list) == FALSE){

    dinucleotide <- ph%>%left_join(bac, by = "dinucleotide")
    dinucleotide <- dinucleotide%>%mutate(log_odds = log2(prop/bac_prop))


    logOddsSummary <- dinucleotide%>%select(dinucleotide, prop, bac_prop, log_odds)%>%unique()%>%
      group_by(dinucleotide)%>%summarise(score = sum(abs(log_odds)))

        dinucleotide <- dinucleotide%>%left_join(logOddsSummary, by = "dinucleotide")


    return(dinucleotide)
  }else{
    print("Error: Files not formated correctly. The 'dinucleotide' column does not exist.")
  }
}




