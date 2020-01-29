#' phageDinucleotideFormat Function
#'
#' @param ph The phage file formatted by phageColumnNames
#' @keywords reformat the phage data
#' @export
#' @examples
#' phageDinucleotideFormat()

phageDinucleotideFormat <- function(ph = phage){

  lookup_list = c("AA", "AC", "GT", "AG", "CC", "TT", "CG", "GG", "GC",
                   "AT", "GA", "TG", "CT", "CA", "TC", "TA")

  test_list <- match(lookup_list, colnames(ph))

  if(anyNA(test_list) == FALSE){
    phageLength <- nrow(ph)

    ##format the codons for plotting
    phageDN <- data.frame(genera = rep(NA, 1), dinucleotide = rep(NA, 1), prop = rep(NA, 1))

    for(i in 76:91){
      genera = ph[,2]
      dinucleotide = rep(colnames(ph)[i], phageLength)
      prop = ph[,i]

      df <- data.frame(genera = genera, dinucleotide = dinucleotide, prop = prop)

      phageDN <- phageDN%>%bind_rows(df)

    }

    phageDN <- phageDN%>%
      mutate(type = rep("phage", nrow(phageDN)))%>%
      filter(!is.na(dinucleotide))
    return(phageDN)
  }else{
    NAindex <- which(is.na(test_list))
    missingAA <- lookup_list[NAindex]
    print("Error: Columns are missing")
    print(missingAA)
  }



}


