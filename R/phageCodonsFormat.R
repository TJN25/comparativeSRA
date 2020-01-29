#' phageCodonsFormat Function
#'
#' @param ph The phage file formatted by phageColumnNames
#' @keywords reformat the phage data
#' @export
#' @examples
#' phageCodonsFormat()

phageCodonsFormat <- function(ph = phage){

  lookup_list = c("CTT", "ATG", "AAG", "AAA", "ATC", "AAC", "ATA", "AGG", "CCT", "ACT", "AGC", "ACA",
                  "AGA", "CAT", "AAT", "ATT", "CTG", "CTA", "CTC", "CAC", "ACG", "CAA", "AGT", "CAG",
                  "CCG", "CCC", "TAT", "GGT", "TGT", "CGA", "CCA", "TCT", "GAT", "CGG", "TTT", "TGC",
                  "GGG", "TAG", "GGA", "TAA", "GGC", "TAC", "TTC", "TCG", "TTA", "TTG", "TCC", "GAA",
                  "TCA", "GCA", "GTA", "GCC", "GTC", "GCG", "GTG", "GAG", "GTT", "GCT", "ACC", "TGA",
                  "GAC", "CGT", "TGG", "CGC")

  test_list <- match(lookup_list, colnames(ph))

  if(anyNA(test_list) == FALSE){
  phageLength <- nrow(ph)

  ##format the codons for plotting
  phageCodons <- data.frame(genera = rep(NA, 1), codon = rep(NA, 1), prop = rep(NA, 1))


  for(i in 12:75){
    genera = ph[,2]
    codon = rep(colnames(ph)[i], phageLength)
    prop = ph[,i]

    df <- data.frame(genera = genera, codon = codon, prop = prop)

    phageCodons <- phageCodons%>%bind_rows(df)

  }

  phageCodons <- phageCodons%>%
    mutate(type = rep("phage", nrow(phageCodons)))%>%
    filter(!is.na(codon))
  return(phageCodons)
  }else{
    NAindex <- which(is.na(test_list))
    missingAA <- lookup_list[NAindex]
    print("Error: Columns are missing")
    print(missingAA)
  }



}


