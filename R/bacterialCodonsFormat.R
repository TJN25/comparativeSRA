#' bacterialCodonsFormat Function
#'
#' @param bac The bacterial file formatted by bacterialColumnNames
#' @keywords reformat the bacterial data
#' @export
#' @examples
#' bacterialCodonsFormat()
bacterialCodonsFormat <- function(bac = bacterial){

  lookup_list = c("CTT", "ATG", "AAG", "AAA", "ATC", "AAC", "ATA", "AGG", "CCT", "ACT", "AGC", "ACA",
                  "AGA", "CAT", "AAT", "ATT", "CTG", "CTA", "CTC", "CAC", "ACG", "CAA", "AGT", "CAG",
                  "CCG", "CCC", "TAT", "GGT", "TGT", "CGA", "CCA", "TCT", "GAT", "CGG", "TTT", "TGC",
                  "GGG", "TAG", "GGA", "TAA", "GGC", "TAC", "TTC", "TCG", "TTA", "TTG", "TCC", "GAA",
                  "TCA", "GCA", "GTA", "GCC", "GTC", "GCG", "GTG", "GAG", "GTT", "GCT", "ACC", "TGA",
                  "GAC", "CGT", "TGG", "CGC")

  test_list <- match(lookup_list, colnames(bac))

  if(anyNA(test_list) == FALSE){
    bacterialLength <- nrow(bac)
    bacCodons <- data.frame(genera = rep(NA, 1), codon = rep(NA, 1), prop = rep(NA, 1))
    for(i in 5:68){
      genera = bac[,2]
      codon = rep(colnames(bac)[i], bacterialLength)
      prop = bac[,i]

      df <- data.frame(genera = genera, codon = codon, prop = prop)

      bacCodons <- bacCodons%>%bind_rows(df)
    }
    bacCodons <- bacCodons%>%
      mutate(type = rep("bacteria", nrow(bacCodons)))%>%
      filter(!is.na(codon))
    bacSummary <- bacCodons%>%group_by(codon)%>%summarise(bac_prop = mean(prop))

    return(bacSummary)
  }else{
    NAindex <- which(is.na(test_list))
    missingAA <- lookup_list[NAindex]
    print("Error: Columns are missing")
    print(missingAA)
  }



}
