#' bacterialDinucleotideFormat Function
#'
#' @param bac The bacterial file formatted by bacterialColumnNames
#' @keywords reformat the bacterial data
#' @export
#' @examples
#' bacterialDinucleotideFormat()


bacterialDinucleotideFormat <- function(bac = bacterial){

  lookup_list = c("AA", "AC", "GT", "AG", "CC", "TT", "CG", "GG", "GC",
                  "AT", "GA", "TG", "CT", "CA", "TC", "TA")

  test_list <- match(lookup_list, colnames(bac))

  if(anyNA(test_list) == FALSE){
    bacterialLength <- nrow(bac)
    bacDinucleotide <- data.frame(genera = rep(NA, 1), dinucleotide = rep(NA, 1), prop = rep(NA, 1))
    for(i in 69:84){
      genera = bac[,2]
      dinucleotide = rep(colnames(bac)[i], bacterialLength)
      prop = bac[,i]

      df <- data.frame(genera = genera, dinucleotide = dinucleotide, prop = prop)

      bacDinucleotide <- bacDinucleotide%>%bind_rows(df)
    }
    bacDinucleotide <- bacDinucleotide%>%
      mutate(type = rep("bacteria", nrow(bacDinucleotide)))%>%
      filter(!is.na(dinucleotide))
    bacDinucleotide <- bacDinucleotide%>%group_by(dinucleotide)%>%summarise(bac_prop = mean(prop))

    return(bacDinucleotide)
  }else{
    NAindex <- which(is.na(test_list))
    missingAA <- lookup_list[NAindex]
    print("Error: Columns are missing")
    print(missingAA)
  }



}
