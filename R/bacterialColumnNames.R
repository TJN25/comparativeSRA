#' bacterialColumnNames Function
#'
#' @param bac The input bacterial file
#' @keywords reformat the bacterial data
#' @export
#' @examples
#' bacterialColumnNames()
bacterialColumnNames <- function(bac = bacterial){

  lookup_list = c("CTT", "ATG", "AAG", "AAA", "ATC", "AAC", "ATA", "AGG", "CCT", "ACT", "AGC", "ACA",
                  "AGA", "CAT", "AAT", "ATT", "CTG", "CTA", "CTC", "CAC", "ACG", "CAA", "AGT", "CAG",
                  "CCG", "CCC", "TAT", "GGT", "TGT", "CGA", "CCA", "TCT", "GAT", "CGG", "TTT", "TGC",
                  "GGG", "TAG", "GGA", "TAA", "GGC", "TAC", "TTC", "TCG", "TTA", "TTG", "TCC", "GAA",
                  "TCA", "GCA", "GTA", "GCC", "GTC", "GCG", "GTG", "GAG", "GTT", "GCT", "ACC", "TGA",
                  "GAC", "CGT", "TGG", "CGC", "AA", "AC", "GT", "AG", "CC", "TT", "CG", "GG", "GC",
                  "AT", "GA", "TG", "CT", "CA", "TC", "TA")

  colname_index = seq(5, ncol(bac)-1, by=2)
  colname_list = as.character(bac[1,colname_index])

  test_list <- match(lookup_list, colname_list)
  if(anyNA(test_list) == FALSE){
    bac <- bac[,-colname_index]
    colnames(bac) <- c("file_path", "genome", "tRNA_gc", "CDS_gc", colname_list)
    bac <- bac[,!is.na(colnames(bac))]
    bac <- bac%>%mutate(dn_sum = sum(AA:TA))
    for(i in 1:nrow(bac)){
    bac[i,85] <- sum(bac[i,69:84])
    }
    for(i in 69:84){
      bac[,i] <- bac[,i]/bac[85]
    }
    for(i in 1:nrow(bac)){
      bac[i,85] <- sum(bac[i,69:84])
    }
    return(bac)
  }else{
    NAindex <- which(is.na(test_list))
    missingAA <- lookup_list[NAindex]
    print("Error: Columns are missing")
    print(missingAA)
  }

}
