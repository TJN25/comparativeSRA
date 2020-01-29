#' phageColumnNames Function
#'
#' @param bac The input phage file
#' @keywords reformat the phage data
#' @export
#' @examples
#' phageColumnNames()


phageColumnNames <- function(ph = phage){

  lookup_list = c("CTT", "ATG", "AAG", "AAA", "ATC", "AAC", "ATA", "AGG", "CCT", "ACT", "AGC", "ACA",
                  "AGA", "CAT", "AAT", "ATT", "CTG", "CTA", "CTC", "CAC", "ACG", "CAA", "AGT", "CAG",
                  "CCG", "CCC", "TAT", "GGT", "TGT", "CGA", "CCA", "TCT", "GAT", "CGG", "TTT", "TGC",
                  "GGG", "TAG", "GGA", "TAA", "GGC", "TAC", "TTC", "TCG", "TTA", "TTG", "TCC", "GAA",
                  "TCA", "GCA", "GTA", "GCC", "GTC", "GCG", "GTG", "GAG", "GTT", "GCT", "ACC", "TGA",
                  "GAC", "CGT", "TGG", "CGC", "AA", "AC", "GT", "AG", "CC", "TT", "CG", "GG", "GC",
                  "AT", "GA", "TG", "CT", "CA", "TC", "TA")


  colname_index = seq(3, ncol(ph)-1, by=2)
  colname_list = as.character(ph[1,colname_index])


  test_list <- match(lookup_list, colname_list)
  if(anyNA(test_list) == FALSE){

    ph <- ph[,-c(colname_index, ncol(ph))]
    colnames(ph) <- c("genome_info", "genome_gc", colname_list)
    ph <- ph%>%separate(col = genome_info, into = c("ID", "Genome", "Host"), sep = "\\|", extra = "merge", remove = F)
    ph <- ph%>%separate(col = Genome, into = c("genome_genera", "genome_species", "genome_extra"), sep = " ", extra = "merge", remove = F)
    ph <- ph%>%separate(col = Host, into = c("host_genera", "host_species", "host_extra"), sep = " ", extra = "merge", remove = F)

    ph <- ph[,!is.na(colnames(ph))]
    ph <- ph%>%mutate(dn_sum = sum(AA:TA))
    for(i in 1:nrow(ph)){
      ph[i,92] <- sum(ph[i,76:91])
    }
    for(i in 76:91){
      ph[,i] <- ph[,i]/ph[92]
    }
    for(i in 1:nrow(ph)){
      ph[i,92] <- sum(ph[i,76:91])
    }

    return(ph)

  }else{
    NAindex <- which(is.na(test_list))
    missingAA <- lookup_list[NAindex]
    print("Error: Columns are missing")
    print(missingAA)
  }
}

