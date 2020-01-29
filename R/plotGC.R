#' plotGC Function
#'
#' @param ph phage dataset
#' @param tr bacterial tRNAs
#' @param cd bacterial CDS
#' @param output plot or write the graph
#' @keywords reformat the phage data
#' @export
#' @examples
#' plotGC()


plotGC <- function(ph = phage, tr = tRNA, cd = CDS, output = "plot"){
  col_count = ncol(tr)

  if(col_count == 2){
    p <- ggplot() +
      geom_density(data = ph, aes(genome_gc), colour="red") +
      geom_density(data = tr, aes(V2), colour="blue") +
      geom_density(data = cd, aes(V2), colour="green")
  }else{
    p <- ggplot() +
      geom_density(data = ph, aes(genome_gc), colour="red") +
      geom_density(data = tr, aes(V3), colour="blue") +
      geom_density(data = cd, aes(V2), colour="green")
  }

  if(output == "plot"){
    p
  }else if(output == "write"){

    print("Saving png.")
    ggsave(filename = paste(input_out_name, "_gc.png", sep = ""), plot = p, device = "png", width = 15, height = 10)
  }else{
    print("Error: Select 'plot' or 'write' for the output option" )
  }
}




