#' plotDinucleotides Function
#'
#' @param cc The dinucleotide dataset
#' @param output plot or write the graph
#' @keywords reformat the phage data
#' @export
#' @examples
#' plotDinucleotides()

plotDinucleotides <- function(cc = dinucleotides, output = "plot"){
  dinucleotide_score <- cc%>%select(score, dinucleotide)%>%unique()


  p <- ggplot(data = cc,aes(x = dinucleotide, y = log_odds)) +
    geom_boxplot( outlier.shape = NA) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+
    #geom_text(data = dinucleotide_score, aes(label = round(score), x = dinucleotide, y = (2.5 + score/50)))

  if(output == "plot"){
    p
  }else if(output == "write"){

    print("Saving png.")
    ggsave(filename = paste(input_out_name, "_dinucleotides.png", sep = ""), plot = p, device = "png", width = 15, height = 10)
  }else{
    print("Error: Select 'plot' or 'write' for the output option" )
  }
}
