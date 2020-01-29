#' plotCodons Function
#'
#' @param cc The codon dataset
#' @param output plot or write the graph
#' @keywords reformat the phage data
#' @export
#' @examples
#' plotCodons()

plotCodons <- function(cc = codons, output = "plot"){
  codons_score <- cc%>%select(score, codon)%>%unique()


  p <- ggplot(data = cc,aes(x = codon, y = log_odds)) +
    geom_boxplot( outlier.shape = NA) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))# +
    #geom_text(data = codons_score, aes(label = round(score), x = codon, y = (2.5 + score/50)))

  if(output == "plot"){
    p
  }else if(output == "write"){

    print("Saving png.")
    ggsave(filename = paste(input_out_name, "_codons.png", sep = ""), plot = p, device = "png", width = 15, height = 10)
  }else{
    print("Error: Select 'plot' or 'write' for the output option" )
  }
}
