#' mauvePlotColourChange Function
#'
#' @param bbone the genoPlotR mauve file
#' @keywords change colors on the bbone file
#' @export
#' @examples
#' mauvePlotColourChange()

mauvePlotColourChange <- function(bbone){
  for(g in 1:length(bbone$dna_segs)){
    lengths <- abs(bbone$dna_segs[[g]]$end - bbone$dna_segs[[g]]$start)
    lenDf <- data.frame(lengths = lengths)
    lenDf <- lenDf %>% mutate(row.num = row_number()) %>% arrange(lengths)

    lenDf$colour <- ""

    for(i in 1:nrow(lenDf)){
      lenDf$colour[i] <- topo.colors(length(lengths))[i]
    }
    lenDf <- lenDf %>% arrange(row.num)
    bbone$dna_segs[[g]]$col <- lenDf$colour
    }

  return(bbone)

}
