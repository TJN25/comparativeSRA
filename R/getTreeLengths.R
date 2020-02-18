#' getTreeLengths Function
#'
#' @param treeLengths list containing tree lengths
#' @param treeDat data frame with node numbers and lengths
#' @param name name of branch for extracting from list
#' @keywords gets tip.label lengths
#' @export
#' @examples
#' getTreeLengths()

getTreeLengths <- function(treeLengths, treeDat, name){
  positions <- treeLengths[[name]]$positions
  next_row <- treeDat %>% filter(X2 == positions[1])
  next_row <- as.character(next_row[1,])
  treeLengths[[name]]$positions <- next_row[1:2]
  treeLengths[[name]]$distances <- c(treeLengths[[name]]$distances, as.numeric(next_row[3]))
  if(treeLengths[[name]]$positions[1] %in% treeDat$X2){
    treeLengths <- getTreeLengths(treeLengths, treeDat, name)
  }
  return(treeLengths)
}
