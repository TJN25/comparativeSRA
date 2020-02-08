#' interset_all Function
#'
#' @param set_val The file reference positons
#' @keywords takes mutiple set values and finds the interrsect of all of them
#' @export
#' @examples
#' mergeSRA()


interset_all <- function(set_val){
  if(is_empty(set_val)){
    return("0")
  }
  setList <- list()
  for(j in 1:length(set_val)){
    value <- set_val[j]
    setList[j] <- strsplit(as.character(value), "-")
  }

  set_val <- paste(Reduce(intersect, setList), collapse = "-")
  return(set_val)
}
