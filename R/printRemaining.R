#' printRemaining Function
#'
#' @param i the position in the loop
#' @param length the length of the loop
#' @param increment how often to print the output
#' @keywords print the remaining number of iterations in a loop
#' @export
#' @examples
#' printRemaining()


printRemaining <- function(i, length, increment = 5){

  increment_val <- round(length/100*increment)
  if(increment_val != 0){

  output <- paste(round(i/length*100), "%, ", sep = "")
  if(i == 1){
    cat(paste("Starting loop with ", formatC(length, big.mark = ","), " iterations:\n", sep = ""))
  }
  if(i %% increment_val == 0){
    cat(output)
  }
  }else{
    if(i == 1){
      cat(paste("Starting loop with ", formatC(length, big.mark = ","), " iterations:\n", sep = ""))
    }
    output <- paste(round(i/length*100), "%, ", sep = "")
    cat(output)

  }
  if(i == length){
    cat("Finished.\n")
  }
}
