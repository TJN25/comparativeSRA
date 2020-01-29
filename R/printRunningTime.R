#' printRunningTime Function
#'
#' @param i the position in the loop
#' @param length the length of the loop
#' @param increment how often to print the output
#' @keywords print the remaining number of iterations in a loop
#' @export
#' @examples
#' printRunningTime()


printRunningTime <- function(runningTime, type = "The function"){
  if(runningTime[3] < 60){
    cat(paste(type," ran for ",  ceiling(as.numeric(runningTime[3]) %% 60), " seconds\n", sep = ""))

  }else{
  cat(paste(type," ran for ", round(as.numeric(runningTime[3])/60), " minutes and ", round(as.numeric(runningTime[3]) %% 60), " seconds\n", sep = ""))
  }
}
