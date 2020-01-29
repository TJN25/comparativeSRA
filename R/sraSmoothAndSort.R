#' sraSmoothandsort Function
#'
#' @param dat The plot file to smooth and sort
#' @param col.num 1 for forward strand and 2 for reverse strand
#' @param binwidth the width to smooth over
#' @param moveby the jumps to make for each smoothing bin
#' @keywords smooth out the plot file, and find ranges where expression is occuring.
#' @export
#' @examples
#' sraSmoothandsort()


sraSmoothandsort <- function(dat, col.num = 1, binwidth = 25, moveby = 5){

  if(col.num == 1){
    TS <- zoo::zoo(dat$V1)
  }else if(col.num == 2){
    TS <- zoo::zoo(dat$V2)
  }else{
    stop("Please select column 1 or 2 with col.num")

  }
tmp <- zoo::rollapply(TS, width = binwidth, by = moveby, FUN = mean, align = "left") ##these values might need to be changed

tmp2 <- as.data.frame(tmp)


df <- data.frame(start = 0, end = 0, median.val = 0, min.value = 0, max.value = 0)
start_i <- 1
start_range <- 1
for(i in 2:nrow(tmp2)){

  if(tmp2[i,1] >= 5){
    if(tmp2[i - 1,1] < 5){
      start_i <- i
      start_range <- row.names(tmp2)[i]
    }
  }else{
    if(tmp2[i - 1,1] >= 5){
      end_range <- row.names(tmp2)[i -1]
      medianValue = median(tmp2[c(start_i:i),1])
      minValue = min(tmp2[c(start_i:i),1])
      maxValue = max(tmp2[c(start_i:i),1])
      df <- rbind(df, data.frame(start = start_range, end = end_range, median.val = medianValue, min.value = minValue, max.value = maxValue))
    }
  }


}

df <- df%>%mutate(length = as.numeric(end) - as.numeric(start))%>%filter(length != 0)
return(df)
}


