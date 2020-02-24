#' collapseAlignment Function
#'
#' @param ref The file reference positons
#' @keywords Takes a reference file with two sequences and merges any rows that contain an ongoing alignment
#' @export
#' @examples
#' collapseAlignment()

collapseAlignment <- function(ref){
  ##being a new start point
  new.seq = T
  dat <- ref[0,]
  i <- 0

  i <- i + 1
  ##loop through the rows and combine the ones that can be
  for(i in 1:nrow(ref)){

    ##if there is no current sequence start a new one
    if(new.seq == T){
      startA <- ref[i,1]
      startB <- ref[i,4]
      endA <- ref[i,2]
      endB <- ref[i,5]
      strandA <- ref[i,3]
      strandB <- ref[i,6]
      new.seq <- F
      ##if this is the final row then write it to the file
      if(i == nrow(ref)){
        tmp <- data.frame(start.a = startA, end.a = endA, strand.a = strandA,
                          start.b = startB, end.b = endB, strand.b = strandB)
        suppressWarnings(dat <- dat%>%bind_rows(tmp))
      }
    }else{

      ##check that the two sequences are the same length
      if((ref[i,2] - ref[i,1]) ==(ref[i,5] - ref[i,4])){

        ##check that the gap between merged sequences is 1
        if((ref[i,1] - endA) == 1 && (ref[i,4] - endB) == 1 ){
          endA <- ref[i,2]
          endB <- ref[i,5]

          ##if this is the final row then write it to the file
          if(i == nrow(ref)){
            tmp <- data.frame(start.a = startA, end.a = endA, strand.a = strandA,
                              start.b = startB, end.b = endB, strand.b = strandB)
            suppressWarnings(dat <- dat%>%bind_rows(tmp))
          }

        }else{
          ##if the gap differs then write the new segment to the file and begin a new start point
          tmp <- data.frame(start.a = startA, end.a = endA, strand.a = strandA,
                            start.b = startB, end.b = endB, strand.b = strandB)
          suppressWarnings(dat <- dat%>%bind_rows(tmp))

          startA <- ref[i,1]
          startB <- ref[i,4]
          endA <- ref[i,2]
          endB <- ref[i,5]
          strandA <- ref[i,3]
          strandB <- ref[i,6]
          new.seq <- F
          if(i == nrow(ref)){
            tmp <- data.frame(start.a = startA, end.a = endA, strand.a = strandA,
                              start.b = startB, end.b = endB, strand.b = strandB)
            suppressWarnings(dat <- dat%>%bind_rows(tmp))
          }
        }
      }else{
        ##if the sequences differ write the new semgent to the file and the current row to a file
        ##and set the new.seq value to TRUE
        tmp <- data.frame(start.a = startA, end.a = endA, strand.a = strandA,
                          start.b = startB, end.b = endB, strand.b = strandB)
        suppressWarnings(dat <- dat%>%bind_rows(tmp))
        new.seq <- T
        dat <- dat%>%bind_rows(ref[i,])
      }

    }





  }

  return(dat)
}
