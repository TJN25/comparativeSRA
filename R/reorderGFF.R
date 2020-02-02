#' reorderGFF Function
#'
#' @param ref The file reference positons
#' @param gff the gff file that needs rearranging
#' @keywords relabels positions for a gff file
#' @export
#' @examples
#' reorderGFF()


reorderGFF <- function(ref, gff, time.it = T, quiet = F, reference.genome = F){
  ref2 <- ref%>%arrange(start.b)

  ptm <- proc.time()

if(reference.genome == F){
  for(i in 1:nrow(gff)){
    if(quiet == F){
    printRemaining(i = i, length = nrow(gff), increment = 5)
    }
    start_val <- gff$start[i]
    end_val <- gff$end[i]
    for(j in 1:nrow(ref2)){
      if(start_val >= ref2$start.b[j] && start_val <= ref2$end.b[j]){
        gff$start[i] <- gff$start[i] + ref2$diff[j]
        gff$end[i] <- gff$end[i] + ref2$diff[j]
        gff$changed[i] <- T
      }
    }

  }
}else{
  for(i in 1:nrow(gff)){
    if(quiet == F){
      printRemaining(i = i, length = nrow(gff), increment = 5)
    }
    start_val <- gff$start[i]
    end_val <- gff$end[i]
    for(j in 1:nrow(ref2)){
      if(start_val >= ref2$start.a[j] && start_val <= ref2$end.a[j]){
        gff$changed[i] <- T
      }
    }

  }
}

  runningTime <- proc.time() - ptm
  if(time.it){
    if(quiet == F){
      printRunningTime(runningTime = runningTime)
    }
  }
  gff <- gff%>%filter(changed == T)
  return(gff)
}



