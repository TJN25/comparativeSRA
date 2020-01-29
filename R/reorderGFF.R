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
    start_val <- gff[i, 3]
    end_val <- gff[i, 4]
    for(j in 1:nrow(ref2)){
      if(start_val >= ref2[j,4] && start_val <= ref2[j,5]){
        gff[i,3] <- gff[i,3] + ref2[j, 7]
        gff[i,4] <- gff[i,4] + ref2[j, 7]
        gff[i, 12] <- T
      }
    }

  }
}else{
  for(i in 1:nrow(gff)){
    if(quiet == F){
      printRemaining(i = i, length = nrow(gff), increment = 5)
    }
    start_val <- gff[i, 3]
    end_val <- gff[i, 4]
    for(j in 1:nrow(ref2)){
      if(start_val >= ref2[j,1] && start_val <= ref2[j,2]){
        gff[i, 12] <- T
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



