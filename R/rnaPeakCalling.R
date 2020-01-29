#' rnaPeakCalling Function
#'
#' @param dat The plot file to smooth and sort
#' @param col.num 1 for forward strand and 2 for reverse strand
#' @param binwidth the width to smooth over
#' @param moveby the jumps to make for each smoothing bin
#' @keywords smooth out the plot file, and find ranges where expression is occuring.
#' @export
#' @examples
#' rnaPeakCalling()


rnaPeakCalling <- function(dat, col.num = 1, quiet = F, time.it = T, plot_threshold = 10,
                           small_peaks = T, binwidth = 25){

  ##convert to tpm
  moveby <- 1
  ##Load the input values
  headerDat <- tryCatch({
    ##Change this path or put the header file in the working directory
    suppressWarnings(df <- read.table("./rnaPeakCallingHeader.txt", as.is = T))

    df
  }, error =  function(e) {
    df <- data.frame(V1 = c("plot_threshold", "plot_threshold_feature", "feature_separation_distance", "feature_length", "r_squared"),
                     V2 = as.numeric(c("10", "0", "50", "50", "0.5")))
    if(quiet == F){
      cat("Using Default Values\n")

    }
    df
  } )

  # headerDat <- data.frame(V1 = c("plot_threshold", "plot_threshold_feature", "feature_separation_distance", "feature_length", "r_squared"),
  #                                      V2 = as.numeric(c("10", "0", "15", "50", "0.5")))

  #plot_threshold <- as.numeric(headerDat$V2[1])
  plot_threshold_feature <- as.numeric(headerDat$V2[2])
  feature_separation_distance <- as.numeric(headerDat$V2[3])
  feature_length <- as.numeric(headerDat$V2[4])
  r_squared <- as.numeric(headerDat$V2[5])

  if(quiet == F){
    cat(paste("Using values: \n",
              "\tplot_threshold ", " = ", plot_threshold, "\n",
              "\tplot_threshold_feature ", " = ", plot_threshold_feature, "\n",
              "\tfeature_separation_distance ", " = ", feature_separation_distance, "\n",
              "\tfeature_length ", " = ", feature_length, "\n",
              "\tr_squared ", " = ", r_squared, "\n", sep = ""))
  }
  dat <- dat%>%
    mutate(nucleotide = row_number())
  ##smooth the data
  if(binwidth > 1){
  if(quiet == F){
    cat(paste("Smoothing data points:\n",
    "\tbinwidth = ", binwidth,
    "\n", sep = ""))
  }
    ptm1 <- proc.time()
    TSRev <- zoo::zoo(dat$V1)

    TSFwd <- zoo::zoo(dat$V2)

  tmpRev <- zoo::rollapply(TSRev, width = binwidth, by = moveby, FUN = mean, align = "left")
  tmpFwd <- zoo::rollapply(TSFwd, width = binwidth, by = moveby, FUN = mean, align = "left")

  datRev <- as.data.frame(tmpRev)
  datFwd <- as.data.frame(tmpFwd)
  datRev <- datRev%>%
    mutate(nucleotide = row_number())
  datFwd <- datFwd%>%
    mutate(nucleotide = row_number())

  datRev <- datRev%>%rename(V1 = tmpRev)
  datFwd <- datFwd%>%rename(V2 = tmpFwd)

  dat <- datRev%>%left_join(datFwd, by = "nucleotide")%>%select(V1, V2, nucleotide)

  if(time.it == T & quiet == F){
    runningTime <- proc.time() - ptm1
      printRunningTime(runningTime = runningTime, type = "The smoothing")
  }

  }else{
    moveby <- 1
  }
  ##add nucleotide infor to the data frame


  ##store original data and keep an editable version of the data
  datTmp <- dat
  callsDat <- data.frame(start = NA, stop = NA, mean = NA)

  ##loop through the process j times
  ptm <- proc.time()

  ##loop through the data table and identify regions above a threshold value
  peak <- F
  feature_found <- F
  for(i in 1:nrow(datTmp)){

    if(quiet == F){
        printRemaining(i = i, length = nrow(datTmp), increment = 5)
    }


    ##check if the value can be ignored
    if(peak == F){
      if(datTmp[i, col.num] < plot_threshold){
        next
      }
    }


    ##get the value for this nucleotide
    mean_val_short <- datTmp[i,col.num]
    if(is.na(mean_val_short)){
      mean_val_short <- 0
    }

    ##check the value is above the cutoff threshold
    if(mean_val_short >= plot_threshold){
      ##if there is not already a feature open, begin a new feature
      if(peak == F){
        start <- i
        peak <- T
        feature_found <- T
      }

      ##if the value is less than the threshold, close and write the new feature
    }else if(mean_val_short < plot_threshold){

      stop <- i
      df <- data.frame(start = start, stop = stop, mean = mean(datTmp[start:stop,col.num]))
      callsDat <- callsDat%>%bind_rows(df)



      peak <- F
    }

  }
  runningTime <- proc.time() - ptm
  if(quiet == F & time.it == T){
    printRunningTime(runningTime = runningTime, type = "The sra peak calling")
  }


  callsCombined <- data.frame(start = 0, stop = 0, mean.score = 0)
  if(feature_found == T){
    ptm1 <- proc.time()
     ##sort the features by start site so that the iterations of the loop are in the correct order for the next step
    # print(nrow(callsDat))
    # print(callsDat)
    if(nrow(callsDat) > 1){


    callsDat <- callsDat%>%filter(!is.na(start))%>%arrange(start)



    #print("callsDat")
    #print(callsDat[1:10,])


    ##loop through the feature calls and combine the features that were close enough together
    peak <- F
    for(i in 1:nrow(callsDat)){
      # print(i)
      # print(nrow(callsDat))

       if(quiet == F){
        printRemaining(i = i, length = nrow(callsDat), increment = 1)
      }
      #print(nrow(callsDat))
      #print(i)
      ##if there is no current feature open a new feature
      if(peak == F){
        start_val <- callsDat[i,1]
        end_val <- callsDat[i,2]
        mean_val <- callsDat[i,3]

        peak <- T
      }else{

        ##This will write a feature when the loop ends if one is open
        if(i == nrow(callsDat)){
          if(abs(end_val - callsDat[i,1]) < feature_separation_distance + (callsDat[i,2] - callsDat[i,1])/50){
            end_val <- callsDat[i,2]
          }
          if(abs(end_val - start_val) >= feature_length){
            meanScore <- mean(dat[start_val:end_val,col.num])
            df <- data.frame(start = start_val, stop = end_val, mean.score = meanScore) #>>>>
            callsCombined <- callsCombined%>%bind_rows(df)
          }
          peak <- F

        }else{

          ##check if the start of the next feature is within 15 nt of the end of the previous feature. If so, change the end site of the feature
          if(abs(end_val - callsDat[i,1]) < feature_separation_distance){
            end_val <- callsDat[i,2]
          }else{
            ##check the length of the feature is more than 50 nt and write the feature
            if(abs(end_val - start_val) >= feature_length){
              meanScore <- mean(dat[start_val:end_val,col.num])
              df <- data.frame(start = start_val, stop = end_val, mean.score = meanScore)
              callsCombined <- callsCombined%>%bind_rows(df)


              ##find peaks within this feature
              mean_val <- mean(dat[start_val:end_val,col.num])
              featureDat <- dat[start_val:end_val,]
              featureDat[,col.num] <- featureDat[,col.num] - mean_val

              #plot(featureDat$V1 ~ featureDat$nucleotide)
              #lines(y = rep(0, length = length(featureDat$nucleotide)), x =featureDat$nucleotide)
              if(small_peaks == T){
                if(end_val - start_val >= 1000){
                  model <- lm(featureDat$V1 ~ poly(featureDat$nucleotide, 2))
                  fit <- summary(model)

                  if(is.na(fit$r.squared)){
                    next
                  }else  if(fit$r.squared < r_squared){
                    #print("featureDat")
                    #print(featureDat[1:10,])
                    df <- rnaPeakCalling(dat = featureDat, col.num = col.num, quiet = T, time.it = F,
                                         plot_threshold = mean_val, small_peaks = F, binwidth = 1)

                    if(nrow(df) > 1){
                      #print("df")
                      #print(df)
                      #print(featureDat)
                      for(k in 2:nrow(df)){
                        df[k,1] <- featureDat[df[k,1],3]
                        df[k,2] <- featureDat[df[k,2],3]
                      }
                      df <- df%>%select(-feature.length)

                      featureCount <- nrow(df)

                      callsCombined <- callsCombined%>%bind_rows(df)
                    }


                  }
                }
              }
            }
            start_val <- callsDat[i,1]
            end_val <- callsDat[i,2]
            mean_val <- callsDat[i,3]
          }
        }
      }

    }



    callsCombined <- callsCombined%>%filter(!is.na(start))%>%mutate(feature.length = stop - start)
    runningTime <- proc.time() - ptm1
    if(quiet == F & time.it == T){
      printRunningTime(runningTime = runningTime, type = "Combining features")
    }
    }
  }
  cat(paste(nrow(callsCombined) - 1, "features found.\n"))
  return(callsCombined)



}
