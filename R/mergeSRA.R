#' mergeSRA Function
#'
#' @param ref The file reference positons
#' @param gff1 The gff file which is the reference gff
#' @param gff2 the gff file that needs rearranging
#' @param filenum1 number or name for the first gff file
#' @param filenum2 number or name for the second gff file
#' @keywords combines multiple gff files from different genomes
#' @export
#' @examples
#' mergeSRA()


mergeSRA <- function(ncRNAgff, gff1, gff2, time.it = T, quiet = F, filenum1 = "1", filenum2 = "2", initial_data = T, align = T){
  # Setup and tests ---------------------------------------------------------

  error_message <- "Either gff1 and gff2 or ncRNAgff are needed:\n"
  stop_val <- 0
  log_file = ""
  if(missing(gff1)){
  error_message <- paste(error_message, "\tArugment gff1 missing\n", sep = "")
  stop_val <- stop_val + 1
  }else{
    gff1 <- gff1%>%mutate(file_id = filenum1)
  }



  if(missing(gff2)){
    error_message <- paste(error_message, "\tArugment gff2 missing\n", sep = "")
    stop_val <- stop_val + 1
    }else{
      gff2b <- gff2%>%mutate(file_id = filenum2)
    }


  if(stop_val > 0){
    if(missing(ncRNAgff)){
      error_message <- paste(error_message, "\tArugment ncRNAgff missing\n", sep = "")
      stop(error_message)
    }
  }
  if(missing(ncRNAgff)){
    ncRNAgff <- gff2b%>%bind_rows(gff1)%>%unique()
  }else{
    if(quiet == F){
    cat("Using the ncRNAgff dataframe:\n")
    }
  }

  ptm <- proc.time()



  # Data frame setup ----------------


  ncRNAgff <- ncRNAgff%>%arrange(start) %>% arrange(strand)

  ##Data will be written into this format
  mergedDat <- data.frame(sequence = as.character("0"), feature = as.character("0"),
                          start = as.integer("0"), end = as.integer("0"),
                          strand = as.character("0"), file_names = as.character("start_row"),
                          row_numbers = as.character("0"), prop_overlap = as.numeric(0), feature_match = F,
                          number_of_features = as.integer("0"),
                          score = as.character("0"),
                          new_feature = F,
                          number_of_rnaseq_files = as.integer("0"),
                          id1 = as.character("0"),
                          id2 = as.character("0"),
                          set_val_1 = as.character("0"),
                          set_val_2 = as.character("0"),
                          file_id = as.character("0"),
                          stringsAsFactors = F)

  ##loop through the combined gff files and combine features that overlap
  #i <- 278
  current_feature <- F #is there a current feature being written?
  new_feature <- F


  # Main loop ---------------------------------------------------------------


  for(i in 1:(nrow(ncRNAgff) - 1)){

    if(quiet ==F){
    printRemaining(i <- i, length = nrow(ncRNAgff) - 1, increment = 5)
    }

    ##write the feature without checks if it could not be remapped
    ##This only applies if theree is a need to align
    if(ncRNAgff$start[i] < 0 && align == T){

      ##define the start and end of the feature
      start_val <- ncRNAgff$start[i]
      start_i <- i ##this will not change and could be removed but might cause problems if it is still called later
      end_val <- ncRNAgff$end[i]

      feature_matched <- F ##this is always false here as no other features have been compared
      prop_val <- 1 ##prop val is always 1 if there is no feature matched

      #select the current feature from the dataframe
      idRows <- ncRNAgff[start_i:i,]

      ##keep the ids from rows of the dataframe that match the filenum (file name)
      ##The purpose of this is to keep track of each original file and row number for each feature as continued merging is done
      id1_val <- paste(idRows$id[idRows$file_id == filenum1], collapse = "-")
      id2_val <- paste(idRows$id[idRows$file_id == filenum2], collapse = "-")

      ##assign each missing file as an id of filenum_0
      ##this needs to be done for each individual file so that it appears as filenum1_0-filenum2_0 not filenum1-filenum2_0
      ##as this will make tracking which files contributed easier
      if(id1_val == ""){
        id_vals <- unlist(strsplit(filenum1, "-"))

        id_vals <- paste(id_vals, "_0", sep = "")

        id1_val <- paste(id_vals, collapse = "-")
      }

      if(id2_val == ""){
        id_vals <- unlist(strsplit(filenum2, "-"))

        id_vals <- paste(id_vals, "_0", sep = "")

        id2_val <- paste(id_vals, collapse = "-")
      }

      ##take the current set value and if one does not exist then set to 1
      set_val_1 <- ifelse(filenum1 == ncRNAgff$file_id[i], ncRNAgff$set_val[i], 0)
      set_val_2 <- ifelse(filenum2 == ncRNAgff$file_id[i], ncRNAgff$set_val[i], 0)

      ##Number of features
      ##Combine each set of ids and get a unique list of each id available
      id_val <- paste(id1_val, id2_val, sep = "-")
      idList <- unique(unlist(strsplit(id_val, "-")))
      ##loop through each ID and check if it is real or a place holder. All real IDs are 1 and placeholders are 0
      for(j in 1:length(idList)){
        number <- unlist(strsplit(as.character(idList[j]), "_"))
        if(length(number)  < 3){
          idList[j] <- as.numeric(0)
        }else{
          number <- number[3]
          if(number == "0"){
            idList[j] <- as.numeric(0)
          }else{
            idList[j] <- as.numeric(1)
          }
        }
      }
      number_of_features <- sum(as.numeric(idList))

      ##get the set values from old features
      if(initial_data == F){
        idSetValues <- ncRNAgff[start_i:i,c(17:(ncol(ncRNAgff)))]
        for(j in 1:ncol(idSetValues)){
          idSetValues[1,j] <- interset_all(idSetValues[idSetValues[,j]!= "0",j])

        }
        idSetValues <- idSetValues[1,]
      }

      tmp <- data.frame(sequence = ncRNAgff$sequence[i],
                        feature = ncRNAgff$feature[i],
                        start = start_val, end = end_val,
                        strand = ncRNAgff$strand[i],
                        file_names = paste(unique(ncRNAgff$file_id[start_i:i]), collapse = ","),
                        row_numbers = as.character(paste(c(start_i:i), collapse = ",")),
                        prop_overlap = prop_val,
                        feature_match = feature_matched,
                        number_of_features = number_of_features,
                        score = as.character(ncRNAgff$score[i]),
                        new_feature = !(F %in% ncRNAgff$new_feature[start_i:i]),
                        number_of_rnaseq_files = sum(as.integer(ncRNAgff$number_of_rnaseq_files[start_i:i])),
                        id1 = id1_val,
                        id2 = id2_val,
                        set_val_1 = as.character(set_val_1),
                        set_val_2 = as.character(set_val_2),
                        file_id = paste(filenum1, "-", filenum2, sep = ""),
                        stringsAsFactors = F)

      if(initial_data == F){
        tmp <- tmp %>% bind_cols(idSetValues)
      }
      mergedDat <- mergedDat%>%bind_rows(tmp)

  ##this is the end of the loop writing a feature that cannot be mapped
  next
    }



    ##if there is no current feature then set a new start value

    if(current_feature == F){
      start_val <- ncRNAgff$start[i]
      start_i <- i
      end_val <- ncRNAgff$end[i]
    }

    ##check if the existing values are outside the range of the initial values
    if(ncRNAgff$end[i] > end_val){
      end_val <- ncRNAgff$end[i]
    }

    ##check if the current end value overlaps with the next starting value and update the end value if it does
    if(end_val > ncRNAgff$start[i + 1] & ncRNAgff$strand[i] == ncRNAgff$strand[i+1]){
      if(ncRNAgff$end[i + 1] > end_val){
        end_val <- ncRNAgff$end[i + 1]
      }
      ##now looking at an ongoing feature
      current_feature <- T
    }else{ ##if the next feature is not overlapping, start writing out the existing featur

      ##check if the subsequent feature was contained within the first feature
      if(ncRNAgff$end[start_i] < end_val){
        prop_val <- (ncRNAgff$end[start_i] - ncRNAgff$start[i])/(end_val - start_val)
      }else{
        prop_val <- 1
      }

      ##check if there was more than one feature overlapping
      feature_matched <- ifelse(length(unique(ncRNAgff$file_id[start_i:i])) > 1, T, F)


      #select the current feature from the dataframe
      idRows <- ncRNAgff[start_i:i,]

      ##keep the ids from rows of the dataframe that match the filenum (file name)
      ##The purpose of this is to keep track of each original file and row number for each feature as continued merging is done
      id1_val <- paste(idRows$id[idRows$file_id == filenum1], collapse = "-")
      id2_val <- paste(idRows$id[idRows$file_id == filenum2], collapse = "-")

      ##assign each missing file as an id of filenum_0
      ##this needs to be done for each individual file so that it appears as filenum1_0-filenum2_0 not filenum1-filenum2_0
      ##as this will make tracking which files contributed easier
      if(id1_val == ""){
        id_vals <- unlist(strsplit(filenum1, "-"))
        id_vals <- paste(id_vals, "_0", sep = "")
        id1_val <- paste(id_vals, collapse = "-")
      }

      if(id2_val == ""){
        id_vals <- unlist(strsplit(filenum2, "-"))
        id_vals <- paste(id_vals, "_0", sep = "")
        id2_val <- paste(id_vals, collapse = "-")
      }

      ##take the current set value and if one does not exist then set to 1
      set_val_1 <- idRows$set_val[idRows$file_id == filenum1]
      set_val_2 <- idRows$set_val[idRows$file_id == filenum2]

      set_val_1 <- interset_all(set_val_1)
      set_val_2 <- interset_all(set_val_2)


      if(set_val_1 == ""){
        set_val_1 <- "0"
      }
      if(set_val_2 == ""){
        set_val_2 <- "0"
      }


      ##Number of features
      ##Combine each set of ids and get a unique list of each id available
      id_val <- paste(id1_val, id2_val, sep = "-")
      idList <- unique(unlist(strsplit(id_val, "-")))
      ##loop through each ID and check if it is real or a place holder. All real IDs are 1 and placeholders are 0
      for(j in 1:length(idList)){
        number <- unlist(strsplit(as.character(idList[j]), "_"))
        if(length(number)  < 3){
          idList[j] <- as.numeric(0)
        }else{
          number <- number[3]
          if(number == "0"){
            idList[j] <- as.numeric(0)
          }else{
            idList[j] <- as.numeric(1)
          }
        }
      }
      number_of_features <- sum(as.numeric(idList))

      ##get the set values from old features
      if(initial_data == F){
        idSetValues <- ncRNAgff[start_i:i,c(17:(ncol(ncRNAgff)))]
        for(j in 1:ncol(idSetValues)){
          idSetValues[1,j] <- interset_all(idSetValues[idSetValues[,j]!= "0",j])

        }
        idSetValues <- idSetValues[1,]
      }

        tmp <- data.frame(sequence = ncRNAgff$sequence[i],
                        feature = ncRNAgff$feature[i],
                        start = start_val, end = end_val,
                        strand = ncRNAgff$strand[i],
                        file_names = paste(unique(ncRNAgff$file_id[start_i:i]), collapse = ","),
                        row_numbers = as.character(paste(c(start_i:i), collapse = ",")),
                        prop_overlap = prop_val,
                        feature_match = feature_matched,
                        number_of_features = number_of_features,
                        score = as.character(ncRNAgff$score[i]),
                        new_feature = !(F %in% ncRNAgff$new_feature[start_i:i]),
                        number_of_rnaseq_files = sum(as.integer(ncRNAgff$number_of_rnaseq_files[start_i:i])),
                        id1 = id1_val,
                        id2 = id2_val,
                        set_val_1 = as.character(set_val_1),
                        set_val_2 = as.character(set_val_2),
                        file_id = paste(filenum1, "-", filenum2, sep = ""),
                        stringsAsFactors = F)
        if(initial_data == F){
          tmp <- tmp %>% bind_cols(idSetValues)
        }
        mergedDat <- mergedDat%>%bind_rows(tmp)

      ##feature has been added so start again
      current_feature <- F
      new_feature <- F
    }
  }

  # Tidy up and finish ------------------------------------------------------


  runningTime <- proc.time() - ptm
  if(time.it){
    if(quiet == F){
      printRunningTime(runningTime = runningTime)
    }
  }



  mergedDat <- mergedDat%>%filter(number_of_features > 0, file_names != "start_row")
  return(mergedDat)

}
