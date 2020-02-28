#' mergeSRAFast Function
#'
#' @param ref The file reference positons
#' @param gff1 The gff file which is the reference gff
#' @param gff2 the gff file that needs rearranging
#' @param filenum1 number or name for the first gff file
#' @param filenum2 number or name for the second gff file
#' @keywords combines multiple gff files from different genomes
#' @export
#' @examples
#' mergeSRAFast()


mergeSRAFast <- function(ncRNAgff, time.it = T, quiet = F, filenum1 = "1", filenum2 = "2", initial_data = T, align = T, no_merge =F){
  # Setup and tests ---------------------------------------------------------
  test_setup <- F
  if(test_setup == T){
    load(file = "~/bin/r_git/R/mergeSRAData.Rda")
    ncRNAgff <- mergeSRAData[["ncRNAgff"]]
    filenum1 <- mergeSRAData[["filenum1"]]
    filenum2 <- mergeSRAData[["filenum2"]]
    initial_data <- mergeSRAData[["initial_data"]]
    align <- mergeSRAData[["align"]]
    quiet <- T
    time.it <- T
    i <- 0
    test_setup <- T
    no_merge <- F
  }

  ptm <- proc.time()

  # Data frame setup ----------------

  ncRNAgff <- ncRNAgff%>%mutate(change = ifelse(start < end, F, T))%>%
    mutate(start.tmp = end)%>%
    mutate(end.tmp = start)%>%
    mutate(start = ifelse(change == T, start.tmp, start))%>%
    mutate(end = ifelse(change == T, end.tmp, end))%>%
    select(-start.tmp, -end.tmp, -change)

  ncRNAgff <- ncRNAgff%>%arrange(start)

  counter = 10
  for(i in 1:nrow(ncRNAgff)){
    if(ncRNAgff$start[i] < 0 && align == T){
      counter <- counter + 1
    }else if(i == nrow(ncRNAgff)){
      counter = counter + 1
    }else{
    if(ncRNAgff$end[i] < ncRNAgff$start[i + 1]){

      counter <- counter + 1
    }
    }
  }





  ##Data will be written into this format
  mergedDat <- data.frame(sequence = rep(as.character("0"), counter), feature = rep(as.character("0"), counter),
                          start = rep(as.integer("0"), counter), end = rep(as.integer("0"), counter),
                          strand = rep(as.character("0"), counter), file_names = rep(as.character("start_row"), counter),
                          row_numbers = rep(as.character("0"), counter), prop_overlap = rep(as.numeric(0), counter), feature_match = rep(F, counter),
                          number_of_features = rep(as.integer("0"), counter),
                          score = rep(as.character("0"), counter),
                          new_feature = rep(F, counter),
                          number_of_rnaseq_files = rep(as.integer("0"), counter),
                          id1 = rep(as.character("0"), counter),
                          id2 = rep(as.character("0"), counter),
                          set_val_1 = rep(as.character("0"), counter),
                          set_val_2 = rep(as.character("0"), counter),
                          file_id = rep(as.character("0"), counter),
                          stringsAsFactors = F)

  ##loop through the combined gff files and combine features that overlap



  current_feature <- F #is there a current feature being written?
  new_feature <- F
  counter <- 0
  #i <-175
  # Main loop ---------------------------------------------------------------

i <- 1
  for(i in 1:(nrow(ncRNAgff))){
    if(quiet ==F){
      printRemaining(i <- i, length = nrow(ncRNAgff) - 1, increment = 5)
    }

    ##write the feature without checks if it could not be remapped
    ##This only applies if theree is a need to align
    if(ncRNAgff$start[i] < 0 && align == T){

      ##define the start and end of the feature
      start_val <- ncRNAgff$start[i]
      start_i <- i
      end_val <- ncRNAgff$end[i]

      feature_matched <- F ##this is always false here as no other features have been compared
      prop_val <- NA ##prop val is always 1 if there is no feature matched

      #select the current feature from the dataframe
      idRows <- ncRNAgff[i,]

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
      counter <- counter + 1
      mergedDat$sequence[counter] <- ncRNAgff$sequence[i]
      mergedDat$feature[counter] <- ncRNAgff$feature[i]
      mergedDat$start[counter] <- start_val
      mergedDat$end[counter] <- end_val
      mergedDat$strand[counter] <- ncRNAgff$strand[i]
      mergedDat$file_names[counter] <- paste(unique(ncRNAgff$file_id[start_i:i]), collapse = ",")
      mergedDat$row_numbers[counter] <- as.character(paste(c(start_i:i), collapse = ","))
      mergedDat$prop_overlap[counter] <- prop_val
      mergedDat$feature_match[counter] <- feature_matched
      mergedDat$number_of_features[counter] <- number_of_features
      mergedDat$score[counter] <- as.character(ncRNAgff$score[i])
      mergedDat$new_feature[counter] <- !(F %in% ncRNAgff$new_feature[start_i:i])
      mergedDat$number_of_rnaseq_files[counter] <- sum(as.integer(ncRNAgff$number_of_rnaseq_files[start_i:i]))
      mergedDat$id1[counter] <- id1_val
      mergedDat$id2[counter] <- id2_val
      mergedDat$set_val_1[counter] <- as.character(set_val_1)
      mergedDat$set_val_2[counter] <- as.character(set_val_2)
      mergedDat$file_id[counter] <- paste(filenum1, "-", filenum2, sep = "")



      ##this is the end of the loop writing a feature that cannot be mapped
      next
    }



    ##if there is no current feature then set a new start value

    if(current_feature == F){
      start_val <- ncRNAgff$start[i]
      start_i <- i
      end_val <- ncRNAgff$end[i]
    }else{

    ##check if the existing values are outside the range of the initial values
    if(ncRNAgff$end[i] > end_val){
      end_val <- ncRNAgff$end[i]
    }
    }
    ##check if the current end value overlaps with the next starting value and update the end value if it does
    if(i == nrow(ncRNAgff)){

      ##check if the subsequent feature was contained within the first feature
      if(ncRNAgff$end[start_i] < end_val){
        prop_val <- (ncRNAgff$end[start_i] - ncRNAgff$start[start_i])/(ncRNAgff$end[i] - ncRNAgff$start[i])
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



      counter <- counter + 1
      mergedDat$sequence[counter] <- ncRNAgff$sequence[i]
      mergedDat$feature[counter] <- ncRNAgff$feature[i]
      mergedDat$start[counter] <- start_val
      mergedDat$end[counter] <- end_val
      mergedDat$strand[counter] <- ncRNAgff$strand[i]
      mergedDat$file_names[counter] <- paste(unique(ncRNAgff$file_id[start_i:i]), collapse = ",")
      mergedDat$row_numbers[counter] <- as.character(paste(c(start_i:i), collapse = ","))
      mergedDat$prop_overlap[counter] <- prop_val
      mergedDat$feature_match[counter] <- feature_matched
      mergedDat$number_of_features[counter] <- number_of_features
      mergedDat$score[counter] <- as.character(ncRNAgff$score[i])
      mergedDat$new_feature[counter] <- !(F %in% ncRNAgff$new_feature[start_i:i])
      mergedDat$number_of_rnaseq_files[counter] <- sum(as.integer(ncRNAgff$number_of_rnaseq_files[start_i:i]))
      mergedDat$id1[counter] <- id1_val
      mergedDat$id2[counter] <- id2_val
      mergedDat$set_val_1[counter] <- as.character(set_val_1)
      mergedDat$set_val_2[counter] <- as.character(set_val_2)
      mergedDat$file_id[counter] <- paste(filenum1, "-", filenum2, sep = "")

      ##feature has been added so start again
      current_feature <- F
      new_feature <- F

    }else{
      if(end_val > ncRNAgff$start[i + 1]){ #& ncRNAgff$strand[i] == ncRNAgff$strand[i+1]){
        if(ncRNAgff$end[i + 1] > end_val){
          end_val <- ncRNAgff$end[i + 1]
        }
        ##now looking at an ongoing feature
        current_feature <- T
      }else{ ##if the next feature is not overlapping, start writing out the existing featur

        ##check if the subsequent feature was contained within the first feature
        if(ncRNAgff$end[start_i] < end_val){
          prop_val <- (ncRNAgff$end[start_i] - ncRNAgff$start[start_i])/(ncRNAgff$end[i] - ncRNAgff$start[i])
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


        counter <- counter + 1
        mergedDat$sequence[counter] <- ncRNAgff$sequence[i]
        mergedDat$feature[counter] <- ncRNAgff$feature[i]
        mergedDat$start[counter] <- start_val
        mergedDat$end[counter] <- end_val
        mergedDat$strand[counter] <- ncRNAgff$strand[i]
        mergedDat$file_names[counter] <- paste(unique(ncRNAgff$file_id[start_i:i]), collapse = ",")
        mergedDat$row_numbers[counter] <- as.character(paste(c(start_i:i), collapse = ","))
        mergedDat$prop_overlap[counter] <- prop_val
        mergedDat$feature_match[counter] <- feature_matched
        mergedDat$number_of_features[counter] <- number_of_features
        mergedDat$score[counter] <- as.character(ncRNAgff$score[i])
        mergedDat$new_feature[counter] <- !(F %in% ncRNAgff$new_feature[start_i:i])
        mergedDat$number_of_rnaseq_files[counter] <- sum(as.integer(ncRNAgff$number_of_rnaseq_files[start_i:i]))
        mergedDat$id1[counter] <- id1_val
        mergedDat$id2[counter] <- id2_val
        mergedDat$set_val_1[counter] <- as.character(set_val_1)
        mergedDat$set_val_2[counter] <- as.character(set_val_2)
        mergedDat$file_id[counter] <- paste(filenum1, "-", filenum2, sep = "")

        ##feature has been added so start again
        current_feature <- F
        new_feature <- F
      }
    }

  }

  # Tidy up and finish ------------------------------------------------------

  mergedDat <- mergedDat%>%filter(number_of_features > 0, file_names != "start_row")

  mergedDat <- mergedDat%>%mutate(change = ifelse(start < end, F, T))%>%
    mutate(start.tmp = end)%>%
    mutate(end.tmp = start)%>%
    mutate(start = ifelse(change == T, start.tmp, start))%>%
    mutate(end = ifelse(change == T, end.tmp, end))%>%
    select(-start.tmp, -end.tmp, -change)

  mergedDat <- mergedDat%>%filter(!is.na(sequence))
  mergedDat[is.na(mergedDat)] <- 0

  runningTime <- proc.time() - ptm
  if(time.it){
    if(quiet == F){
      printRunningTime(runningTime = runningTime)
    }
  }



  return(mergedDat)

}
