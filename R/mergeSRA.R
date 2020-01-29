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


mergeSRA <- function(ncRNAgff, gff1, gff2, time.it = T, quiet = F, filenum1 = "1", filenum2 = "2", print_log = F){

  ##setup and tests
   error_message <- "Either gff1 and gff2 or ncRNAgff are needed:\n"
  stop_val <- 0
  log_file = ""
  if(missing(gff1)){
  error_message <- paste(error_message, "\tArugment gff1 missing\n", sep = "")
  stop_val <- stop_val + 1
  }else{
    gff1 <- gff1%>%mutate(filenum = filenum1)
  }

if(missing(gff2)){
  error_message <- paste(error_message, "\tArugment gff2 missing\n", sep = "")
  stop_val <- stop_val + 1
  }else{
    gff2b <- gff2%>%mutate(filenum = filenum2)
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


##sort through the features and merge overlapping features
ncRNAgff <- ncRNAgff%>%arrange(start)


mergedDat <- data.frame(sequence = as.character("0"), feature = as.character("0"),
                        start = as.integer("0"), end = as.integer("0"),
                        strand = as.character("0"), file_names = as.character("start_row"),
                        row_numbers = as.character("0"), prop_overlap = as.numeric(0), feature_match = F,
                        number_of_features = as.integer("0"),
                        score = as.character("0"),
                        new_feature = F,
                        number_of_rnaseq_files = as.integer("0"),
                         id = as.character("0"),
                         set_val = as.character("0"),
                        file_id = as.character("0"),
                        stringsAsFactors = F)

##loop through the combined gff files and combine features that overlap
i <- 892
current_feature <- F #is there a current feature being written?
new_feature <- F
for(i in 1:(nrow(ncRNAgff) - 1)){

  if(quiet ==F){
  printRemaining(i <- i, length = nrow(ncRNAgff) - 1, increment = 5)
  }

  ##write the feature without checks if it could not be remapped
  if(ncRNAgff[i,3] < 0){
    start_val <- ncRNAgff[i,3]
    start_i <- i
    end_val <- ncRNAgff[i,4]

    feature_matched <- ifelse(length(unique(ncRNAgff[start_i:i, 14])) > 1, T, F)
    prop_val <- 1
      idRows <- ncRNAgff[start_i:i,]
      id1_val <- idRows[idRows[,14] == filenum1,13]
      id2_val <- idRows[idRows[,14] == filenum2,13]
      if(is_empty(id1_val)){
        id1_val <- paste(filenum1, "_0", sep = "")
      }
      if(is_empty(id2_val)){
        id2_val <- paste(filenum2, "_0", sep = "")
      }
      set_val1 <- ifelse(filenum1 == ncRNAgff$filenum[i], 1, 0)

      set_val2 <- ifelse(filenum2 == ncRNAgff$filenum[i], 1, 0)
    tmp <- data.frame(sequence = ncRNAgff[i,1],
                      feature = ncRNAgff[i,2],
                      start = start_val, end = end_val,
                      strand = ncRNAgff[i,5],
                      file_names = paste(unique(ncRNAgff[start_i:i, 14]), collapse = ","),
                      row_numbers = paste(c(start_i:i), collapse = ","),
                      prop_overlap = prop_val,
                      feature_match = feature_matched,
                      number_of_features = length(start_i:i),
                      score = as.character(ncRNAgff[i,11]),
                      new_feature = !(F %in% ncRNAgff[start_i:i, 9]),
                      number_of_rnaseq_files = sum(as.integer(ncRNAgff[start_i:i, 10])),
                      id = paste(id1_val, "-", id2_val, sep = ""),
                      set_val = paste(unique(c(set_val1, set_val2)), collapse = "-"),
                      file_id = paste(filenum1, "-", filenum2, sep = ""),
                      stringsAsFactors = F)
    mergedDat <- mergedDat%>%bind_rows(tmp)


next
  }



  ##if there is no current feature then set a new start value

  if(current_feature == F){
    start_val <- ncRNAgff[i,3]
    start_i <- i
    end_val <- ncRNAgff[i,4]
  }



  ##set the new end value
  if(ncRNAgff[i, 4] > end_val){
    end_val <- ncRNAgff[i,4]
  }

  ##check if the current end value overlaps with the next starting value and update the end value if it does
  if(end_val > ncRNAgff[i + 1, 3] & ncRNAgff[i,5] == ncRNAgff[i+1, 5]){
    if(ncRNAgff[i + 1, 4] > end_val){
      end_val <- ncRNAgff[i + 1,4]
    }
    current_feature <- T
  }else{

    ##check if the subsequent feature was contained within the first feature
    if(ncRNAgff[start_i, 4] < end_val){
      prop_val <- (ncRNAgff[start_i, 4] - ncRNAgff[i, 3])/(end_val - start_val)
    }else{
      prop_val <- 1
    }
    feature_matched <- ifelse(length(unique(ncRNAgff[start_i:i, 14])) > 1, T, F)


    #if(length(start_i:i) <= 2){
       idRows <- ncRNAgff[start_i:i,]
       id1_val <- paste(idRows[idRows[,14] == filenum1,13], collapse = "-")
       id2_val <- paste(idRows[idRows[,14] == filenum2,13], collapse = "-")
       if(id1_val == ""){
         id1_val <- paste(filenum1, "_0", sep = "")
       }
       if(id2_val == ""){
         id2_val <- paste(filenum2, "_0", sep = "")
       }
      #id1_val <- "1"
      #id2_val <- "2"


      tmp <- data.frame(sequence = ncRNAgff[i,1],
                      feature = ncRNAgff[i,2],
                      start = start_val, end = end_val,
                      strand = ncRNAgff[i,5],
                      file_names = paste(unique(ncRNAgff[start_i:i, 14]), collapse = ","),
                      row_numbers = paste(c(start_i:i), collapse = ","),
                      prop_overlap = prop_val,
                      feature_match = feature_matched,
                      number_of_features = length(start_i:i),
                      score = as.character(ncRNAgff[i,11]),
                      new_feature = !(F %in% ncRNAgff[start_i:i, 9]),
                      number_of_rnaseq_files = sum(as.integer(ncRNAgff[start_i:i, 10])),
                      id = paste(id1_val, "-", id2_val, sep = ""),
                      set_val = paste(ifelse(id1_val == paste(filenum1, "_0", sep = ""), 0, 1), "-", ifelse(id2_val == paste(filenum2, "_0", sep = ""), 0, 1), sep = ""),
                      file_id = paste(filenum1, "-", filenum2, sep = ""),
                      stringsAsFactors = F)
    mergedDat <- mergedDat%>%bind_rows(tmp)
    #}else{
      #log_file <- paste(log_file, start_i, "to", i, "contains too many peaks. There should be one or two. \nThis is limited by the need for an ID for each peak.\n")
    #}
    current_feature <- F
    new_feature <- F
  }
}

runningTime <- proc.time() - ptm
if(time.it){
  if(quiet == F){
    printRunningTime(runningTime = runningTime)
  }
}

if(print_log){
  cat(log_file)
}

mergedDat <- mergedDat%>%filter(number_of_features > 0, file_names != "start_row")
return(mergedDat)

}
