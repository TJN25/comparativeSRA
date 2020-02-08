#' scoreRNAReadDepth Function
#'
#' @param file_path location of all the files to work with
#' @param gff_names location of the gff files relative to the file path
#' @param sra the gff file to use
#' @param allData the file to append the data to
#' @keywords take in a list of calls and gets the read depths from all the RNASeq files
#' @export
#' @examples
#' scoreRNAReadDepth()


scoreRNAReadDepth <- function(file_path, gff_names, sra, allData){
  gffDat <- read.table(paste(file_path, "/", sra, "_new_calls.txt", sep = ""), comment.char = "#", quote = "", sep = "\t", as.is = T, header = T)

  if(dir.exists(paste(file_path, "/plot_files/", sep = ""))){
    rnaFiles <- list.files( paste(file_path, "/plot_files/", sep = ""), pattern = "_ncRNA.plot$")
    rnaData <- read.table(paste(file_path, "/plot_files/", rnaFiles[1], sep = ""), sep = "")

  }else{
    rnaFiles <- list.files( paste(file_path, "/", sep = ""), pattern = "_ncRNA.plot$")
    rnaData <- read.table(paste(file_path, "/", rnaFiles[1], sep = ""), sep = "")
  }

  randomDat <- read.table(paste(file_path, "/random_sequences/shifted/", sra, "_shifted_random_new_calls.txt", sep = ""), comment.char = "#", quote = "", sep = "\t", as.is = T, header = T)


  rnaData <- rnaData %>% mutate(maxVal = ifelse(V1 > V2, V1, V2)) %>% select(maxVal)
  for(i in 2:length(rnaData)){

    if(dir.exists(paste(file_path, "/plot_files/", sep = ""))){
      tmp <- read.table(paste(file_path, "/plot_files/", rnaFiles[2], sep = ""), sep = "")

    }else{
      tmp <- read.table(paste(file_path, "/", rnaFiles[2], sep = ""), sep = "")
    }
    tmp <- tmp %>% mutate(maxVal = ifelse(V1 > V2, V1, V2)) %>% select(maxVal)

    rnaData <- rnaData %>% bind_cols(tmp)
  }

  rnaTotal <- sum(unlist(rnaData[,]))
  rnaTotal <- rnaTotal/1000000

  randomDat <- randomDat %>% mutate(score_2 = 0)
  i <- 1
  for(i in 1:nrow(randomDat)){
    randomDat$score_2[i] <- max(unlist(rnaData[randomDat$start[i]:randomDat$end[i],], use.names = F))/rnaTotal
  }

  for(i in 1:nrow(gffDat)){
    gffDat$score_2[i] <- max(unlist(rnaData[gffDat$start[i]:gffDat$end[i],], use.names = F))/rnaTotal
  }

  d1 <- gffDat %>%filter(new_feature == F, grepl(pattern = "sra_calls", x = file_names)) %>% mutate(group = "1") %>% select(group, score_2, id)
  d2 <- gffDat%>% filter(new_feature == T) %>%  mutate(group = "2") %>% select(group, score_2, id)
  d3 <- randomDat %>% mutate(group = "3")%>% select(group, score_2, id)

  if(missing("allData")){
    allData <- d1 %>% bind_rows(d2) %>% bind_rows(d3)
  }else{
    allData <- allData %>% bind_rows(d1) %>% bind_rows(d2) %>% bind_rows(d3)
  }
  return(allData)
}
