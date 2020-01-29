#' reorderSraPlot Function
#'
#' @param reference The reference sra file
#' @param sra the sra file to rearrange
#' @param alignment the alignment file
#' @keywords Reorder the sra plot file to allow mutplie genomes to be compared
#' @export
#' @examples
#' reorderSraPlot()

reorderSraPlot <- function(sra, alignment, reference, seqA = 1, seqB = 2, quiet = F, time.it = T){
  ptm <- proc.time()
if(!is.numeric(seqA)){
  cat(paste("SeqA changed from ", seqA, " to 1 as ", seqA, " is not numeric\n", sep = ""))
  seqA <- 1
}

  if(!is.numeric(seqB)){
    cat(paste("SeqB changed from ", seqB, " to 2 as ", seqB, " is not numeric\n", sep = ""))
    seqB <- 2
  }


colnames(sra) <- c("a", "b")
sra <- sra%>%mutate(nucleotide = row_number())%>%mutate(nucleotide.adj = NA)

alignment <- alignment[,c((seqA*2 - 1):(seqA*2), (seqB*2 - 1):(seqB*2))]
colnames(alignment) <- c("seq0_leftend", "seq0_rightend", "seq1_leftend", "seq1_rightend")


alignment <- alignment%>%mutate(strand = ifelse(seq0_leftend < 0, "-", "+"))%>%
  mutate(seq0_leftend = abs(seq0_leftend))%>%
  mutate(seq0_rightend = abs(seq0_rightend))%>%
  mutate(diff = seq0_leftend - seq1_leftend)

# i <- 6
for(i in 1:nrow(alignment)){
  #print(i)
  #print(alignment[i,3])
  if(quiet == F){
    printRemaining(i = i, length = nrow(alignment), increment = 5)
  }


  if(alignment[i,1] == 0 && alignment[i,2] == 0){

    sra[alignment[i,3]:alignment[i,4],4] <-  rep(0, length(alignment[i,3]:alignment[i,4]))

  }else if(alignment[i,3] != 0 || alignment[i,4] != 0){
    sra[alignment[i,3]:alignment[i,4],4] <-  sra[alignment[i,3]:alignment[i,4],3] + alignment[i,6]
   }#else{
  #   len = length(alignment[i,1]:alignment[i,2])
  #   tmp <- data.frame(a = rep(0, len), b = rep(0, len), nucleotide = rep(0, len), nucleotide.adj = alignment[i,1]:alignment[i,2])
  #   sra <- sra%>%bind_rows(tmp)
  # }

}


sra <- sra%>%arrange(nucleotide.adj)
sra <- sra%>%filter(nucleotide.adj != 0)

genomeFill <- data.frame(nucleotide.adj = c(1:nrow(reference)))

sraAll <- sra%>%full_join(genomeFill)
sraAll <- sraAll%>%arrange(nucleotide.adj)%>%
  mutate(a = ifelse(is.na(a), 0, a))%>%
  mutate(b = ifelse(is.na(b), 0, b))
sraAll <- sraAll%>%group_by(nucleotide.adj)%>%arrange(-a)%>%arrange(-b)%>%mutate(order.num = row_number())
sraAll <- sraAll%>%filter(order.num == 1)%>%arrange(nucleotide.adj)%>%ungroup()

runningTime <- proc.time() - ptm
if(time.it){
  if(quiet == F){
    printRunningTime(runningTime = runningTime)
  }
}

return(sraAll)
}
