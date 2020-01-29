#' removeCDSregions Function
#'
#' @param plotDat The file with the plot values
#' @param gffDat The gff file that contains the CDS locations
#' @keywords Convert a cmscan output file to a gff formatted file. This will not include the top header info.
#' @export
#' @examples
#' removeCDSregions()

removeCDSregions <- function(plotDat, gffDat, buffer = 150, quiet = F, stranded = F, time.it = F){
ptm <- proc.time()
  if(colnames(gffDat)[1] == "V1"){
    if(ncol(gffDat) == 9){
      colnames(gffDat) <- c("sequence", "source", "feature", "start", "end", "score", "strand", "phase", "Atrribute")
    }else{
      stop("GFF file does not appear to contain the correct number of rows")
    }
  }

  if(is.na(match('feature', colnames(gffDat)))){
    stop("GFF file does not appear to contain a feature column")
  }

  if('CDS' %in% gffDat$feature == F){
    stop("No CDS regions to remove")
  }

  gffDat <- gffDat%>%
    filter(feature == "CDS")

  gffDat <- gffDat%>%
    mutate(feature.length = end - start + buffer*2)

  featureTotal <- sum(gffDat$feature.length)

  plotDat <- plotDat%>%
    mutate(nucleotide = row_number())


  codingRegions <- data.frame(nucleotide = NA, strand = NA)
  if(quiet == F){
    cat("Adjusting the plot values\n")
  }


  for(i in 1:nrow(gffDat)){
    if(quiet == F){
      printRemaining(i = i, length = nrow(gffDat))
    }
    codingRegions <- codingRegions%>%bind_rows(data.frame(nucleotide = c(((gffDat[i,4]) - buffer):((gffDat[i,5]) + buffer)),
                                                          strand = gffDat[i,7]))
  }

  codingRegions <- unique(codingRegions)

  codingRegions <- codingRegions%>%mutate(keep = F)

if(stranded == T){
  plotDatRev <- plotDat%>%select(V1, nucleotide)
  plotDatFwd <- plotDat%>%select(V2, nucleotide)

  codingRegionsFwd <- codingRegions%>%filter(strand == "+")
  codingRegionsRev <- codingRegions%>%filter(strand == "-")

  plotDatRev <- plotDatRev%>%
    left_join(codingRegionsRev)%>%
    mutate(keep = ifelse(is.na(keep), T, keep))

  plotDatFwd <- plotDatFwd%>%
    left_join(codingRegionsFwd)%>%
    mutate(keep = ifelse(is.na(keep), T, keep))

  plotDatRev <- plotDatRev%>%
    mutate(V1 = ifelse(keep == F & strand == "-", 0, V1))

  plotDatFwd <- plotDatFwd%>%
    mutate(V2 = ifelse(keep == F& strand == "+", 0, V2))
  plotDatRev <- plotDatRev%>%select(V1, nucleotide)
  plotDatFwd <- plotDatFwd%>%select(V2, nucleotide)

  plotDat <- plotDatFwd%>%full_join(plotDatRev)%>%select(V1,V2,nucleotide)

}else{
  codingRegionsUnstranded <- codingRegions%>%select(nucleotide, keep)%>%unique()

  plotDat <- plotDat%>%
    left_join(codingRegionsUnstranded)%>%
    mutate(keep = ifelse(is.na(keep), T, keep))

  plotDat <- plotDat%>%
    mutate(V1 = ifelse(keep == F, 0, V1))%>%
    mutate(V2 = ifelse(keep == F, 0, V2))%>%
    select(-keep)

}

  if(quiet == F & time.it == T){
    runningTime <- proc.time() - ptm
    printRunningTime(runningTime = runningTime, type = "removeCDSregions")
  }


  return(plotDat)
}
