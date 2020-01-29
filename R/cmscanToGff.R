#' cmscanToGff Function
#'
#' @param rfamRes The cmscan output file
#' @keywords Convert a cmscan output file to a gff formatted file. This will not include the top header info.
#' @export
#' @examples
#' cmscanToGff()


cmscanToGff <- function(rfamRes){
  colnames(rfamRes) <- c("idx", "target.name", "accession1", "query.name", "accession2", "clan.name", "mdl", "mdl.from",  "mdl.to", "seq.from",   "seq.to",
                         "strand", "trunc", "pass",   "gc",  "bias",  "score",   "E.value", "inc", "olp", "anyidx", "afrct1", "afrct2", "winidx", "wfrct1", "wfrct2", "description.of.target")

  rfamRes <- rfamRes%>%
    mutate(seq.from2 = seq.from)%>%
    mutate(seq.from = ifelse(seq.from > seq.to, seq.to, seq.from))%>%
    mutate(seq.to = ifelse(seq.from2 > seq.to, seq.from2, seq.to))

  gff <- data.frame(seqname = rfamRes$query.name,
                    source = rep("rfam", nrow(rfamRes)),
                    feature = rep("ncRNA", nrow(rfamRes)),
                    start = rfamRes$seq.from,
                    end = rfamRes$seq.to,
                    score = rfamRes$score,
                    strand = rfamRes$strand,
                    frame = rep(".", nrow(rfamRes)),
                    attribute = rfamRes$target.name)
  return(gff)
}

