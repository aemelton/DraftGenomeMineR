### AE Melton, 2020

GetScaffolds <- function(BlastHitFile, genome){

  # Extract scaffolds from draft genome and make a fasta file
  cl.filt.unique <- read.csv(file = BlastHitFile)
  header <- NULL
  seq <- NULL
  #S.start <- NULL
  #S.end <- NULL
  if(nrow(cl.filt.unique) != 1){
    for(i in 1:nrow(cl.filt.unique)){
  
      header[i] <- as.character(cl.filt.unique$SubjectID[i])
      seq[i] <- genome[c(match(paste(">", header[i], sep=''), genome)+1)]
      }
    } else {
      header <- as.character(cl.filt.unique$SubjectID)
      seq <- genome[c(match(paste(">", header, sep=''), genome)+1)]
    }
  x <- dplyr::tibble(name = header, seq = seq)
  return(x)
}
