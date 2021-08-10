### AE Melton, 2021
# Mine transcriptomes with trinotate annotations
ExtractTranscriptsOfInterest <- function(FASTA, TranscriptsOfInterest, output.filename){
  # Load the transcriptome. Must be one-line fasta
  transcriptome <- readLines(con = FASTA, encoding = "us-ascii")
  head(transcriptome)
  
  # Read in the trinotate spreadsheet
  df <- read.csv(TranscriptsOfInterest, header = T)
  #head(annotation)
  
  # Now, umm, if a line has a keyword, like, say for example... Aquaporin! Extract the scaffold.
  # We will need to make a list of scaffold IDs to extract.
  header <- NULL
  seq <- NULL
  
  for(i in 1:nrow(df)){
    
    #Add DNA sequence adapted to strand
    #Extract seq from FASTA file
    #header[i] <- paste0(">", cl.filt.unique$SubjectID[i])
    header[i] <- as.character(df$transcript_id[i])
    seq[i] <- transcriptome[c(grep(paste(">", df$transcript_id[i], sep=''), transcriptome)+1, fixed = T)]
    
    #S.start[i] <- cl.filt.unique$S.start[i]
    #S.end[i] <- cl.filt.unique$S.end[i]
    
    #Package sequence and extract start and end
    #header[i] <- paste0(">", cl$SubjectID[i])
    #seq[i] <- paste(strsplit(seqRaw, split='')[[1]][as.numeric(start):as.numeric(end)], collapse = "")
  }
  x <- dplyr::tibble(name = header, seq = seq) # This will assemble the headers and sequences into an object that
  # can be written into a FASTA format file
  writeFasta(data = x, filename = output.filename)
}