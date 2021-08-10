### AE Melton, 2021
# Mine transcriptomes with trinotate annotations
#setwd("~/Dropbox/BSU_Research/Transcriptomics_HudsonAlpha/Assembled_Transcriptomes/171_Leaves/")

FindTranscriptsOfInterest <- function(FASTA, Trinotation, KeyWord, OutputFile){
# Load the transcriptome. Must be one-line fasta
  transcriptome <- readLines(con = FASTA, encoding = "us-ascii")
#head(transcriptome)

# Read in the trinotate spreadsheet
  annotation <- read.csv(Trinotation, sep = "\t", header = T)
#head(annotation)

# Now, umm, if a line has a keyword, like, say for example... Aquaporin! Extract the scaffold.
# We will need to make a list of scaffold IDs to extract.
  tmp <- NULL
  transcripts2extract <- NULL

  for(i in 1:ncol(annotation)){
    tmp <- annotation[c(grepl(pattern = KeyWord, x = annotation[,i], ignore.case = T)),]
    transcripts2extract <- rbind(transcripts2extract, tmp)
    }
#head(transcripts2extract)
#nrow(transcripts2extract)
  out <- unique(transcripts2extract)
#nrow(out)
  write.csv(x = out, OutputFile, row.names = F)
}