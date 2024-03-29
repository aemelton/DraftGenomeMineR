### AE Melton, 2020
# Extract promoter sequences for analyses.

GetPromoterSequences <- function(orfs.report,
                                 blast.out,
                                 pattern.to.keep,
                                 scaffold.fasta,
                                 promoter.csv.file.out,
                                 promoter.sequence.fasta,
                                 promoters.folder,
                                 perc.ident.threshold){
  
### Get just the orfs that match the gene  
orf.blast.out <- read.csv(blast.out)
orf.blast.out.filt <- filter(orf.blast.out, orf.blast.out$Perc.Ident >= perc.ident.threshold)
orf.blast.out.filt.keep <- str_remove_all(string = orf.blast.out.filt$SubjectID, pattern = pattern.to.keep)
#

#
csv.full <- read.csv(file = orfs.report, sep = ",")
csv <- csv.full[csv.full$ORFID == orf.blast.out.filt.keep,]
csvsub <- csv.full %>%
  filter(csv.full$ORFID %in% orf.blast.out.filt.keep)
csv <- csvsub[which(as.numeric(csvsub$start) >= 1500),]

#

#Create data frame and populate it
Promoters <- data.frame(matrix(ncol=6, nrow=nrow(csv)))
colnames(Promoters) <- c("ORFid", "ScaffoldID","Strand", "Start","End", "Sequence")

Promoters$ORFid <- csv$ORFID
Promoters$End <- csv$start
Promoters$ScaffoldID <- as.vector(csv$scaffoldID)
Promoters$Strand <- as.vector(csv$strand)
Promoters$Start <- as.numeric(Promoters$End) - 1500
Promoters$Start[Promoters$start < 0] <- 1

#Read FASTA file (line by line)
scaffold <- readLines(scaffold.fasta)

#Extract sequences to be mined for promoter sequences using PALACE
setwd(promoters.folder)
for(i in 1:nrow(Promoters)){
  #Add DNA sequence adapted to strand
  #Extract seq from FASTA file
  seqRaw <- scaffold[c(grep(paste0(">", Promoters$ScaffoldID[i]), scaffold)+1)]
  
  #Package sequence and extract start and end
  if(Promoters$Strand[i] == "+"){
    Promoters$Sequence[i] <- paste(strsplit(seqRaw, split='')[[1]][as.numeric(Promoters$Start[i]):as.numeric(Promoters$End[i])], collapse = "")
  }
  if(Promoters$Strand[i] == "-"){
    revComp <- as.vector(reverseComplement(DNAStringSet(seqRaw)))
    Promoters$Sequence[i] <- paste(strsplit(revComp, split='')[[1]][as.numeric(Promoters$Start[i]):as.numeric(Promoters$End[i])], collapse = "")
  }
}
#

#
write.csv(Promoters, promoter.csv.file.out, row.names = F, quote = F)
#

#
csv <- Promoters
seqType <- "DNA"
FASTA <- NULL
for(i in 1:nrow(csv)){
  if(seqType == "DNA"){
    DNA <- paste(paste(">", csv$ScaffoldID[i], "_", csv$Gene_hypothesis[i], sep=''), csv$Sequence[i], sep='\n')
    FASTA <- rbind(FASTA, DNA) 
  }
}
#

#
write.table(FASTA, file = promoter.sequence.fasta, col.names = F, row.names = F, quote = F)
#

}
