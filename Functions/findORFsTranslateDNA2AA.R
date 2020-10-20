# Sven Buerki
#A function to findORFs in scaffold file (FASTA format) and produce AA sequences as well as report
findORFsTranslateDNA2AA <- function(scaffold, scaffoldID){
  #Extract scaffold sequence (and convert to right format)
  seqRaw <- scaffold[c(grep(paste(">", scaffoldID, sep=''), scaffold)+1)] 
  seqs <- DNAStringSet(seqRaw)
  
  #positive strands 
  pos <- ORFik::findORFs(seqs, startCodon = "ATG", minimumLength = 40) 
  #negative strands (DNAStringSet only if character) 
  
  neg <- ORFik::findORFs(reverseComplement(DNAStringSet(seqs)), startCodon = "ATG", minimumLength = 40) 
  pos <- relist(c(GRanges(pos,strand = "+"),GRanges(neg,strand = "-")),skeleton = merge(pos,neg))
  
  #Process output
  pos <- as.data.frame(pos)
  pos$scaffoldID <- rep(scaffoldID, nrow(pos))
  pos$ORFID <- paste("ORF_", seq(from=1, to=nrow(pos)), sep='')
  
  #Extract and write ORFs from seq and produce AA sequence
  for(i in 1:nrow(pos)){
    if(as.character(pos$strand[i]) == "+"){
      extSeq <- as.DNAbin(DNAStringSet(paste(strsplit(seqRaw, split='')[[1]][c(pos$start[i]:pos$end[i])], collapse = "")))
      AAseq <- paste(paste(">", pos$scaffoldID[i], "_", pos$ORFID[i], sep=''), paste(as.character(trans(extSeq, code = 1, codonstart = 1))[[1]], collapse = ''), sep='\n')
      write.table(AAseq, file = paste("AA_ORFs/", pos$scaffoldID[i], "_", pos$ORFID[i], ".fa", sep=''), col.names = F, row.names = F, quote = F)
    }
    if(as.character(pos$strand[i]) == "-"){
      revComp <- reverseComplement(DNAStringSet(seqs))
      
      extSeq <- as.DNAbin(DNAStringSet(paste(strsplit(as.vector(revComp), split='')[[1]][c(pos$start[i]:pos$end[i])], collapse = "")))
      AAseq <- paste(paste(">", pos$scaffoldID[i], "_", pos$ORFID[i], sep=''), paste(as.character(trans(extSeq, code = 1, codonstart = 1))[[1]], collapse = ''), sep='\n')
      write.table(AAseq, file = paste("AA_ORFs/", pos$scaffoldID[i], "_", pos$ORFID[i], ".fa", sep=''), col.names = F, row.names = F, quote = F)
    }
  }
  
  #Merge all ORFs for BLAST analysis
  system(paste("cat AA_ORFs/", pos$scaffoldID[i], "* > AA_ORFs/", pos$scaffoldID[i], "_ORFs.fa", sep=''))
  system(paste("rm AA_ORFs/", pos$scaffoldID[i], "_ORF_*", sep=""))
  #Save pos file
  write.table(pos, file = paste("ORFs_report/",pos$scaffoldID[i], "_", "ORFs.csv", sep=''), col.names = T, row.names = F, quote = F)
} 
