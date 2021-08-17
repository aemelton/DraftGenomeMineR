### AE Melton, 2020

AnnotateORFsV1DEFUNCT <- function(AA.BlastDB.folder = "~/Dropbox/Genome_PlayGround/AA_BlastDB/",
                         AA.ORF.folder = "~/Dropbox/Genome_PlayGround/AA_ORFs/",
                         annotated.genes.file = "~/Dropbox/Genome_PlayGround/FASTAs/PIP1_3.fa",
                         AA.FASTA.out.folder = "~/Dropbox/Genome_PlayGround/AA_FASTA/",
                         BlastDB.type = "prot",
                         blast.type = "blastp"){
  
#setwd("AA_BlastDB/")
setwd(AA.ORF.folder)
orf.files <- list.files()
file.copy(orf.files, AA.BlastDB.folder)
### All of this will need to be in a loop to loop over each scaffold, generate a db for each, and annotate the ORFs
#genes.seq <- readLines(con = "Scaffold18599_ORFs.fa")
setwd(AA.BlastDB.folder)
for(i in 1:length(orf.files)){
  makeblastdb(file = orf.files[i], dbtype = BlastDB.type)
}
#

setwd(AA.ORF.folder)
orf.files <- list.files()

for(i in 1:length(orf.files)){
  setwd(AA.BlastDB.folder)
  blast.db.path <- orf.files[i]
  annotated.fasta <- readAAStringSet(filepath = annotated.genes.file)
  bl <- blast(db = blast.db.path, type = blast.type)
  cl <- predict(bl, annotated.fasta[1,])
  
  if(nrow(cl) > 0){
    
    blast.csv.filename <- paste0(orf.files[i], "_BlastOut.csv")
    write.csv(x = cl, file = blast.csv.filename)
    
    # Extract and Annotate ORFs
    genes.seq <- readLines(con = orf.files[i])
    header <- NULL
    seq <- NULL
    
    setwd(AA.FASTA.out.folder)
    header <- as.character(cl$QueryID)
    seq <- genes.seq[c(grep(paste(">", cl$QueryID[1], sep=''), genes.seq)+1)] 
    filename <- paste0(header, "_", as.character(cl$SubjectID[1]), ".fasta")
    x <- dplyr::data_frame(name = header, seq = seq)
    writeFasta(data = x, filename = filename)
    
    #for(q in 1:nrow(cl)){
    #  setwd(AA.FASTA.out.folder)
    #  header <- as.character(cl$QueryID[q])
    #  seq <- genes.seq[c(grep(paste(">", cl$SubjectID[q], sep=''), genes.seq)+1)] 
    #  filename <- paste0(header, "_", as.character(cl$SubjectID[q]), ".fasta")
    #  
    #  x <- dplyr::data_frame(name = header, seq = seq)
    #  writeFasta(data = x, filename = filename)
    #  }
    }else{
      next
    }
  }
}