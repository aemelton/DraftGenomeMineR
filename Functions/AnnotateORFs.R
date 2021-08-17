### AE Melton, 2021
# Annotate ORFs

AnnotateORFs <- function(AA.ORF.folder, AA.BlastDB.folder, annotated.genes.file, AA.FASTA.out.folder, BlastDB.type, blast.type){
# Make data base of genes of interest to annotate ORFs
setwd(AA.ORF.folder)
orf.files <- list.files()
#BlastDB.type <- "prot"
#blast.type <- "blastp"
#setwd("AA_BlastDB/")
file.copy(orf.files, AA.BlastDB.folder)
### All of this will need to be in a loop to loop over each scaffold, generate a db for each, and annotate the ORFs
setwd(AA.BlastDB.folder)
for(i in 1:length(orf.files)){
  makeblastdb(file = orf.files[i], dbtype = BlastDB.type)
}

setwd(AA.ORF.folder)
orf.files <- list.files()

for(i in 1:length(orf.files)){
  setwd(AA.BlastDB.folder)
  blast.db.path <- orf.files[i]
  annotated.fasta <- readAAStringSet(filepath = annotated.genes.file)
  bl <- blast(db = blast.db.path, type = blast.type)
  cl <- predict(bl, annotated.fasta) #annotated.fasta[1,]
  blast.csv.filename <- paste0(orf.files[i], "_BlastOut.csv")
  write.csv(x = cl, file = blast.csv.filename, row.names = F)
}  
#if(nrow(cl) > 0){



# Extract and Annotate ORFs
genes.seq <- readLines(con = orf.files[i])
header <- NULL
seq <- NULL

setwd(AA.FASTA.out.folder)
for(i in 1:nrow(cl)){
  header <- as.character(cl$QueryID[i])
  seq <- genes.seq[c(grep(paste(">", cl$QueryID[i], sep=''), genes.seq)+1)] 
  filename <- paste0(header, "_", as.character(cl$SubjectID[i]), ".fasta")
  x <- dplyr::tibble(name = header, seq = seq)
  writeFasta(data = x, filename = filename)
  }
}
