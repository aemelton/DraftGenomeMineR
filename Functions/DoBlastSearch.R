### AE Melton, 2020

DoBlastSearch <- function(query.file.path,
                          genome.file.name,
                          genome.path,
                          blast.db.path,
                          AA.BlastDB.folder,
                          AA.ORF.folder,
                          min.e = 0.0005,
                          query.type,
                          blast.type,
                          make.BlastDB = T,
                          BlastDB.type){
  # Read in the draft genome to be mined
  genome <- readLines(con = genome.path)
  head(genome)
  #
  
  # Do you need to make a new blast database?
  if(make.BlastDB == TRUE){
    setwd("BlastDBs/")
    makeblastdb(file = genome.file.name, dbtype = BlastDB.type)
  }
  #
  
  # Read in the query
  setwd(project.folder)
  if (query.type == "DNA") {
    query <- readDNAStringSet(filepath = query.file.path,
                              format = "fasta")
  } else {
    query <- readAAStringSet(filepath = query.file.path,
                             format = "fasta")
  }
  #
  
  # Set up the blast
  bl <- blast(db = blast.db.path, type = blast.type)
  #
  
  # Do the blast query
  cl <- predict(bl, query)
  cl
  #
  
  # Filter out hits to just have unique scaffolds to extract from draft genome
  cl.filt <- subset(x = cl, E <= min.e)
  cl.filt
  cl.filt.unique <- cl.filt[!duplicated(cl.filt[,c('QueryID')]),]
  cl.filt.unique
  nrow(cl)
  nrow(cl.filt)
  nrow(cl.filt.unique)
  write.csv(x = cl.filt.unique, file = "Unique_Filtered_Blast_Hit_Info.csv", row.names = F)
  #
}
