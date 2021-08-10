### AE Melton, 2021
# ASsmeble a protein sequence from translated ORFs

GetProteinSequence <- function(tmp.header, outfile.name){
  
  AA.orf.files <- list.files()
  
  AA.orfs <- sapply(AA.orf.files, readLines)
  
  orfs <- AA.orfs[c(grep(">", AA.orfs)+1)] 
  
  protein <- paste(orfs, collapse="")
  protein <- str_replace_all(protein, "[[:punct:]]", "")
  
  header <- tmp.header
  x <- dplyr::tibble(name = header, seq = protein)
  setwd(protein.fasta.folder)
  writeFasta(data = x, filename = outfile.name)
  
}

