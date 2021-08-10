### AE Melton, 2021
# Subset a draft genome assembly

SubsetDraftGenomeAssembly <- function(genome, sample.size, output.file){
  
headers <- genome[c(grep(">", genome))]
head(headers)

TO.GET <- sample(headers, sample.size)

header <- NULL
seq <- NULL

for(i in 1:length(TO.GET)){

  header[i] <- as.character(TO.GET[i])
  seq[i] <- genome[c(match(TO.GET[i], genome)+1, fixed = T)]

  }
x <- dplyr::tibble(name = header, seq = seq)
writeFasta(data = x, filename = output.file)
}
###