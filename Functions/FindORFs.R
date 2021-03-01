### AE Melton, 2020

FindORFs <- function(OutputFasta, Minimum.Length){
  
# Find ORFs in scaffolds
scaffold <- readLines(OutputFasta)
scaffoldID <- grep(pattern = "^>", x = scaffold, value = T)
scaffoldID <- gsub(pattern = ">", replacement = "", x = scaffoldID)
tryCatch(
  {
    for(i in 1:length(scaffoldID)){
      findORFsTranslateDNA2AA(scaffold = scaffold, scaffoldID = scaffoldID[i], MinLen = Minimum.Length)
    }
  })
}
