## AE Melton, 2021
# Modified code from scripts by SB. This script generates a function to look through New Place output files and the CSV with promoter data to generate counts of RDEs per scaffold.

GetPromoterElementsPerScaffold <- function(Promoters.Folder,
                                           Promoters,
                                           New.Place.Folder,
                                           Output.Folder,
                                           Output.Filename){
setwd(Promoters.Folder)
Promoters <- read.csv(Promoters)
PLACEfiles <- list.files(path = New.Place.Folder, pattern='.txt', full.names = T)

OUT <- NULL
for(i in 1:length(PLACEfiles)){
ScaffoldID <- gsub(pattern = "_New_PLACE.txt", replacement = "", x = PLACEfiles[i])
ScaffoldID <- gsub(pattern = "New_PLACE_txt/", replacement = "", x = ScaffoldID)
x <- readLines(PLACEfiles[i])
RDE <- sapply(strsplit(x[3:length(x)], " "),"[[",1)
OUT <- as.data.frame(table(RDE))
OUT <- data.frame(OUT, rep(ScaffoldID, nrow(OUT)))
colnames(OUT) <- c("RDE", "Count", "ScaffoldID")
setwd(Output.Folder)
file.name <- paste0(ScaffoldID, "_RDE_Counts.csv")
write.csv(x = OUT, file = file.name, row.names = F)
setwd(Promoters.Folder)
}

setwd(Output.Folder)
csv.files <- list.files()
csv <- NULL
df <- data.frame(RDE = character(), ScaffoldID = character())
for(i in 1:length(csv.files)){
  csv <- read.csv(file = csv.files[i])
  newID <- csv$ScaffoldID[1]
  csv.boop <- csv[,1:2]
  colnames(csv.boop) <- c("RDE", as.character(newID))
  df <- full_join(df, csv.boop, 
                  by = "RDE")    #rbind(df, csv)
}

df <- df[-2] # It makes an empty scaffoldID column. KICK IT OUT.
df[is.na(df)] <- 0
write.csv(x = df, file = Output.Filename, row.names = F)
}    
