### AEM, 2020
#

#
library(tidyr)
library(dplyr)
library(stringr)
#library(gsubfn)
#

#
setwd("~/Dropbox/BSU_Research/Aquaporin/OutFiles_AEM/Promoters/")
#

###~~~
#What mechanisms regulate aquaporin gene expressions?
###~~~

#
place_seq <- readLines(con = "New_PLACE_data_site/place_seq.txt")
#

#
start.pattern <- "^ID   "
end.pattern <- "//"
#

#
place_seq.split <- split(place_seq, "//", sep = "//")
place_seq.split
#

#
pattern <- "^ID   "
RDE.ID <- subset(place_seq, grepl(pattern, place_seq) )
RDE.ID <- gsub("^.{0,5}", "", RDE.ID)
length(RDE.ID)
#

#
pattern <- "^KW   "
Categories <- subset(place_seq, grepl(pattern, place_seq) )
Categories <- gsub("^.{0,5}", "", Categories)
length(Categories)
#

#
cat <- data.frame(RDE.ID, Categories)
head(cat)
cat
write.csv(x = cat, file = "RDE_and_Categories.csv", row.names = F)
#

#
cat <- read.csv(file = "RDE_and_Categories.csv", header = T)
head(cat)
#cat.split <- separate(data = cat,
#                      col = "Categories",
#                      into = rep("Category", 21),
#                      sep = ";")
cat.split <- read.csv("OutFiles_AEM/Promoters/RDE_and_Categories_split.csv")
cat.split
#

#
plant_parts <- c("flower", "nodule", "leaf",
                 "xylem", "phloem", "root",
                 "shoot", "seed", "endosperm",
                 "wood", "bud", "fruit",
                 "cotyledon", "meristem", "guard cell",
                 "tuber", "embryo", "stem",
                 "vascular", "hair", "trichome",
                 "pollen", "primordia", "anther",
                 "hypocotyl")
#

#
light <- c("light", "chloroplast", "circadian",
                "chlorophyll", "phytochrome", "dark",
                "plastid", "etiolation", "photo")

ABA <- c("ABA", "abscisic", "Abscisic")

water_stress <- c("drought", "water", "dehydration", "osmotic")

temperature_stress <- c("cold", "heat", "temperature", "freezing")

stress_hormone <- c("ethylene", "auxin", "cytokinin", "gibberellin")

#storage_proteins <- c("glutenin", "legumin", "calmodulin",
#                     "phaseolin", "napin")

#Molecular.Function.cats <- c("")

growth <- c("sugar", "elongation", "cell")

other_stress <- c("xenobiotic", "wound", "hypoxic", "oxygen", "anaerobic", "salt")
#

# Filter out RDEs that don't filter into a listed category
df <- data.frame(RDE.ID = character(469),
                 light = character(469),
                 ABA = character(469),
                 water_stress = character(469),
                 temperature_stress = character(469), 
                 stress_hormone = character(469),
                 #storage_proteins = character(469),
                 growth = character(469), 
                 other_stress = character(469),
                 plant_parts = character(469))
df$RDE.ID <- cat$RDE.ID
#

# Function by alexis_jar on stackoverflow (https://stackoverflow.com/questions/19747384/create-new-column-in-dataframe-based-on-partial-string-matching-other-column)
ff = function(x, patterns, replacements = patterns, fill = 0, ...)
{
  stopifnot(length(patterns) == length(replacements))
  
  ans = rep_len(as.character(fill), length(x))    
  empty = seq_along(x)
  
  for(i in seq_along(patterns)) {
    greps = grepl(patterns[[i]], x[empty], ...)
    ans[empty[greps]] = replacements[[i]]  
    empty = empty[!greps]
  }
  
  return(ans)
}
#
### LOOPS!
###############################################################################
#
for(j in 1:nrow(cat)){
  #
  #
  for(i in 1:length(plant_parts)){
    df$plant_parts <- ff(x = cat$Categories, patterns = plant_parts)
  }
  for(i in 1:length(light)){
    df$light <- ff(x = cat$Categories, patterns = light)
  }
  for(i in 1:length(ABA)){
    df$ABA <- ff(x = cat$Categories, patterns = ABA)
  }
  for(i in 1:length(water_stress)){
    df$water_stress <- ff(x = cat$Categories, patterns = water_stress)
  }
  for(i in 1:length(temperature_stress)){
    df$temperature_stress <- ff(x = cat$Categories, patterns = temperature_stress)
  }  
  for(i in 1:length(stress_hormone)){
    df$stress_hormone <- ff(x = cat$Categories, patterns = stress_hormone)
  }
  #for(i in 1:length(storage_proteins)){
  #  df$storage_proteins <- ff(x = cat$Categories, patterns = storage_proteins)
  #}
  for(i in 1:length(growth)){
    df$growth <- ff(x = cat$Categories, patterns = growth)
  }
  for(i in 1:length(other_stress)){
    df$other_stress <- ff(x = cat$Categories, patterns = other_stress)  
  }
}
df
write.csv(x = df, file = "RDE_and_Categories_V6.csv")
#

#
df.bin <- df
for(r in 1:nrow(df.bin)){
  #
  #
  for(c in 2:ncol(df)){
    if(df.bin[r,c] != 0)
    df.bin[r,c] <- 1
  }
}
df.bin.trim <- df.bin[,1:9]
write.csv(x = df.bin.trim, file = "RDE_and_PlatPart_COUNTS_V6.csv")
#

#
tmp <- NULL
CatMat <- data.frame(RDE.ID = character(), Categories = character())
for(i in 1:nrow(df.bin)){
  tmp$RDE.ID <- df.bin[i,1]
  tmp$Categories <- sum(as.numeric(df.bin[i,2:8]))
  CatMat <- rbind(CatMat, tmp)
}
CatMat
max(CatMat$Categories)
#

#
PivotTable <- read.csv("RDE_analyses_sagebrush_Aquaporins_w_Total.csv")
PivotTable <- PivotTable[-c(2:3)]
#

#
#Genes <- colnames(PivotTable[,2:45])
df.bin <- df.bin[,1:8]
Category <- colnames(df.bin[-1])
#columns <- c("Category", Genes)
big.df <- data.frame(Category, PivotTable[1:7,2:45])
big.df[big.df == 1] <- 0
big.df
#

#
PivotTable # Tells us whether a specific gene contains a specific RDE
df.bin # Tells us whether a specific RDE falls into a specific category
big.df # Will tell us how many RDEs for a given gene fall into a given category
#

#
PivotTable$RDE <- gsub(x = PivotTable$RDE, pattern = "_.*", replacement = "")
df.bin.trim <- df.bin[(df.bin$RDE.ID %in% PivotTable$RDE),]
nrow(df.bin.trim)
#

#
tmp <- big.df
for(i in 2:ncol(PivotTable)){
  for(j in 1:nrow(PivotTable)){
    if(PivotTable[j,i] == 1){
      for(q in 2:ncol(df.bin.trim)){
        if(df.bin.trim[j,q] == 1){
          g <- q-1
          tmp[g,i] <- tmp[g,i] + 1
        }
      }
    }
  }
}
#

#

#
tmp
for(i in 1:nrow(tmp)){
tmp[i,46] <- sum(tmp[i,2:45])
}
tmp
write.csv(x = tmp, file = "Promoter_Category_Count_By_Scaffold_v5.csv", row.names = F)
#
###############################################################################
#

#
PivotTable <- read.csv("OutFiles_AEM/Promoters/RDE_analyses_sagebrush_Aquaporins_w_Total.csv")
PivotTable
PivotTable.sort <- PivotTable[order(-PivotTable$Total),]
PivotTable.sort
#

#
PivotTable.sort$RDE <- str_remove_all(PivotTable.sort$RDE, "_.*")
PivotTable.sort$Category <- FinalCat$Category_Final[match(PivotTable.sort$RDE, FinalCat$RDE.ID)]
write.csv(x = PivotTable.sort, file = "RDE_analyses_sagebrush_Aquaporins_w_Total_Partial_category_V3_7_28.csv", row.names = F)
#

#
PivotTable.sort %>%
count(Category)
#

# Cheating - pulling in all the one Sven had categorized
Sven_Table <- read.csv("../Promoters/RDE_analyses_sagebrush_Aquaporins_reduced_Promoters.csv")
Sven_Table
PivotTable$Category <- Sven_Table$Category_FINAL[match(PivotTable$RDE, Sven_Table$RDE)]
PivotTable
write.csv(x = PivotTable, file = "Promoters/RDE_analyses_sagebrush_Aquaporins_w_Total_Partial_category_V2_7_28.csv", row.names = F)
#

############ Fail
#
df <- as.data.frame(Categories)
Split.Categories <- separate_rows(df$Categories, sep = ";")
df
#