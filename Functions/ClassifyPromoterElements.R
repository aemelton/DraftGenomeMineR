### AE Melton, 2021
# Group New Place RDEs into categorical bins

ClassifyPromoterElements <- function(New.Place.DB, out.file){

### First, read in the New Place database 
#
place_seq <- readLines(con = New.Place.DB)
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
write.csv(x = cat, file = "RDE_and_Categories_one_row_per_RDE.csv", row.names = F)
#

#
#cat <- read.csv(file = "RDE_and_Categories.csv", header = T)


#
plant_parts <- c("flower", "nodule", "leaf",
                 "xylem", "phloem", "root",
                 "shoot", "seed", "endosperm",
                 "wood", "bud", "fruit",
                 "cotyledon", "meristem", "guard cell",
                 "tuber", "embryo", "stem",
                 "vascular", "hair", "trichome",
                 "pollen", "primordia", "anther",
                 "hypocotyl", "mesophyll")
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

df.bin <- df
for(r in 1:nrow(df.bin)){
  for(c in 2:ncol(df)){
    if(df.bin[r,c] != 0)
      df.bin[r,c] <- 1
  }
}

write.csv(x = df.bin, file = out.file, row.names = F)
}
