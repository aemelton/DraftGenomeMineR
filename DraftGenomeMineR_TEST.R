### AE Melton, 2020

#################################################################################################################
# Module 1: Set up the environment for genome mining
#################################################################################################################

# Load all the libraries!
source("~/Dropbox/R_Packages/DraftGenomeMineR/RequiredLibraries.R") # This script contains a list of packages that will be installed, if needed, and load the libraries
#

# Directory for where everything shall live; all folders will be in this "master" directory
project.folder <- "~/Dropbox/Genome_PlayGround/"
setwd(project.folder)
#

#################################################################################################################
# Module 2: Load data and perform a blast search
#################################################################################################################


# The following are variables that are listed in a function to perform a blast search
# The function requires:
# fasta for draft genome, a blast database, a fasta to use as a query, specify whether the query is AA or DNA sequence, and an output fasta file name

#
setwd(project.folder)
query.file.path <- "FASTAs/PIP1_3.fa"
genome.file.name <- "Draft_Genome_Assembly.fasta"
genome.path <- "FASTAs/Draft_Genome_Assembly.fasta"
blast.db.path <- "BlastDBs/Draft_Genome_Assembly.fasta"
AA.BlastDB.folder <- "~/Dropbox/Genome_PlayGround/AA_BlastDB/"
AA.ORF.folder <- "~/Dropbox/Genome_PlayGround/AA_ORFs/"
max.e <- 0.00001
min.bit <- 50
perc.ident <- 50.000
query.type <- "AA"
blast.type <- "tblastn"
make.BlastDB <- T
BlastDB.type <- "nucl"
#

# Read in the draft genome to be mined. readLines will read in the fasta file line by line. 
# Be aware of the return characters in your text files, as different operating systems may read these differently.
genome <- readLines(con = genome.path)
head(genome) # Print the top 6 lines of the fasta. There should be no spaces in what is printed. Each header should be
# on one line, followed by its entire scaffold on the next.
#

# Do you need to make a new blast database?
if(make.BlastDB == TRUE){
setwd("BlastDBs/")
makeblastdb(file = genome.file.name, dbtype = BlastDB.type)
}
#

# Read in the query. Specify whether it is a DNA (or RNA; these will be read the same) or amino acid sequence.
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

# Filter out hits to just have unique scaffolds to extract from draft genome (no need to extract the same scaffold 
# multiple times if it has multiple hits)
cl.filt <- subset(x = cl, Perc.Ident >= 75) # Several ways to filter: percent identity, E-value, scaffold ID... SubjectID =="LBNU01006105.1"cl.filt
cl.filt.unique <- cl.filt[!duplicated(cl.filt[,c('SubjectID')]),]
cl.filt.unique
nrow(cl)
nrow(cl.filt)
nrow(cl.filt.unique)
write.csv(x = cl.filt.unique, file = "Unique_Filtered_Blast_Hit_Info.csv", row.names = F)
#

##########
DoBlastSearch(query.file.path <- "FASTAs/Dreb19.fna",
              genome.file.name <- "final.contigs.fa",
              genome.path <- "FASTAs/final.contigs.fa",
              blast.db.path <- "BlastDBs/final.contigs.fa",
              AA.BlastDB.folder <- "~/Dropbox/Genome_PlayGround/AA_BlastDB/",
              AA.ORF.folder <- "~/Dropbox/Genome_PlayGround/AA_ORFs/",
              min.e <- 0.00e+00,
              query.type <- "DNA",
              blast.type <- "blastn",
              make.BlastDB <- T,
              BlastDB.type <- "nucl")
##########

#################################################################################################################
# Module 3: Extract scaffolds of interest identified in blast search
#################################################################################################################

# This chunk of code uses the output of the previous chunk, the blast search, to find and extract scaffolds of interest

# Extract scaffolds from draft genome and make a fasta file
cl.filt.unique <- read.csv(file = "Unique_Filtered_Blast_Hit_Info.csv") # Output of Chunk1
header <- NULL # Generate empty objects to store the headers and sequences for the scaffolds of interest
seq <- NULL
#S.start <- NULL
#S.end <- NULL
for(i in 1:nrow(cl.filt.unique)){

  #Add DNA sequence adapted to strand
  #Extract seq from FASTA file
  #header[i] <- paste0(">", cl.filt.unique$SubjectID[i])
  header[i] <- as.character(cl.filt.unique$SubjectID[i])
  seq[i] <- genome[c(match(paste(">", cl.filt.unique$SubjectID[i], sep=''), genome)+1, fixed = T)]
  
  #S.start[i] <- cl.filt.unique$S.start[i]
  #S.end[i] <- cl.filt.unique$S.end[i]
  
  #Package sequence and extract start and end
  #header[i] <- paste0(">", cl$SubjectID[i])
  #seq[i] <- paste(strsplit(seqRaw, split='')[[1]][as.numeric(start):as.numeric(end)], collapse = "")
}
x <- dplyr::tibble(name = header, seq = seq) # This will assemble the headers and sequences into an object that
# can be written into a FASTA format file
x
#

###
x <- GetScaffolds(BlastHitFile = "Unique_Filtered_Blast_Hit_Info.csv", genome = genome)
x
writeFasta(data = x, filename = "~/Dropbox/Genome_PlayGround/Output_FASTAs/OUT.fa")
###

#################################################################################################################
# Module 4: Identify ORFs in the extracted scaffolds
#################################################################################################################

# Find ORFs in scaffolds; This is pretty memory intense. It will write out a lot of fasta
# files - one for each ORF. Larger scafolds may not be able to annotated with this on
# computers without a lot of free hard drive space

#scaffold <- readLines("Output_FASTAs/Scaffold128070.fasta")
#scaffoldID <- grep(pattern = "^>", x = scaffold, value = T)
#scaffoldID <- gsub(pattern = ">", replacement = "", x = scaffoldID)
#tryCatch(
#  {
#  for(i in 1:length(scaffoldID)){
#    findORFsTranslateDNA2AA(scaffold = scaffold, scaffoldID = scaffoldID[i])
#  }
#  }, warning = function(w) {
#    warning-handler-code
#  }, error = function(e) {
#    error-handler-code
#  }, finally = {
#    cleanup-code
#  })
#

scaffold <- readLines("Output_FASTAs/OUT.fa")
scaffoldID <- grep(pattern = "^>", x = scaffold, value = T)
scaffoldID <- gsub(pattern = ">", replacement = "", x = scaffoldID)
tryCatch(
  {
    for(i in 1:length(scaffoldID)){
      findORFsTranslateDNA2AA(scaffold = scaffold, scaffoldID = scaffoldID[i])
    }
  })

###
FindORFs(OutputFasta = "Output_FASTAs/OUT.fa", Minimum.Length = 40)
###
#################################################################################################################
# Module 5: Annotate ORFs and write out a fasta containing the amino acid sequence for each scaffold
#################################################################################################################

# Make data base of genes of interest to annotate ORFs
setwd(AA.ORF.folder)
orf.files <- list.files()
BlastDB.type <- "prot"
blast.type <- "blastp"
#setwd("AA_BlastDB/")
file.copy(orf.files, AA.BlastDB.folder)
### All of this will need to be in a loop to loop over each scaffold, generate a db for each, and annotate the ORFs
#genes.seq <- readLines(con = "Scaffold18599_ORFs.fa")
setwd(AA.BlastDB.folder)
for(i in 1:length(orf.files)){
  makeblastdb(file = orf.files[i], dbtype = BlastDB.type)
}


# Annotate ORFs of interest and write them to their own fasta
annotated.genes.file <- "~/Dropbox/Genome_PlayGround/FASTAs/PIP1_3.fa"
AA.FASTA.out.folder <- "~/Dropbox/Genome_PlayGround/AA_FASTA/"
BlastDB.type <- "prot"
blast.type <- "blastp"

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
    #for(q in 1:nrow(cl)){
    #  setwd(AA.FASTA.out.folder)
    #  header <- as.character(cl$QueryID[q])
    #  seq <- genes.seq[c(grep(paste(">", cl$SubjectID[q], sep=''), genes.seq)+1)] 
    #  filename <- paste0(header, "_", as.character(cl$SubjectID[q]), ".fasta")
    #  
    #  x <- dplyr::data_frame(name = header, seq = seq)
    #  writeFasta(data = x, filename = filename)
    #  }
#  }else{
#    next
#  }
#}
#
##########
setwd(project.folder)
AnnotateORFs(AA.BlastDB.folder <- "~/Dropbox/Genome_PlayGround/AA_BlastDB/",
             AA.ORF.folder <- "~/Dropbox/Genome_PlayGround/AA_ORFs/",
             BlastDB.type = "prot",
             blast.type = "blastp",
             annotated.genes.file <- "~/Dropbox/Genome_PlayGround/FASTAs/PIP1_3.fa",
             AA.FASTA.out.folder <- "~/Dropbox/Genome_PlayGround/AA_FASTA/")
##########
##################################################################################################################
# Module 6: Assemble proteins from AA_ORF files
#################################################################################################################



##################################################################################################################
# Module 7: Use the blast results to extract promoter sequences for analyses
#################################################################################################################

# This function will need some love. Currently, there is one filtering step and one string modifcation step that are example specific.
setwd(project.folder)
GetPromoterSequences(orfs.report = "ORFs_report/Scaffold151535_ORFs.csv",
                     blast.out = "AA_BlastDB/Scaffold151535_ORFs.fa_BlastOut.csv",
                     scaffold.fasta = "Output_FASTAs/OUT.fa",
                     promoter.csv.file.out = "PROMOTER_OUT_TEST.csv",
                     promoter.sequence.fasta = "PROMOTER_OUT_TEST.fa",
                     pattern.to.keep = "Scaffold151535_",
                     promoters.folder = "Promoters/",
                     perc.ident.threshold = 50)

##################################################################################################################
# Module X: Align assembled proteins using MAFFT
#################################################################################################################

DoMafftAlignment(align.method = "--auto",
in.file = "~/Dropbox/BSU_Research/Aquaporin/OutFiles_AEM/FASTA/AQP_AEM_trim.fa",
out.file = "~/Desktop/Test.fa")

##################################################################################################################
# Module X: Reconstruct phylogeny using RAxML
#################################################################################################################

setwd("~/Desktop/")
DoRAxMLReconstruction(algorithm = "a", 
                      in.file = "AQP_AEM_trim.fa", 
                      out.file = "raxTest", 
                      model = "PROTGAMMAGTR", 
                      parsimony.random.seed = 12345, 
                      bootstrap.random.seed = 12345,
                      BS.rep.count = 100)

##################################################################################################################
# Module X: THE NEXT SECTIONS ARE UNDERCONSTRUCTION. PLEASE MAINTAIN A SAFE SPEED AND BE AWARE OF ROAD WORK SIGNS FOR THE NEXT 5 MILES.
#################################################################################################################


######################

ClassifyPromoterElements(New.Place.DB = "~/Dropbox/BSU_Research/Aquaporin/OutFiles_AEM/Promoters/New_PLACE_data_site/place_seq.txt",
                         out.file = "~/Desktop/testPromoEleClass.csv")

######################
setwd(project.folder)
GetPromoterElementsPerScaffold(Promoters.Folder <- "~/Dropbox/Genome_PlayGround/Promoters/",
                               Promoters <- "~/Desktop/PROMOTER_OUT_TEST.csv",
                               New.Place.Folder <- "~/Dropbox/Genome_PlayGround/New_PLACE_txt",
                               Output.Folder <- "~/Desktop/Promoter_CSV_TEST/",
                               Output.Filename <- "~/Desktop/RDE_per_Scaffold_Counts.csv")

######################

GetORFMap(scaffold.file = "Output_FASTAs/Scaffold128070.fasta",
          ORF.data.file = "ORFs_report/Scaffold128070_ORFs.csv",
          user.sep = " ",
          out.file = "~/Desktop/TEST.pdf")

######################