# AE Melton, 2020
# Original script by Sven Buerki

#
#install.packages("ape")
#install.packages("bit64")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.11")
#BiocManager::install("ORFik")
#install.packages("readr")
#

#
#library(ape)
#library(ORFik)
#library(readr)
#source("Scripts_SB/findORFsTranslateDNA2AA.R")
#source("Scripts_AEM/processFile.R")
#

#
setwd(project.folder)

#Load csv file
csv <- read.csv("G1_out_subset.csv")

#Subset according to prediction of full gene
#csv <- subset(csv, csv$Prediction_full_gene_in_scaffold == "Possible")

#Make list of scaffold to run function in loop
scaffoldID <- as.vector(csv$SubjectID)

#Read FASTA file (line by line)
#scaffold <- readr::read_lines("Output_FASTAs/G1AQPScaffolds.fa")
scaffold <- processFile("Output_FASTAs/G1AQPScaffolds.fa")
head(scaffold)

for(i in 1:length(scaffoldID)){
  findORFsTranslateDNA2AA(scaffold = scaffold, scaffoldID = scaffoldID[i], MinLen = 40)
}

####
#Merge data on ORFs coding for Aquaporin proteins (confirmed by BLASTp)
####
#List all files with ORFs and info if ORF codes for Aquaporin (=1)
ORFfiles <- list.files(path = "ORFs_report", pattern = ".csv", full.names = T)
BLASTORFs <- list.files(path = "../ORFs/BLAST_AA", pattern = '.csv', full.names = T)

OUT <- NULL
for(i in 1:length(ORFfiles)){
  tmp <- read.csv(ORFfiles[i])
  #Subset to only include ORFs coding for target gene
  tmp <- subset(tmp, tmp$Aquaporin == "1")
  if(nrow(tmp) > 0){
    #Add a col in tmp to store GenBank acc# from BLAST search
    tmp$BLAST_top_hit <- rep("NA", nrow(tmp))
    
    #Build a vector to fetch info from BLAST report (on ORFs)
    BLASTID2search <- paste(as.vector(tmp$scaffoldID), "_", as.vector(tmp$ORFID), sep='')
    #Find BLAST report for Scaffold and extract info for ORFID
    BLAST <- read.csv(BLASTORFs[grep(paste("BLAST_AA/", tmp$scaffoldID[1], "_", sep=''), BLASTORFs)], header=F)
    BLASTred <- subset(BLAST, BLAST$V1 %in% BLASTID2search)
    for(j in 1:length(BLASTID2search)){
      BLASTtmp <- subset(BLASTred, BLASTred == BLASTID2search[j])
      tmp$BLAST_top_hit[j] <- paste(as.vector(BLASTtmp$V2), collapse = ":")
    }  
  }else{
    #If there are no ORFs underpinning an Aquaporin gene then go to next file (but don't break the loop)
    next
  }
  OUT <- rbind(OUT, data.frame(tmp))
}
FINAL <- as.data.frame(OUT[,c(8,9,4,5,7,11)])
write.table(FINAL, "Scaffolds_ORFs_AA.csv", row.names = F)

###
#Edit csv to add strand, N ORFs and Charset
###
#Load csv file
csv <- read.csv("../ORFs/all_BLAST_PCR_scaffolds_predictions_23Jan2020.csv")

#Add cols to csv
csv$strand <- rep("NA", nrow(csv))
csv$N_ORF <- rep("NA", nrow(csv))
csv$ORFID <- rep("NA", nrow(csv))
csv$charset_ORFs <- rep("NA", nrow(csv))
csv$BLAST_ORFs <- rep("NA", nrow(csv))
csv$SeqDNAORF <- rep("NA", nrow(csv))
csv$SeqAAORF <- rep("NA", nrow(csv))

#Add info to csv based on FINAL
scaffoldORF <- unique(FINAL$scaffoldID)

for(i in 1:length(scaffoldORF)){
  #which row in csv is equal to scaffoldORF
  rwocsv <- match(scaffoldORF[i], csv$ScaffoldID)
  
  #Subset FINAL to scaffoldORF
  tmpScaff <- subset(FINAL, FINAL$scaffoldID == scaffoldORF[i])
  
  #Add strand (majority of ORFs)
  strandScaff <-  names(sort(table(tmpScaff$strand), decreasing = T)[1])
  csv$strand[rwocsv] <- strandScaff
  
  #Count N ORFs using strandScaff to further subset tmpScaff
  tmpScaff <- subset(tmpScaff, tmpScaff$strand == strandScaff)
  #Sort according to start
  tmpScaff <- tmpScaff[order(tmpScaff$start),]
  csv$N_ORF[rwocsv] <- nrow(tmpScaff)
  
  #Add ORFs IDs
  csv$ORFID[rwocsv] <- paste(tmpScaff$ORFID, collapse = "; ")
  
  #Add charset
  csv$charset_ORFs[rwocsv] <- paste(paste(tmpScaff$start, tmpScaff$end, sep='-'), collapse = "; ")
  
  #Add BLAST
  csv$BLAST_ORFs[rwocsv] <- paste(tmpScaff$BLAST_top_hit, collapse = ":")
  
  #Add DNA sequence adapted to strand
  #Extract seq from FASTA file
  seqRaw <- scaffold[c(grep(paste(">", scaffoldORF[i], sep=''), scaffold)+1)] 
  
  #Package sequence and extract start and end
  if(strandScaff == "+"){
    csv$SeqDNAORF[rwocsv] <- paste(strsplit(seqRaw, split='')[[1]][c(min(tmpScaff$start):max(tmpScaff$end))], collapse = "")
  }
  if(strandScaff == "-"){
    revComp <- as.vector(reverseComplement(DNAStringSet(seqRaw)))
    csv$SeqDNAORF[rwocsv] <- paste(strsplit(revComp, split='')[[1]][c(min(tmpScaff$start):max(tmpScaff$end))], collapse = "")
  }
  
  ##Generate AA seq for each ORF, but merge them
  AAseq <- NULL
  for(j in 1:nrow(tmpScaff)){
    if(as.character(tmpScaff$strand[j]) == "+"){
      extSeq <- as.DNAbin(DNAStringSet(paste(strsplit(seqRaw, split='')[[1]][c(tmpScaff$start[j]:tmpScaff$end[j])], collapse = "")))
      AAORF <- paste(as.character(trans(extSeq, code = 1, codonstart = 1))[[1]], collapse = '')
      AAseq <- c(AAseq, AAORF)
    }
    if(as.character(tmpScaff$strand[j]) == "-"){
      revComp <- as.vector(reverseComplement(DNAStringSet(seqRaw)))
      
      extSeq <- as.DNAbin(DNAStringSet(paste(strsplit(as.vector(revComp), split='')[[1]][c(tmpScaff$start[j]:tmpScaff$end[j])], collapse = "")))
      AAORF <- paste(as.character(trans(extSeq, code = 1, codonstart = 1))[[1]], collapse = '')
      AAseq <- c(AAseq, AAORF)
    }
    
  }
  csv$SeqAAORF[rwocsv] <- paste(AAseq, collapse = '')
}
write.csv(csv, file="all_BLAST_PCR_scaffolds_ORFs_26Jan20.csv", row.names = F, quote = F)

###
#Create FASTA files: AA or DNA
###
#rm scaffolds without Aquaporin genes
csv <- csv[-grep("NA", csv$N_ORF),]
seqType <- "DNA"
FASTA <- NULL
for(i in 1:nrow(csv)){
  if(seqType == "AA"){
    AA <- paste(paste(">", csv$ScaffoldID[i], "_", csv$Gene_hypothesis[i], sep=''), csv$SeqAAORF[i], sep='\n')
    FASTA <- rbind(FASTA, AA)  
  }
  if(seqType == "DNA"){
    DNA <- paste(paste(">", csv$ScaffoldID[i], "_", csv$Gene_hypothesis[i], sep=''), csv$SeqDNAORF[i], sep='\n')
    FASTA <- rbind(FASTA, DNA) 
  }
}

#write.table(FASTA, file="Aquaporin_AA_seq.fa", col.names = F, row.names = F, quote = F)
write.table(FASTA, file="Aquaporin_DNA_seq.fa", col.names = F, row.names = F, quote = F)

###
#Build map of scafflod with ORFs
###
#Read FASTA file (line by line)
scaffold <- readLines("Scaffolds_BLAST.fasta")
#CSV master file
csv <- read.csv("all_BLAST_PCR_scaffolds_ORFs_NPA_TMH_PSORT_30Jan20.csv")
#List all files with ORFs and info if ORF codes for Aquaporin (=1)
ORFfiles <- list.files(path = "ORFs_report", pattern = ".csv", full.names = T)

for(i in 1:length(ORFfiles)){
  print(i)
  ORF <- read.csv(ORFfiles[i])
  strandSeq <- as.vector(csv$strand[match(ORF$scaffoldID[1], as.vector(csv$ScaffoldID))])
  ORFstrand <- subset(ORF, ORF$strand == strandSeq)
  aquaORF <- which(ORFstrand$Aquaporin == "1")
  seq <- strsplit(scaffold[grep(paste(">", ORF$scaffoldID[1], sep=''), scaffold)+1], split='')
  
  #Create plot
  if(length(aquaORF > 0)){
    pdf(paste("ORFs_map/", ORF$scaffoldID[1], "_ORFs_annotated.pdf", sep=''))
    ORFAqua <- ORFstrand[aquaORF,]
    ORFstrand <- ORFstrand[-aquaORF,]
    plot(x=1, y=1, xlim=c(0,length(seq[[1]])), ylim=c(0,2), type='n', bty="n", axes=F, xlab = "", ylab='')
    #Title
    text(x=5, y=2, paste(ORF$scaffoldID[1], " (strand: ", strandSeq, ")", sep=''), adj=0, cex=.8)
    #Create sequence
    segments(x0=0, x1=length(seq[[1]]), y0=1, y1=1, col='black', lwd=3)
    #Add ORFs
    rect(xleft=ORFstrand$start, xright=ORFstrand$end, ybottom=0.75, ytop=1.25, col='grey')
    rect(xleft=ORFAqua$start, xright=ORFAqua$end, ybottom=0.75, ytop=1.25, col='blue')
    text(x=(ORFAqua$start + ORFAqua$end)/2, y=0.7, paste(ORFAqua$ORFID, " (", ORFAqua$start, ":", ORFAqua$end, ")", sep=''), srt=90, col='blue', adj=1, cex=0.4)
    text(x=(ORFstrand$start + ORFstrand$end)/2, y=1.3, paste(ORFstrand$ORFID, " (", ORFstrand$start, ":", ORFstrand$end, ")", sep=''), srt=90, col='black', adj=0, cex=0.4)
    #Add x axis
    axis(side = 1)
    mtext("Sequence (bp)", pos = c(0,0.5), side=1, line=2, cex.lab=0.6,las=1)
    dev.off()
  }else{
    next
  }
  
}

#####
#Check for NPA motifs (expect 2) in AA sequences
#####
csv$N_NPA <- rep("NA", nrow(csv))
csv$NPA_motif <- rep("NA", nrow(csv))
csv$StartNPA <- rep("NA", nrow(csv))
  
for(i in 1:length(csv$SeqAAORF)){
  if(is.na(csv$SeqAAORF[i]) == FALSE){
    #Match to find NP and locations along sequence
    tmp <- matchPattern(pattern = "NP", AAString(csv$SeqAAORF[i])) 
    tmp <- as.data.frame(ranges(tmp))
    
    if(nrow(tmp) > 0){
      #N NP
      csv$N_NPA[i] <- nrow(tmp)
      #Infer motifs
      csv$NPA_motif[i] <- paste(paste(rep("NP", nrow(tmp)), strsplit(as.vector(csv$SeqAAORF[i]), split='')[[1]][as.numeric(tmp[,2]+1)], sep = ''), collapse = '/')
      #Start motifs
      csv$StartNPA[i] <- paste(tmp[,1], collapse = '/')
    }else{
      tmp <- matchPattern(pattern = "PA", AAString(csv$SeqAAORF[i]))
      tmp <- as.data.frame(ranges(tmp))
      if(nrow(tmp) > 0){
        #N NP
        csv$N_NPA[i] <- nrow(tmp)
        #Infer motifs
        csv$NPA_motif[i] <- paste(paste(strsplit(as.vector(csv$SeqAAORF[i]), split='')[[1]][as.numeric(tmp[,1]-1)], rep("PA", nrow(tmp)), sep = ''), collapse = '/')
        #Start motifs
        csv$StartNPA[i] <- paste(tmp[,1]-1, collapse = '/')
      }else{
        next
      }
    }
  }else{
    next
  }
}
write.csv(csv, file='all_BLAST_PCR_scaffolds_ORFs_NPA_30Jan20.csv', row.names = F, quote = F)
###
#Download Aquaporin AA seq from csv
###

#Create a list with unique GenBank acc# to be used to download seq
AllAABLAST <- paste(csv$BLAST_ORFs, collapse = ":")
AllAABLAST <- sort(unique(strsplit(AllAABLAST, split=':')[[1]]))
write.csv(AllAABLAST, file = "BLAST_seq_AA.txt", col.names = F, row.names = F, quote = F)

####
#Process output of RAxML analyses
####

#ART and Sagebrush
FASTA <- read.FASTA("CLEAN_Sequences_aquaporin_sagebrush/Aquaporin_ARTA_Sagebrush_AA.fasta")
OLDNAMES <- names(FASTA)

RenameTip <- matrix(ncol=3, nrow=length(FASTA))
RenameTip[,1] <- OLDNAMES
RenameTip[,2] <- OLDNAMES
RenameTip[,2] <- gsub("Scaffold", "S", RenameTip[,2])
RenameTip[,3] <- rep('blue', nrow(RenameTip))


toFind <- grep("_ARATH", OLDNAMES)
RenameTip[toFind,3] <- "black"

ART <- OLDNAMES[toFind]
ART <- sapply(strsplit(ART,split="_ARATH") , "[[", 1)
RenameTip[toFind,1] <- paste(ART, "_ARATH", sep='')
RenameTip[toFind,2] <- paste("ARATH_", sapply(strsplit(ART,split="[|]") , "[[", 3), sep='')

#Rename tips

tr <- read.tree("CLEAN_Sequences_aquaporin_sagebrush/RAxML_AA/RAxML_bipartitions_AA_ART_sagebrush.tre")

tr$tip.label <- RenameTip[match(tr$tip.label, RenameTip[,1]),2]
tr <- root.phylo(tr, outgroup = tr$tip.label[match("ARATH_TIP41", tr$tip.label)])
tr <- ladderize(tr)

pdf("RAxML_Aquaporin.pdf")
#Fan
plot(tr, type='fan', align.tip.label=T, cex=.4, use.edge.length = T, tip.color = RenameTip[match(tr$tip.label, RenameTip[,2]),3], no.margin  = T, font=2, x.lim = c(-13,13), y.lim=c(-14,14))
#Unrooted (which is more correct)
#plot(tr, type='unrooted', align.tip.label=T, cex=.4, use.edge.length = T, tip.color = RenameTip[match(tr$tip.label, RenameTip[,2]),3], font=2)
dev.off()
write.tree(tr, file='RAxML_renamed_ART_sagebrush.tre')
plotTree(tr,type="fan",fsize=0.4,lwd=1,
         ftype="i")

#Read in all FASTA file and process to rename tips in tree
FASTAall <- read.FASTA("CLEAN_Sequences_aquaporin_sagebrush/Aquaporin_all.fa")
OLDNAMES <- names(FASTAall)
RenameTip <- matrix(ncol=4, nrow=length(FASTAall))
colnames(RenameTip) <- c("Acc#", "Prot", "Species", "All")
RenameTip[,4] <- OLDNAMES

#Arabidopsis
toFindAt <- grep("_ARATH", OLDNAMES)
RenameTip[toFindAt,3] <- "AraTh"
ART <- OLDNAMES[toFindAt]
ART <- sapply(strsplit(ART,split="_ARATH") , "[[", 1)
RenameTip[toFindAt,1] <- paste(ART, "_ARATH", sep='')
RenameTip[toFindAt,2] <- sapply(strsplit(ART, split="[|]"), "[[", 3)

#Sagebrush
toFindSage <- grep("Scaffold", OLDNAMES)
RenameTip[toFindSage,3] <- "At"
SAGE <- OLDNAMES[toFindSage]
RenameTip[toFindSage,1] <- SAGE
RenameTip[toFindSage,2] <- sapply(strsplit(SAGE,split="_") , "[[", 2)

#Asteraceae samples
toFindAst <- setdiff(seq(from=1, to=nrow(RenameTip), by=1), c(toFindAt, toFindSage))
AST <- OLDNAMES[toFindAst]
RenameTip[toFindAst,1] <- sapply(strsplit(AST, split=' '), "[[",1)
RenameTip[toFindAst,3] <- gsub("]","", sapply(strsplit(AST, split='\\['), "[[",2))

write.csv(RenameTip, "CLEAN_Sequences_aquaporin_sagebrush/RAxML_AA/Rename_all.csv", row.names = F, quote = F)

#Rename tips in all tree
RenameTip <- read.csv("CLEAN_Sequences_aquaporin_sagebrush/RAxML_AA/Rename_all.csv")
RenameTip$Newnames <- paste(RenameTip$Prot, RenameTip$Species, RenameTip$Acc., sep='_')
RenameTip$Color <- rep("black", nrow(RenameTip))
RenameTip$Color[which(RenameTip$Species == "At")] <- "blue"

tr <- read.tree("CLEAN_Sequences_aquaporin_sagebrush/RAxML_AA/RAxML_bipartitions_AA_all.tre")
tr <- ladderize(tr)
toFindtips <- match(tr$tip.label, RenameTip[,1])
tr$tip.label <- RenameTip[toFindtips,7]
tr$node.label[which(as.numeric(tr$node.label) < 50)] <- ""
write.tree(tr, file = "CLEAN_Sequences_aquaporin_sagebrush/RAxML_AA/renamed_RAxML_bipartitions_AA_all.tre")
colEdge <- rep("black", nrow(tr$edge))
colEdge[match(grep("At", tr$tip.label), tr$edge[,2])] <- "blue"
pdf("RAxML_Aquaporin_all_trial.pdf")
par(mfrow=c(1,2))
plot(tr, align.tip.label=F, cex=.1, use.edge.length = F, font=2, no.margin = T, label.offset = 0, show.tip.label = F, edge.color = colEdge, tip.color = RenameTip[toFindtips,8], show.node.label = F, edge.width=0.5)

species.order <- tr$tip.label[tr$edge[tr$edge[,2] <= length(tr$tip.label),2]]

plot(x=0, y=10, xlim=c(0, 10), ylim=c(0, length(tr$tip.label)), type='n', axes = F, xlab="", ylab="")

#Add all At prot sequences
#text(y=1+grep("At", species.order), x=0, labels=species.order[grep("At", species.order)], cex=.4, col='blue', adj=0)

#SIP
segments(x0=-0.2,y0=403, x1=-0.2, y1=434)
text(x=1, y=(403+434)/2, labels="SIP", srt=0, cex=, adj=1)

#TIP
segments(x0=-0.2,y0=400, x1=-0.2, y1=273)
text(x=1, y=(400+273)/2, labels="TIP", srt=0, cex=, adj=1)

#NIP
segments(x0=-0.2,y0=270, x1=-0.2, y1=157)
text(x=1, y=(270+157)/2, labels="NIP", srt=0, cex=, adj=1)

#PIP
segments(x0=-0.2,y0=155, x1=-0.2, y1=1)
text(x=1, y=(155+1)/2, labels="PIP", srt=0, cex=, adj=1)

dev.off()

#####
#Figure with phylogenetic relationship and subcellular location of protein
#####
require(ape)
tr <- read.tree(file = "CLEAN_Sequences_aquaporin_sagebrush/RAxML_AA/renamed_RAxML_bipartitions_AA_all.tre")
#Load csv file
csv <- read.csv("all_BLAST_PCR_scaffolds_ORFs_NPA_TMH_PSORT_30Jan20.csv")
 
#Subset tree to only include artemisia and Arabidopsis accessions
toFind <- grep(paste(c("At","Aa","AraTh"),collapse="|"), tr$tip.label, value=TRUE)
trSage <- keep.tip(tr, toFind)

#Color branches
colEdge <- rep("black", nrow(trSage$edge))
colEdge[match(grep("At_", trSage$tip.label), trSage$edge[,2])] <- "blue"
colEdge[match(grep("Aa_", trSage$tip.label), trSage$edge[,2])] <- "green"

#Subset csv to tips
csvTree <- subset(csv, csv$IdTree %in% trSage$tip.label[grep("At_", trSage$tip.label)])
#Extract subcell localization
csvTree$loc <- sapply(strsplit(as.vector(csvTree$WoLF_PSORT), split='[:]'), "[[", 1)
csvTree$loc <- gsub("cysk", "cyto", csvTree$loc)

csvTree$locCol <- csvTree$loc
csvTree$locCol[which(csvTree$locCol == "chlo")] <- "green"
csvTree$locCol[which(csvTree$locCol == "cyto")] <- "blue"
csvTree$locCol[which(csvTree$locCol == "E.R.")] <- "red"
csvTree$locCol[which(csvTree$locCol == "extr")] <- "pink"
csvTree$locCol[which(csvTree$locCol == "mito")] <- "yellow"
csvTree$locCol[which(csvTree$locCol == "nucl")] <- "black"
csvTree$locCol[which(csvTree$locCol == "plas")] <- "grey"
csvTree$locCol[which(csvTree$locCol == "vacu")] <- "brown"
csvTree$locCol[which(csvTree$locCol == "NA")] <- ""

trSage$tip.label <- gsub("_At", "", trSage$tip.label)


#Plot of tree and subcellular localization
pdf("Tree_Aquaporin_sagebrush_annua.pdf")
plot(trSage, use.edge.length=F, align.tip.label = F, show.node.label = F, edge.color = colEdge, cex=.3, font=1, no.margin = F, label.offset = 5)
tiplabels(pch=16, col=csvTree$locCol, cex=1.4)
nodelabels(trSage$node.label, adj=c(1, 1), frame = 'none', cex=.4)
dev.off()
#SIP
segments(x0=16.5,y0=length(trSage$tip.label)+1, x1=16.5, y1=length(trSage$tip.label)-2)
text(x=16.7, y=length(trSage$tip.label)-1, labels="SIP", cex=.6, srt=90, adj=0.5)
#TIP
segments(x0=16.5,y0=length(trSage$tip.label)-3, x1=16.5, y1=length(trSage$tip.label)-18)
text(x=16.7, y=length(trSage$tip.label)-10, labels="TIP", cex=.6, srt=90, adj=0.5)
#NIP
segments(x0=16.5,y0=length(trSage$tip.label)-19, x1=16.5, y1=length(trSage$tip.label)-29)
text(x=16.7, y=length(trSage$tip.label)-24, labels="NIP", cex=.6, srt=90, adj=0.5)

#PIP
segments(x0=16.5,y0=length(trSage$tip.label)-30, x1=16.5, y1=1)
text(x=16.7, y=length(trSage$tip.label)-37, labels="PIP", cex=.6, srt=90, adj=0.5)

# Add a legend
legend("bottomleft", legend=unique(csvTree$loc)[1:8],
       pch=16, col=unique(csvTree$locCol)[1:8], cex=0.7)
dev.off()

#NIP
segments(x0=-0.2,y0=270, x1=-0.2, y1=157)
text(x=1, y=(270+157)/2, labels="NIP", srt=0, cex=, adj=1)

#PIP
segments(x0=-0.2,y0=155, x1=-0.2, y1=1)
text(x=1, y=(155+1)/2, labels="PIP", srt=0, cex=, adj=1)




plot(, cex=.1)
###########
##Functions
###########
#A function to findORFs in scaffold file (FASTA format) and produce AA sequences as well as report
findORFsTranslateDNA2AA <- function(scaffold, scaffoldID){
  #Extract scaffold sequence (and convert to right format)
  seqRaw <- scaffold[c(grep(paste(">", scaffoldID, sep=''), scaffold)+1)] 
  seqs <- DNAStringSet(seqRaw)
  
  #positive strands 
  pos <-findORFs(seqs, startCodon = "ATG", minimumLength = 40) 
  #negative strands (DNAStringSet only if character) 
  
  neg <-findORFs(reverseComplement(DNAStringSet(seqs)), startCodon = "ATG", minimumLength = 40) 
  pos <- relist(c(GRanges(pos,strand = "+"),GRanges(neg,strand = "-")),skeleton = merge(pos,neg))
  
  #Process output
  pos <- as.data.frame(pos)
  pos$scaffoldID <- rep(scaffoldID, nrow(pos))
  pos$ORFID <- paste("ORF_", seq(from=1, to=nrow(pos)), sep='')
  
  #Extract and write ORFs from seq and produce AA sequence
  for(i in 1:nrow(pos)){
    if(as.character(pos$strand[i]) == "+"){
      extSeq <- as.DNAbin(DNAStringSet(paste(strsplit(seqRaw, split='')[[1]][c(pos$start[i]:pos$end[i])], collapse = "")))
      AAseq <- paste(paste(">", pos$scaffoldID[i], "_", pos$ORFID[i], sep=''), paste(as.character(trans(extSeq, code = 1, codonstart = 1))[[1]], collapse = ''), sep='\n')
      write.table(AAseq, file = paste("AA_ORFs/", pos$scaffoldID[i], "_", pos$ORFID[i], ".fa", sep=''), col.names = F, row.names = F, quote = F)
    }
    if(as.character(pos$strand[i]) == "-"){
      revComp <- reverseComplement(DNAStringSet(seqs))
      
      extSeq <- as.DNAbin(DNAStringSet(paste(strsplit(as.vector(revComp), split='')[[1]][c(pos$start[i]:pos$end[i])], collapse = "")))
      AAseq <- paste(paste(">", pos$scaffoldID[i], "_", pos$ORFID[i], sep=''), paste(as.character(trans(extSeq, code = 1, codonstart = 1))[[1]], collapse = ''), sep='\n')
      write.table(AAseq, file = paste("AA_ORFs/", pos$scaffoldID[i], "_", pos$ORFID[i], ".fa", sep=''), col.names = F, row.names = F, quote = F)
    }
  }
  
  #Merge all ORFs for BLAST analysis
  system(paste("cat AA_ORFs/", pos$scaffoldID[i], "* > AA_ORFs/", pos$scaffoldID[i], "_ORFs.fa", sep=''))
  system(paste("rm AA_ORFs/", pos$scaffoldID[i], "_ORF_*", sep=""))
  #Save pos file
  write.table(pos, file = paste("ORFs_report/",pos$scaffoldID[i], "_", "ORFs.csv", sep=''), col.names = T, row.names = F, quote = F)
} 
