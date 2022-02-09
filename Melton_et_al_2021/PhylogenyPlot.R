### AE Melton, 2020
# Generate figure showing NPA motif evolution in AQPs

#
#install.packages("ape")
#install.packages("phytools")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ggtree")
#install.packages('treeio')
#install.packages("evobiR")
#

#
library(ape)
library(ggtree)
library(phytools)
library(dplyr)
library(evobiR)
#

############################# 8/19
label.info <- read.csv(file = "OutFiles_AEM/Rename_all_AEM_V2.csv")
label.info
tree.id <- label.info$TreeID
Prot.name <- label.info$Prot
Species <- label.info$Species
df <- data.frame(tree.id, Prot.name, Species)
df <- subset(x = df, Species == "At" | Species == "Aa" | Species == "AraTh")
df

for(i in 1:nrow(df)){
        if (df$Species[i] == "AraTh") {
                df$Symbol[i] <- 0
        }
        if (df$Species[i] == "Aa") {
                df$Symbol[i] <- 2
        }
        if (df$Species[i] == "At") {
                df$Symbol[i] <- 1
        }
}
df

for(i in 1:nrow(df)){
        df$Subfamily[i] <- gsub(x = df$Prot.name[i], pattern = "[0-9].*", replacement = "")
}
df

for(i in 1:nrow(df)){
        if (df$Subfamily[i] == "NIP") {
                df$Col[i] <- "blue"
        }
        if (df$Subfamily[i] == "PIP") {
                df$Col[i] <- "yellow"
        }
        if (df$Subfamily[i] == "SIP") {
                df$Col[i] <- "pink"
        }
        if (df$Subfamily[i] == "TIP") {
                df$Col[i] <- "green"
        }
        if (df$Subfamily[i] == "Unk") {
                df$Col[i] <- "grey"
        }
}
df

write.csv(x = df, file = "OutFiles_AEM/RAxML_AEM/table_for_making_plot_V2.csv")
###

####
setwd("~/Dropbox/BSU_Research/Aquaporin/")
tree <- read.tree(file = "OutFiles_AEM/RAxML_AEM/AQP_RaxML_GUI_V5_same_as_V4_but_removed_bad_genes/RAxML_bestTree.AQP_aln_V3_EDITS.tre")
csv <- read.csv(file = "OutFiles_AEM/RAxML_AEM/table_for_making_plot_V3.csv")
csv
dummy <- tree # Let's make a duplicate file, 'cause I know I am going to mess this up.
#

#
# 
dummy$tip.label
csv$tree.id

plot(dummy, show.tip.label = FALSE)
tiplabels(pch = csv$Symbol, col = csv$Col) # OK! This worked, but they aren't in the same order (see above). dundundun....

#csv.trim <- csv[(dummy$tip.label %in% csv$tree.id),]
dummy <- midpoint.root(tree = dummy)
dummy.order <- ReorderData(tree = dummy, data = csv, taxa.names = 1)

df <- as.data.frame(dummy.order)

dummy$tip.label
dummy.order$tree.id

pdf("~/Dropbox/Manuscripts/BSU/AQP/AQP_PHYLO_GGTREE.pdf")
p <- ggtree(dummy, layout = "circular")
p <- p %<+% df + geom_tippoint(aes(color = Subfamily, shape = Species)) + theme(legend.position = c(0.1,0.25))
plot(p)
dev.off()
