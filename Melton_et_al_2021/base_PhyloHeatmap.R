### AE Melton, 2020
# Phylo + heatmap plot

#
library(ape)
library(ggtree)
library(phytools)
library(dplyr)
library(evobiR)
#

#
setwd("~/Dropbox/BSU_Research/Aquaporin/OutFiles_AEM/")
#

#
tree <- read.tree(file = "RAxML_AEM/AQP_RaxML_GUI_V5_same_as_V4_but_removed_bad_genes/RAxML_bestTree.AQP_aln_V3_EDITS.tre")
csv <- read.csv("Promoters/Promoter_Cat_Count_by_Scaff_EDITS.csv", row.names = 1, header = T)[,-30]
head(csv)
csv <- as.data.frame(t(csv))
nrow(csv)
rownames(csv) <- gsub(replacement = "-", x = rownames(csv), pattern = "\\.")
scaff.names <- rownames(csv)
dumdum <-  keep.tip(phy = tree, tip = scaff.names)
dummy <- midpoint.root(tree = dumdum)
dummy.order <- ReorderData(tree = dumdum, data = csv)
#

setwd("~/Desktop/")
#
dummy$tip.label <- " "
par(mfrow=c(1,2))
plot.phylo(x = dumdum, show.tip.label = T, align.tip.label = T)
heatmap(as.matrix(dummy.order), margins = c(10,10), col = cm.colors(256), scale = "row")

#