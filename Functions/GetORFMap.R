### AE Melton, 2021
# Generate a map of ORFs for a given gene along a scaffold

GetORFMap <- function(scaffold.file, ORF.data.file, user.sep, out.file){

scaffold <- readLines(con = scaffold.file)
scaffold.length <- nchar(scaffold[2])
ORF.data <- read.csv(ORF.data.file, sep = user.sep)
ORF.data.sort <- ORF.data[order(ORF.data$start),]

#
Not.ORF.width <- NULL
Not.ORF.start <- NULL
first.NON.ORF.end <- ORF.data.sort$start[1] - 1
first.NON.ORF.start <- 1
first.NON.ORF.width <- ORF.data.sort$start[1] - 1
first.NON.ORF <- data.frame(Category = "Not.ORF", start = first.NON.ORF.start, width = first.NON.ORF.width)

for(i in 1:(nrow(ORF.data.sort) - 1)){
Not.ORF.width[i] <- (ORF.data.sort$start[i+1] - 1) - ORF.data.sort$end[i]
Not.ORF.start[i] <- ORF.data.sort$end[i] + 1
}

Not.ORF <- NULL
Not.ORF.df <- NULL
Not.ORF$Category <- rep(x = "Not.ORF", times = length(Not.ORF.width))
Not.ORF$start <- Not.ORF.start
Not.ORF$width <- Not.ORF.width
Not.ORF.df <- data.frame(Category = Not.ORF$Category, start = Not.ORF$start, width = Not.ORF$width)
Not.ORF.df <- rbind(Not.ORF.df, first.NON.ORF)
#Not.ORF.df
#

#
ORF <- NULL
ORF$Category <- rep(x = "ORF", times = nrow(ORF.data.sort))
ORF$start <- ORF.data.sort$start
ORF$width <- ORF.data.sort$width
ORF.df <- data.frame(Category = ORF$Category, start = ORF$start, width = ORF$width)
#

#
big.df <- NULL
big.df <- rbind(Not.ORF.df, ORF.df)
#

#
big.df.sort <- big.df[order(big.df$start),]
#

#
pdf(out.file)
barplot(as.matrix(big.df.sort$width), horiz = T, col = c("red","black"))
dev.off()

}
