### AE Melton, 2021
# Do some MAFFT alignments ... in R!!!

DoMafftAlignment <- function(align.method = "--auto", in.file, out.file){

# Paste the peices together and run mafft
mafft.command <- paste("mafft", align.method, in.file, " > ", out.file)
system(mafft.command)
#
}
