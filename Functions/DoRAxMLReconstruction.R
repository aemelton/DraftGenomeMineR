### AE Melton, 2021
# Do RAxML phylogenetic reconstruction

DoRAxMLReconstruction <- function(algorithm = "a", in.file, out.file, model, parsimony.random.seed = 12345, bootstrap.random.seed = 12345, BS.rep.count = 1000){
  
# Paste the peices together and run mafft
raxml.command <- paste("raxmlHPC-SSE3",
                       "-f",
                       algorithm,
                       "-s",
                       in.file,
                       "-n",
                       out.file,
                       "-m",
                       model,
                       "-p",
                       parsimony.random.seed,
                       "-x",
                       bootstrap.random.seed,
                       "-#",
                       BS.rep.count)
system(raxml.command)
}
#