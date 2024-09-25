library(ape)
library(phyloch)

wd = "C:/SkimmingLoci_reduced/4_coal_sims/seq-gen"
output.dir = "C:/SkimmingLoci_reduced/4_coal_sims"
setwd(wd)

##########################################
### import simulate seqs
##########################################

list.files(pattern=".fas") -> files


## align

for (i in 1:length(files)) {
  files[i] -> f0
  sub(".fas", ".phy", f0) -> f1
  call <- paste0("muscle -in ", f0, " -out temp.fas")
  system(call)
  read.dna("temp.fas", "fasta") -> a1
  write.phy(a1, file=f1)
  unlink("temp.fas")
  cat("\r",i)
}

##########################################
###  Run RAxML
##########################################

boot = 10
Th = 4

names(aligns) -> files
sub(".phy", "", files) -> labs

for (i in 1:length(aligns)) {
  cat("\r", i)
  call <- paste("raxmlHPC-PTHREADS -f a -m GTRGAMMA -p 12345 -x 12345 -T", Th, "-N", boot, "-s", files[i], "-n", labs[i])
  system(call, show.output.on.console = F)
  unlink(list.files(pattern=".reduced"))
  unlink(list.files(pattern="BranchLabels"))
  unlink(list.files(pattern="RAxML_info"))
  unlink(list.files(pattern="RAxML_bestTree"))
  unlink(files[i])
}


##########################################
### Merge trees and export single file
##########################################

list.files(pattern="RAxML_bipartitions") -> files
sapply(files, read.tree, simplify = F) -> trees
class(trees) <- "multiPhylo"

setwd(output.dir)

file = "sim_gene.trees"

write.tree(trees, file)


### get annotated species tree 

call <- "raxmlHPC -f b -t speciestree.tre -z sim_gene.trees -m GTRGAMMA -n ERR"
system(call)
