library(ape)


################################
### dirs
################################

wd = "C:/SkimmingLoci_reduced/3_quartet_concordance/phyparts"
sp.dir = "C:/SkimmingLoci_reduced/3_quartet_concordance/"
genes.dir  = "C:/SkimmingLoci_reduced/1_gene_trees"
genes.out = "C:/SkimmingLoci_reduced/3_quartet_concordance/phyparts/genetrees"

################################
### trees
################################

setwd(sp.dir)

read.tree("speciestree_geneboot_rooted.tre") -> tree

setwd(genes.dir)
read.csv("gene_names.csv") -> labs
read.tree("gene_trees.trees") -> genes
labs[,1] -> labs

################################
### Export
################################

setwd(wd)

write.tree(tree, "speciestree.tre")

setwd(genes.out)

for (i in 1:length(genes)) {
  genes[i] -> t0
  write.tree(t0, file=paste(labs[i], ".tre"))
}

################################
### call
################################

### java -Xmx16g -jar ~/phyparts/phyparts.jar -a 1 -v -d /mnt/c/SkimmingLoci_reduced/3_quartet_concordance/phyparts/genetrees -m /mnt/c/SkimmingLoci_reduced/3_quartet_concordance/phyparts/speciestree.tre -o phyparts
