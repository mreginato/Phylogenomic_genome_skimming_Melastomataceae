library(ape)
library(phyloch)
library(phytools)
library(skimmingLociR)
library(phangorn)
library(xlsx)
source("F:/R Scripts/_Functions/supportColors.R")

###################################################
### Dirs
###################################################

astral.dir = "C:/SkimmingLoci_reduced/2_Astral"
singles.dir = "C:/SkimmingLoci_reduced/1_gene_trees"

setwd(astral.dir)

###################################################
### Astral - All
###################################################

run.astral = T

jar.folder = "C:/Programas/Astral 5.7.8"

### Conc trees

setwd(singles.dir)

read.tree("gene_trees_rooted.trees") -> trees

setwd(astral.dir)
write.tree(trees, "gene_trees.trees")

### Run Astral

setwd(jar.folder)
file.copy(paste(astral.dir, "/gene_trees.trees", sep=""), paste(jar.folder, "/gene_trees.trees", sep=""), overwrite = T)
call <- paste("java -jar astral.5.7.8.jar -i gene_trees.trees -o speciestree.tre")
call2 <- paste("java -jar astral.5.7.8.jar -i gene_trees.trees -o speciestree_geneboot.tre --gene-only")
call
call2

if (run.astral) {
  shell(call)
  file.copy(paste(jar.folder, "/speciestree.tre", sep=""), paste(astral.dir, "/speciestree.tre", sep=""), overwrite = T)
  shell(call2)
  file.copy(paste(jar.folder, "/speciestree_geneboot.tre", sep=""), paste(astral.dir, "/speciestree_geneboot.tre", sep=""), overwrite = T)
}

### Tree

setwd(astral.dir)
readNexus("speciestree_rooted.tre", "raxml") -> astral.reg
readNexus("speciestree_geneboot_rooted.tre", "raxml") -> astral.gene

RF.dist(astral.reg, astral.gene)

ladderize(astral.gene) -> astral
plot(astral, cex=0.5)

as.numeric(astral$node.label) -> astral$node.label

### plot

cuts <- c("[50,80]", "[81,95]", "[96,100]")

dev.off()

fixNodes(astral) -> astral

plot(astral, edge.width = 1.5, edge.color = "gray80", cex=0.6, use.edge.length=F)

### Pdf

cairo_pdf("Astral_gene_support.pdf", height = 7, width=12, onefile = T)

astral -> t1
supportColors(as.numeric(t1$node.label), breaks = c(50, 80, 95, 100)) -> cols
c(1,2) -> miss
(Ntip(t1)++1):(Ntip(t1)++Nnode(t1)) -> nodes
cols$support -> c0
c0[-miss] -> c0
nodes[-miss] -> nodes
plot(astral, edge.width = 2, edge.color = "gray80", cex=0.4)
nodelabels(text=rep("", length(c0)), node=nodes, frame="circ", bg=c0, cex=0.8)
title("Astral - raw")

legend("bottomleft", fill=, legend=cuts, title="Support", pch=21, pt.bg=cols$cols, box.col=NA, bg=NA, pt.cex = 2)
dev.off()





