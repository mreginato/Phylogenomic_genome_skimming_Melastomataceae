library(ape)
library(phyloch)
library(phangorn)

wd = "C:/SkimmingLoci_reduced/4_coal_sims"
setwd(wd)

##########################################
### import trees
##########################################

read.tree("speciestree.tre") -> sptree
read.tree("concatenate.phy.raxml.bestTree_rooted.tre") -> raxml

ladderize(sptree) -> sptree
ladderize(raxml) -> raxml

fixNodes(sptree) -> sptree
fixNodes(raxml) -> raxml

### check topology

RF.dist(sptree, raxml)
sptree$tip.label == raxml$tip.label

### flag polytomies in astral

pol =  0.001

which(sptree$edge.length <= pol) -> polys

##########################################
### get theta
##########################################

raxml$edge.length/sptree$edge.length -> t
#(raxml$edge.length - sptree$edge.length)/(sptree$edge.length+raxml$edge.length) -> t
#t[polys] <- 0
round(range(t),2)

range(raxml$edge.length)
range(sptree$edge.length)

sptree -> tree -> tree.log
boxplot(t)

log(t+1) -> t2
tree.log$edge.length <- t2
tree$edge.length <- NULL

par(mfrow=c(1,3))
plot(sptree, show.tip.label=F)
title("Coalescent units")
plot(raxml, show.tip.label=F)
title("Mutation unit")
plot(tree.log, show.tip.label=F)
title("Theta")

### theta node

c((Ntip(tree)+1):(Ntip(tree)+Nnode(tree))) -> nodes

sapply(nodes, FUN=function(x)(which.edge(tree, x))) -> edges
#t[edges] -> node.thetas
t2[edges] -> node.thetas
tree$node.label <- node.thetas
tree$node.label[1] <- 0

write.tree(tree.log, "theta_tree_log.tre")
write.tree(tree, "theta_tree_nodes.tre")
