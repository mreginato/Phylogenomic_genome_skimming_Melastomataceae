library(ape)
library(phyloch)
library(phyloU)
library(phytools)

setwd("C:/SkimmingLoci_reduced/3_quartet_concordance")
getwd() -> wd
paste(wd, "/phyparts", sep="") -> out.dir

setwd("C:/SkimmingLoci_reduced/")

read.csv("Metadata_classification_tips.csv") -> meta
rownames(meta) <- meta$old.terminal
meta[which(meta$tribe != "CAP"),] -> meta

setwd(wd)

#######################################
### Phyparts output
#######################################

### nodelabels
### n supporting, % -->  c("Support", "Conflict (most common)", "Conflict (others)", "Neutral")

setwd(out.dir)

total_genes = 684

sptreefile = "speciestree.tre"
keyfile = "out.node.key"
concon_treefile = "out.concon.tre"
pdistfile = "phyparts_dist.csv"
piesfile = "phyparts_pies.csv"

read.tree(concon_treefile)[[1]] -> concon
read.tree(sptreefile) -> sptree
ladderize(concon) -> concon
ladderize(sptree) -> sptree

sptree$tip.label == concon$tip.label
sptree$node.label <- concon$node.label
sptree -> concon

read.table(keyfile) -> keyfile
colnames(keyfile) <- c("node", "taxa")

read.csv(pdistfile) -> pdist
read.csv(piesfile) -> ppies

setwd(wd)

### Get ape node

keyfile$node.ape <- NA

for (i in 1:nrow(keyfile)) {
  keyfile$taxa[i] -> t0
  gsub("\\)[0-9]*", " ", t0) -> t0
  gsub("* .[0-9]", " ", t0) -> t0
  gsub("* [0-9] *", " ", t0) -> t0
  gsub("(", "", t0, fixed=T) -> t0
  gsub(")", "", t0, fixed=T) -> t0
  strsplit(t0, ",")[[1]] -> t0
  trimws(t0) -> t0
  unlist(lapply(strsplit(t0, " "), "[", 1)) -> t0
  getMRCA(concon, t0) -> keyfile$node.ape[i]
}

keyfile[,c(1,3)]
keyfile[order(keyfile$node.ape),c(1,3)] -> keyfile

ppies[match(keyfile$node, ppies$node),] -> ppies
pdist[match(keyfile$node, pdist$node),] -> pdist

rbind(rep(NA,4), ppies[,-1]) -> piedat

### Export tree with number of genes and piedata

concon$node.label -> n
apply(piedat, MARGIN=1, paste, collapse="|") -> x
paste(n, x, sep="|") -> nodelabs

concon$node.label <- nodelabs

write.tree(concon, "phyparts_concon_labels.tre")
write.csv(keyfile, "phyparts_ape_nodes_trans.csv", row.names = F)


#######################################
### Plot - All
#######################################

setwd(wd)

# Blue: Support the shown topology
# Green: Conflict with the shown topology (most common conflicting bipartion)
# Red: Conflict with the shown topology (all other supported conflicting bipartitions)
# Gray: Have no support for conflicting bipartion

## Numbers:
## Above = support

cols <- c("blue", "red", "lightcoral", "gray90")
leg.txt <- c("Support", "Conflict (most common)", "Conflict (others)", "Neutral")


pdf("phyparts_conconplot.pdf", width=8.5, height=11)
par(mar=c(0,0,0,0))
plot(concon, edge.width = 2, label.offset = 0.1, cex=0.4)
title("Phyparts")
#nodelabels(concon$node.label, frame="none", adj=c(1.7,1.7), cex=0.8)
#nodelabels(round(as.numeric(trees[[i]]$node.label)/total_genes,2)*100, frame="none", adj=c(2,2), cex=0.8)
nodelabels(pie=piedat, piecol = cols, cex=0.5)
legend("bottomleft", legend=leg.txt, fill=cols, box.col=NA)
dev.off()

system("open phyparts_conconplot.pdf")


###################################################
### Tribes summary
###################################################

concon -> tree
match(tree$tip.label, meta$terminal) -> x
which(is.na(x)) -> drop
drop.tip(tree, drop) -> tree

tree$tip.label -> tips
meta[which(!duplicated(meta$terminal)),] -> meta
which(is.na(match(tips, meta$terminal))) 
rownames(meta) <- meta$terminal
meta[tree$tip.label,] -> meta

sort(unique(meta$tribe)) -> tribes

tribes.dat <- data.frame(tribe=tribes, node=NA, internal.nodes=NA)

### get tribes nodes

for (i in 1:length(tribes)) {
  tribes[i] -> t0
  meta$terminal[which(meta$tribe == t0)] -> tips0
  if (length(tips0) > 1) {
    if (is.monophyletic(tree, tips0) == F) {
      stop("not monophyletic")
    }
    getMRCA(tree, tips0) -> node0
    node0 -> tribes.dat$node[i]
    descendants(tree, node0, type="both") -> i0
    paste(i0, collapse=",") -> tribes.dat$internal.nodes[i]
  } else {
    match(tips0, tree$tip.label) -> n0
    n0 -> tribes.dat$node[i]
    n0 -> tribes.dat$internal.nodes[i]
  }
  
}
tribes.dat

### figure backbone

(Ntip(tree)+1):(Ntip(tree)+Nnode(tree)) -> nodes

tribes.dat[,1:2] -> nodes.dat
nodes -> backbone
backbone[which(is.na(match(backbone, tribes.dat$node)))] -> backbone
as.numeric(unlist(strsplit(tribes.dat$internal.nodes, ","))) -> internal.nodes
backbone[which(is.na(match(backbone, internal.nodes)))] -> backbone

plot(tree, use.edge.length=F)
nodelabels(node=backbone)
nodes -> types
types[which(is.na(match(nodes, internal.nodes))==F)] <- "Internal nodes"
types[which(is.na(match(nodes, backbone))==F)] <- "Backbone"
types[match(nodes.dat$node[-17], nodes)] <- nodes.dat$tribe[-17]

tree$node.label -> piedat
do.call(rbind, strsplit(piedat, "|", fixed=T)) -> piedat
mode(piedat) <- "numeric"
colnames(piedat) <- c( "concord_absolute", "concord", "most_conflict", "other_conflict", "the_rest")

data.frame(type=types, piedat) -> summary.out

write.csv(summary.out, "phyparts_summary_with_node_info.csv", row.names=F)

###################################################
### Plot (tribes)
###################################################

which(duplicated(meta$tribe)) -> duplis
drop.tip(tree, duplis) -> dtree
ladderize(dtree) -> dtree
fixNodes(dtree) -> dtree
checkPhylo(dtree, meta)$data -> dclades
table(meta$tribe) -> clades.n
dclades$n[match(names(clades.n), dclades$tribe)] <- clades.n
dtree$tip.label == rownames(dclades)

rep(0.2, nrow(dclades)) -> depths


trans <- data.frame(tip.label=rownames(dclades), clade.label=dclades$tribe,
                    N=rep(5,nrow(dclades)), depth=depths)

trans.2 <- data.frame(tip.label=rownames(dclades), clade.label=dclades$tribe,
                      N=dclades$n, depth=depths)

dtree -> dtree2
dtree2$edge.length[] <- 1
chronopl(dtree2, 1) -> dtree2

phylo.toBackbone(dtree2, trans) -> btree
phylo.toBackbone(dtree, trans.2) -> btree.2
plot(btree)
plot(btree.2)

btree$node.label -> piedat
do.call(rbind, strsplit(piedat, "|", fixed = T)) -> piedat
mode(piedat) <- "numeric"
piedat[,-1] -> piedat

### plot

cols <- c("blue", "red", "lightcoral", "gray90")
leg.txt <- c("Support", "Conflict (most common)", "Conflict (others)", "Neutral")

plot(btree, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(pie=piedat, piecol = cols, cex=0.5)
title("Phyparts")
legend("bottomleft", legend=leg.txt, fill=cols, box.col=NA, title="% of genes")

plot(btree.2, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(pie=piedat, piecol = cols, cex=0.5)
title("Phyparts")
legend("bottomleft", legend=leg.txt, fill=cols, box.col=NA, title="% of genes")





cairo_pdf("phyparts_tribes.pdf", width=6, height=6, onefile=T)

plot(btree, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(pie=piedat, piecol = cols, cex=0.5)
title("Phyparts")

plot(btree.2, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(pie=piedat, piecol = cols, cex=0.5)
title("Phyparts")
plot(1, xlab="", ylab="", axes=F, col="white")
legend("bottomleft", legend=leg.txt, fill=cols, box.col=NA, title="% of genes")

dev.off()

cairo_pdf("phyparts_tribes_narrow.pdf", width=4.5, height=6, onefile=T)
plot(btree, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(pie=piedat, piecol = cols, cex=0.5)
title("Phyparts")

plot(btree.2, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(pie=piedat, piecol = cols, cex=0.5)
title("Phyparts")
plot(1, xlab="", ylab="", axes=F, col="white")
legend("bottomleft", legend=leg.txt, fill=cols, box.col=NA, title="% of genes")
dev.off()
