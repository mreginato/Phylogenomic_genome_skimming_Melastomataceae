library(ape)
library(phyloch)
library(phyloU)
library(phytools)
source("F:/R Scripts/_Functions/supportColors.R")

wd = "C:/SkimmingLoci_reduced/4_coal_sims"
genetrees.dir = "C:/SkimmingLoci_reduced/1_gene_trees"
phyparts.dir = "C:/SkimmingLoci_reduced/3_quartet_concordance/phyparts"


setwd("C:/SkimmingLoci_reduced/")

read.csv("Metadata_classification_tips.csv") -> meta
rownames(meta) <- meta$old.terminal
meta[which(meta$tribe != "CAP"),] -> meta

setwd(wd)

#######################################
### Gene trees
#######################################

setwd(genetrees.dir)

read.tree("gene_trees_rooted.trees") -> genes
read.csv("gene_names.csv")[,1] -> genes.labs
names(genes) <- genes.labs

### get rates

gene.rates <- vector(length=length(genes))
names(gene.rates) <- names(genes)

for (i in 1:length(genes)) {
  genes[[i]] -> t0
  sum(t0$edge.length)/Nedge(t0) -> gene.rates[i]
}

#######################################
### Phyparts output
#######################################


setwd(phyparts.dir)

sptreefile = "speciestree.tre"
keyfile = "phyparts_ape_nodes_trans.csv"
read.tree(sptreefile) -> tree
ladderize(tree) -> tree

read.csv(keyfile) -> nodes.key
rbind(c(999, 275), nodes.key) -> nodes.key

### import concor genes

list.files(pattern="concord.node") -> concon.files
list.files(pattern="conflict.node") -> confli.files

nodes.key -> summary.out
summary.out$concord_rate <- NA
summary.out$conflict_rate <- NA

for (i in 1:length(concon.files)) {
  concon.files[i] -> f0
  as.numeric(strsplit(f0, ".", fixed=T)[[1]][4]) -> node0
  read.table(f0)[,1] -> genes0
  unlist(lapply(strsplit(genes0, "/"), FUN=function(x)(x[[length(x)]]))) -> genes0
  sub(".tre$", "", genes0) -> genes0
  match(genes0, names(gene.rates)) -> x
  median(gene.rates[x]) -> r0
  summary.out$concord_rate[match(node0, summary.out$node)] <- r0
}


for (i in 1:length(confli.files)) {
  confli.files[i] -> f0
  as.numeric(strsplit(f0, ".", fixed=T)[[1]][4]) -> node0
  read.table(f0)[,1] -> genes0
  unlist(lapply(strsplit(genes0, "/"), FUN=function(x)(x[[length(x)]]))) -> genes0
  sub(".tre$", "", genes0) -> genes0
  match(genes0, names(gene.rates)) -> x
  median(gene.rates[x]) -> r0
  summary.out$conflict_rate[match(node0, summary.out$node)] <- r0
}

summary.out$concord_rate/summary.out$conflict_rate -> summary.out$rate_ratio
summary.out


#######################################
### Plot - All
#######################################

tree$node.label <- summary.out$rate_ratio

setwd(wd)

range(summary.out$rate_ratio, na.rm = T)

cols <- rev(c("red", "lightcoral", "skyblue", "darkblue"))
qc.breaks <- c(0,0.6, 1.2, 1.8, 2.4)
supportColors(tree$node.label, breaks = qc.breaks, colors = cols) -> qc.cols
qc.cols$support -> cols.nodes
#cols.nodes[is.na(cols.nodes)] <- "snow"

plot(tree, cex=0.5, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=cols.nodes, cex=0.5)
title("Median rate ratio (concord/discord)")
legend("bottomleft", legend=qc.cols$cuts, fill=qc.cols$cols)



pdf("Rates_concordance-all.pdf", width=11, height=14)
plot(tree, edge.width = 2, label.offset = 0.1, cex=0.4)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=cols.nodes, cex=0.5)
title("Median rate ratio (concord/discord)")
legend("bottomleft", legend=qc.cols$cuts, fill=qc.cols$cols)
dev.off()

tree -> dtree
round(dtree$node.label,2) -> dtree$node.label
dtree$node.label[is.na(dtree$node.label)] <- 0

write.tree(dtree, "Rates_concordance_median_rate.tre")


###################################################
### Tribes summary
###################################################

apply(round(summary.out[,3:5],2), MARGIN=1, paste, collapse="|") -> tree$node.label

tree$node.label

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

do.call(rbind, strsplit(tree$node.label, "|", fixed=T)) -> x
colnames(x) <- c("concord_rate", "conflict_rate", "rate_ratio")

data.frame(type=types, x) -> summary.out

write.csv(summary.out, "Rates_concordance_summary_with_node_info.csv", row.names=F)


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


### plot

as.numeric(unlist(lapply(strsplit(btree$node.label, "|", fixed=T), "[", 3))) -> btree$node.label

cols <- rev(c("red", "lightcoral", "skyblue", "darkblue"))
qc.breaks <- c(0,0.6, 1.2, 1.8, 2.4)
supportColors(btree$node.label, breaks = qc.breaks, colors = cols) -> qc.cols
qc.cols$support -> cols.nodes
#cols.nodes[is.na(cols.nodes)] <- "snow"

plot(btree, cex=0.5, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=cols.nodes, cex=0.5)
title("Median rate ratio (concord/discord)")
legend("bottomleft", legend=qc.cols$cuts, fill=qc.cols$cols)

plot(btree.2, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=cols.nodes, cex=0.5)
title("Median rate ratio (concord/discord)")
legend("bottomleft", legend=qc.cols$cuts, fill=qc.cols$cols)




cairo_pdf("Rates_concordance_tribes.pdf", width=6, height=6, onefile=T)

plot(btree, cex=0.5, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=cols.nodes, cex=0.5)
title("Median rate ratio (concord/discord)")

plot(btree.2, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=cols.nodes, cex=0.5)
title("Median rate ratio (concord/discord)")

plot(1, xlab="", ylab="", axes=F, col="white")
legend("bottomleft", legend=qc.cols$cuts, fill=qc.cols$cols)

dev.off()

cairo_pdf("Rates_concordance_tribes_narrow.pdf", width=4.5, height=6, onefile=T)
plot(btree, cex=0.5, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=cols.nodes, cex=0.5)
title("Median rate ratio (concord/discord)")

plot(btree.2, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=cols.nodes, cex=0.5)
title("Median rate ratio (concord/discord)")

plot(1, xlab="", ylab="", axes=F, col="white")
legend("bottomleft", legend=qc.cols$cuts, fill=qc.cols$cols)
dev.off()
