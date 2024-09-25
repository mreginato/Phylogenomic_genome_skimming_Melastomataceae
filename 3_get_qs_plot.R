library(ape)
library(phyloch)
library(scales)
library(phytools)
library(phyloU)
source("F:/R Scripts/_Functions/supportColors.R")

###################################################
### Dirs
###################################################

wd = "C:/SkimmingLoci_reduced/3_quartet_concordance"

setwd(wd)

setwd("C:/SkimmingLoci_reduced/")

read.csv("Metadata_classification_tips.csv") -> meta

setwd(wd)

###################################################
### Import concon
###################################################

#read.tree("Astral.concon.tre")[[1]] -> gn.supporting
#read.tree("Astral.concon.tre")[[2]] -> gn.conflicting
#gn.supporting$node.label[1] <- NA
#gn.conflicting$node.label[1] <- NA
#as.numeric(gn.supporting$node.label) -> gn.supporting$node.label
#as.numeric(gn.conflicting$node.label) -> gn.conflicting$node.label

#fixNodes(gn.supporting) -> gn.supporting
#fixNodes(gn.conflicting) -> gn.conflicting

###################################################
### Import files
###################################################

qc <- read.tree("qs_astral.labeled.tre.qc")
qd <- read.tree("qs_astral.labeled.tre.qd")
qi <- read.tree("qs_astral.labeled.tre.qi")

ladderize(qc) -> qc
ladderize(qd) -> qd
ladderize(qi) -> qi

### remove ougroup

in.node = 285

extract.clade(qc,in.node) -> qc
extract.clade(qd,in.node) -> qd
extract.clade(qi,in.node) -> qi

fixNodes(qc) -> qc
fixNodes(qi) -> qi
fixNodes(qd) -> qd

# process node labels of above three labeled trees

# qc tree
qc$node.label <- as.numeric(gsub("qc=","",qc$node.label))
# qd tree
qd$node.label <- as.numeric(gsub("qd=","",qd$node.label))
# qi tree
qi$node.label <- as.numeric(gsub("qi=","",qi$node.label))

# add a customized label for internode or inter-branch, i.e., qc/qd/qI

qc -> tree
score_raw = paste(qc$node.label,"/",qd$node.label,"/",qi$node.label,sep="")
score_raw = gsub("NA/NA/NA","",score_raw)
tree$node.label <- score_raw


#tree$tip.label == gn.conflicting$tip.label

write.tree(tree, "speciestree_qc_annnotated.tre")

###################################################
### Plot (all)
###################################################

range(qc$node.label, na.rm=T)
range(qd$node.label, na.rm=T)
range(qi$node.label, na.rm=T)

cols <- c("red", "lightcoral", "skyblue", "darkblue")
cols2 <- c("red", "skyblue", "darkblue")
qc.breaks <- c(-1,-0.2,0,0.2,1)
qd.breaks <- qi.breaks <- c(0,0.3,0.7,1)


supportColors(qc$node.label, breaks = qc.breaks, colors = cols) -> qc.cols
supportColors(qi$node.label, breaks = qi.breaks, colors = cols2) -> qi.cols
supportColors(qd$node.label, breaks = qd.breaks, colors = cols2) -> qd.cols

plot(tree, cex=0.5, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=qc.cols$support, cex=0.5)
title("Quartet Concordance")
legend("bottomleft", legend=qc.cols$cuts, fill=qc.cols$cols)

plot(tree, cex=0.5, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qi.cols$support)), frame="circ", bg=qi.cols$support, cex=0.5)
title("Quartet Informativeness")
legend("bottomleft", legend=qi.cols$cuts, fill=qi.cols$cols)

plot(tree, cex=0.5, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qd.cols$support)), frame="circ", bg=qd.cols$support, cex=0.5)
title("Quartet Diferential")
legend("bottomleft", legend=qd.cols$cuts, fill=qd.cols$cols)

pdf("Quartet_concordance-all_tips.pdf", width=11, height=14)
plot(tree, cex=0.3, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=T)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=qc.cols$support, cex=0.5)
title("Quartet Concordance")
legend("bottomleft", legend=qc.cols$cuts, fill=qc.cols$cols)

plot(tree, cex=0.3, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=T)
nodelabels(text=rep("", length(qi.cols$support)), frame="circ", bg=qi.cols$support, cex=0.5)
title("Quartet Informativeness")
legend("bottomleft", legend=qi.cols$cuts, fill=qi.cols$cols)

plot(tree, cex=0.3, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=T)
nodelabels(text=rep("", length(qd.cols$support)), frame="circ", bg=qd.cols$support, cex=0.5)
title("Quartet Diferential")
legend("bottomleft", legend=qd.cols$cuts, fill=qd.cols$cols)

dev.off()

###################################################
### Plot (summary)
###################################################

fixNodes(tree) -> tree

(Ntip(tree)+1):(Ntip(tree)+Nnode(tree)) -> nodes


data.frame(node=nodes, qc=qc$node.label, qd=qd$node.label, qi=qi$node.label) -> sum
sum -> sum.c
sum.c$qc[] <- "Strong support"
sum.c$qc[sum$qc <= 0.2] <- "Weak support"
sum.c$qc[sum$qc <= -0.2] <- "Counter suppport"

sum.c$qd[] <- "Equal support"
sum.c$qd[sum$qd <= 0.3] <- "Alternative topology"

sum.c$qi[] <- "Informative"
sum.c$qi[sum$qi <= 0.3] <- "Lack information"

paste(sum.c$qc, sum.c$qd, sum.c$qi, sep="|") -> sum.c$summary
paste(sum.c$qc, sum.c$qi, sep="|") -> sum.c$summary.2
head(sum.c)
table(sum.c$summary)
table(sum.c$summary.2)

#which(sum.c$summary.2 == "Alternative topology|Informative") -> hard.incongruence
which(sum.c$summary.2 == "Counter suppport|Informative") -> hard.incongruence
grep("Weak", sum.c$summary) -> x
grep("Alternative", sum.c$summary) -> y
grep("Informative", sum.c$summary) -> z
intersect(intersect(x,y),z) -> hard.2
c(hard.incongruence,hard.2) -> hard.incongruence

sum.c$node[hard.incongruence] -> hard.incongruence

plot(tree, cex=0.5, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=qc.cols$support, cex=0.5)
nodelabels(text="*", node=hard.incongruence, frame="none", cex=2, adj=1.5)
title("Quartet Concordance")

###################################################
### Tribes summary
###################################################

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

data.frame(type=types, sum, sum.c) -> summary.out

write.csv(summary.out, "qs_summary_with_node_info.csv", row.names=F)

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
#dtree2$edge.length[] <- 1
#chronopl(dtree2, 1) -> dtree2

phylo.toBackbone(dtree2, trans) -> btree
phylo.toBackbone(dtree, trans.2) -> btree.2
plot(btree)
plot(btree.2)

### fix node labels

as.numeric(unlist(lapply(strsplit(btree$node.label, "/"), "[", 1))) -> btree$qc
as.numeric(unlist(lapply(strsplit(btree$node.label, "/"), "[", 2))) -> btree$qd
as.numeric(unlist(lapply(strsplit(btree$node.label, "/"), "[", 3))) -> btree$qi

### plot


cols <- c("red", "lightcoral", "skyblue", "darkblue")
cols2 <- c("red", "skyblue", "darkblue")
qc.breaks <- c(-1,-0.2,0,0.2,1)
qd.breaks <- qi.breaks <- c(0,0.3,0.7,1)

match(trans$clade.label, summary.out$type) -> x
summary.out$qc[x] -> qc.tips

supportColors(qc.tips, breaks = qc.breaks, colors = cols) -> qc.tips.cols
supportColors(btree$qc, breaks = qc.breaks, colors = cols) -> qc.cols
supportColors(btree$qi, breaks = qi.breaks, colors = cols2) -> qi.cols
supportColors(btree$qd, breaks = qd.breaks, colors = cols2) -> qd.cols

plot(btree, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=qc.cols$support, cex=0.5)
title("Quartet Concordance")
legend("bottomleft", legend=qc.cols$cuts, fill=qc.cols$cols)

plot(btree, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qi.cols$support)), frame="circ", bg=qi.cols$support, cex=0.5)
title("Quartet Informativeness")
legend("bottomleft", legend=qi.cols$cuts, fill=qi.cols$cols)

plot(btree, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qd.cols$support)), frame="circ", bg=qd.cols$support, cex=0.5)
title("Quartet Diferential")
legend("bottomleft", legend=qd.cols$cuts, fill=qd.cols$cols)


## final

which(btree$qi < 0.2) -> lack
(Ntip(btree)+1):(Ntip(btree)+Nnode(btree)) -> nodes.b
nodes.b[lack] -> lack

cairo_pdf("Quartet_concordance.pdf", width=6, height=6, onefile=T)

plot(btree, cex=0.8, use.edge.length=F, col="gray90", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=qc.cols$support, cex=0.7)
tiplabels(text=rep("", length(qc.tips.cols$support)), frame="circ", bg=qc.tips.cols$support, cex=0.7, offset=-0.4)
#nodelabels(text="*", node=lack, frame="none", cex=2)
nodelabels(btree$node.label, frame="none", adj=-0.2, cex=0.5)

plot(btree.2, cex=0.8, use.edge.length=F, col="gray90", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=qc.cols$support, cex=0.7)
tiplabels(text=rep("", length(qc.tips.cols$support)), frame="circ", bg=qc.tips.cols$support, cex=0.7, offset=-0.4)
#nodelabels(text="*", node=lack, frame="none", cex=2)
nodelabels(btree$node.label, frame="none", adj=-0.2, cex=0.5)

title("Quartet Concordance")
plot(1, axes=F)
#legend("center", legend=c(qc.cols$cuts, "* lack information"), fill=c(qc.cols$cols, NA))
legend("center", legend=qc.cols$cuts, fill=qc.cols$cols)

dev.off()

cairo_pdf("Quartet_concordance_narrow.pdf", width=4.5, height=6, onefile=T)

plot(btree, cex=0.8, use.edge.length=F, col="gray90", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=qc.cols$support, cex=0.7)
#tiplabels(text=rep("", length(qc.tips.cols$support)), frame="circ", bg=qc.tips.cols$support, cex=0.7, adj=2)
#nodelabels(text="*", node=lack, frame="none", cex=2)
#nodelabels(btree$node.label, frame="none", adj=-0.2, cex=0.5)
title("Quartet Concordance")

plot(btree.2, cex=0.8, use.edge.length=F, col="gray90", edge.width=2, show.tip.label=F)
nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=qc.cols$support, cex=0.7)
#tiplabels(text=rep("", length(qc.tips.cols$support)), frame="circ", bg=qc.tips.cols$support, cex=0.7, adj=2)
#nodelabels(text="*", node=lack, frame="none", cex=2)
#nodelabels(btree$node.label, frame="none", adj=-0.2, cex=0.5)
title("Quartet Concordance")

plot(1, axes=F)
legend("center", legend=c(qc.cols$cuts, "* lack information"), fill=c(qc.cols$cols, NA))

dev.off()
