library(ape)
library(phyloch)
library(phangorn)
library(phytools)
library(phyloU)
library(scales)
library(xlsx)
source("F:/R Scripts/_Functions/supportColors.R")


##########################################
###  Dirs
##########################################


wd = "C:/SkimmingLoci_reduced/7_plots"
trees.dir = "C:/SkimmingLoci_reduced/6_trees"

### meta

setwd("C:/SkimmingLoci_reduced/")

read.csv("Metadata_classification_tips.csv") -> meta
meta[-which(duplicated(meta$terminal)),] -> meta
rownames(meta) <- meta$terminal

setwd(wd)


##########################################
###  Import trees
##########################################

setwd(trees.dir)

readNexus("concatenate_iqtree_rooted.tre", format="raxml") -> conc
readNexus("speciestree_geneboot_rooted.tre", format="raxml") -> astral
readNexus("RAxML_plastome_rooted.tre", format="raxml") -> plastome
readNexus("RAxML_mitochondrial_rooted.tre", format="raxml") -> mito
readNexus("RAxML_ribosomal_rooted.tre", format="raxml") -> ribo

setwd(wd)

### list all

list(Mitochondrial=mito, Plastome=plastome, Ribosomal=ribo, Concatenate=conc, Astral=astral) -> trees

astral$tip.label -> tips

for (i in 1:length(trees)) {
  trees[[i]] -> t0
  which(is.na(match(t0$tip.label, tips))) -> drop
  drop.tip(t0, drop) -> t0
  ladderize(t0) -> t0
  fixNodes(t0) -> t0
  round(as.numeric(t0$node.label),3) -> t0$node.label
  t0 -> trees[[i]]
}

unlist(lapply(trees, Ntip))

lapply(trees, ladderize) -> trees

##########################################
###  Plot all
##########################################

cuts <- c("0-50", "51-80", "81-95", "96-100")

try(dev.off())

plot(trees$Astral, edge.width = 1.5, edge.color = "gray80", cex=0.6, use.edge.length=F)

### Pdf

cairo_pdf("Trees_support_all.pdf", height = 12.5, width=8, onefile = T)

for (i in 1:length(trees)) {
  trees[[i]] -> t1
  supportColors(as.numeric(t1$node.label), breaks = c(0, 50, 80, 95, 100)) -> cols
  (Ntip(t1)++1):(Ntip(t1)++Nnode(t1)) -> nodes
  cols$support -> c0
  #c0[which(is.na(c0)==F)] -> c0
  #nodes[which(is.na(c0)==F)] -> nodes
  plot(t1, edge.width = 2, edge.color = "gray80", cex=0.4)
  nodelabels(text=rep("", length(c0))[-1], node=nodes[-1], frame="circ", bg=c0[-1], cex=0.6)
  title(names(trees)[i])
  legend("bottomleft", legend=cuts, title="Support", pch=21, pt.bg=cols$cols, box.col=NA, bg=NA, pt.cex = 2)
}

plot(1, xlab="", ylab="", col="white", axes=F)

dev.off()

system("open Trees_support_all.pdf")


##########################################
###  Tribes - backbone
##########################################

backs <- vector("list", length=length(trees))
names(backs) <- names(trees)
scales.length <- backs
tribes.boot <- backs
dtrees <- backs

for (i in 1:length(trees)) {
  trees[[i]] -> tree
  
  ## backbone phylo
  checkPhylo(tree, meta) -> c0
  c0$data -> meta0
  which(duplicated(meta0$tribe)) -> duplis
  drop.tip(tree, duplis) -> dtree
  ladderize(dtree) -> dtree
  fixNodes(dtree) -> dtree
  checkPhylo(dtree, meta0)$data -> dclades
  table(meta0$tribe) -> clades.n
  dclades$n[match(names(clades.n), dclades$tribe)] <- clades.n
  dtree$tip.label == rownames(dclades)
  plot(dtree) -> x
  x$x.lim[2]*0.02 -> size0
  rep(size0, nrow(dclades)) -> depths
  trans <- data.frame(tip.label=rownames(dclades), clade.label=dclades$tribe,
                      N=dclades$n, depth=depths)
  phylo.toBackbone(dtree, trans) -> btree
  btree -> backs[[i]]
  round(size0,3) -> scales.length[[i]]
  dtree -> dtree2
  dtree2$tip.label <- trans$clade.label
  dtree2 -> dtrees[[i]]
  ### tribes support
  meta0$tribe[match(dtree$tip.label, meta0$terminal)] -> tribes0
  support0 <- vector()
  for (k in 1:length(tribes0)) {
    tribes0[k] -> t0
    meta0$terminal[which(meta0$tribe == t0)] -> tips0
    if (length(tips0) > 1) {
      getMRCA(tree, tips0) -> node0
      tree$node.label[node0-Ntip(tree)] -> s0
    } else {
      s0 <- NA
    }
    c(support0,s0) -> support0
  }
  support0 -> tribes.boot[[i]]
  
}


##########################################
###  Tribes - Plot
##########################################

### include nodelabels all

cuts <- c("0-50", "51-80", "81-95", "96-100")

cairo_pdf("Trees_support_all-tribes.pdf", width=5.5, height=7, onefile=T)
for (i in 1:length(trees)) {
  backs[[i]] -> btree
  trees[[i]] -> tree
  tribes.boot[[i]] -> tribes.boot0
  scales.length[[i]] -> scale0
  c((Ntip(btree)+1):(Ntip(btree)+Nnode(btree))) -> nodes
  which(is.na(tribes.boot0)==F) -> keep.tips
  
  supportColors(as.numeric(btree$node.label), breaks = c(0, 50, 80, 95, 100)) -> cols
  cols$support -> c0
  supportColors(tribes.boot0, breaks = c(0, 50, 80, 95, 100))$support -> t0
  
  plot(btree, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=1.7, show.tip.label=F)
  nodelabels(text=rep("", length(c0))[-1], node=nodes[-1], frame="circ", bg=c0[-1], cex=0.6)
  tiplabels(tip=keep.tips, pch=21, bg=t0[keep.tips], cex=1, offset=-0.9*scale0)
  add.scale.bar(length=scale0, lwd=2, cex=0.7)
  title(names(trees)[i])
}
plot(1, xlab="", ylab="", axes=F, col="white")
legend("center", legend=cuts, title="Support", pch=21, pt.bg=cols$cols, box.col=NA, bg=NA, pt.cex = 2)  
dev.off()

system("open Trees_support_all-tribes.pdf")


cairo_pdf("Trees_support_all-tribes_wider.pdf", width=6.5, height=7, onefile=T)
for (i in 1:length(trees)) {
  backs[[i]] -> btree
  trees[[i]] -> tree
  tribes.boot[[i]] -> tribes.boot0
  scales.length[[i]] -> scale0
  c((Ntip(btree)+1):(Ntip(btree)+Nnode(btree))) -> nodes
  which(is.na(tribes.boot0)==F) -> keep.tips
  
  supportColors(as.numeric(btree$node.label), breaks = c(0, 50, 80, 95, 100)) -> cols
  cols$support -> c0
  supportColors(tribes.boot0, breaks = c(0, 50, 80, 95, 100))$support -> t0
  
  plot(btree, cex=1.1, use.edge.length=F, edge.color="gray", edge.width=1.7, show.tip.label=F)
  nodelabels(text=rep("", length(c0))[-1], node=nodes[-1], frame="circ", bg=c0[-1], cex=0.6)
  tiplabels(tip=keep.tips, pch=21, bg=t0[keep.tips], cex=1, offset=-0.9*scale0)
  add.scale.bar(length=scale0, lwd=2, cex=0.7)
  title(names(trees)[i])
}
plot(1, xlab="", ylab="", axes=F, col="white")
legend("center", legend=cuts, title="Support", pch=21, pt.bg=cols$cols, box.col=NA, bg=NA, pt.cex = 2)  
dev.off()

system("open Trees_support_all-tribes_wider.pdf")


##########################################
###  Density plot
##########################################

table(unlist(lapply(dtrees, FUN=function(x)(x$tip.label)))) -> x
names(which(x != max(x))) -> drop
lapply(dtrees, drop.tip, drop) -> dtrees
for (i in 1:length(dtrees)) {
  dtrees[[i]] -> t0
  t0$edge.length[] <- 1
  chronopl(t0, 1) -> t0
  t0 -> dtrees[[i]]
}
rev(dtrees) -> dtrees
class(dtrees) <- "multiPhylo"
dtrees

densityTree(dtrees, use.edge.length = F, compute.consensus = F, fix.depth = T)

write.nexus(dtrees, file="density_trees.tre")
write.nexus(dtrees[[1]], file="density_trees_major.tre")


##########################################
###  Bootstrap support
##########################################

names(trees)
names(tribes.boot)
names(backs)

matrix(ncol=length(trees), nrow=3) -> boot.t
rownames(boot.t) <- c("All nodes", "Backbone", "Tribes")
colnames(boot.t) <- names(trees)

lapply(trees, FUN=function(x)(x$node.label)) -> all.boots
round(unlist(lapply(all.boots, mean, na.rm=T)),1) -> x
lapply(all.boots, range, na.rm=T) -> y
paste0("[", unlist(lapply(y, paste, collapse="-")), "]") -> y
paste(x,y) -> boot.t[1,]


lapply(backs, FUN=function(x)(x$node.label)) -> back.boots
round(unlist(lapply(back.boots, mean, na.rm=T)),1) -> x
lapply(back.boots, range, na.rm=T) -> y
paste0("[", unlist(lapply(y, paste, collapse="-")), "]") -> y
paste(x,y) -> boot.t[2,]

round(unlist(lapply(tribes.boot, mean, na.rm=T)),1) -> x
lapply(tribes.boot, range, na.rm=T) -> y
paste0("[", unlist(lapply(y, paste, collapse="-")), "]") -> y
paste(x,y) -> boot.t[3,]


write.xlsx2(boot.t, "bootstrap_summary_trees.xlsx")


