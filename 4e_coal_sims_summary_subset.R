library(ape)
library(phyloch)
library(phangorn)
library(phytools)
library(relaimpo)
library(phyloU)
library(scales)
library(lattice)
source("F:/R Scripts/_Functions/supportColors.R")


##########################################
###  Dirs
##########################################


wd = "C:/SkimmingLoci_reduced/4_coal_sims"
qc.dir = "C:/SkimmingLoci_reduced/3_quartet_concordance"
boot.dir = "C:/SkimmingLoci_reduced/2_astral"
wd.plot = "C:/SkimmingLoci_reduced/7_plots"

### meta

setwd("C:/SkimmingLoci_reduced/")

read.csv("Metadata_classification_tips.csv") -> meta

setwd(wd)


##########################################
###  Import trees
##########################################

setwd(qc.dir)

### qc

read.tree("speciestree_qc_annnotated.tre") -> qc -> qd -> qi
as.numeric(unlist(lapply(strsplit(qc$node.label, "/"), "[", 1))) -> qc$node.label
as.numeric(unlist(lapply(strsplit(qd$node.label, "/"), "[", 2))) -> qd$node.label
as.numeric(unlist(lapply(strsplit(qi$node.label, "/"), "[", 3))) -> qi$node.label
ladderize(qc) -> qc
fixNodes(qc) -> qc
ladderize(qd) -> qd
fixNodes(qd) -> qd
ladderize(qi) -> qi
fixNodes(qi) -> qi

qd$node.label[is.na(qd$node.label)] <- 0
qc$tip.label -> tips

### phyparts

read.tree("phyparts_concon_labels.tre") -> phy
strsplit(phy$node.label, "|", fixed=T) -> x
do.call(rbind, x) -> x
phy -> phy.concon -> phy.confli -> phy.unin -> phy.ratio
phy.concon$node.label <- round(as.numeric(x[,2]),2)
phy.confli$node.label <- round(as.numeric(x[,3]),2)+round(as.numeric(x[,4]),2)
phy.unin$node.label <- round(as.numeric(x[,5]),2)

round(phy.concon$node.label/100,2) -> phy.concon$node.label
round(phy.confli$node.label/100,2) -> phy.confli$node.label
round(phy.unin$node.label/100,2) -> phy.unin$node.label

### make it deviation from 0.5
phy.concon$node.label - phy.confli$node.label -> x
abs(x) -> x
(x-1)*-1 -> x
x -> phy.ratio$node.label

setwd(boot.dir)

readNexus("speciestree_geneboot_rooted.tre", "raxml") -> boot

setwd(wd)

### error (invert -> the higher the value, the higher the error)
readNexus("RAxML_bipartitions.ERR_rooted.tre", "raxml") -> error
as.numeric(error$node.label) -> error$node.label
error$node.label -> x
(x-100)*-1 -> error$node.label
error$node.label[1:2] <- 0

### theta
read.tree("theta_tree_nodes.tre") -> theta

### introgression
read.tree("unbalancedTriplet.perc-MR.tre") -> intro

### rate

read.tree("Rates_concordance_median_rate.tre") -> rate
as.numeric(rate$node.label) -> rate$node.label
range(rate$node.label)

### tree plot (brlen)

read.tree("concatenate.phy.raxml.bestTree_rooted.tre") -> brtree

### list all

list(brtree,boot,qc,qd,qi,phy.ratio, phy.unin, error, theta, intro, rate) -> trees
length(trees)
names(trees) <- c("plot", "bootstraps", "qc", "qd", "qi", "conc_confli_ratio", "uninformative", "error", "theta", "reticulation", "rate")

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

trees[[1]] -> tree
trees[-1] -> trees

setwd(wd.plot)

##########################################
###  Tribes
##########################################

c((Ntip(tree)++1):(Nnode(tree)+Ntip(tree))) -> nodes

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

plot(tree)
plot(tree, use.edge.length=F)
nodelabels(node=backbone)
nodes -> types
types[which(is.na(match(nodes, internal.nodes))==F)] <- "Internal nodes"
types[which(is.na(match(nodes, backbone))==F)] <- "Backbone"
types[match(nodes.dat$node[-17], nodes)] <- nodes.dat$tribe[-17]

data.frame(matrix(nrow=length(nodes), ncol=length(trees)+2)) -> summary.out
colnames(summary.out) <- c("node", "type", names(trees))
summary.out$node <- nodes
summary.out$type <- types

for (i in 1:length(trees)) {
  trees[[i]] -> t0
  t0$node.label -> x
  summary.out[,i+2] <- x
}

summary.out$error[1] <- 0

### make 0 - 1

apply(summary.out[,-c(1:2)], MARGIN=2, FUN=range, na.rm=T)

summary.out$bootstraps/100 -> summary.out$bootstraps
round(summary.out$qc, 2) -> summary.out$qc
summary.out$error/100 -> summary.out$error

summary.out$rate -1 -> x
x[x == 0] <- NA
x/max(x, na.rm = T) -> x
range(x, na.rm = T) 
x -> summary.out$rate

head(summary.out)
colnames(summary.out)
colnames(summary.out) <- c("Node", "Type", "Gene Bootstrap", "QC", "QD", "QI", "ILS (gene trees)", 
                          "Uninformative genes", "Inference error", "ILS (theta)",  "RI", "Rate ratio")

summary.out[,-c(5:7)] -> summary.out
write.csv(summary.out, "summary_coal_sims-all_nodes.csv", row.names=F)


tree -> node.tree
c((Ntip(tree)+1):(Nnode(tree)+Ntip(tree))) -> node.tree$node.label
write.tree(node.tree, "summary_coal_sims-all_nodes_key.tre")

## theta outliers

summary.out$`ILS (theta)`-> x
boxplot(x)
x[x > 0.6] <- 0.6
x -> summary.out$`ILS (theta)`

### scale 0 - 1

summary.out -> summary.out.raw

apply(summary.out[,-c(1:2)], MARGIN=2, range, na.rm=T)

summary.out$`Uninformative genes` -> x
round(x/max(x),2) -> summary.out$`Uninformative genes`

summary.out$`Inference error` -> x
round(x/max(x, na.rm = T), 2) -> summary.out$`Inference error`

summary.out$RI -> x
x/max(x, na.rm = T) -> summary.out$RI

summary.out$`ILS (theta)`-> x
x/max(x) -> summary.out$`ILS (theta)`

summary.out$`Rate ratio` -> x
(x+1)/2 -> summary.out$`Rate ratio`


write.csv(summary.out, "summary_coal_sims-all_nodes-scaled.csv", row.names=F)

## export all trees

trees[c(1,2,6:10)] -> trees

for (i in 1:length(trees)) {
  write.tree(trees[[i]], file=paste0("summary_coal_sims-", names(trees)[i], ".tre"))
}

##########################################
###  Check response corr
##########################################

### subset variables

colnames(summary.out)

cols.res = c(3,4)
cols.pred = c(5:ncol(summary.out)) ## boot
colnames(summary.out)[cols.res]
colnames(summary.out)[cols.pred]
dat.res = summary.out[,cols.res]
dat.pre = summary.out[,cols.pred]
colnames(dat.res)
colnames(dat.pre)

### correlogram (sreponse)

res.mat <- matrix(nrow=ncol(dat.res), ncol=ncol(dat.res))
colnames(res.mat) <- rownames(res.mat) <- colnames(dat.res)

for (i in 1:ncol(dat.res)) {
  dat.res[,i] -> d1
  for (k in 1:ncol(dat.res)) {
    dat.res[,k] -> d2
    cor.test(d1,d2) -> cor0
    cor0$estimate -> res.mat[i,k]
  }
}

abs(res.mat) -> res.mat
diag(res.mat) <- NA
levelplot(res.mat)
res.mat -> res.mat.r

### correlogram (predictors)

res.mat <- matrix(nrow=ncol(dat.pre), ncol=ncol(dat.pre))
colnames(res.mat) <- rownames(res.mat) <- colnames(dat.pre)

for (i in 1:ncol(dat.pre)) {
  dat.pre[,i] -> d1
  for (k in 1:ncol(dat.pre)) {
    dat.pre[,k] -> d2
    cor.test(d1,d2) -> cor0
    cor0$estimate -> res.mat[i,k]
  }
}

abs(res.mat) -> res.mat
diag(res.mat) <- NA
res.mat -> res.mat.p

### correlogram (all)

cbind(dat.res,dat.pre) -> dat.all

res.mat <- matrix(nrow=ncol(dat.all), ncol=ncol(dat.all))
colnames(res.mat) <- rownames(res.mat) <- colnames(dat.all)
res.mat.pv <- res.mat

for (i in 1:ncol(dat.all)) {
  dat.all[,i] -> d1
  for (k in 1:ncol(dat.all)) {
    dat.all[,k] -> d2
    cor.test(d1,d2) -> cor0
    cor0$estimate -> res.mat[i,k]
    cor0$p.value -> res.mat.pv[i,k]
  }
}

abs(res.mat) -> res.mat
diag(res.mat.pv) <- diag(res.mat) <- NA
res.mat -> res.mat.all

write.csv(res.mat.all, "correlogram_pearson.csv")
write.csv(res.mat.pv, "correlogram_pearson_pvalues.csv")


### levelplot


levelplot(res.mat.r,border="black", col.regions = c("white", rev(heat.colors(20))), cuts=10, pretty=F, colorkey=T, xlab="", ylab="")
levelplot(res.mat.p,border="black", col.regions = c("white", rev(heat.colors(20))), cuts=10, pretty=F, colorkey=T, xlab="", ylab="")
levelplot(res.mat.all,border="black", col.regions = c("white", rev(heat.colors(20))), cuts=10, pretty=F, colorkey=T, xlab="", ylab="")


cairo_pdf(filename="summary_coal_sims-correlograms.pdf", width=6, height=6, onefile = T)

levelplot(res.mat.r,border="black", col.regions = c("white", rev(heat.colors(20))), cuts=10, pretty=F, colorkey=T, xlab="", ylab="")
levelplot(res.mat.p,border="black", col.regions = c("white", rev(heat.colors(20))), cuts=10, pretty=F, colorkey=T, xlab="", ylab="")
levelplot(res.mat.all,border="black", col.regions = c("white", rev(heat.colors(20))), cuts=10, pretty=F, colorkey=T, xlab="", ylab="")

dev.off()


##########################################
###  Var imp
##########################################

colnames(dat.res)

data.frame(Boot=dat.res[,1], dat.pre) -> dat
colnames(dat)

dat[which(summary.out$Type == "Backbone"),] -> dat.back
dat[which(summary.out$Type == "Internal nodes"),] -> dat.int
dat[intersect(which(summary.out$Type != "Backbone"), which(summary.out$Type != "Internal nodes")),] -> dat.tribes

calc.relimp(dat, type = "lmg",rela=F) -> imp.all
imp.all
calc.relimp(dat.back, type = "lmg",rela=F) -> imp.back
imp.back
calc.relimp(dat.tribes, type = "lmg",rela=F) -> imp.tribes
imp.tribes
calc.relimp(dat.int, type = "lmg",rela=F) -> imp.int
imp.int

bootimpo.all <- boot.relimp(dat, b = 100, type = c("lmg"), rank = TRUE, diff = TRUE, rela = F)
booteval.relimp(bootimpo.all,lev=0.9,nodiff=TRUE) -> boot.all

## Plot

### node key
col.nodes <- c("red4", "blue4", "yellow3", "gray30")
alpha(col.nodes, 0.7) -> col.nodes

summary.out$Type -> node.types
node.types[which(is.na(match(node.types, tribes))==F)] <- "Tribes (MRCA)"
node.types[node.types == "Internal nodes"] <- "Tribes (internal nodes)"
as.factor(node.types) -> node.types
cex=c(1.5,1,1.5)

summary.out$Node[match(tribes, summary.out$Type)] -> nodes.tribes
names(nodes.tribes) <- tribes
grep("Rupestr", tree$tip.label) -> nodes.tribes[17]
nodes.tribes

col.nodes[node.types] -> col.nodes.plot
cex[node.types] -> cex.plot

tree -> ptree
ptree$tip.label <- rep("aaa", Ntip(tree))
chronopl(ptree, 1) -> ptree

names(nodes.tribes) -> x
max(nchar(x)) -> m
for (i in 1:length(x)) {
  x[i] -> x0
  (m-nchar(x0))+1 -> r0
  paste0(x0, paste(rep(" ", r0), collapse="")) -> x[i]
}
names(nodes.tribes) <- x

#plot(tree, show.tip.label=F, use.edge.length=T, edge.width=2, edge.color="gray90", x.lim=0.9)
plot(tree, show.tip.label=F, use.edge.length=F, edge.width=2, edge.color="gray90", x.lim=2000)
for (i in 1:length(nodes.tribes)) {
  cladelabels(text=names(nodes.tribes)[i], node=nodes.tribes[i], orientation = "horizontal", offset=0.1)
}

nodelabels(pch=21, bg=col.nodes.plot, cex=cex.plot)
legend("bottomleft", legend=levels(node.types), pch=21, pt.bg = col.nodes, box.col=NA, bg=NA, pt.cex = cex)

chronos(tree, 1) -> tree.p
fixNodes(tree.p) -> tree.p

setwd(wd.plot)

pdf("Var_imp_tree_nodes.pdf", width=10, height=10)
plot(tree.p, show.tip.label=F, use.edge.length=T, edge.width=2, edge.color="gray90", x.lim=3)
for (i in 1:length(nodes.tribes)) {
  cladelabels(text=names(nodes.tribes)[i], node=nodes.tribes[i], orientation = "horizontal", offset=0.1, cex=0.5, wing.length = 0)
}
nodelabels(pch=21, bg=col.nodes.plot, cex=cex.plot)
plot(tree, show.tip.label=F, use.edge.length=F, edge.width=2, edge.color="gray90", x.lim=2000)
for (i in 1:length(nodes.tribes)) {
  cladelabels(text=names(nodes.tribes)[i], node=nodes.tribes[i], orientation = "horizontal", offset=0.1, cex=0.5, wing.length = 0)
}
nodelabels(pch=21, bg=col.nodes.plot, cex=cex.plot)

plot(tree, show.tip.label=F, use.edge.length=F, edge.width=2, edge.color="gray90", x.lim=2000, direction="left")
for (i in 1:length(nodes.tribes)) {
  cladelabels(text=names(nodes.tribes)[i], node=nodes.tribes[i], orientation = "horizontal", offset=100, cex=0.5, wing.length = 0)
}
nodelabels(pch=21, bg=col.nodes.plot, cex=cex.plot)



plot(1, xlab="", ylab="", axes=F, col="white")
legend("center", legend=levels(node.types), pch=21, pt.bg = col.nodes, box.col=NA, bg=NA, pt.cex = cex)
dev.off()

pdf("Var_imp_tree_nodes_fan.pdf", width=8, height=6)
plot(tree, show.tip.label=F, use.edge.length=F, edge.width=2, edge.color="gray90", type="fan")
for (i in 1:length(nodes.tribes)) {
  arc.cladelabels(text=names(nodes.tribes)[i], node=nodes.tribes[i], orientation = "horizontal", stretch = 1.5, ln.offset=1.03, lab.offse=1.04, cex=0.6, wing.length = 0.5, mark.node=F)
}
nodelabels(pch=21, bg=col.nodes.plot, cex=cex.plot*.8)
plot(1, xlab="", ylab="", axes=F, col="white")
legend("center", legend=levels(node.types), pch=21, pt.bg = col.nodes, box.col=NA, bg=NA, pt.cex = cex)
dev.off()

### barplot

range(imp.back@lmg)
ylim = c(0,0.6)

plot(boot.all)

### colocar boxplot of bootstrap do lado


pdf("Var_imp_barplots.pdf", width=4, height=10)

par(mfrow=c(4,1), mar=c(4,10,2,2))

barplot(imp.all@lmg, xlab="% of response variance", col=col.nodes[4], xlim = ylim, horiz=T, las=2)
#barplot(imp.all@lmg, ylab="% of response variance", col=col.nodes[4], ylim = ylim)
title("All nodes - Relative importance for QC")
round(imp.all@R2*100,1) -> x
legend("topright", legend=paste0("R2 = ", x, " %"), box.col = NA, cex=1.1)

barplot(imp.back@lmg, xlab="% of response variance", col=col.nodes[1], xlim = ylim, horiz=T, las=2)
#barplot(imp.back@lmg, ylab="% of response variance", col=col.nodes[1], ylim = ylim)
title("Backbone - Relative importance for QC")
round(imp.back@R2*100,1) -> x
legend("topright", legend=paste0("R2 = ", x, " %"), box.col = NA, cex=1.1)

barplot(imp.tribes@lmg, xlab="% of response variance", col=col.nodes[3], xlim = ylim, horiz=T, las=2)
#barplot(imp.tribes@lmg, ylab="% of response variance", col=col.nodes[2], ylim = ylim)
title("Tribes (MRCA) - Relative importance for QC")
round(imp.tribes@R2*100,1) -> x
legend("topright", legend=paste0("R2 = ", x, " %"), box.col = NA, cex=1.1)

barplot(imp.int@lmg, xlab="% of response variance", col=col.nodes[2], xlim = ylim, horiz=T, las=2)
#barplot(imp.int@lmg, ylab="% of response variance", col=col.nodes[3], ylim = ylim)
title("Tribes (internal nodes) - Relative importance for QC")
round(imp.int@R2*100,1) -> x
legend("topright", legend=paste0("R2 = ", x, " %"), box.col = NA, cex=1.1)

dev.off()

system("open Var_imp_barplots.pdf")

pdf("Var_imp_barplots_boxplots.pdf", width=2, height=10)
par(mfrow=c(4,1), mar=c(4,10,2,2))

boxplot(dat$Boot, col=col.nodes[4], ylim=c(0,1), axes=F)
axis(4)
boxplot(dat.back$Boot, col=col.nodes[1], ylim=c(0,1), axes=F)
axis(4)
boxplot(dat.tribes$Boot, col=col.nodes[3], ylim=c(0,1), axes=F)
axis(4)
boxplot(dat.int$Boot, col=col.nodes[2], ylim=c(0,1), axes=F)
axis(4)

dev.off()

system("open Var_imp_barplots_boxplots.pdf")


pdf("Var_imp_barplots_biplots_back.pdf", width=4, height=4)
#par(mfrow=c(1,1), mar=c(4,10,2,2))

plot(dat.back$Boot ~ dat.back$Uninformative.genes, axes=F, cex=1.5, pch=21, bg=col.nodes[1], ylab="Gene Bootstrap", xlab="Uninformative genes (%)")
axis(1); axis(2)
abline(lm(dat.back$Boot ~ dat.back$Uninformative.genes), lty="dashed", lwd=2, col="gray30")

plot(dat.back$Boot ~ dat.back$Inference.error, axes=F, cex=1.5, pch=21, bg=col.nodes[1], ylab="Gene Bootstrap", xlab="Inference error")
axis(1); axis(2)
abline(lm(dat.back$Boot ~ dat.back$Inference.error), lty="dashed", lwd=2, col="gray30")

plot(dat.back$Boot ~ dat.back$ILS..theta., axes=F, cex=1.5, pch=21, bg=col.nodes[1], ylab="Gene Bootstrap", xlab="ILS (theta)")
axis(1); axis(2)
abline(lm(dat.back$Boot ~ dat.back$ILS..theta.), lty="dashed", lwd=2, col="gray30")

plot(dat.back$Boot ~ dat.back$RI, axes=F, cex=1.5, pch=21, bg=col.nodes[1], ylab="Gene Bootstrap", xlab="RI (reticulation index)")
axis(1); axis(2)
abline(lm(dat.back$Boot ~ dat.back$RI), lty="dashed", lwd=2, col="gray30")

plot(dat.back$Boot ~ dat.back$Rate.ratio, axes=F, cex=1.5, pch=21, bg=col.nodes[1], ylab="Gene Bootstrap", xlab="Rate ratio")
axis(1); axis(2)
abline(lm(dat.back$Boot ~ dat.back$Rate.ratio), lty="dashed", lwd=2, col="gray30")

dev.off()

system("open Var_imp_barplots_biplots_back.pdf")


sink("Var_imp_out.txt")
cat("All", fill=T)
imp.all
cat("Backbone", fill=T)
imp.back
cat("Tribes (mrca)", fill=T)
imp.tribes
cat("Internal", fill=T)
imp.int
sink()



##########################################
###  Plot - Trees
##########################################


### include nodelabels all

summary.out$Node == summary.out.raw$Node

#apply(summary.out[,-c(1:2)], 1, paste, collapse="|") -> node.labs
apply(summary.out.raw[,-c(1:2)], 1, paste, collapse="|") -> node.labs
tree$node.label <- node.labs

## backbone phylo

which(duplicated(meta$tribe)) -> duplis
drop.tip(tree, duplis) -> dtree
ladderize(dtree) -> dtree
fixNodes(dtree) -> dtree
checkPhylo(dtree, meta)$data -> dclades
table(meta$tribe) -> clades.n
dclades$n[match(names(clades.n), dclades$tribe)] <- clades.n
dtree$tip.label == rownames(dclades)

rep(0.02, nrow(dclades)) -> depths
rep(0.2, nrow(dclades)) -> depths.2

trans <- data.frame(tip.label=rownames(dclades), clade.label=dclades$tribe,
                      N=dclades$n, depth=depths)

trans.2 <- data.frame(tip.label=rownames(dclades), clade.label=dclades$tribe,
                    N=rep(5,nrow(dclades)), depth=depths.2)



dtree -> dtree2
#dtree2$edge.length[] <- 1
chronopl(dtree2, 1) -> dtree2

phylo.toBackbone(dtree, trans) -> btree
phylo.toBackbone(dtree2, trans.2) -> btree.2
plot(btree)
plot(btree.2)

### fix node labels

colnames(summary.out)[-c(1:2)]

do.call(rbind, strsplit(btree$node.label, "|", fixed=T)) -> pie.mat
colnames(pie.mat) <- colnames(summary.out)[-c(1:2)]

mode(pie.mat) <- "numeric"

colnames(pie.mat)
#-1*(pie.mat[,9]-1) -> pie.mat[,9]

### plot

apply(pie.mat, 2, range)

cols1 <- c("red", "lightcoral", "skyblue", "darkblue")
cols2 <- rev(cols1) ## QD, ILS, unin genes, inf error, theta, RI, rate ratio
breaks1 <- c(0, 0.6, 0.9, 0.95, 1) ## boot
breaks2 <- c(-1,-0.2,0,0.2,1) ## QC
breaks3 <- c(0,0.3,0.6, 0.9,1) ## QD, QI, ILS, unin genes
breaks4 <- c(0,0.1,0.3,0.5,0.7) ## inf error
breaks5 <- c(0,0.05,0.11,0.15,0.2) ## theta
breaks6 <- c(-1,-0.5,0,0.5,1) ## rate ratio
breaks7 <- c(0,0.15,0.3,0.45,0.6) ## inf error



colnames(pie.mat)

all.breaks <- list(breaks1, breaks2, breaks3, breaks4, breaks5, breaks7, breaks6)
all.cols <- list(cols1, cols1, cols2, cols2, cols2, cols2, cols2)
ncol(pie.mat)
length(all.breaks)
length(all.cols)

match(trans$clade.label, summary.out$Type) -> tips.x
summary.out.raw -> summary.out

cairo_pdf("Var_imp_trees-all.pdf", width=4, height=6, onefile=T)
for (i in 1:ncol(pie.mat)) {
  as.numeric(pie.mat[,i]) -> d0
  colnames(pie.mat)[i] -> lab0
  summary.out[tips.x,lab0] -> qc.tips
  all.breaks[[i]] -> b0
  all.cols[[i]] -> cols
  supportColors(qc.tips, breaks = b0, colors = cols) -> qc.tips.cols
  supportColors(d0, breaks = b0, colors = cols) -> qc.cols
  which(is.na(qc.tips.cols$support)) -> miss.tips
  c(1:Ntip(btree))[-miss.tips] -> tips0
  
  
  plot(btree, cex=0.8, use.edge.length=F, edge.color="gray", edge.width=2, show.tip.label=F)
  nodelabels(text=rep("", length(qc.cols$support)), frame="circ", bg=qc.cols$support, cex=0.5)
  tiplabels(tip = tips0, text=rep("", length(tips0)), frame="circ", bg=qc.tips.cols$support[-miss.tips], cex=0.5, offset = -0.018)
  title(lab0)
  #legend("bottomleft", legend=qc.cols$cuts, pch=21, pt.bg=qc.cols$cols, pt.cex=1.5, box.col=NA, bg=NA, title=lab0, title.cex = 1)
  legend("bottomright", legend=qc.cols$cuts, pch=21, pt.bg=qc.cols$cols, pt.cex=1.2, box.col=NA, bg=NA)
  
}
dev.off()

system("open Var_imp_trees-all.pdf")


