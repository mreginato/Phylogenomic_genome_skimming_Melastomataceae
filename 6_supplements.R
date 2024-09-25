library(ape)
library(phyloch)
library(phytools)
library(xlsx)

##########################################
###  Dirs
##########################################


wd = "C:/SkimmingLoci_reduced/"
wd.plot = "C:/SkimmingLoci_reduced/7_plots"
wd.trees = "C:/SkimmingLoci_reduced/6_trees"

setwd(wd)

read.csv("Metadata_classification_tips.csv") -> meta
meta[-which(duplicated(meta$terminal)),] -> meta
rownames(meta) <- meta$terminal

read.xlsx2("Skimming_samples.xlsx", 1) -> vouchers

setwd(wd)

### format meta

meta[,c(4,3,1)] -> meta
colnames(meta) <- c("Tribe", "Genus", "Terminal")

strsplit(meta$Terminal, "-") -> x
sub("_", " ", unlist(lapply(x, "[", 1))) -> meta$Species
sub("_", " ", unlist(lapply(x, "[", 2))) -> meta$Code
meta[,-3] -> meta
meta$Source <- ""
meta$Voucher <- ""

match(vouchers$NextGen_Cod, meta$Code) -> x
vouchers[which(is.na(x)==F),] -> z
vouchers$Species[which(is.na(x))]
match(z$NextGen_Code, meta$Code) -> x
meta$Voucher[x] <- z$Voucher

##########################################
### Trees
##########################################

setwd(wd.trees)

list.files(pattern=".tre") -> files
sapply(files, readNexus, simplify = F) -> trees

names(trees) <- c("Concatenate", "Mitochondrial", "Plastidial", "Ribosomal", "Astral")

unlist(lapply(trees, FUN=function(x)(x$tip.label))) -> all.tips

which(is.na(match(all.tips, rownames(meta))))
which(is.na(match(rownames(meta), all.tips)))

### export singe tree files

trees -> trees.s
class(trees.s) <- "multiPhylo"

writeNexus(trees.s, file="Trees-all.tre")
write.tree(trees.s, file="Trees-all.tre")

### vouchers

matrix(ncol=5, nrow=nrow(meta)) -> mat
colnames(mat) <- names(trees)

for (i in 1:ncol(mat)) {
  trees[[i]] -> t0
  match(t0$tip.label, rownames(meta)) -> x
  mat[x,i] <- "x"
}

data.frame(meta, mat) -> meta

setwd(wd)

write.xlsx2(meta, file="Metadata_classification_tips-vouchers.xlsx", row.names = F)

##########################################
### Metrics
##########################################

setwd(wd.plot)

list.files(pattern="summary_coal_sims-") -> files
files[grep("*.tre", files)] -> files

sapply(files, read.tree) -> trees
class(trees) <- "multiPhylo"
sub(".tre", "", unlist(lapply(strsplit(files, "-"), "[", 2))) -> x
names(trees) <- x
writeNexus(trees, "All-metrics.tre")
