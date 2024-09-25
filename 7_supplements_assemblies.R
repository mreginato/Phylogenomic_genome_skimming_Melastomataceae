library(ape)
library(phyloch)
library(skimmingLociR)
library(xlsx)

##########################################
###  Dirs
##########################################


wd = "C:/SkimmingLoci_reduced/"
wd.trees = "C:/SkimmingLoci_reduced/6_trees"
wd.aligns = "C:/SkimmingLoci_reduced/0_Plastome_ribo_mito"
wd.conc = "C:/SkimmingLoci_reduced/0_concatenate"

setwd(wd)

read.csv("Metadata_classification_tips.csv") -> meta
meta[-which(duplicated(meta$terminal)),] -> meta
rownames(meta) <- meta$terminal

setwd(wd)


##########################################
###  Aligns - stats
##########################################

### conc

setwd(wd.conc)

read.dna("concatenate.phy") -> conc
read.csv("concatenate_map.csv", row.names=1) -> conc.map

rownames(conc) -> x
which(is.na(match(x, rownames(meta))))

n <- vector()
stats <- vector("list", length=nrow(conc.map))

for (i in 1:nrow(conc.map)) {
  conc[,conc.map[i,1]:conc.map[i,2]] -> a0
  delete.empty.cells(a0) -> a0
  nrow(a0) -> n0
  alignStats(a0) -> s0
  c(N=n0,s0) -> stats[[i]]
  cat("\r",i)
}

do.call(rbind, stats) -> conc.stats
rownames(conc.map) -> rownames(conc.stats)
alignStats(conc) -> conc.s
c(N=nrow(conc), conc.s) -> conc.s
rbind(conc.stats, conc.s) -> conc.stats

setwd(wd.trees)

colnames(conc.stats)[1] <- "Terminals"
write.xlsx2(conc.stats, "Stats-aligns-nuclear_low_copy.xlsx", row.names=T)

### others

setwd(wd.aligns)

### mito

read.dna("mitochondrial.phy") -> align -> mito
read.csv("mitochondrial_map.csv", row.names=1) -> align.map -> mito.map
n <- vector()
stats <- vector("list", length=nrow(align.map))

for (i in 1:nrow(align.map)) {
  align[,align.map[i,1]:align.map[i,2]] -> a0
  delete.empty.cells(a0) -> a0
  nrow(a0) -> n0
  alignStats(a0) -> s0
  c(N=n0,s0) -> stats[[i]]
  cat("\r",i)
}

do.call(rbind, stats) ->stats
sub(".phy$", "", rownames(align.map)) -> rownames(stats)
alignStats(align) -> align.s
c(N=nrow(align), align.s) -> aligns.s
rbind(stats,aligns.s) -> stats
stats -> mito.stats

### ribo

read.dna("ribosomal.phy") -> align -> ribo
read.csv("ribosomal_map.csv", row.names=1) -> align.map -> ribo.map
n <- vector()
stats <- vector("list", length=nrow(align.map))

for (i in 1:nrow(align.map)) {
  align[,align.map[i,1]:align.map[i,2]] -> a0
  delete.empty.cells(a0) -> a0
  nrow(a0) -> n0
  alignStats(a0) -> s0
  c(N=n0,s0) -> stats[[i]]
  cat("\r",i)
}

do.call(rbind, stats) ->stats
sub(".phy$", "", rownames(align.map)) -> rownames(stats)
alignStats(align) -> align.s
c(N=nrow(align), align.s) -> aligns.s
rbind(stats,aligns.s) -> stats
stats -> ribo.stats

### plastome

read.dna("plastome.phy") -> align -> plastome
read.csv("plastome_map.csv", row.names=1) -> align.map -> plastome.map
n <- vector()
stats <- vector("list", length=nrow(align.map))

for (i in 1:nrow(align.map)) {
  align[,align.map[i,1]:align.map[i,2]] -> a0
  delete.empty.cells(a0) -> a0
  nrow(a0) -> n0
  alignStats(a0) -> s0
  c(N=n0,s0) -> stats[[i]]
  cat("\r",i)
}

do.call(rbind, stats) ->stats
sub(".phy$", "", rownames(align.map)) -> rownames(stats)
alignStats(align) -> align.s
c(N=nrow(align), align.s) -> aligns.s
rbind(stats,aligns.s) -> stats
stats -> plastome.stats

setwd(wd.trees)

colnames(mito.stats)[1] <- "Terminals"
colnames(ribo.stats)[1] <- "Terminals"
colnames(plastome.stats)[1] <- "Terminals"

write.xlsx2(mito.stats, "Stats-aligns-mitochondrial.xlsx", row.names=T)
write.xlsx2(ribo.stats, "Stats-aligns-ribosomal.xlsx", row.names=T)
write.xlsx2(plastome.stats, "Stats-aligns-plastome.xlsx", row.names=T)


### mean range

keep.cols <- c(1,3:6)
mito.stats[,keep.cols] -> mito0
mode(mito0) <- "numeric"
data.frame(Data="Mitochondrial", mito0) -> mito0
mito0 -> mito.s

ribo.stats[,keep.cols] -> mito0
mode(mito0) <- "numeric"
data.frame(Data="Ribosomal", mito0) -> mito0
mito0 -> ribo.s

plastome.stats[,keep.cols] -> mito0
mode(mito0) <- "numeric"
data.frame(Data="Plastidial", mito0) -> mito0
mito0 -> plastome.s

conc.stats[,keep.cols] -> mito0
mode(mito0) <- "numeric"
data.frame(Data="Nuclear", mito0) -> mito0
mito0 -> conc.s

rownames(conc.s) <- NULL
rownames(ribo.s) <- NULL
rownames(mito.s) <- NULL
rownames(plastome.s) <- NULL


total.terminals = c(mito.s[nrow(mito.s),2], conc.s[nrow(conc.s),2], plastome.s[nrow(plastome.s),2], ribo.s[nrow(ribo.s),2])

conc.s[-nrow(conc.s),] -> conc.s
plastome.s[-nrow(plastome.s),] -> plastome.s
mito.s[-nrow(mito.s),] -> mito.s
ribo.s[-nrow(ribo.s),] -> ribo.s

rbind(conc.s, ribo.s, mito.s, plastome.s) -> aligns.s

aggregate(aligns.s$Terminals, by=list(aligns.s$Data), FUN=length) -> total.loci

aggregate(aligns.s$Terminals, by=list(aligns.s$Data), FUN=median) -> med.terminals
aggregate(aligns.s$Terminals, by=list(aligns.s$Data), FUN=range) -> range.terminals
apply(range.terminals$x, 1, paste, collapse="-") -> range.terminals

aggregate(aligns.s$Aligned_bp, by=list(aligns.s$Data), FUN=median) -> med.bp
aggregate(aligns.s$Aligned_bp, by=list(aligns.s$Data), FUN=sum) -> sum.bp

aggregate(aligns.s$Variable, by=list(aligns.s$Data), FUN=median) -> med.var
aggregate(aligns.s$Variable, by=list(aligns.s$Data), FUN=sum) -> sum.var

aggregate(aligns.s$PIS, by=list(aligns.s$Data), FUN=median) -> med.pis
aggregate(aligns.s$PIS, by=list(aligns.s$Data), FUN=sum) -> sum.pis

aggregate(aligns.s$Missing.data, by=list(aligns.s$Data), FUN=median) -> med.miss
aggregate(aligns.s$Missing.data, by=list(aligns.s$Data), FUN=range) -> range.miss
apply(range.miss$x, 1, paste, collapse="-") -> range.miss


data.frame(Loci_n=total.loci$x, Terminals_total=total.terminals, Terminals_loci_median=med.terminals$x,
           Terminals_loci_range=range.terminals, Bp_loci_total=sum.bp$x,
           Bp_loci_median=med.bp$x, Variable_loci_total=sum.var$x,
           Variable_loci_median=med.var$x, PIS_loci_total=sum.pis$x,
           PIS_loci_median=med.pis$x, Missing_loci_med=med.miss$x, 
           Missing_loci_range=range.miss) -> a.stats
rownames(a.stats) <- med.var$Group.1

write.xlsx2(a.stats, "Stats-aligns-All.xlsx")



##########################################
###  Samples - stats
##########################################

rownames(meta) -> spp

samples <- vector("list", length=length(spp))
names(samples) <- spp

for (k in 1:length(spp)) {
  spp[k] -> sp0
  ## conc
  conc.map -> map
  conc -> align
  match(sp0, rownames(align)) -> r0
  if (is.na(r0) == F) {
    stats0 <- vector(length=nrow(map))
    stats0[] <- NA
    for (i in 1:nrow(map)) {
      align[r0,conc.map[i,1]:conc.map[i,2]] -> a0
      delete.empty.cells(a0, quiet = T) -> a0
      if (nrow(a0) > 0) {
        alignStats(a0) -> s0
        s0[2] -> s0
        as.numeric(s0) -> stats0[i]
      }
      cat("\r",i)
    }
    length(which(is.na(stats0)==F)) -> loci.n
    range(stats0, na.rm = T) -> range.n
    median(stats0, na.rm=T) -> median.n
    sum(stats0, na.rm = T) -> total.n
    c(loci_n=loci.n, total_bp=total.n, median_bp_loci=median.n, range_bp_loci=paste(range.n, collapse="-")) -> stats0
  } else {
    c(loci_n=0, total_bp=0, median_bp_loci=0, range_bp_loci=0) -> stats0
  }
  stats0 -> conc0
  
  ## plastome
  plastome.map -> map
  plastome -> align
  match(sp0, rownames(align)) -> r0
  if (is.na(r0) == F) {
    stats0 <- vector(length=nrow(map))
    stats0[] <- NA
    for (i in 1:nrow(map)) {
      align[r0,map[i,1]:map[i,2]] -> a0
      delete.empty.cells(a0, quiet = T) -> a0
      if (nrow(a0) > 0) {
        alignStats(a0) -> s0
        s0[2] -> s0
        as.numeric(s0) -> stats0[i]
      }
      cat("\r",i)
    }
    length(which(is.na(stats0)==F)) -> loci.n
    range(stats0, na.rm = T) -> range.n
    median(stats0, na.rm=T) -> median.n
    sum(stats0, na.rm = T) -> total.n
    c(loci_n=loci.n, total_bp=total.n, median_bp_loci=median.n, range_bp_loci=paste(range.n, collapse="-")) -> stats0
  } else {
    c(loci_n=0, total_bp=0, median_bp_loci=0, range_bp_loci=0) -> stats0
  }
  stats0 -> plastome0
  
  
  ## mito
  mito.map -> map
  mito -> align
  match(sp0, rownames(align)) -> r0
  if (is.na(r0) == F) {
    stats0 <- vector(length=nrow(map))
    stats0[] <- NA
    for (i in 1:nrow(map)) {
      align[r0,map[i,1]:map[i,2]] -> a0
      delete.empty.cells(a0, quiet = T) -> a0
      if (nrow(a0) > 0) {
        alignStats(a0) -> s0
        s0[2] -> s0
        as.numeric(s0) -> stats0[i]
      }
      cat("\r",i)
    }
    length(which(is.na(stats0)==F)) -> loci.n
    range(stats0, na.rm = T) -> range.n
    median(stats0, na.rm=T) -> median.n
    sum(stats0, na.rm = T) -> total.n
    c(loci_n=loci.n, total_bp=total.n, median_bp_loci=median.n, range_bp_loci=paste(range.n, collapse="-")) -> stats0
  } else {
    c(loci_n=0, total_bp=0, median_bp_loci=0, range_bp_loci=0) -> stats0
  }
  stats0 -> mito0
  
  ## ribo
  ribo.map -> map
  ribo -> align
  match(sp0, rownames(align)) -> r0
  if (is.na(r0) == F) {
    stats0 <- vector(length=nrow(map))
    stats0[] <- NA
    for (i in 1:nrow(map)) {
      align[r0,map[i,1]:map[i,2]] -> a0
      delete.empty.cells(a0, quiet = T) -> a0
      if (nrow(a0) > 0) {
        alignStats(a0) -> s0
        s0[2] -> s0
        as.numeric(s0) -> stats0[i]
      }
      cat("\r",i)
    }
    length(which(is.na(stats0)==F)) -> loci.n
    range(stats0, na.rm = T) -> range.n
    median(stats0, na.rm=T) -> median.n
    sum(stats0, na.rm = T) -> total.n
    c(loci_n=loci.n, total_bp=total.n, median_bp_loci=median.n, range_bp_loci=paste(range.n, collapse="-")) -> stats0
  } else {
    c(loci_n=0, total_bp=0, median_bp_loci=0, range_bp_loci=0) -> stats0
  }
  stats0 -> ribo0
  
  ###
  rbind(conc0, plastome0, mito0, ribo0) -> t0 
  rownames(t0) <- c("Nuclear_low_copy", "Plastidial", "Mitochondrial", "Ribosomal")
  t(t0) -> t0
  data.frame(Terminal=sp0, t0) -> t0
  t0 -> samples[[k]]
  cat(k, fill=T)

}


do.call(rbind, samples) -> samples
write.xlsx2(samples, "stats-samples.xlsx", row.names=F)

### 1 per row

### order: loci_n, total_bp, median_bp

sort(unique(samples$Terminal)) -> spp

matrix(ncol=4, nrow=length(spp)) -> samples.s
data.frame(samples.s, row.names=spp) -> samples.s
colnames(samples.s) <- colnames(samples)[-1]

for (i in 1:length(spp)) {
  spp[i] -> sp0
  samples[which(samples$Terminal == sp0),-1] -> d0
  paste0("loci_n=", d0[1,]) -> l0
  paste0("total_bp=", d0[2,]) -> l1
  paste0("median_bp_loci=", d0[3,]) -> l2
  paste(l0,l1,l2, sep="; ") -> samples.s[i,]
}

as.matrix(samples.s) -> samples.s
samples.s[samples.s == "loci_n=0; total_bp=0; median_bp_loci=0"] <- "loci_n=0"


write.xlsx2(samples.s, "stats-samples-1row.xlsx", row.names=T)

### summary

### order: loci_n, total_bp, median_bp

unique(samples$Terminal) -> spp

matrix(ncol=4, nrow=length(spp)) -> samples.s
data.frame(samples.s, row.names=spp) -> samples.s
colnames(samples.s) <- colnames(samples)[-1]

match(spp, samples$Terminal) -> l1
l1+1 -> l2
l1+2 -> l3

as.matrix(samples[l1,-1]) -> loci.n
as.matrix(samples[l2,-1]) -> total.bp
as.matrix(samples[l3,-1]) -> median.bp
mode(loci.n) <- "numeric"
mode(total.bp) <- "numeric"
mode(median.bp) <- "numeric"
loci.n[loci.n == 0] <- NA
total.bp[total.bp == 0] <- NA
median.bp[median.bp == 0] <- NA

rownames(loci.n) <- rownames(total.bp) <- rownames(median.bp) <- spp
loci.n[spp,] -> loci.n
total.bp[spp,] -> total.bp

c("Punica_granatum-gb", "Alzatea_verticillata-AM235484") -> rem
loci.n[-match(rem, spp),] -> loci.n

apply(loci.n, 2, median, na.rm=T)
apply(loci.n, 2, range, na.rm=T)
rownames(loci.n)[apply(loci.n, 2, which.min)]

total.bp[-match(rem, spp),] -> total.bp

apply(total.bp, 2, median, na.rm=T)
apply(total.bp, 2, range, na.rm=T)
rownames(loci.n)[apply(total.bp, 2, which.min)]
