library(ape)
library(skimmingLociR)
library(outliers)

#############################################
### Dirs
#############################################

wd = "C:/SkimmingLoci_all"
assembly.dir = "C:/SkimmingLoci_all/2_Assembly"
new.align.dir = paste(assembly.dir, "/7_concatenated_filtered", sep="")
dir.create(new.align.dir)

setwd(assembly.dir)

read.table("stats.loci.final.txt", row.names=1, sep="\t") -> loci.s
read.table("stats.samples.loci.txt", row.names=1, sep="\t") -> samples.s

setwd(wd)

head(loci.s)
head(samples.s)

#############################################
### Aligns
#############################################

setwd(paste(assembly.dir, "/5_aligns_final",sep=""))

list.files(pattern=".fas") -> files
sub(".fas", "", files) -> labels

sapply(files, read.dna, "fasta", simplify=F) -> aligns
names(aligns) <- labels

setwd(wd)


#############################################
### Flag libraries
#############################################

samples.f <- vector(length=nrow(samples.s))
names(samples.f) <- rownames(samples.s)
samples.f[] <- "keep"

################################
### By number of loci

outlier(samples.s$loci_n, opposite =F) -> out.loci
out.loci

n.thresh = 30
#n.thresh = out.loci

hist(samples.s$loci_n)
abline(v=n.thresh, lty="dashed", col="red")

length(which(samples.s$loci_n <= n.thresh))
rownames(samples.s)[which(samples.s$loci_n <= n.thresh)]
samples.f[which(samples.s$loci_n <= n.thresh)] <- "remove"

################################
### By median coverage

outlier(samples.s$coverage_median, opposite = T) -> out.cov
out.cov

# n.thresh = 0.05
c.thresh = out.cov

hist(samples.s$coverage_median)
abline(v=c.thresh, lty="dashed", col="red")

rownames(samples.s)[which(samples.s$coverage_median <= c.thresh)]
samples.f[which(samples.s$coverage_median <= c.thresh)] <- "remove"


################################
### with too many ambiguities

outlier(samples.s$ambiguities_median, opposite = F) -> out.a
out.a

a.thresh = 0.04
#a.thresh = out.a

hist(samples.s$ambiguities_median)
abline(v=a.thresh, lty="dashed", col="red")

rownames(samples.s)[which(samples.s$ambiguities_median >= a.thresh)]
#samples.s[which(samples.s$ambiguities_median >= a.thresh)] <- "remove"


################################

data.frame(samples.f)

#############################################
### Flag loci
#############################################

loci.f <- vector(length=length(labels))
names(loci.f) <- labels
loci.f[] <- "keep"

################################
### By length

outlier(loci.s$aligned_bp, opposite = T) -> out.length
out.length

l.thresh = 450 
#l.thresh = out.length

hist(loci.s$aligned_bp, breaks=100)
abline(v=l.thresh, lty="dashed", col="red")

length(which(loci.s$aligned_bp <= l.thresh))
nrow(loci.s)
loci.f[which(loci.s$aligned_bp <= l.thresh)] <- "remove"


################################
### with too many ambiguities

loci.s$sites_with_ambiguities / loci.s$aligned_bp -> loci.s$ambiguities_rel

outlier(loci.s$ambiguities_rel) -> out.a
out.a

#a.thresh = out.a
a.thresh = 0.6

hist(loci.s$ambiguities_rel, breaks=100)
abline(v=a.thresh, lty="dashed", col="red")

length(which(loci.s$ambiguities_rel >= a.thresh))
loci.f[which(loci.s$ambiguities_rel >= a.thresh)] <- "remove"

################################
### coverage sd

outlier(loci.s$coverage_sd, opposite = T) -> out.c
out.c

c.thresh = out.c
c.thresh = 0.39

hist(loci.s$coverage_sd, breaks=100)
abline(v=c.thresh, lty="dashed", col="red")

length(which(loci.s$coverage_sd >= c.thresh))
loci.f[which(loci.s$coverage_sd >= c.thresh)] <- "remove"


#############################################
### New Alignments
#############################################

### Remove loci

aligns[which(loci.f == "keep")] -> aligns.f

length(aligns.f) / length(aligns) ## percent of aligns kept

### Remove libraries

names(which(samples.f == "keep")) -> keep.s
lapply(aligns.f, FUN=function(x)(x[which(is.na(match(rownames(x), keep.s))==F),])) -> aligns.f

### Trim (optional)
### This is an example of how one alignment looked before trimming

image(aligns.f[[1]])

lapply(aligns.f, trimAligns, min.missing=0.5) -> aligns.f

### Let's get rid of putative empty cells across the entire alignment

lapply(aligns.f, trimAligns, min.missing = 1, edges.only = F) -> aligns.f

## and after trimming

image(aligns.f[[1]])

#############################################
### New stats
#############################################


lapply(aligns.f, alignStats, include.amb=F) -> aligns.stats
do.call(rbind, aligns.stats) -> stats.n

#############################################
### Export
#############################################

### Export stats

setwd(assembly.dir)

write.table(stats.n, "loci.stats.filtered.txt", sep="\t", quote=F, na="")

setwd(new.align.dir)

### Export new align

fastConc(aligns.f, fill.with.gaps = T, map = T) -> conc

writeRAxML(conc$align, file="concatenated.phy", map = conc$map)

setwd(wd)
