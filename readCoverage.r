library(ggplot2)  
#Adapted from: https://www.biostars.org/p/5165/

setwd("/users/PAS0107/osu6702/project/annotation/alignment")

cov=read.table("barcode01.merged.clean..sorted.bam.stats.coverage", sep="\t")  
cov[1,]  
ggplot(cov, aes(x=V3, y=V4)) + geom_bar(stat='identity') + xlab('coverage') + ylab('count')

##OR

library(Rsamtools)
bamcoverage <- function (bamfile) {
  # read in the bam file
  bam <- scanBam(bamfile)[[1]] # the result comes in nested lists
  # filter reads without match position
  ind <- ! is.na(bam$pos)
  ## remove non-matches, they are not relevant to us
  bam <- lapply(bam, function(x) x[ind])
  ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
  ## names of the bam data frame:
  ## "qname"  "flag"   "rname"  "strand" "pos"    "qwidth"
  ## "mapq"   "cigar"  "mrnm"   "mpos"   "isize"  "seq"    "qual"
  ## construc: genomic ranges object containing all reads
  ranges <- GRanges(seqnames=Rle(bam$rname), ranges=ranges, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname )
  ## returns a coverage for each reference sequence (aka. chromosome) in the bam file
  return (mean(coverage(ranges)))      
}

bamcoverage("barcode01.merged.clean..sorted.bam")
bamcoverage("barcode02.merged.clean..sorted.bam")
bamcoverage("barcode03.merged.clean..sorted.bam")
bamcoverage("barcode04.merged.clean..sorted.bam")
bamcoverage("unclassified.merged.clean..sorted.bam")
