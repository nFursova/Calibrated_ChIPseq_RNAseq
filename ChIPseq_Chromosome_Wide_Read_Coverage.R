#!/usr/bin/Rscript

require('GenomicRanges')
require('Rsamtools')
require('GenomicAlignments')
require('rtracklayer')
require('BSgenome.Mmusculus.UCSC.mm10')
require('ggplot2')
require('reshape2')

args=commandArgs(trailingOnly=TRUE)

file=args[1]
name=args[2]

bins250k <- tileGenome(seqinfo(Mmusculus), tilewidth=250000, cut.last.tile.in.chrom=TRUE) # you can use any size of tiles 

chromatinExpress<- function (r1,name){
# function returns sum of the average ChIP-seq signal based on bigwig intervals in 250 kb bins
# input: bw files; name for output file 
# output: RData file ("name.RData") saved in the working directory with GRanges object containing expression values for 250 kb bins
 
r1L <- import(r1, as="RleList")
value4<-seqlevels(r1L)
bins250k<-keepSeqlevels(bins250k, value4,pruning.mode="coarse")
binnedAverage <- function(bins, numvar, mcolname)
{
    stopifnot(is(bins, "GRanges"))
    stopifnot(is(numvar, "RleList"))
    stopifnot(identical(seqlevels(bins), names(numvar)))
    bins_per_chrom <- split(ranges(bins), seqnames(bins))
    means_list <- lapply(names(numvar),
        function(seqname) {
            views <- Views(numvar[[seqname]],
                           bins_per_chrom[[seqname]])
            viewMeans(views)
        })
    new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
    mcols(bins)[[mcolname]] <- new_mcol
    bins
}

bins <- binnedAverage(bins250k, r1L, "binned_var1")
binsTotal<-bins250k
binsTotal$score<-bins$binned_var1

save(binsTotal, file=paste(name, ".RData",sep=""))

}

chromatinExpress(file, name)

