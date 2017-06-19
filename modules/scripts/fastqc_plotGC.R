#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))

plotGC_f <- function(gc_in, gc_out, gc_thumb_out) {
    data <- read.table(gc_in, sep="\t", header=TRUE, check.names=F )
    #print(data)

    x <- data.frame(data)
    colnames(x) <- c('GC','Count')
    gc <- as.vector(x['GC'])
    count <- as.vector(x['Count'])

    #full
    png(gc_out, height=8,width=8, unit='in', res=300)
    plot(x, type='l', col=rainbow(1)[1])
    #mark the 42% mark
    #abline(v=42, col="steelblue")
    abline(v=42, col="black")

    #NOTE: if you do this before  next plot it will create an Rplot.pdf file!!
    #DON'T do this!!
    #junk <- dev.off() 
    
    #thumbnail
    png(gc_thumb_out, height=85,width=85, unit='px')
    par(mar=c(0,0,0,0))
    plot(x, type='l',col=rainbow(1)[1],bty='n', xaxs='i', yaxs='i', ann=FALSE, xaxt='n', yaxt='n', bty='n')
    #mark the 42% mark
    abline(v=42, col="black")
    junk <- dev.off()
    
}
args <- commandArgs( trailingOnly = TRUE )
arg_in = args[1]
arg_full_out = args[2]
arg_thumb_out = args[3]

plotGC_f(arg_in, arg_full_out, arg_thumb_out)

