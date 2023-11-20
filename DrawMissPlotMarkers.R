#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
library(ggplot2)

if (length(args)==0) {
  stop("One argument must be supplied (input file)", call.=FALSE)
}
fileName = args[1]

lmiss <- read.table(paste(fileName,".lmiss",sep=""), header=T)
try(lmiss[is.na(lmiss$F_MISS),]$F_MISS <- 0, silent=TRUE)
try(lmiss[lmiss$F_MISS<=0.0001,]$F_MISS <- 0.00011)

pdf(paste0(fileName, "_lmiss.pdf"), width=13/2.54, height=13/2.54, pointsize=12)
ggplot(data = lmiss, aes(x = log10(F_MISS))) + geom_histogram(bins = 50) +
  theme_bw() +
  scale_x_continuous(breaks = c(seq(-4, 0, 1)), labels = c(0.0001, 0.001, 0.01, 0.1, 1), limits = c(-5, 0)) +
  geom_vline(aes(xintercept = log10(0.10)), linetype=2) +
  xlab("Proportion of missing data") +
  ylab("Number of SNPs")
dev.off()
