#!/usr/bin/env Rscript

library("ggplot2")

args = commandArgs(trailingOnly=TRUE)
fileName = args[1]

if (length(args)==0) {
  stop("One argument must be supplied (input file)", call.=FALSE)
}

imiss=read.table(paste(fileName,".imiss",sep=""),h=T)
het=read.table(paste(fileName,".het",sep=""),h=T)
imiss$logF_MISS = log10(imiss$F_MISS)
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.
df <- merge(het, imiss, by = "IID")

pdf(paste0(fileName, "_imiss-vs-het.pdf"), width=13/2.54, height=13/2.54, pointsize=12)
ggplot(data = df, aes(x = logF_MISS, y = meanHet)) + geom_point() +
  geom_hline(aes(yintercept = mean(het$meanHet, na.rm=TRUE)-(6*sd(het$meanHet, na.rm=TRUE))), linetype=2) +
  geom_hline(aes(yintercept = mean(het$meanHet, na.rm=TRUE)+(6*sd(het$meanHet, na.rm=TRUE))), linetype=2) +
  geom_vline(aes(xintercept = log10(0.10)), linetype=2) +
  scale_x_continuous(breaks = c(seq(-5, 0, 1)), labels = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1), limits = c(-5, 0)) +
  ylim(c(0, 1)) +
  xlab("Proportion of missing genotypes") +
  ylab("Heterozygosity rate")
dev.off()


MISS <- data.frame(imiss[imiss$F_MISS > 0.10, c("FID","IID")])
HET  <- data.frame(het[het$meanHet > (mean(het$meanHet, na.rm=TRUE)+(6*sd(het$meanHet, na.rm=TRUE))) | het$meanHet < (mean(het$meanHet)-(6*sd(het$meanHet))), c("FID","IID")])

names(MISS) <- c("ID")
names(HET)  <- c("ID")

EXCLUDE <- unique(rbind(MISS,HET))

if (!all(is.na(EXCLUDE))) {
  write.table(EXCLUDE, file="IndividualsToExlcudeMissHet.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, na="")
}
