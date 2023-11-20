args = commandArgs(trailingOnly=TRUE)
pedfile = args[1]

print(pedfile)
fr <- read.table(pedfile, sep="\t", skip=10, header=FALSE)
fr <- fr[order(fr$V2),]

write.table(fr, "OrderedPedFile.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
