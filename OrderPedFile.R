args = commandArgs(trailingOnly=TRUE)
pedfile = args[1]

fr <- read.table(pedfile, sep="\t", skip=9, header=TRUE)
fr <- fr[order(fr$Sample.ID),]

write.table(fr, "OrderedPedFile.txt", sep="\t", quote=FALSE, row.names=FALSE)