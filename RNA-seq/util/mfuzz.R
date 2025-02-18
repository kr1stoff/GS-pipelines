#!/usr/bin env Rscript

args <- commandArgs(trailingOnly = TRUE)

input   <- args[1]
samples <- args[2]
out_dir   <- args[3]
type      <- args[4]



library(Mfuzz)

data <- read.table(input, header=T, sep="\t", row.names=1, comment.char="", check.names=F)
sample <- unlist(strsplit(samples, split = ","))

sample


if (type == "transcript") {
    data <- data[, -1]
}


data <- data[, sample]

count <- data.matrix(data)

eset <- new("ExpressionSet",exprs = count)
eset <- filter.std(eset,min.std=0)
eset <- standardise(eset)
c <- 6
m <- mestimate(eset)
cl <- mfuzz(eset, c = c, m = m)


table <- paste(out_dir, "/gene.cluster.xls", sep="" )

write.table(cl$cluster, table, sep = "\t",  quote = F)

pdf_file <- paste(out_dir, "/gene.cluster.pdf", sep="" )

pdf(pdf_file)
mfuzz.plot(
eset,
cl,
mfrow=c(2,3),
new.window= FALSE)
dev.off()


