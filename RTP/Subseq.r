options(warn=-1)
suppressWarnings(suppressMessages(library("Biostrings", quietly = TRUE)))

args <- commandArgs(trailingOnly = TRUE)

file_in <- args[1]
file_out <- args[2]
split_point <- as.numeric(args[3])
downstream <- as.numeric(args[4])

fa <- readDNAStringSet(file_in)
seqid <- names(fa)
len <- width(fa)

start <- rep(1, length(seqid))
end <- rep(split_point - 1, length(seqid))

if (downstream > 0) {
	start <- rep(split_point, length(seqid))
	end <- len
} 

subseq <- unlist(lapply(1:length(seqid), FUN=function(x){fa[[x]][seq(start[x],end[x])]}))
subseq <- DNAStringSet(subseq)
names(subseq) <- seqid

writeXStringSet(subseq, file_out, format="fasta", width = 2 * max(len))