options(warn=-1)
suppressWarnings(suppressMessages(library("Biostrings", quietly = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
file_in <- args[1] #"test_set_rdt.fna"
bed_in <- args[2] #"test_set_rdt.bed"
file_out <- args[3] #"test_set_rdt_pe.fna"
split_in <- as.numeric(args[4]) #75
class <- args[5] #"Independent"
train <- as.numeric(args[6])
upto <- 325

fa <- readDNAStringSet(file_in)
seqid <- names(fa)
len <- width(fa)

df_split <- data.frame()
if (file.exists(bed_in)) {
	df <- read.delim2(bed_in, header=FALSE)
	df_split <- data.frame(ID = as.character(df$V1), SPLIT = as.numeric(df$V2))
} else {
	df_split <- data.frame(ID = seqid, SPLIT = rep(split_in, length(seqid)))
}

df_seq <- data.frame(ID = seqid, LENGTH = len)
df_seq <- merge(df_seq, df_split, by = c("ID"), all.x=TRUE)
df_seq$SPLIT[is.na(df_seq$SPLIT)] <- 75

start <- as.numeric(df_seq$SPLIT) - as.numeric(df_seq$SPLIT) + 1
end <- as.numeric(df_seq$SPLIT) + (upto - as.numeric(df_seq$SPLIT))

b1 <- as.numeric(start > 0)
b2 <- as.numeric(end <= as.numeric(df_seq$LENGTH))
b3 <- b1 * b2
df_seq_sel <- df_seq[which(b3 > 0),]

fa_sel <- fa[which(names(fa) %in% as.character(df_seq_sel$ID))]

if (train > 0) {
	names(fa_sel) <- paste(names(fa_sel), class, sep="|")
}

start <- as.numeric(df_seq_sel$SPLIT) - as.numeric(df_seq_sel$SPLIT) + 1
end <- as.numeric(df_seq_sel$SPLIT) + (upto - as.numeric(df_seq_sel$SPLIT))

subseq <- unlist(lapply(1:length(names(fa_sel)), FUN=function(x){fa_sel[[x]][seq(start[x],end[x])]}))
subseq <- DNAStringSet(subseq)
names(subseq) <- names(fa_sel)

#fa_sel
writeXStringSet(subseq, file_out, format="fasta", width = 2 * max(width(subseq)))
