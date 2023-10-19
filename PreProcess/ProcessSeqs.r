options(warn=-1)
suppressWarnings(suppressMessages(library("Biostrings", quietly = TRUE)))

args <- commandArgs(trailingOnly = TRUE)

file_in <- args[1] #"../test_data/test_set_rit.fna"
outdir <- args[2] #"./test_out_bkgd_rit_2"
genome_in <- args[3] #"../test_data/GCF_000750555.fna"
training <- as.numeric(args[4]) #1
always_add_rand_seqs <- as.numeric(args[5]) #0
mer <- as.numeric(args[6])
min_seq_len <- 40
rand_seq_probs_str <- "0.25,0.25,0.25,0.25"
rand_seq_probs <- as.numeric(strsplit(as.character(rand_seq_probs_str), ",")[[1]])[1:4]
basename_no_ext <- basename(gsub("\\.\\w+$", "", file_in))

randDNAStr <- function(len, probs) {
	seq <- paste(apply(rmultinom(len, size = 1, prob = probs), 2, FUN=function(x){c("A", "C", "G", "T")[x>0]}), collapse="")
	return(seq)
}

randDNAStrGenome <- function(len, gen) {
	nseq <- length(gen)
	ind <- sample(1:nseq,1)
	glen <- width(gen[ind])
	parts_ind <- 1:floor(glen/mer)
	#len_parts <- length(parts)
	parts_ind_samp <- sample(parts_ind, ceiling(len/mer) + 1, replace=FALSE)
	take_inds <- do.call(c, lapply(parts_ind_samp, FUN=function(x){(mer * (x-1) + 1):(mer * (x-1) + mer)}))
	seq <- gen[[ind]][take_inds[1:len]]
	#seq <- gen[[ind]][sample(1:width(gen[ind]),len,replace=FALSE)]
	return(seq)
}

fa <- readDNAStringSet(file_in)
seqid <- names(fa)

genome_seq <- c()
if (file.exists(genome_in)) {
	genome_seq <- readDNAStringSet(genome_in)
}

ntoks <- length(strsplit(as.character(seqid[1]), "\\|")[[1]])

class_name <- rep(NA, length(seqid))

if ((ntoks > 1) && (training > 0)) {
	class_name <- unlist(lapply(seqid, FUN = function(x){toks <- strsplit(as.character(x), "\\|")[[1]]; return(toks[ntoks])}))
}

score <- rep(NA, length(seqid))
if ((ntoks > 2) && (training > 0)) {
	score <- unlist(lapply(seqid, FUN = function(x){toks <- strsplit(as.character(x), "\\|")[[1]]; return(toks[ntoks - 1])}))
}

id_full_seq <- unlist(lapply(seqid, FUN = function(x){toks <- strsplit(as.character(x), "\\|")[[1]]; return(toks[1])}))

seqid_short <- paste("SEQ_", (1:length(fa)), sep="")
names(fa) <- seqid_short

chk <- as.numeric(sum(! is.na(class_name)) != 0) * (as.numeric(length(unique(class_name)) < 2) + as.numeric(always_add_rand_seqs > 0))
chk <- chk * training

combined_fa <- fa

if (chk > 0) {
	rnd_seqs <- c()
	if (length(genome_seq) > 0) {
		rnd_seqs <- unlist(lapply(1:length(fa), FUN=function(x){randDNAStrGenome(nchar(toString(fa[x])),genome_seq)}))
	} else {
		rnd_seqs <- unlist(lapply(1:length(fa), FUN=function(x){randDNAStr(nchar(toString(fa[x])),rand_seq_probs)}))
	}
	rnd_seqs <- DNAStringSet(rnd_seqs)
	
	names(rnd_seqs) <- paste("SEQ_", (1:length(fa))+length(fa), sep="")
	combined_fa <- c(fa, rnd_seqs)
	seqid_short <- c(seqid_short, names(rnd_seqs))
	
	if (length(unique(class_name)) < 2) {
		the_class <- class_name[1]
		class_name <- c(class_name, rep(paste("non-", the_class, sep=""), length(fa)))
		if (as.numeric(sum(! is.na(score)) == 0)) {
			seqid <- c(seqid, paste("RAND_SEQ_",1:length(fa),"|non-",the_class,sep=""))
		} else {
			seqid <- c(seqid, paste("RAND_SEQ_",1:length(fa),"|0|non-",the_class,sep=""))
		}
	} else {
		class_name <- c(class_name, rep("Other", length(fa)))
		if (as.numeric(sum(! is.na(score)) == 0)) {
			seqid <- c(seqid, paste("RAND_SEQ_",1:length(fa),"|Other",sep=""))
		} else {
			seqid <- c(seqid, paste("RAND_SEQ_",1:length(fa),"|0|Other",sep=""))
		}
	}
	
	id_full_seq <- c(id_full_seq, paste("RAND_SEQ_",1:length(fa),sep=""))
	
	if (as.numeric(sum(! is.na(score)) == 0)) {
		score <- c(score, rep(NA, length(fa)))
	} else {
		score <- c(score, rep(0, length(fa)))
	}
}

seq_lens <- width(combined_fa)
has_n <- rep(0, length(combined_fa))
too_short <- rep(0, length(combined_fa))
has_n[grep("N", toString(combined_fa))] <- 1
too_short[seq_lens < min_seq_len] <- 1

state <- rep("accept", length(combined_fa))
state[has_n > 0] <- "reject: contains N"
state[too_short > 0] <- "reject: sequence too short"

seqid <- gsub(",", ";", seqid)
id_full_seq <- gsub(",", ";", id_full_seq)

df_lookup <- data.frame(Names_seq = as.character(names(combined_fa)), 
                        ID_full_seq = as.character(id_full_seq), 
                        Status = as.character(state), 
                        Header_content = as.character(seqid),
                        Seq_length = as.numeric(seq_lens), 
                        Seq_score = as.numeric(score), 
                        Term_class = as.character(class_name))

combined_fa_sel <- combined_fa[as.character(df_lookup$Status) == "accept"]

bkgd_samp <- c()

if (length(genome_seq) > 0) {
	bkgd_samp <- unlist(lapply(1:length(combined_fa_sel), FUN=function(x){randDNAStrGenome(nchar(toString(combined_fa_sel[x])),genome_seq)}))
} else {
	bkgd_samp <- unlist(lapply(1:length(combined_fa_sel), FUN=function(x){randDNAStr(nchar(toString(combined_fa_sel[x])),rand_seq_probs)}))
}
bkgd_samp <- DNAStringSet(bkgd_samp)
names(bkgd_samp) <- paste("BKGD_", seq(1, length(combined_fa_sel)), sep = "")

if (dir.exists(outdir)) {
	system(paste("rm -rf ", outdir, sep=""))
}
system(paste("mkdir ", outdir, sep=""))

writeXStringSet(combined_fa_sel, paste(outdir, "/", basename_no_ext, ".proc.fna", sep=""), format="fasta", width = 2 * max(seq_lens))
writeXStringSet(bkgd_samp, paste(outdir, "/", basename_no_ext, ".bkgd.fna", sep=""), format="fasta", width = 2 * max(seq_lens))
write.table(df_lookup, file=paste(outdir, "/", basename_no_ext, ".lookup.csv", sep=""), sep=",", row.names=FALSE, quote=FALSE)
#write.csv(df_lookup, file=paste(outdir, "/", basename_no_ext, ".lookup.csv", sep=""), row.names=FALSE)

#readDNAStringSet("../../OperonExtract_bed/ecoli_complete_operon_clusters_rdt_rit_new/GCF_000027125/input_genome.op.flank3p.fna")
#readDNAStringSet("../test_data/test_set_rdt.fna")