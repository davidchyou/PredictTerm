options(warn=-1)
suppressWarnings(suppressMessages(library("dplyr", quietly = TRUE)))

args <- commandArgs(trailingOnly = TRUE)

gff_in <- args[1] #"test_set_rdt.proc-geneMode-rnie.gff"
fa_path <- args[2] #"../test_out/test_set_rdt.proc.fna"
csv_out <- args[3] #"test_set_rdt.proc.fna.rnie.csv"

df_from_RL <- function(vec) {
	proc_col <- function(x, i) {
		vi <- unlist(lapply(1:length(x), FUN=function(t){strsplit(x[t], "\t")[[1]][i]}))
		if (i == 1) {
			vi <- as.character(vi)
		} else {
			vi <- as.numeric(vi)
		}
		dftemp <- data.frame(V = vi)
		names(dftemp) <- paste("V", i, sep="")
		return(dftemp)
	} 
	dfout <- do.call(cbind, lapply(1:3, function(x){proc_col(vec, x)}))
	return(dfout)
}

rl_gff <- system(paste("awk '{printf \"%s\\t%s\\t%s\\n\", $1, gensub(/.*;E_value=(\\S+);Bit.*/, \"\\\\1\", $9),$6}' ", gff_in, sep=""), intern = TRUE)
#rl_gff_1 <- system(paste("awk '{printf \"%s\\t%s\\t%s\\n\", $1, gensub(/.*;E_value=\\S+;Bits=\\S+;sumBits=(\\S+);Note.*/, \"\\\\1\", $9),$6}' ", gff_in, sep=""), intern = TRUE)
df_gff_0 <- data.frame(Names_seq = character(), RNIE_evalue = numeric(), RNIE_bit_score = numeric())
#df_gff_1 <- data.frame(Names_seq = character(), RNIE_evalue = numeric(), RNIE_sum_bit_score = numeric())
if (length(rl_gff) > 0) {
	df_gff_0 <- df_from_RL(rl_gff)
	#df_gff_1 <- df_from_RL(rl_gff_1)
	names(df_gff_0) <- c("Names_seq","RNIE_evalue","RNIE_bit_score")
	#names(df_gff_1) <- c("Names_seq","RNIE_evalue","RNIE_sum_bit_score")
	df_gff_0$Names_seq <- as.character(df_gff_0$Names_seq)
	df_gff_0$RNIE_evalue <- as.numeric(df_gff_0$RNIE_evalue)
	df_gff_0$RNIE_bit_score <- as.numeric(df_gff_0$RNIE_bit_score)
	
	df_gff_0 <- df_gff_0 %>% group_by(Names_seq) %>% slice(which.min(RNIE_evalue))
	df_gff_0 <- as.data.frame(df_gff_0)
	df_gff_0$Names_seq <- as.character(df_gff_0$Names_seq)
	df_gff_0$RNIE_evalue <- as.numeric(df_gff_0$RNIE_evalue)
	df_gff_0$RNIE_bit_score <- as.numeric(df_gff_0$RNIE_bit_score)
}

rl_head <- system(paste("grep '^>' ", fa_path," | sed 's/^>//g'", sep=""), intern = TRUE)
df_base <- data.frame(Names_seq = rl_head)
df_base$Names_seq <- as.character(df_base$Names_seq)

df_gff <- merge(df_base, df_gff_0, by = c("Names_seq"), all.x=TRUE)
df_gff$RNIE_no_return <- rep(0, dim(df_gff)[1])
df_gff$RNIE_no_return[is.na(df_gff$RNIE_evalue)] <- 1
df_gff$RNIE_evalue[is.na(df_gff$RNIE_evalue)] <- -1 * log(1 - (1 - 1e-6))
df_gff$RNIE_bit_score[is.na(df_gff$RNIE_bit_score)] <- 0

write.table(df_gff, file=csv_out, sep=",", row.names=FALSE, quote=FALSE)