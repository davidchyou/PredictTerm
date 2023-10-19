options(warn=-1)
suppressWarnings(suppressMessages(library("dplyr", quietly = TRUE)))

args <- commandArgs(trailingOnly = TRUE)

gff_in <- args[1]
fa_path <- args[2]
csv_out <- args[3]

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
	dfout <- do.call(cbind, lapply(1:2, function(x){proc_col(vec, x)}))
	return(dfout)
}

rl_gff <- system(paste("awk -F, 'NR>1 {printf \"%s\\t%s\\n\", $1,$6}' ", gff_in, sep=""), intern = TRUE)

df_gff_0 <- data.frame(Names_seq = character(), RTP_score = numeric())

if (length(rl_gff) > 0) {
	df_gff_0 <- df_from_RL(rl_gff)
	names(df_gff_0) <- c("Names_seq","RTP_score")
	df_gff_0$Names_seq <- as.character(df_gff_0$Names_seq)
	df_gff_0$RTP_score <- as.numeric(df_gff_0$RTP_score)
	
	df_gff_0 <- df_gff_0 %>% group_by(Names_seq) %>% slice(which.max(RTP_score))
	df_gff_0 <- as.data.frame(df_gff_0)
	df_gff_0$Names_seq <- as.character(df_gff_0$Names_seq)
	df_gff_0$RTP_score <- as.numeric(df_gff_0$RTP_score)
}

rl_head <- system(paste("grep '^>' ", fa_path," | sed 's/^>//g'", sep=""), intern = TRUE)
df_base <- data.frame(Names_seq = rl_head)
df_base$Names_seq <- as.character(df_base$Names_seq)

df_gff <- merge(df_base, df_gff_0, by = c("Names_seq"), all.x=TRUE)
df_gff$RTP_no_return <- rep(0, dim(df_gff)[1])
df_gff$RTP_no_return[is.na(df_gff$RTP_score)] <- 1
df_gff$RTP_score[is.na(df_gff$RTP_score)] <- 0

write.table(df_gff, file=csv_out, sep=",", row.names=FALSE, quote=FALSE)