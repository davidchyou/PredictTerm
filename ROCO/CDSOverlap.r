options(warn=-1)
suppressWarnings(suppressMessages(library("Biostrings", quietly = TRUE)))
suppressWarnings(suppressMessages(library("dplyr", quietly = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
fa_path <- args[1] #"/mnt/SSD/brownlab/RDT_RIT_2nd_wave_paper/RDT_RIT_Output_paper_2/test_TermClassify_1mer/RDT/PreProcess/input_seq.proc.fna"
genome <- args[2] #"/mnt/SSD/brownlab/RDT_RIT_2nd_wave_paper/test_data_paper/GCF_000750555.fna"
gff <- args[3] #"/mnt/SSD/brownlab/RDT_RIT_2nd_wave_paper/test_data_paper/GCF_000750555.cds.gff"
blast <- args[4] #"/mnt/SSD/brownlab/RDT_RIT_2nd_wave_paper/ROCO/ROCO/input_seq_blast.bed"
conv_flag <- args[5] #"/mnt/SSD/brownlab/RDT_RIT_2nd_wave_paper/ROCO/ROCO/conv_div.txt"
mode <- args[6]
out <- args[7]

rl_head <- system(paste("grep '^>' ", fa_path," | sed 's/^>//g'", sep=""), intern = TRUE)
df_gff <- read.delim2(gff, header=FALSE)
df_blast <- read.delim2(blast, header=FALSE)
df_conv <- read.delim2(conv_flag, header=FALSE)
fa_genomes <- readDNAStringSet(genome)
ids <- names(fa_genomes)
ids <- gsub("(\\S+)\\s+.*", "\\1", ids)
glens <- width(fa_genomes)

cds_accession <- as.character(df_gff$V1)
cds_start <- as.numeric(df_gff$V4)
cds_end <- as.numeric(df_gff$V5)
cds_ori <- rep(0, length(cds_accession))
cds_ori[as.character(df_gff$V7) == "+"] <- 1
df_gff_simple <- data.frame(ACCESSION=cds_accession, START=cds_start, END=cds_end, ORI=cds_ori)

blast_qid <- as.character(df_blast$V1)
blast_accession <- as.character(df_blast$V2)
blast_start <- as.numeric(df_blast$V9)
blast_end <- as.numeric(df_blast$V10)
blast_ori <- rep(0, length(blast_qid))
blast_ori[as.character(df_blast$V13) == "plus"] <- 1
temp <- blast_start[blast_start > blast_end]
blast_start[blast_start > blast_end] <- blast_end[blast_start > blast_end]
blast_end[blast_start > blast_end] <- temp

buffer_cds_lst <- lapply(1:length(ids),FUN=function(x){rep(0,glens[x])})
#buffer_cds_lst_plus <- lapply(1:length(ids),FUN=function(x){rep(0,glens[x])})
#buffer_cds_lst_minus <- lapply(1:length(ids),FUN=function(x){rep(0,glens[x])})
names(buffer_cds_lst) <- ids
#names(buffer_cds_lst_plus) <- ids
#names(buffer_cds_lst_minus) <- ids

for (i in 1:(length(cds_accession))) {
	nbase <- length(buffer_cds_lst[[cds_accession[i]]])
	if (nbase > 0) {
		start <-  cds_start[i]
		end <- cds_end[i]
		
		if (end > nbase) {
			end <- nbase
		}
		
		if (start > nbase) {
			start <- nbase
		}
		buffer_cds_lst[[cds_accession[i]]][start:end] <- 1
		
		#t_ori <- cds_ori[i]
		#if (t_ori > 0) {
		#	buffer_cds_lst_plus[[cds_accession[i]]][start:end] <- 1
		#} else {
		#	buffer_cds_lst_minus[[cds_accession[i]]][start:end] <- 1
		#}
		
	}
	#print(i)
}

prcov <- rep(0, length(blast_qid))
prcov_samp <- rep(0, length(blast_qid))
#conv_flag <- rep(0, length(blast_qid))
for (i in 1:(length(blast_qid))) {
	#bf_plus <- buffer_cds_lst_plus[[blast_accession[i]]]
	#bf_minus <- buffer_cds_lst_minus[[blast_accession[i]]]
	#bf_combined <- bf_minus + (2 * bf_plus)
	
	nbase <- length(buffer_cds_lst[[blast_accession[i]]])
	if (nbase > 0) {
		start <-  blast_start[i]
		end <- blast_end[i]
		ori <- blast_ori[i]
		
		if (end > nbase) {
			end <- nbase
		}
		
		if (start > nbase) {
			start <- nbase
		}
		
		if (end < 1) {
			end <- 1
		}
		
		if (start < 1) {
			start <- 1
		}
		
		#chk <- round((start + end) / 2)
		#if (bf_combined[chk] > 0) {
		#	conv_flag[i] <- as.numeric((ori + 1) != bf_combined[chk])
		#} else {
		#	vlhs <- bf_combined[1:(chk-1)]
		#	vrhs <- bf_combined[-(1:chk)]
		#	vlhs <- vlhs[vlhs > 0]
		#	vrhs <- vrhs[vrhs > 0]
		#	vlhs_last <- vlhs_last[length(vlhs_last)]
		#	vrhs_first <- vrhs[1]
		#	if (ori > 0) {
		#		conv_flag[i] <- as.numeric((ori + 1) != vrhs_first)
		#	} else {
		#		conv_flag[i] <- as.numeric((ori + 1) != vlhs_last)
		#	}
		#}
		
		prcov[i] <- mean(buffer_cds_lst[[blast_accession[i]]][start:end])
		
		start <- sample(1:nbase,1,replace=TRUE)
		end <- start + 325 - 1
		prcov_samp[i] <- mean(buffer_cds_lst[[blast_accession[i]]][start:end])
	}
	#print(i)
}

df <- df_gff_simple[order(df_gff_simple$ACCESSION, df_gff_simple$START),]
df <- df %>% group_by(ACCESSION) %>% mutate(ORI_PREV=lag(ORI), ORI_NEXT=lead(ORI), START_NEXT=lead(START), END_PREV=lag(END))
df <- as.data.frame(df)

is_tandem <- rep(1, dim(df)[1])
is_tandem[df$ORI == 0] <- (as.numeric(df$ORI == df$ORI_PREV))[df$ORI == 0]
is_tandem[df$ORI == 1] <- (as.numeric(df$ORI == df$ORI_NEXT))[df$ORI == 1]
is_tandem[is.na(is_tandem)] <- 1
is_conv <- 1 - is_tandem

df_out_temp_conv <- df_conv
names(df_out_temp_conv) <- c("Names_seq", "Conv")
df_out_temp_conv$Names_seq <- as.character(df_out_temp_conv$Names_seq)

df_out_temp_prcov <- data.frame(Names_seq=blast_qid,CDS_Cover=prcov)
df_out <- data.frame(Names_seq=rl_head)
df_out <- merge(df_out, df_out_temp_conv, by = c("Names_seq"), all.x=TRUE)
df_out <- merge(df_out, df_out_temp_prcov, by = c("Names_seq"), all.x=TRUE)

n_na_conv <- sum(as.numeric(is.na(df_out$Conv)))
n_na_cover <- sum(as.numeric(is.na(df_out$CDS_Cover)))

if (n_na_conv > 0) {
	if (mode == "RDT") {
		df_out[["Conv"]][is.na(df_out[["Conv"]])] <- rbinom(n_na_conv, 1, 0.12)
	} else if (mode == "IT") {
		df_out[["Conv"]][is.na(df_out[["Conv"]])] <- rbinom(n_na_conv, 1, 0.5)
	} else {
		df_out[["Conv"]][is.na(df_out[["Conv"]])] <- sample(is_conv, n_na_conv, replace=TRUE)
	}
}

if (n_na_cover > 0) {
	df_out[["CDS_Cover"]][is.na(df_out[["CDS_Cover"]])] <- sample(prcov_samp, n_na_cover, replace=TRUE)
}

write.csv(df_out, out, row.names=FALSE)