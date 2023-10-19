options(warn=-1)
args <- commandArgs(trailingOnly = TRUE)

file_in <- args[1] #"../testTermPredict_retrain/TermClassifyMainProcess/Pred/TermClassify/prob_rdt_rit_combined.csv"
bed_in <- args[2] #"../testTermPredict_retrain/MakeWindows/input_seq.window.bed"
file_out <- args[3] #"result.bed"
th_rdt <- as.numeric(args[4])
th_it <- as.numeric(args[5])
th_disc <- as.numeric(args[6])
ref_bed <- args[7]
revcom <- as.numeric(args[8])

if (th_rdt < 0 || th_rdt > 1) {
	th_rdt <- 0.5
}

if (th_it < 0 || th_it > 1) {
	th_it <- 0.5
}

if (th_disc < 0 || th_disc > 1) {
	th_disc <- 0.05
}

df_ttcsv <- read.csv(file_in)
df_bed <- read.delim2(bed_in, header=FALSE)
df_ttcsv_simple <- data.frame(Window_ID = as.character(df_ttcsv$Header_content), 
                              Prob_Independent = as.numeric(df_ttcsv$Prob_Independent),
                              Prob_rho_dependent = as.numeric(df_ttcsv$Prob_rho_dependent))
df_bed_simple <- data.frame(Window_ID = as.character(df_bed$V4), 
                            Header_content = as.character(df_bed$V1), 
                            Start = as.numeric(df_bed$V2), 
                            End = as.numeric(df_bed$V3),
                            Ori = as.character(df_bed$V6)) 
df_combined <- merge(df_ttcsv_simple, df_bed_simple, by=c("Window_ID"), all.x=TRUE)
s_it <- df_combined$Prob_Independent
s_rdt <- df_combined$Prob_rho_dependent
tt_class <- rep("Unclassified", dim(df_combined)[1])
tt_class[(as.numeric(s_it > th_it) + as.numeric(s_rdt > th_rdt)) > 0] <- "IT+RDT"
tt_class[(as.numeric(s_it > th_it) * as.numeric((s_it - s_rdt) > th_disc)) > 0] <- "IT"
tt_class[(as.numeric(s_rdt > th_rdt) * as.numeric((s_rdt - s_it) > th_disc)) > 0] <- "RDT"
score_str <- paste(tt_class,s_it, s_rdt, sep=",")
df_bed_out <- data.frame(V1 = as.character(df_combined$Header_content), 
                         V2 = as.numeric(df_combined$Start), 
                         V3 = as.numeric(df_combined$End), 
                         V4 = as.character(df_combined$Window_ID),
                         V5 = score_str, V6=as.character(df_combined$Ori))

if (file.exists(ref_bed)) {
	df_ref <- read.delim2(ref_bed, header=FALSE)
	df_ref$V1 <- as.character(df_ref$V1)
	df_ref$V4 <- as.character(df_ref$V4)
	df_ref$V6 <- as.character(df_ref$V6)
	df_ref$V2 <- as.numeric(df_ref$V2)
	df_ref$V3 <- as.numeric(df_ref$V3)
	df_ref$V5 <- as.numeric(df_ref$V5)
	names(df_ref) <- c("V1a", "V2a", "V3a", "V1", "V5a", "V6a")
	
	if (revcom > 0) {
		ori_new <- rep("+", dim(df_ref)[1])
		ori_new[as.character(df_ref$V6a) == "+"] <- "-"
		df_ref$V6a <- ori_new
	}
	
	df_bed_merge <- merge(df_bed_out, df_ref, by = c("V1"), all.x=TRUE)
	g_chr_id <- as.character(df_bed_merge$V1a)
	g_seq_id <- as.character(df_bed_merge$V1)
	g_window_id <- as.character(df_bed_merge$V4)
	g_score_str <- as.character(df_bed_merge$V5)
	g_ori <- as.character(df_bed_merge$V6a)
	g_ori[is.na(g_ori)] <- "+"
	g_start_plus <- df_bed_merge$V2 + df_bed_merge$V2a - 1
	g_start_minus <- df_bed_merge$V3a - df_bed_merge$V3
	g_end_plus <- df_bed_merge$V3 + df_bed_merge$V2a
	g_end_minus <- df_bed_merge$V3a - df_bed_merge$V2 + 1
	g_start <- g_start_plus
	g_start[g_ori == "-"] <- g_start_minus[g_ori == "-"]
	g_end <- g_end_plus
	g_end[g_ori == "-"] <- g_end_minus[g_ori == "-"]
	g_start[is.na(g_start)] <- (df_bed_merge$V2)[is.na(g_start)]
	g_end[is.na(g_end)] <- (df_bed_merge$V3)[is.na(g_end)]
	g_chr_id[is.na(g_chr_id)] <- g_seq_id[is.na(g_chr_id)]
	df_bed_out <- data.frame(V1 = g_chr_id, V2 = g_start, 
                         	 V3 = g_end, V4 = g_window_id,
                         	 V5 = g_score_str, V6 = g_ori)
}

df_bed_out <- df_bed_out[order(as.numeric(df_bed_out$V2)),]
write.table(df_bed_out, file=paste(file_out, ".bed", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

out_s_it <- as.numeric(unlist(lapply(as.character(df_bed_out$V5), FUN=function(x){strsplit(x, ",")[[1]][2]})))
out_s_rdt <- as.numeric(unlist(lapply(as.character(df_bed_out$V5), FUN=function(x){strsplit(x, ",")[[1]][3]})))

df_bed_out_it <- df_bed_out
df_bed_out_it$V5 <- out_s_it
write.table(df_bed_out_it, file=paste(file_out, ".ritscore.bed", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

df_bed_out_rdt <- df_bed_out
df_bed_out_rdt$V5 <- out_s_rdt
write.table(df_bed_out_rdt, file=paste(file_out, ".rdtscore.bed", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)




