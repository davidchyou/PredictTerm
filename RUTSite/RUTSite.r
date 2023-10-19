options(warn=-1)
suppressWarnings(suppressMessages(library("stringr", quietly = TRUE)))
suppressWarnings(suppressMessages(library("zoo", quietly = TRUE)))
suppressWarnings(suppressMessages(library("dplyr", quietly = TRUE)))
suppressWarnings(suppressMessages(library("Biostrings", quietly = TRUE)))
suppressWarnings(suppressMessages(library("randomForest", quietly = TRUE)))

var_excl_train <- c("Names_seq", "Header_content")
var_excl_predict <- c("Names_seq", "Header_content", "Term_class")
var_excl_outdata <- c("Names_seq", "Header_content", "Term_class")
val_excl_bkgd_adj <- c("Names_seq","RNIE_evalue","RNIE_bit_score","RNIE_no_return", "Header_content", "Term_class")
expr_all <- "Term_class~."

bubble_geometry <- function(faseq, seqname, te) {
	start <- ifelse(length(faseq) - te > 100, te, length(faseq) - 100 + 1)
	start <- ifelse(start > 0, start, 1)
	seqin <- faseq[seq(start, length(faseq))]
	
	pc_nt <- letterFrequencyInSlidingView(seqin, 50, c("A", "C", "G", "T")) * (100 / 50)
	diff_c_g <- as.numeric(pc_nt %*% c(0,1,-1,0))
	pc_a <- as.numeric(pc_nt %*% c(1,0,0,0))
	pc_c <- as.numeric(pc_nt %*% c(0,1,0,0))
	pc_g <- as.numeric(pc_nt %*% c(0,0,1,0))
	pc_t <- as.numeric(pc_nt %*% c(0,0,0,1))
	b_c_gr_g <- as.numeric(diff_c_g > 0)
	rr <- rle(b_c_gr_g)
	offset <- cumsum(rr$lengths)
	nn <- length(rr$value)
	blk_id <- rep(1:nn, rr$lengths)
	
	df_area <- data.frame(BLK_ID = blk_id, 
	                      VAL1 = diff_c_g * b_c_gr_g,
	                      PC_A = pc_a * b_c_gr_g,
	                      PC_C = pc_c * b_c_gr_g,
	                      PC_G = pc_g * b_c_gr_g,
	                      PC_T = pc_t * b_c_gr_g)
	                      
	df_area_agg <- df_area %>% group_by(BLK_ID) %>% summarize(S_AREA = sum(VAL1),
	                                                          MAX_CG_DIFF = max(VAL1),
	                                                          MEAN_CG_DIFF = mean(VAL1),
	                                                          MAX_PC_C = max(PC_C),
	                                                          MIN_PC_G = min(PC_G),
	                                                          MEAN_PC_A = mean(PC_A),
	                                                          MEAN_PC_C = mean(PC_C),
	                                                          MEAN_PC_G = mean(PC_G),
	                                                          MEAN_PC_T = mean(PC_T))
	df_area_agg <- as.data.frame(df_area_agg)

	df_sum <- data.frame(C_GR_G = rr$value, 
	                     START = offset - rr$lengths + 1, 
	                     END = offset, 
	                     LENGTH = rr$lengths, 
	                     S_AREA = df_area_agg$S_AREA,
	                     MAX_CG_DIFF = df_area_agg$MAX_CG_DIFF,
	                     MEAN_CG_DIFF = df_area_agg$MEAN_CG_DIFF,
	                     MAX_PC_C = df_area_agg$MAX_PC_C,
	                     MIN_PC_G = df_area_agg$MIN_PC_G,
	                     MEAN_PC_A = df_area_agg$MEAN_PC_A,
	                     MEAN_PC_C = df_area_agg$MEAN_PC_C,
	                     MEAN_PC_G = df_area_agg$MEAN_PC_G,
	                     MEAN_PC_T = df_area_agg$MEAN_PC_T)	
	df_sum <- df_sum[df_sum$C_GR_G > 0,]
	
	nrow <- dim(df_sum)[1]
	seqlen <- length(seqin)
	ind_longest <- ifelse(nrow > 0,which.max(df_sum$LENGTH)[1],0)
	start_longest <- ifelse(nrow > 0,as.numeric(df_sum$START)[ind_longest],0)
	end_longest <- ifelse(nrow > 0,as.numeric(df_sum$END)[ind_longest],0)
	
	max_bubble_length <- ifelse(nrow > 0,max(df_sum$LENGTH),0)
	max_bubble_area <- ifelse(nrow > 0,max(df_sum$S_AREA),0)
	sum_bubble_length <- ifelse(nrow > 0,sum(df_sum$LENGTH),0)
	sum_bubble_area <- ifelse(nrow > 0,sum(df_sum$S_AREA),0)
	norm_sum_bubble_length <- ifelse(nrow > 0,sum_bubble_length / seqlen,0)
	norm_sum_bubble_area <- ifelse(nrow > 0,sum_bubble_area / seqlen,0)
	bubble_area_longest <- ifelse(nrow > 0,as.numeric(df_sum$S_AREA)[ind_longest],0)
	max_cg_diff_longest <- ifelse(nrow > 0,as.numeric(df_sum$MAX_CG_DIFF)[ind_longest],0)
	mean_cg_diff_longest <- ifelse(nrow > 0,as.numeric(df_sum$MEAN_CG_DIFF)[ind_longest],0)
	max_pc_c_longest <- ifelse(nrow > 0,as.numeric(df_sum$MAX_PC_C)[ind_longest],0)
	min_pc_g_longest <- ifelse(nrow > 0,as.numeric(df_sum$MIN_PC_G)[ind_longest],0)
	mean_pc_a_longest <- ifelse(nrow > 0,as.numeric(df_sum$MEAN_PC_A)[ind_longest],0)
	mean_pc_c_longest <- ifelse(nrow > 0,as.numeric(df_sum$MEAN_PC_C)[ind_longest],0)
	mean_pc_g_longest <- ifelse(nrow > 0,as.numeric(df_sum$MEAN_PC_G)[ind_longest],0)
	mean_pc_t_longest <- ifelse(nrow > 0,as.numeric(df_sum$MEAN_PC_T)[ind_longest],0)
	
	subseq <- seqin
	if (nrow > 0) {subseq <- seqin[start_longest:end_longest]}
	count_yc_dimer_rep_1 <- ifelse(nrow > 0,1000*length(str_extract_all(subseq, "(C|T)C[ACGTN]{9,13}")[[1]]) / max_bubble_length,0)
	count_yc_dimer_rep_2 <- ifelse(nrow > 0,1000*length(str_extract_all(subseq, "(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}")[[1]]) / max_bubble_length,0)
	count_yc_dimer_rep_3 <- ifelse(nrow > 0,1000*length(str_extract_all(subseq, "(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}")[[1]]) / max_bubble_length,0)
	count_yc_dimer_rep_4 <- ifelse(nrow > 0,1000*length(str_extract_all(subseq, "(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}")[[1]]) / max_bubble_length,0)
	count_yc_dimer_rep_5 <- ifelse(nrow > 0,1000*length(str_extract_all(subseq, "(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}")[[1]]) / max_bubble_length,0)
	count_yc_dimer_rep_6 <- ifelse(nrow > 0,1000*length(str_extract_all(subseq, "(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}")[[1]]) / max_bubble_length,0)
	
	count_yc_dimer_rep_1_full <- 1000*length(str_extract_all(seqin, "(C|T)C[ACGTN]{9,13}")[[1]]) / seqlen
	count_yc_dimer_rep_2_full <- 1000*length(str_extract_all(seqin, "(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}")[[1]]) / seqlen
	count_yc_dimer_rep_3_full <- 1000*length(str_extract_all(seqin, "(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}")[[1]]) / seqlen
	count_yc_dimer_rep_4_full <- 1000*length(str_extract_all(seqin, "(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}")[[1]]) / seqlen
	count_yc_dimer_rep_5_full <- 1000*length(str_extract_all(seqin, "(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}")[[1]]) / seqlen
	count_yc_dimer_rep_6_full <- 1000*length(str_extract_all(seqin, "(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}(C|T)C[ACGTN]{9,13}")[[1]]) / seqlen
	
	df1 <- as.data.frame(100*t(oligonucleotideFrequency(seqin, 1, as.prob=TRUE)))
	names(df1) <- paste("Pc", names(df1), sep = "")
	df2 <- as.data.frame(100*t(oligonucleotideFrequency(seqin, 2, as.prob=TRUE)))
	names(df2) <- paste("Pc", names(df2), sep = "")
	df3 <- as.data.frame(100*t(oligonucleotideFrequency(seqin, 3, as.prob=TRUE)))
	names(df3) <- paste("Pc", names(df3), sep = "")
	
	df_out <- data.frame(Names_seq = c(seqname),
	                     MAX_BUBBLE_LENGTH = c(max_bubble_length),
	                     MAX_BUBBLE_AREA = c(max_bubble_area),
	                     SUM_BUBBLE_LENGTH = c(sum_bubble_length),
	                     SUM_BUBBLE_AREA = c(sum_bubble_area),
	                     NORM_SUM_BUBBLE_LENGTH = c(norm_sum_bubble_length),
	                     NORM_SUM_BUBBLE_AREA = c(norm_sum_bubble_area),
	                     MAX_CG_DIFF_LONGEST_BUBBLE = max_cg_diff_longest,
	                     MEAN_CG_DIFF_LONGEST_BUBBLE = mean_cg_diff_longest,
	                     BUBBLE_AREA_LONGEST_BUBBLE = bubble_area_longest,
	                     YC_DIMER_REP_1_LONGEST_BUBBLE = count_yc_dimer_rep_1,
	                     YC_DIMER_REP_2_LONGEST_BUBBLE = count_yc_dimer_rep_2,
	                     YC_DIMER_REP_3_LONGEST_BUBBLE = count_yc_dimer_rep_3,
	                     YC_DIMER_REP_4_LONGEST_BUBBLE = count_yc_dimer_rep_4,
	                     YC_DIMER_REP_5_LONGEST_BUBBLE = count_yc_dimer_rep_5,
	                     YC_DIMER_REP_6_LONGEST_BUBBLE = count_yc_dimer_rep_6,
	                     YC_DIMER_REP_1_FULL = count_yc_dimer_rep_1_full,
	                     YC_DIMER_REP_2_FULL = count_yc_dimer_rep_2_full,
	                     YC_DIMER_REP_3_FULL = count_yc_dimer_rep_3_full,
	                     YC_DIMER_REP_4_FULL = count_yc_dimer_rep_4_full,
	                     YC_DIMER_REP_5_FULL = count_yc_dimer_rep_5_full,
	                     YC_DIMER_REP_6_FULL = count_yc_dimer_rep_6_full,
	                     MAX_PcC_LONGEST_BUBBLE = max_pc_c_longest,
	                     MIN_PcG_LONGEST_BUBBLE = min_pc_g_longest,
	                     MEAN_PcA_LONGEST_BUBBLE = mean_pc_a_longest,
	                     MEAN_PcC_LONGEST_BUBBLE = mean_pc_c_longest,
	                     MEAN_PcG_LONGEST_BUBBLE = mean_pc_g_longest,
	                     MEAN_PcT_LONGEST_BUBBLE = mean_pc_t_longest)
	df_out <- cbind(df_out, df1, df2, df3)
	return(df_out)
}

bubble_geometry_data <- function(fa, te) {
	seqid <- names(fa)
	df <- as.data.frame(bind_rows(lapply(1:length(fa), FUN=function(x){bubble_geometry(fa[[x]], seqid[x], te)})))
	#df <- do.call(rbind, lapply(1:length(fa), FUN=function(x){bubble_geometry(fa[[x]], seqid[x], te)}))
	return(df)
}

intr_data <- function(rnie_gff, fa, te) {
	seqid <- names(fa) 
	
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
		dfout <- as.data.frame(bind_cols(lapply(1:3, function(x){proc_col(vec, x)})))
		#dfout <- do.call(cbind, lapply(1:3, function(x){proc_col(vec, x)}))
		return(dfout)
	}
	
	start <- rep(te - 30 + 1, length(seqid))
	end <- rep(te, length(seqid))

	subseq <- unlist(lapply(1:length(seqid), FUN=function(x){fa[[x]][seq(start[x],end[x])]}))
	subseq <- DNAStringSet(subseq)
	names(subseq) <- seqid

	df1 <- as.data.frame(oligonucleotideFrequency(subseq, 1, as.prob=TRUE))
	names(df1) <- c("PcA_MY_1", "PcC_MY_1", "PcG_MY_1", "PcT_MY_1")
	df3 <- as.data.frame(oligonucleotideFrequency(subseq, 2, as.prob=TRUE))
	names(df3) <- paste("Pc", names(df3), "_MY_1", sep="") 
	df5 <- as.data.frame(oligonucleotideFrequency(subseq, 3, as.prob=TRUE))
	names(df5) <- paste("Pc", names(df5), "_MY_1", sep="")

	start <- rep(te + 1, length(seqid))
	end <- rep(te + 50, length(seqid))

	subseq <- unlist(lapply(1:length(seqid), FUN=function(x){fa[[x]][seq(start[x],end[x])]}))
	subseq <- DNAStringSet(subseq)
	names(subseq) <- seqid

	df2 <- as.data.frame(oligonucleotideFrequency(subseq, 1, as.prob=TRUE))
	names(df2) <- c("PcA_MY_2", "PcC_MY_2", "PcG_MY_2", "PcT_MY_2")
	df4 <- as.data.frame(oligonucleotideFrequency(subseq, 2, as.prob=TRUE))
	names(df4) <- paste("Pc", names(df4), "_MY_2", sep="") 
	df6 <- as.data.frame(oligonucleotideFrequency(subseq, 3, as.prob=TRUE))
	names(df6) <- paste("Pc", names(df6), "_MY_2", sep="")
	
	df7 <- as.data.frame(bind_cols(data.frame(Names_seq = as.character(seqid)), df1, df2, df3, df4, df5, df6))
	#df7 <- cbind(data.frame(Names_seq = as.character(seqid)), df1, df2, df3, df4, df5, df6)
	df8 <- data.frame(Names_seq = as.character(seqid), RNIE_evalue = -1 * log(1 - (1 - 1e-6)), RNIE_bit_score = 0, RNIE_no_return = 1)
	
	if (file.exists(rnie_gff)) {
		rl_gff_default <- c()
		rl_gff <- tryCatch(expr = {system(paste("awk '{printf \"%s\\t%s\\t%s\\n\", $1, gensub(/.*;E_value=(\\S+);Bit.*/, \"\\\\1\", $9),$6}' ", rnie_gff, sep=""), intern = TRUE)}, 
		                   error = function(error_message){return(rl_gff_default)})
		
		df_gff_default <- data.frame(Names_seq = character(), RNIE_evalue = numeric(), RNIE_bit_score = numeric())                   
		df_gff_0 <- df_gff_default
		if (length(rl_gff) > 0) {
			df_gff_0 <- tryCatch(expr = {df_from_RL(rl_gff)}, 
			                     error = function(error_message){return(df_gff_default)})
			if (dim(df_gff_0)[1] > 0) {
				names(df_gff_0) <- c("Names_seq","RNIE_evalue","RNIE_bit_score")
				df_gff_0$Names_seq <- as.character(df_gff_0$Names_seq)
				df_gff_0$RNIE_evalue <- as.numeric(df_gff_0$RNIE_evalue)
				df_gff_0$RNIE_bit_score <- as.numeric(df_gff_0$RNIE_bit_score)
			
				df_gff_0 <- df_gff_0 %>% group_by(Names_seq) %>% dplyr::slice(which.min(RNIE_evalue))
				df_gff_0 <- as.data.frame(df_gff_0)
				df_gff_0$Names_seq <- as.character(df_gff_0$Names_seq)
				df_gff_0$RNIE_evalue <- as.numeric(df_gff_0$RNIE_evalue)
				df_gff_0$RNIE_bit_score <- as.numeric(df_gff_0$RNIE_bit_score)
			}
		}
		
		df_base <- data.frame(Names_seq = as.character(seqid))
		df_gff <- merge(df_base, df_gff_0, by = c("Names_seq"), all.x=TRUE)
		df_gff$RNIE_no_return <- rep(0, dim(df_gff)[1])
		df_gff$RNIE_no_return[is.na(df_gff$RNIE_evalue)] <- 1
		df_gff$RNIE_evalue[is.na(df_gff$RNIE_evalue)] <- -1 * log(1 - (1 - 1e-6))
		df_gff$RNIE_bit_score[is.na(df_gff$RNIE_bit_score)] <- 0
		
		df8 <- df_gff
	}
	
	df_out <- merge(df7, df8, by = c("Names_seq"), all.x=TRUE)
	return(df_out)
}

do_background_adj <- function(df_in, df0_in) {
	ind_all <- 1:(dim(df_in)[2])
	all_vars <- names(df_in)
	ind_all_vars_adj <- ind_all[-1*c(which(all_vars %in% val_excl_bkgd_adj), grep("USER_", all_vars))]
	
	df_bkgd_adj <- df0_in[,ind_all_vars_adj]
	df_bkgd_adj_agg <- apply(df_bkgd_adj, 2, mean)
	df_bkgd_adj_agg_rep <- as.data.frame(bind_cols(lapply(names(df_bkgd_adj_agg), FUN=function(x){df<-data.frame(X=rep(df_bkgd_adj_agg[x],dim(df_in)[1]));names(df)<-c(x);return(df)})))
	#df_bkgd_adj_agg_rep <- do.call(cbind, lapply(names(df_bkgd_adj_agg), FUN=function(x){df<-data.frame(X=rep(df_bkgd_adj_agg[x],dim(df_in)[1]));names(df)<-c(x);return(df)}))
	
	df_orig <- df_in
	psc <- 1e-6
	df_orig[,ind_all_vars_adj] <- log((df_orig[,ind_all_vars_adj] + psc) / (df_bkgd_adj_agg_rep + psc), 2)
	names(df_orig)[ind_all_vars_adj] <- paste("LOG2_OE_", names(df_orig)[ind_all_vars_adj], sep = "")
	
	return(df_orig)
}

getROC_AUC = function(probs, true_Y){
    probsSort <- sort(probs, decreasing = TRUE, index.return = TRUE)
    val <- unlist(probsSort$x)
    idx <- unlist(probsSort$ix)  

    roc_y <- true_Y[idx];
    stack_x <- cumsum(roc_y == 0)/sum(roc_y == 0)
    stack_y <- cumsum(roc_y == 1)/sum(roc_y == 1)    

    auc <- sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
    return(list(stack_x=stack_x, stack_y=stack_y, auc=auc))
}

rf_lhocv <- function(data) {
	inds <- 1:(dim(data)[1])
	ind_excl_train <- which(names(data) %in% var_excl_train)
	ind_excl_predict <- which(names(data) %in% var_excl_predict)
	seqid <- as.character(data$Names_seq)
	term_class <- as.character(data$Term_class)
	
	term_class_unique <- unique(term_class)
	
	ind_samp <- c()
	for (i in 1:length(term_class_unique)) {
		class <- term_class_unique[i]
		ind_class <- inds[which(term_class == class)]
		n_class <- length(ind_class)
		n_samp <- max(c(1, floor(n_class)/2))
		ind_class_samp <- sample(ind_class, n_samp, replace=FALSE)
		ind_samp <- c(ind_samp, ind_class_samp)
	}
	
	rf1 <- randomForest(formula = as.formula(expr_all), data=data[-ind_samp,-ind_excl_train], ntree=500, importance=TRUE)
	df_prob1 <- as.data.frame(predict(rf1,data[ind_samp, -ind_excl_predict],type="prob"))
	df_prob1 <- cbind(df_prob1, data.frame(Names_seq = as.character(seqid)[ind_samp], Term_class = as.character(term_class)[ind_samp]))
	
	rf2 <- randomForest(formula = as.formula(expr_all), data=data[ind_samp,-ind_excl_train], ntree=500, importance=TRUE)
	df_prob2 <- as.data.frame(predict(rf2,data[-ind_samp, -ind_excl_predict],type="prob"))
	df_prob2 <- cbind(df_prob2, data.frame(Names_seq = as.character(seqid)[-ind_samp], Term_class = as.character(term_class)[-ind_samp]))
	
	df_prob <- rbind(df_prob1, df_prob2)
	return(df_prob)
}

roc_data_one_class <- function(roc_data, class_name) {
	true_values <- ifelse(roc_data$Term_class==class_name,1,0)
	aList <- getROC_AUC(roc_data[[class_name]], true_values)
	
	xx <- unlist(aList$stack_x) #perf@x.values[[1]]
	yy <- unlist(aList$stack_y) #perf@y.values[[1]]
	auc <- unlist(aList$auc) #aperf@y.values[[1]]
	
	lbl <- paste(as.character(class_name), " (AUC=", round(auc,4), ")", sep = "")
	lblrep <- rep(lbl, length(xx))
	
	df <- data.frame(FPR = xx, TPR = yy, LBL = lblrep)
	return(df)
}

rf_train <- function(data) {
	pdata <- data
	ind_excl_train <- which(names(pdata) %in% var_excl_train)
	rf <- randomForest(formula = as.formula(expr_all), data=pdata[,-ind_excl_train], ntree=500, importance=TRUE)
	
	prediction_for_roc_curve <- rf_lhocv(pdata)

	classes <- levels(pdata$Term_class)
	roc_data <- do.call(rbind, lapply(classes, FUN=function(x){roc_data_one_class(prediction_for_roc_curve,x)})) 
	
	prediction_for_roc_curve <- as.data.frame(prediction_for_roc_curve)
	names(prediction_for_roc_curve) <- c(as.character(levels(pdata$Term_class)), "Names_seq", "Term_class")
	
	return(list(model=rf, roc=roc_data, validation=prediction_for_roc_curve))
}

rf_predict <- function(data_in, model_in) {
	pdata <- data_in
	rf <- model_in
	
	ind_excl_predict <- which(names(pdata) %in% var_excl_predict)
	ind_excl_outdata <- which(names(pdata) %in% var_excl_outdata)
	pred <- predict(rf, pdata[,-ind_excl_predict],type="prob")
	pred <- as.data.frame(pred)
	names(pred) <- paste("Prob_", names(pred), sep = "")
	
	dfout <- cbind(data.frame(Names_seq = pdata[,1]), pred, pdata[,-ind_excl_outdata])
	
	return(dfout)
}

handle_user_data <- function(path, main_df) {
	df <- read.csv(path, stringsAsFactors = FALSE)
	names(df)[1] <- "Names_seq"
	names(df)[-1] <- paste("USER_", names(df)[-1], sep="")
	df_all <- merge(main_df, df, by = c("Names_seq"), all.x=TRUE)
	
	ind_user <- grep("^USER_", names(df_all))
	for (i in 1:length(ind_user)) {
		v <- df_all[,ind_user[i]]
		v <- as.numeric(v)
		v[is.na(v)] <- 0
		#if (class(v) == "character") {
		#	v[is.na(v)] <- "NA"
		#} else if (class(v) == "factor") {
		#	v <- as.character(v)
		#	v <- as.numeric(v)
		#	v[is.na(v)] <- 0	
		#} else {
		#	v <- as.numeric(v)
		#	v[is.na(v)] <- 0
		#}
		df_all[,ind_user[i]] <- v
	}
	return(df_all)
}

args <- commandArgs(trailingOnly = TRUE)
seq_path <- args[1] #"../PreProcess/test_PreProcess_pred_te/test_all_ds_te.proc.fna"
lookup_path <- args[2] #"../PreProcess/test_PreProcess_pred_te/test_all_ds_te.lookup.csv"
rnie_gff <- args[3] #"../RNIE/test_RNIE_pred_te/test_all_ds_te.proc.fna.work-geneMode-rnie.gff"
model_in <- args[4] #"test_RUTSite_rit_train/rf_model.rds"
outdir <- args[5] #"test_RUTSite_rit_pred"
train <- as.numeric(args[6]) #0
data_only <- as.numeric(args[7]) #0
te <- as.numeric(args[8]) #75
bkgd_path <- args[9] 
user_data_path <- args[10]

model_out <- paste(outdir, "/rf_model.rds", sep="")
roc_data_out <- paste(outdir, "/roc_curve.csv", sep="")
combined_raw_data_out <- paste(outdir, "/combined_raw_data.csv", sep="")
combined_raw_data_adj_out <- paste(outdir, "/combined_raw_data_adj.csv", sep="")
validation_data_out_combined <- paste(outdir, "/combined_validation.csv", sep="")
predict_data_out_combined <- paste(outdir, "/combined_prediction.csv", sep="")

if (dir.exists(outdir)) {
	system(paste("rm -rf ", outdir, sep=""))
}
system(paste("mkdir ", outdir, sep=""))

fa <- readDNAStringSet(seq_path)

df1 <- bubble_geometry_data(fa, te)
df2 <- intr_data(rnie_gff, fa, te)

df <- merge(df1, df2, by = c("Names_seq"))
df_cpy <- df

df_lookup <- data.frame()
df_lookup_sel_sub <- data.frame()
if (file.exists(lookup_path)) {
	df_lookup <- read.csv(lookup_path)
	df_lookup_sel <- df_lookup[as.character(df_lookup$Status) == "accept",]
	df_lookup_sel_sub <- df_lookup_sel[,which(names(df_lookup_sel) %in% c("Names_seq", "Header_content", "Term_class"))]
} else {
	seqid <- names(fa)
	class_name <- unlist(lapply(seqid, FUN = function(x){toks <- strsplit(as.character(x), "\\|")[[1]];
		                                             	 ntoks <- length(toks); 
		                                             	 return(toks[ntoks])}))
	header_content <- seqid
	df_lookup_sel_sub <- data.frame(Names_seq = as.character(seqid), Header_content = as.character(header_content), Term_class = as.character(class_name))
	
}
df <- merge(df, df_lookup_sel_sub, by = c("Names_seq"))
if (file.exists(user_data_path)) {
	df <- handle_user_data(user_data_path, df)
}
write.csv(df, combined_raw_data_out, row.names=FALSE)

df_class <- as.data.frame(table(as.character(df_lookup_sel_sub$Term_class)))
n_class <- dim(df_class)[1]
min_data_per_class <- 1
df$Term_class <- as.factor(df$Term_class)

if (dim(df_class)[2] > 1) {
	min_data_per_class <- min(as.numeric(df_class[,2]))
}

df_bkgd <- data.frame()
df_adj <- data.frame()
if (file.exists(bkgd_path)) {
	fa_bkgd <- readDNAStringSet(bkgd_path)
	seqid_bkgd <- names(fa_bkgd)
	df_lookup_sel_sub_bkgd <- data.frame(Names_seq = as.character(seqid_bkgd), Header_content = "NA", Term_class = "NA")
	df1_bkgd <- bubble_geometry_data(fa_bkgd, te)
	df2_bkgd <- intr_data(rnie_gff, fa_bkgd, te)
	df_bkgd <- merge(df1_bkgd, df2_bkgd, by = c("Names_seq"))
	df_bkgd <- merge(df_bkgd, df_lookup_sel_sub_bkgd, by = c("Names_seq"))
	df_adj <- do_background_adj(df, df_bkgd)
	if (file.exists(user_data_path)) {
		df_adj <- handle_user_data(user_data_path, df_adj)
	}
	write.csv(df_adj, combined_raw_data_adj_out, row.names=FALSE)
	df_adj$Term_class <- as.factor(df_adj$Term_class)	
}

rf_model <- list()
df_combined_validation <- data.frame()
df_pred_combined <- data.frame()
if (data_only < 1) {
	if (train > 0 && n_class <= 10 && n_class > 1 && min_data_per_class >= 5) {
		obj <- list()
		if (file.exists(bkgd_path)) {
			obj <- rf_train(df_adj)
		} else {
			obj <- rf_train(df)
		}
		df_combined_validation <- merge(obj$validation, df, by = c("Names_seq", "Term_class"))
		names(df_combined_validation) <- gsub("\\W", "_", names(df_combined_validation))
		saveRDS(obj, model_out)
		write.csv(obj$roc, roc_data_out, row.names=FALSE)
		write.csv(df_combined_validation, validation_data_out_combined, row.names=FALSE)
	} else {
		if (file.exists(model_in)) {
			obj <- readRDS(model_in)
			n_classes <- length(obj$model$classes)
			df_pred <- data.frame()
			if (file.exists(bkgd_path)) {
				df_pred <- rf_predict(df_adj, obj$model)
			} else {
				df_pred <- rf_predict(df, obj$model)
			}
			df_pred_combined <- merge(df_pred[,seq(1,n_classes+1)], df_cpy, by = c("Names_seq"))
			df_pred_combined <- merge(df_pred_combined, df_lookup_sel_sub[,1:2], by = c("Names_seq"))
			names(df_pred_combined) <- gsub("\\W", "_", names(df_pred_combined))
			write.csv(df_pred_combined, predict_data_out_combined, row.names=FALSE)
		}
	}
}
