read_data <- function(path, id) {
	pdata <- read.csv(path)
	pdata$Names_seq <- as.character(pdata$Names_seq)
	pdata$Names_seq <- paste(pdata$Names_seq, id, sep = "|")
	return(pdata)
}

args <- commandArgs(trailingOnly = TRUE)
infiles_str <- args[1] #"NA"
csv_in <- args[2] #"rit_group_1mer_bkgd_te.csv"
csv_out <- args[3] #"test_CAR_rit_bkgd"

infiles <- strsplit(infiles_str, ",")[[1]]

paths_tbl <- read.csv(csv_in, header=FALSE)
infiles <- unique(c(infiles, as.character(paths_tbl$V1)))
infiles <- infiles[file.exists(infiles)]
             
df_overall <- data.frame()
df_overall <- tryCatch(expr = {do.call(rbind, lapply(1:length(infiles), FUN=function(x){read_data(infiles[x], x)}))}, error = function(error_message){}) 

write.csv(df_overall, csv_out, row.names=FALSE)