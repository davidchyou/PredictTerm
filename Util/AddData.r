args <- commandArgs(trailingOnly = TRUE)

master_file <- args[1]
tthp_file <- args[2]
rtp_file <- args[3]
out_file <- args[4]

df_master <- read.csv(master_file)
df_tthp <- read.csv(tthp_file)
df_rtp <- read.csv(rtp_file)

df_master$Names_seq <- as.character(df_master$Names_seq)
df_tthp$Names_seq <- as.character(df_tthp$Names_seq)
df_tthp$TTHP_score <- as.numeric(df_tthp$TTHP_score)
df_tthp$TTHP_no_return <- as.numeric(df_tthp$TTHP_no_return)
df_rtp$Names_seq <- as.character(df_rtp$Names_seq)
df_rtp$RTP_score <- as.numeric(df_rtp$RTP_score)
df_rtp$RTP_no_return <- as.numeric(df_rtp$RTP_no_return)

df_master <- merge(df_master, df_tthp, by = c("Names_seq"))
df_master <- merge(df_master, df_rtp, by = c("Names_seq"))

write.csv(df_master, out_file, row.names=FALSE)