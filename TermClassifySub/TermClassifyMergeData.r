args <- commandArgs(trailingOnly = TRUE)

info_file <- args[1] #"../test_out_predict_rdt/path_tbl.csv"
rdt_pred_file <- args[2] #"../test_out_predict_rdt/RFPredict/prediction.csv"
rit_pred_file <- args[3] #"../test_out_predict_rit/RFPredict/prediction.csv"
file_out <- args[4] #"prob_rdt_rit_combined.csv"
rdt_name <- args[5] #"rho_dependent"
rit_name <- args[6] #"Independent"

df_rdt <- read.csv(rdt_pred_file)
df_rit <- read.csv(rit_pred_file)
df_lookup <- read.csv(info_file)

rdt_colname <- paste("Prob_", rdt_name, sep="")
rit_colname <- paste("Prob_", rit_name, sep="")

df_lookup <- df_lookup[,1:4]
pr_rdt <- df_rdt[[rdt_colname]]
pr_rit <- df_rit[[rit_colname]]

pdata <- data.frame(Names_seq = as.character(df_rit$Names_seq), 
                    Prob_Independent = pr_rit, 
                    Prob_rho_dependent = pr_rdt)

pdata <- merge(df_lookup, pdata, by = c("Names_seq"))
write.csv(pdata, file_out, row.names=FALSE)