###############################################################################################################
### varscan2annovar.R > to convert varscan VCF to annovar input file
###############################################################################################################

# setting arguments from command line
args = commandArgs(trailingOn=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  args[2] = "out.txt"
}

# Creating annovar dataframe
df_varscan <- read.delim(args[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
df_annovar <- cbind(df_varscan[,1:2], df_varscan[,2:ncol(df_varscan)]) #add a duplication of position
colnames(df_annovar) <- c("chrom", "start", "end", "ref", "alt", "normal_reads1", "normal_reads2", "normal_var_freq", 
                          "normal_gt", "tumor_reads1", "tumor_reads2", "tumor_var_freq", "tumor_gt", "somatic_status", 
                          "variant_p_value", "somatic_p_value", rep("info", 8)) #adjust colnames

# Tunning indels 
for(i in 1:nrow(df_annovar)){
  # deletions
  if(grepl("^-", df_annovar$alt[i])){
    df_annovar$ref[i] <- gsub("-", "", df_annovar$alt[i]) #add sequence to ref
    df_annovar$alt[i] <- "-" #delete sequence from alt (it's a deletion)
    df_annovar$end[i] <- df_annovar$start[i] + nchar(df_annovar$ref[i]) #sum deletion length (to indicate size of deletion)
    df_annovar$start[i] <- df_annovar$start[i] + 1 #add 1 because deletion begins one base after ref
  }
  # insertions
  if(grepl("^\\+", df_annovar$alt[i])){
    df_annovar$ref[i] <- "-"
    df_annovar$alt[i] <- gsub("\\+", "", df_annovar$alt[i])
  }
}

write.table(df_annovar, args[2], quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

