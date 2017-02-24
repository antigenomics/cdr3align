library(reshape2) 
library(gplots)
library(plyr) 
library(MASS)

folder = "seqdata/HomoSapiens/indels/"
outfolder = "seqdata/HomoSapiens/indels/"
inpname = 'indels_indels_'
info_from_indels <- function(folder, inpname) {
df <- read.csv(paste(folder, inpname, "in.csv", sep =""), stringsAsFactors = F) 
df <- melt(df)
df$type <- "inner"
df.1 <- read.csv(paste(folder, inpname,"out.csv", sep =""), stringsAsFactors = F) 
df.1 <- melt(df.1)
df.1$type <- "outer"
df <- rbind(df, df.1)
colnames(df) <- c("to", "from", "count", "type") 
print(df)
write.table(df, file = paste(folder, inpname, "all", sep = ""))
}
df <- read.csv(paste(folder, inpname, "in", sep =""), sep = "\t", stringsAsFactors = F) 
df <- melt(df)
df$type <- "inner"
df.1 <- read.csv(paste(folder, inpname, "out", sep =""), sep = "\t", stringsAsFactors = F) 
df.1 <- melt(df.1)
df.1$type <- "outer"
df.1$X18 <- NULL
df.1$X19 <- NULL
print(df.1)
df <- rbind(df, df.1)
colnames(df) <- c("len", "pos", "#indels", "type")
write.table(df, file = paste(folder, inpname, "all", sep = ""))