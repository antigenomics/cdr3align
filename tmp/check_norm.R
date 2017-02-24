library(pROC)
library(car)

#folder <- "seqdata/HomoSapiens/scores/"
#inpname <- "pairscores_3_0_0_3.txt"
#df <- read.csv(paste(folder, inpname, sep =""), sep = '\t', stringsAsFactors = F)
#colnames(df) <- c("sameantigen", "score", "pair")
#qqPlot(df$score)
#shapiro.test(df$score)
#testroc <- roc(sameantigen ~ score, df, plot=TRUE)
#plot(testroc)
read_file_for_ROC <- function(folder, sub, del, ins, tot, name) {
  suffix <- paste(sub, del, ins, tot, sep="_")
  inpname <- paste(paste(name, suffix, sep="_"), '.txt', sep='')
  df <- read.csv(paste(folder, inpname, sep =""), sep = '\t', stringsAsFactors = F)
  colnames(df) <- c("sameantigen", "score", "pair")
  df$dataset <- suffix
  return(df)
}

folder <- "seqdata/HomoSapiens/scores/"
for (s in 1:1) {
  df <- read_file_for_ROC(folder, s, 0, 0, s, 'pairscores')
  facet_wrap()
  dfroc <- roc(sameantigen ~ score, df)
  plot(dfroc)
}

