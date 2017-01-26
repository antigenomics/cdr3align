#!/usr/bin/Rscript

args <- commandArgs(trailingOnly=T)
#args = c("HomoSapiens", 0, 5)
names(args) <- c("species", "min.score", "min.cdr3.count")

require(dplyr)

df <- read.table("vdjdb.slim.txt", header=T, sep = "\t") %>%
      filter(species == args["species"] &
         vdjdb.score >= as.integer(args["min.score"])) %>%
      distinct(cdr3, antigen.epitope, v.end, j.start, gene) %>%
      group_by(antigen.epitope) %>%
      mutate(cdr3.count = length(cdr3)) %>%
      filter(cdr3.count >= as.integer(args["min.cdr3.count"])) %>%
      select(cdr3, v.end, j.start, gene, antigen.epitope) %>% 
      droplevels()

write.table(
   df,
   "vdjdb.txt",
   quote = F, sep = "\t", row.names = F, col.names = F
)

print(summary(df))
print(paste("Unique antigens:", length(unique(df$antigen.epitope))))
print(summary(df$antigen.epitope, maxsum = 9999))
