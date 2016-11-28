#!/usr/bin/Rscript

args <- commandArgs(trailingOnly=T)
names(args) <- c("species", "gene", "mhc", 
   "min.score", "min.cdr3.count")

require(dplyr)

df <- read.table("vdjdb.slim.txt", header=T, sep = "\t") %>%
      filter(species == args["species"] &
         gene == args["gene"] &
         mhc.class == args["mhc"] &
         vdjdb.score >= as.integer(args["min.score"])) %>%
      distinct(cdr3, antigen.epitope) %>%
      group_by(antigen.epitope) %>%
      mutate(cdr3.count = length(cdr3)) %>%
      filter(cdr3.count >= as.integer(args["min.cdr3.count"])) %>%
      select(cdr3, antigen.epitope) %>% 
      droplevels()

write.table(
   df,
   "vdjdb.txt",
   quote = F, sep = "\t", row.names = F, col.names = F
)

print(summary(df))
print(paste("Unique antigens:", length(unique(df$antigen.epitope))))
print(summary(df$antigen.epitope, maxsum = 9999))
