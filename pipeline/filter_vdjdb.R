#!/usr/bin/Rscript
set.seed(51103)

args = commandArgs(trailingOnly=T)
#args = c("HomoSapiens", "TRB", 0, 5)
names(args) = c("species", "gene", "min.score", "min.cdr3.count")

require(dplyr)

df = read.table("vdjdb.slim.txt", header=T, sep = "\t") %>%
      filter(species == args["species"] & gene == args["gene"] &
         vdjdb.score >= as.integer(args["min.score"])) %>%
      distinct(cdr3, antigen.epitope, gene, v.end, j.start) %>%
      group_by(antigen.epitope) %>%
      mutate(cdr3.count = length(cdr3)) %>%
      filter(cdr3.count >= as.integer(args["min.cdr3.count"])) %>%
      dplyr::select(cdr3, v.end, j.start, gene, antigen.epitope) %>% 
      droplevels()

print(summary(df))
print(paste("Unique antigens:", length(unique(df$antigen.epitope))))
print(summary(df$antigen.epitope, maxsum = 9999))

df = as.data.frame(df)
df$cdr3 = as.character(df$cdr3)
df$gene = as.character(df$gene)
df$antigen.epitope = as.character(df$antigen.epitope)

df.2 = read.table("naive.pool.txt", header=T, sep = "\t", stringsAsFactors=F) %>%
        filter(species == args["species"] & gene == args["gene"]) %>%
        dplyr::select(cdr3, v.end, j.start, gene, antigen.epitope) %>%
        sample_n(100000)

df.2 = as.data.frame(df.2)

write.table(
   rbind(df, df.2),
   "vdjdb.txt",
   quote = F, sep = "\t", row.names = F, col.names = F
)

