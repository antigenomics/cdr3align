library(reshape2) 
library(gplots)
library(plyr) 
library(MASS)

info_from_matrices <- function(inpname, diagonal=FALSE) {
folder = "seqdata/HomoSapiens/"
outfolder = "seqdata/HomoSapiens/pictures/"
inpname = inpname

df <- read.csv(paste(folder, "msubstitution_matrix_", inpname, ".csv", sep =""), stringsAsFactors = F) 
df <- melt(df)
df$type <- "outer"
df.1 <- read.csv(paste(folder, "substitution_matrix_", inpname, ".csv", sep =""), stringsAsFactors = F) 
df.1 <- melt(df.1)

df.1$type <- "inner"
df <- rbind(df, df.1)
colnames(df) <- c("from", "to", "count", "type") 
df <- subset(df, from != "C" & to != "C")

df$count <- df$count + 1
# Pii
df <- ddply(df, .(from, type), transform, P = count / sum(count))
# Pij
df <- ddply(df, .(from, type), transform, P = ifelse(from == to, P,
                                                     (1 - P[which(from == to)]) * count / sum(count[which(from != to)])))
df.inv <- data.frame(from = df$to, to = df$from, type = df$type, count = df$count, P = df$P) 
df <- merge(df, df.inv, by = c("from", "to", "type", "count"))
df$P <- mapply(function(x, y) min(x, y), df$P.x, df$P.y)                                 
df.1 <- ddply(df, .(from, to), summarize, Q = log10(P[which(type=="inner")] / P[which(type=="outer")]))
df.m <- dcast(df.1, from ~ to)

rownames(df.m) <- df.m$from 
df.m$from <- NULL
df.m <- data.matrix(df.m)
df.m <- df.m[order(rownames(df.m)), order(colnames(df.m))]
#png(paste(outfolder, inpname, '_heatmap.png', sep=''))
heatmap.2(df.m, scale="none", symm = T, trace = "none",
          #hclustfun = function(x) hclust(dist(x), method = "ward"), 
          col=colorpanel(100, "#2166ac", "#f7f7f7", "#b2182b"))
#dev.off()

colors <- c(rep("blue", 8), rep("red", 6), rep("yellow", 6))
names(colors) <- strsplit("I V L F C M A W G T S Y P H N D Q E K R", " ")[[1]]
df.mds <- isoMDS(as.dist(-df.m + 2), k = 2)
#png(paste(outfolder, inpname, '_mds.png', sep=''))
plot(df.mds$points, type = "n")
text(df.mds$points, labels = rownames(df.m), col = colors[rownames(df.m)])
#dev.off()

if (diagonal == TRUE) {
df.diagonal <- ddply (df.1, .(from, to), transform, Q = ifelse(from==to, Q, 0))
conname = "_diagonal_1"
df.mdiagonal <- dcast(df.diagonal, from ~ to)

rownames(df.mdiagonal) <- df.mdiagonal$from 
df.mdiagonal$from <- NULL
df.mdiagonal <- data.matrix(df.mdiagonal)
df.mdiagonal <- df.mdiagonal[order(rownames(df.mdiagonal)), order(colnames(df.mdiagonal))]
#png(paste(outfolder, inpname, conname, '_heatmap.png', sep=''))
heatmap.2(df.mdiagonal, scale="none", symm = T, trace = "none",
          #hclustfun = function(x) hclust(dist(x), method = "ward"), 
          col=colorpanel(100, "#2166ac", "#f7f7f7", "#b2182b"))
#dev.off()

colors <- c(rep("blue", 8), rep("red", 6), rep("yellow", 6))
names(colors) <- strsplit("I V L F C M A W G T S Y P H N D Q E K R", " ")[[1]]
df.mdsdiagonal <- isoMDS(as.dist(-df.mdiagonal + 2), k = 2)
#png(paste(outfolder, inpname, conname, '_mds.png', sep=''))
plot(df.mdsdiagonal$points, type = "n")
text(df.mdsdiagonal$points, labels = rownames(df.mdiagonal), col = colors[rownames(df.mdiagonal)])
#dev.off()
}
}