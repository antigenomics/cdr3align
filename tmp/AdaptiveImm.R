library(igraph)
library(RColorBrewer)
library(plyr)
args<-commandArgs(TRUE)
inpname<-args[1]
print(inpname)
nodes <- read.csv("trb_sequences_teach_antigens.csv", header=T, as.is=T)
ag_count <- ddply(nodes, .(antigen), summarize, count = length(id))
nodes <- subset(nodes, !grepl("_", antigen) & antigen %in% subset(ag_count, count > 5)$antigen)

ags <- unique(nodes$antigen)
ag.colors <- data.frame(antigen = ags, 
                        color = colorRampPalette(brewer.pal(12, "Paired"))(length(ags)),
                        stringsAsFactors = F)

nodes <- merge(nodes, ag.colors)

node.colors <- nodes$color
names(node.colors) <- nodes$id

links <- read.csv(paste("matrices_from_alignment/",inpname,".txt", sep=''), header=T, row.names=1)
links2 <- subset(links, row.names(links) %in% nodes$id)
net <- graph_from_incidence_matrix(links2)
V(net)$color <- node.colors[V(net)$name]
E(net)$width <- 0.2
V(net)$size <- log2(degree(net)+1)
V(net)$label.color <- "black"
V(net)$label <- NA
pdf(paste(inpname,'.pdf', sep=''))
layout <- layout_with_fr(net)
plot(net, layout=layout)
dev.off()