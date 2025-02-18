args = commandArgs(trailingOnly=T)

if(length(args) != 2){

        cat(
        "============== do network topology ===========
        Usage:
        Rscript igraph.network.R ppi outputpath
                parameters ->
                          ppi: [file -> always interactions.xls];
                       outputpath: [path -> path for output]; \n")
        options("show.error.messages" = F)
        stop()

}


input <- normalizePath(args[1])
outputpath <- args[2]

library(igraph)
x <- read.table(input, sep = "\t", header = T)
net <- graph_from_data_frame(d = x,directed = F)


# network
network_attr <- data.frame(vcount(net), ecount(net), edge_density(net), transitivity(net, type="global"))
colnames(network_attr) <- c("vertices", "edges", "density", "clustering coefficient")

out_file <- paste(outputpath, "network.attr.xls", sep = "/")
write.table(network_attr, out_file, sep = "\t", row.names = F, quote = F)

# nodes

node_attr <- data.frame(
    degree(net, mode="all"),
    betweenness(net),
    closeness(net),
    transitivity(net, type="local")
)
colnames(node_attr) <- c("degree", "betweenness centrality", "closeness centrality", "clustering coefficient")
out_file <- paste(outputpath, "node.attr.xls", sep = "/")
write.table(node_attr, out_file, sep = "\t", quote = F, col.names = NA)


