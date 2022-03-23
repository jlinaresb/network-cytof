require(igraph)
autoNameForGraph <- function(inGraph) {
  vlabels <- get.vertex.attribute(inGraph, "name")
  if (is.null(vlabels)) {
    vlabels <- get.vertex.attribute(inGraph, "label")
    if (is.null(vlabels)) {
      vlabels <- paste("BF", 1:vcount(inGraph), sep="")
    } 
    outGraph <- set.vertex.attribute(inGraph, "name", value=vlabels); 
  } else {
    outGraph <- inGraph
  }
  return (outGraph)
}

#
# Hand this function the graph, the name of the output file, and an edge label.
#

igraphToSif <- function(inGraph, outfile="output.sif", edgeLabel="label") {
  
  sink(outfile)
  
  singletons <- as.list(get.vertex.attribute(inGraph, "name"))
  edgeList <- get.edgelist(inGraph, names=FALSE)
  nodeNames <- get.vertex.attribute(inGraph, "name")
  numE <- ecount(inGraph)
  
  for (i in 1:numE) {
    node1 <- edgeList[i,1]
    node2 <- edgeList[i,2]
    singletons <- singletons[which(singletons != nodeNames[node1])]
    singletons <- singletons[which(singletons != nodeNames[node2])]
    cat(nodeNames[node1], "\t", edgeLabel, "\t", nodeNames[node2], "\n")
  }
  
  for (single in singletons) {
    cat(single, "\n")
  }
  
  sink()
}

#
# Test snippet
#

# igraphToSif(graph, "~/tmp/myGraph.sif", "myLabel")
