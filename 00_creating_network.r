# Convert to Gene symbols
require(AnnotationDbi)
require(OmnipathR)
require(org.Hs.eg.db)
require(dplyr)
require(igraph)
require(ggplot2)
require(ggraph)
require(dbparser)
require(XML)


# Initialise OmniPath database
# ===
# Download protein-protein interactions
interactions = import_omnipath_interactions(resources=c('SignaLink3', 'PhosphoSite', 'SIGNOR')) %>% as_tibble()

# Convert to igraph objects:
OPI_g = interaction_graph(interactions = interactions)

# Proteins of Interest
# ===
poi = c('STAT1', 'STAT3', 'STAT4', 'RPS6', 'MAPK3',
 'MAPK1', 'MAPKAPK2', 'AKT1', 'RELA', 'CREB1', 'PTPN11')
estimulations = c('TLR9', 'TLR7')


# Quality control
# ===
poi = poi[which(poi %in% interactions$target_genesymbol == TRUE)]
estimulations = estimulations[which(estimulations %in% interactions$source_genesymbol == TRUE)]


# Extract paths
# ===
collected_path_nodes = list()
for (i in seq_along(estimulations)) {
  paths = shortest_paths(OPI_g, 
                         from = estimulations[i],
                         to = poi,
                         output = 'vpath')
  path_nodes = lapply(paths$vpath, names) %>% unlist() %>% unique()
  collected_path_nodes[[i]] = path_nodes
}
collected_path_nodes = unlist(collected_path_nodes) %>% unique()  
nodes = c(poi, estimulations, collected_path_nodes) %>% unique()

# Build network
# ===
network = induced_subgraph(graph = OPI_g, vids = nodes)

# Annotate network
# ===
# Vertix
V(network)$node_type = ifelse(
  V(network)$name %in% estimulations, "estimulated",
  ifelse(
    V(network)$name %in% poi, "measured", "intermediated"))

# Edges
E(network)$direction = ifelse(
  E(network)$is_inhibition == 1, 'inhibition', 'activation')


# Plotting
# ===
ggraph(
  network,
  layout = "lgl",
  area = vcount(network)^2.3,
  repulserad = vcount(network)^1.2,
  coolexp = 1.1
) +
  geom_edge_link(
    aes(
      start_cap = label_rect(node1.name),
      end_cap = label_rect(node2.name)),
    arrow = arrow(length = unit(4, 'mm')
    ),
    edge_width = .5,
    edge_alpha = .2
  ) +
  geom_node_point() +
  geom_node_label(aes(label = name, color = node_type)) +
  scale_color_discrete(
    guide = guide_legend(title = 'Node type')
  ) +
  theme_bw() +
  xlab("") +
  ylab("") +
  ggtitle("Induced network")


# Export to cytoscape
# ===
require(RCy3)
createNetworkFromIgraph(
  network,
  title = "From igraph")


# Convert to data.frame
# ===
df = as.data.frame(get.edgelist(network))
df$direction = ifelse(get.edge.attribute(network)$direction == 'activation', 1, -1)
names(df) = c('source', 'target', 'direction')
df = df[, c(1, 3, 2)]

write.table(df,
            quote = F,
            col.names = F,
            row.names = F,
            sep = '\t',
            file = '~/projects/networks-cytof/networks/network_v3.sif')
