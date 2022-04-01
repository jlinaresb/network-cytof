# Convert to Gene symbols
require(AnnotationDbi)
require(OmnipathR)
require(org.Hs.eg.db)
require(dplyr)
require(igraph)
require(ggplot2)
require(ggraph)
# require(dbparser)
# require(XML)

# Required packages
# ===
require(AnnotationDbi)
require(OmnipathR)
require(org.Hs.eg.db)
require(dplyr)
require(igraph)
require(ggplot2)
require(ggraph)

annot = read.csv('~/projects/networks-cytof/annotations/annotation_markers.csv', header = F, skip = 1)
annot = annot[, -1]
names(annot) = annot[1,]
annot = annot[-1, ]
row.names(annot) = 1:nrow(annot)
head(annot)

# Load interactions
# ===
interactions = import_omnipath_interactions(resources = 'SIGNOR') %>% as_tibble()
OPI_g = interaction_graph(interactions = interactions)

# Quality control
# ===
poi = annot[which(annot$annot == 'measured'), ]$`gene symbol`
estimulations = annot[which(annot$annot == 'estimulated'), ]$`gene symbol`

poi = poi[which(poi %in% interactions$target_genesymbol == TRUE)]
estimulations = estimulations[which(estimulations %in% interactions$source_genesymbol == TRUE)]

setdiff(annot$`gene symbol`, poi)

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
# ggraph(
#   network,
#   layout = "lgl",
#   area = vcount(network)^2.3,
#   repulserad = vcount(network)^1.2,
#   coolexp = 1.1
# ) +
#   geom_edge_link(
#     aes(
#       start_cap = label_rect(node1.name),
#       end_cap = label_rect(node2.name)),
#     arrow = arrow(length = unit(4, 'mm')
#     ),
#     edge_width = .5,
#     edge_alpha = .2
#   ) +
#   geom_node_point() +
#   geom_node_label(aes(label = name, color = node_type)) +
#   scale_color_discrete(
#     guide = guide_legend(title = 'Node type')
#   ) +
#   theme_bw() +
#   xlab("") +
#   ylab("") +
#   ggtitle("Induced network")


# Export to cytoscape
# ===
# require(RCy3)
# createNetworkFromIgraph(
#   network,
#   title = "all")


# Convert to data.frame
# ===
df = as.data.frame(get.edgelist(network))
df$direction = ifelse(get.edge.attribute(network)$direction == 'activation', 1, -1)
names(df) = c('source', 'target', 'direction')
df = df[, c(1, 3, 2)]

# Filter
# ===
curation = get.edge.attribute(network)$curation_effort
df$curation = curation
# df = df[which(df$curation > 1),]


df = subset(df, select = -c(curation))


# Change names to proteins
# ====
# In poi's
v = unique(c(df$source, df$target))
poi_g = intersect(annot$`gene symbol`, v)
poi_p = annot[match(intersect(annot$`gene symbol`, v), annot$`gene symbol`),]$marker
names(poi_p) = poi_g


# In intermediates
intermediate_g = setdiff(v, annot$`gene symbol`)
complexes_g = intermediate_g[grep('_', intermediate_g)]

intermediate_g = setdiff(intermediate_g, complexes_g)
intermediate_p = intermediate_g
names(intermediate_p) = intermediate_g

# In complexes
if (length(complexes_g > 0)) {
  cplx = import_omnipath_complexes(resources = 'SIGNOR')
  complexes_p = cplx[match(complexes_g, cplx$components_genesymbols),]$name
  names(complexes_p) = complexes_g
}

v.names = c(poi_p, intermediate_p, complexes_p)

# 

s.names = v.names[match(unique(df$source), names(v.names))]
t.names = v.names[match(unique(df$target), names(v.names))]
df$source = factor(df$source, labels = s.names)
df$target = factor(df$target, labels = t.names)


# Save table to .sif
# =====
write.table(df,
            quote = F,
            col.names = F,
            row.names = F,
            sep = '\t',
            file = '~/projects/networks-cytof/networks/network_v5.sif')


# Read network with cellNOptR
# ===
require(CellNOptR)
model = readSIF('~/projects/networks-cytof/networks/network_v5.sif')
plotModel(model, 
          stimuli = 'TLR9',
          signals = poi_p)



# Complex model example
# plotModel(readSIF('~/R/x86_64-pc-linux-gnu-library/4.1/CellNOptR/DREAMmodel/LiverPKNDREAM.sif'))
