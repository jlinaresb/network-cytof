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
interactions = import_omnipath_interactions() %>% as_tibble()

# Convert to igraph objects:
OPI_g = interaction_graph(interactions = interactions )


# Proteins of Interest
# ===
POI = tibble(protein = c('STAT1', 'STAT3', 'STAT4', 'RPS6', 'MAPK3',
 'MAPK1', 'MAPKAPK2', 'AKT1', 'RELA', 'CREB1', 'PTPN11'))
estimulations = c('TLR9', 'TLR7')

# Drug targets
# ===
## parse data from XML and save it to memory
# read_drugbank_xml_db("~/projects/networks/data/drugbank_all_full_database.xml/full database.xml")
# dbparser::drugs(save_csv = T, csv_path = '~/projects/networks/data/drugbank_all_full_database.xml/')
# drugs = dbparser::run_all_parsers()


# Quality control
# ===
POI = POI %>% mutate(in_OP = protein %in% interactions$target_genesymbol)

# Build network
# ===
collected_path_nodes = list()
for (i in seq_along(estimulations)) {
  paths = shortest_paths(OPI_g, 
                         from = estimulations[i],
                         to = POI$protein,
                         output = 'vpath')
  path_nodes = lapply(paths$vpath, names) %>% unlist() %>% unique()
  collected_path_nodes[[i]] = path_nodes
}
collected_path_nodes = unlist(collected_path_nodes) %>% unique()  
nodes = c(POI$protein, estimulations, collected_path_nodes) %>% unique()

network = induced_subgraph(graph = OPI_g, vids = nodes)

# Annotate network
# ===
V(network)$node_type = ifelse(
  V(network)$name %in% POI$protein, "POI",
  ifelse(
    V(network)$name %in% estimulations, "direct estimulation", "intermediate node"))

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

