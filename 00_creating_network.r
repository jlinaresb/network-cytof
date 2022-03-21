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
POI = tibble(protein = c('STAT1', 'STAT3', 'STAT4', 'RPS6', 'MAPK3', 'MAPK1', 'MAPKAPK2', 'AKT', 'NFKB1', 'CREB', 'PTPN11'))
# poi = mapIds(org.Hs.eg.db,
#              keys = proteins,
#              column = 'UNIPROT',
#              keytype = 'ALIAS')


# Drug targets
# ===
## parse data from XML and save it to memory
read_drugbank_xml_db("~/projects/networks/data/drugbank_all_full_database.xml/full database.xml")


dbparser::drugs(save_csv = T, csv_path = '~/projects/networks/data/drugbank_all_full_database.xml/')

drugs = dbparser::run_all_parsers()



# Quality control
# ===
POI = POI %>% mutate(in_OP = protein %in% interactions$target_genesymbol)

# Build network
# ===
network = induced_subgraph(graph = OPI_g, vids = c(1:10))

# Plotting
# ===
ggraph(
  network,
  layout = "lgl",
  area = vcount(network)^2.3,
  repulserad = vcount(network)^1.2,
  coolexp = 1.1
) +
  geom_node_point() +
  theme_bw() +
  xlab("") +
  ylab("") +
  ggtitle("Induced network")

