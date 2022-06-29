require(CellNOptR)

setwd('~/git/network-cytof/')

midas=CNOlist('02_preprocess_cytoff/data/GEN2.2_norm.csv')
pknmodel=readSIF('01_create_network/data/network.sif')

res = readRDS('03_run_CellNOptR/data/results.RDS')


indices = indexFinder(midas, pknmodel, verbose = T)
model = compressModel(pknmodel, indices)

# model = preprocessing(midas, model,
#                       expansion = T, compression = F, cutNONC = F,
#                       verbose = T)

setwd('~/git/network-cytof/04_cytoscape/data/')
writeScaffold(
  modelComprExpanded = model,
  optimResT1 = res$OptT1,
  optimResT2 = res$OptT5,
  modelOriginal = pknmodel, 
  CNOlist = midas
)
