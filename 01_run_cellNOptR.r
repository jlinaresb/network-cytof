# Run CellNOptR
# ====
library(CellNOptR)
library(igraph)

# Load data
# ===
# pknmodel.x = readSIF()
# cnolist.x = CNOlist()

data("ToyModel", package="CellNOptR")
data("CNOlistToy", package="CellNOptR")

pknmodel = ToyModel
cnolist = CNOlist(CNOlistToy)

checkSignals(CNOlistToy,pknmodel)
plotModel(model = pknmodel, CNOlist = cnolist)

# Preprocessing
# ===
indices = indexFinder(CNOlistToy, pknmodel, verbose = T)
NCNOindices = findNONC(pknmodel, indices, verbose=TRUE)
NCNOcut = cutNONC(pknmodel, NCNOindices)
indicesNCNOcut = indexFinder(CNOlistToy, NCNOcut)


# Compressing the model
# ===
