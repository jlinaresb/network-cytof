# Run CellNOptR
# ====
# Example from: https://bioconductor.org/packages/release/bioc/vignettes/CellNOptR/inst/doc/CellNOptR-vignette.pdf

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
# indices = indexFinder(CNOlistToy, pknmodel, verbose = T)
# NCNOindices = findNONC(pknmodel, indices, verbose=TRUE)
# NCNOcut = cutNONC(pknmodel, NCNOindices)
# indicesNCNOcut = indexFinder(CNOlistToy, NCNOcut)

# Compressing the model
# NCNOcutComp = compressModel(NCNOcut,indicesNCNOcut)
# indicesNCNOcutComp = indexFinder(CNOlistToy, NCNOcutComp)

# Expanding the gates
# model = expandGates(NCNOcutComp, maxInputsPerGate = 3)

model = preprocessing(CNOlistToy, pknmodel, expansion = T,
                      compression = T, cutNONC = T, verbose = T)

# Training the model
# ===
initBstring = rep(1,length(model$reacID))
opt = gaBinaryT1(CNOlist = CNOlistToy, model = model,
                 initBstring = initBstring, verbose = T)


# Plot results
# ===
cutAndPlot(model=model, bStrings=list(opt$bString),
           CNOlist=CNOlistToy,plotPDF=F)

plotModel(model, CNOlistToy, bString = opt$bString)
bs = mapBack(model, pknmodel, opt$bString)
plotModel(pknmodel, CNOlistToy, bs, compressed=model$speciesCompressed)
