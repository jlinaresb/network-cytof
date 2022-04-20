# Run CellNOptR
# ====
# Example from: https://bioconductor.org/packages/release/bioc/vignettes/CellNOptR/inst/doc/CellNOptR-vignette.pdf

library(CellNOptR)
library(igraph)

# Load data
# ===
pknmodel = readSIF('~/projects/networks-cytof/networks/network_v5_modified.sif')
plotModel(pknmodel)

cno = readMIDAS('~/projects/networks-cytof/data/MIDAS/GEN2.2.csv')
cno = makeCNOlist(cno, subfield = F)

checkSignals(cno, pknmodel)
plotModel(model = pknmodel, CNOlist = cno)

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

model = preprocessing(cno, pknmodel, expansion = T,
                      compression = T, cutNONC = T, verbose = T)

# Training the model
# ===
initBstring = rep(1,length(model$reacID))
opt = gaBinaryTN(CNOlist = CNOlistToy, model = model,
                 initBstring = initBstring, verbose = T)

?gaBinaryTN
# Plot results
# ===
cutAndPlot(model=model, bStrings=list(opt$bString),
           CNOlist=CNOlistToy,plotPDF=F)

plotModel(model, CNOlistToy, bString = opt$bString)
bs = mapBack(model, pknmodel, opt$bString)
plotModel(pknmodel, CNOlistToy, bs, compressed=model$speciesCompressed)
