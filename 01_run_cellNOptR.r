# Run CellNOptR
# ====
# Example from: https://bioconductor.org/packages/release/bioc/vignettes/CellNOptR/inst/doc/CellNOptR-vignette.pdf

library(CellNOptR)
library(igraph)

# Load data
# ===
pknmodel = readSIF('~/projects/networks-cytof/networks/network_v5_modified.sif')
plotModel(pknmodel)

cno = readMIDAS('~/projects/networks-cytof/data/MIDAS/GEN2.2_norm.csv')
cno = makeCNOlist(cno, subfield = F)

cno

checkSignals(cno, pknmodel)
plotModel(model = model, CNOlist = cno)

# Preprocessing
# ===
model = preprocessing(cno, pknmodel,
                      expansion = F, compression = T, cutNONC = F, 
                      verbose = T)






# Training the model
# ===
initBstring = rep(1,length(model$reacID))

cutAndPlot(cno, model, list(initBstring))



opt = gaBinaryTN(CNOlist = CNOlistToy,
                 model = model,
                 bStrings = initBstring,
                 verbose = T)

?gaBinaryTN
# Plot results
# ===
cutAndPlot(model=model, bStrings=list(opt$bString),
           CNOlist=CNOlistToy,plotPDF=F)

plotModel(model, CNOlistToy, bString = opt$bString)
bs = mapBack(model, pknmodel, opt$bString)
plotModel(pknmodel, CNOlistToy, bs, compressed=model$speciesCompressed)
