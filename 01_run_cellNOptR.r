# Run CellNOptR
# ====
# Example from: https://bioconductor.org/packages/release/bioc/vignettes/CellNOptR/inst/doc/CellNOptR-vignette.pdf

library(CellNOptR)
library(igraph)

# Load data
# ===
pknmodel = readSIF('~/projects/networks-cytof/networks/network_v5.sif')
plotModel(pknmodel)

cno = readMIDAS('~/projects/networks-cytof/data/MIDAS/GEN2.2_norm.csv')
cno = makeCNOlist(cno, subfield = F)

checkSignals(cno, pknmodel)
plotModel(model = pknmodel, CNOlist = cno)

# Preprocessing
# ===
model = preprocessing(cno, pknmodel,
                      expansion = F, compression = T, cutNONC = F, 
                      verbose = T)

plotModel(model = model, CNOlist = cno)


# HASTA AQUÍ!!!
# ========  


# Training the model
# ===
initBstring = rep(1,length(model$reacID))
cutAndPlot(cno, model, list(initBstring))
opt = gaBinaryTN(CNOlist = cno,
                 model = model,
                 bStrings = list(initBstring),
                 verbose = T)

# Plot results
# ===
cutAndPlot(model=model, 
           bStrings=list(opt$bString),
           CNOlist=cno,
           plotPDF=F)

plotModel(model, cno, bString = opt$bString)
bs = mapBack(model, pknmodel, opt$bString)
plotModel(pknmodel, CNOlistToy, bs, compressed=model$speciesCompressed)
