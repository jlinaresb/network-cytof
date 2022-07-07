# =======
# Run CellNOptR
# =======
library(CellNOptR)
library(NMF)

# Arguments
# ===
setwd('~/git/network-cytof/')
example = 'ours' #ours saez toy
cell.line = 'GEN2.2'    #SIDT1 GEN2.2
expand = T
maxGens = 100

# Load model and midas
# ===
if (example == 'saez') {
  midas=CNOlist('~/git/combiMS/data/phosphos_processed/CH003.csv')
  model=readSIF('~/git/combiMS/files/model/combiMS_PKN_No_Duplication_Activation_sign_PREPROCESSED.sif')
} else if (example == 'ours') {
  midas=CNOlist(paste0('02_preprocess_cytoff/data/', cell.line, '_norm.csv'))
  pknmodel=readSIF('01_create_network/data/network2.sif')
} else if (example == 'toy'){
  data(CNOlistToy2,package="CellNOptR")
  data(ToyModel2,package="CellNOptR")
  model = ToyModel2
  midas = CNOlist(CNOlistToy2)
  rm(CNOlistToy2, ToyModel2)
}

# Compress model
# ===
if (example == 'ours') {
  indices = indexFinder(midas, pknmodel, verbose = T)
  model = compressModel(pknmodel, indices)
}


# Plot Model
# ===
numInteractions=length(model$reacID)
sprintf("*********After compressing without expanding logic gates, the model has %s nodes and %s reactions",
        length(model$namesSpecies),numInteractions)

plotModel(model,
          midas,
          graphvizParams=list(fontsize=35,nodeWidth=1,nodeHeight=1))

# Option to expand model
if (expand == T) {
  model = preprocessing(midas, model,
                        expansion = T, compression = F, cutNONC = F,
                        verbose = T)
  plotModel(model,
            midas,
            graphvizParams=list(fontsize=35,nodeWidth=1,nodeHeight=1))
}


# Cunt and plot model
# ===
initBstring = sample(1, length(model$reacID), replace = T)
cutAndPlot(midas,
           model,
           list(initBstring),
           plotParams=list(maxrow=25, margin=0.1, width=20, height=20)) #simulate


# First time point (at 10 mins)
# ===
Opt = gaBinaryT1(CNOlist=midas,
                 model=model,
                 initBstring=initBstring,
                 stallGenMax=600,
                 maxTime=10*60, 
                 maxGens=maxGens,
                 verbose=T,
                 popSize=100,
                 elitism=2)     

# plot fits
cutAndPlot(model=model,
           bStrings=list(Opt$bString),
           CNOlist=midas, 
           plotPDF=FALSE, 
           plotParams=list(maxrow=25, margin=0.1, width=40, height=20, cex = .5))

# plot objective function
plotFit(optRes=Opt)

# plot network
plotModel(model, midas, Opt$bString)
bestNw=Opt$bString

# compare with PKN
bs = mapBack(model, pknmodel, Opt$bString)
plotModel(pknmodel, midas, bs, compressed = model$speciesCompressed)


# Second time point (at 20 mins)
# ===
OptT2 = gaBinaryTN(CNOlist=midas,
                   model=model,
                   bStrings= list(Opt$bString),
                   stallGenMax=600,
                   maxTime=10*60, 
                   maxGens=maxGens,
                   verbose=T,
                   popSize=100,
                   elitism=2)

#plot fits
cutAndPlot(model=model,
           bStrings=list(Opt$bString, OptT2$bString),
           CNOlist=midas, 
           plotPDF=F,
           plotParams=list(maxrow=25, margin=0.1, width=40, height=20, cex = .5))

# plot objective function
plotFit(optRes=OptT2)

# plot network
# plotModel(model, midas, bString = list(Opt$bString, OptT2$bString))
# bestNw=Opt$bString



# Third time point (at 30 mins)
# ===
OptT3 = gaBinaryTN(CNOlist=midas,
                   model=model,
                   bStrings= list(Opt$bString, OptT2$bString),
                   stallGenMax=600,
                   maxTime=10*60, 
                   maxGens=maxGens,
                   verbose=T,
                   popSize=100,
                   elitism=2,
                   timeIndex = 3)

#plot fits
cutAndPlot(model=model,
           bStrings=list(Opt$bString, OptT2$bString, OptT3$bString),
           CNOlist=midas, 
           plotPDF=F,
           plotParams=list(maxrow=25, margin=0.1, width=20, height=20))

# plot objective function
plotFit(optRes=OptT3)

# compare with PKN

# PKN model
plotModel(pknmodel, midas)
# At time 1 (10)
bs1 = mapBack(model, pknmodel, Opt$bString)
plotModel(pknmodel, midas, bs1, compressed = model$speciesCompressed)
# At time 2 (60)
bs2 = mapBack(model, pknmodel, OptT2$bString)
plotModel(pknmodel, midas, bs2, compressed = model$speciesCompressed)
# At time 3 (120)
bs3 = mapBack(model, pknmodel, OptT3$bString)
plotModel(pknmodel, midas, bs3, compressed = model$speciesCompressed)



# interactionsTN = buildBitString(list(Opt$bString,OptT2$bString,
#                                      OptT3$bString))
# 
# 
# plotModel(model,
#           midas,
#           bString=interactionsTN$bs,
#           indexIntegr=interactionsTN$bsTimes)



#   
# calculate average times an interaction survives in the all networks within relTol
# there is a best network expressed as interaction vector and there is the same for all networks within a relTol 
sprintf("*********the number of networks within relTols is %s", dim(Opt$stringsTol)[1])
numEdges=dim(Opt$stringsTol)[2]           
cat(sprintf("*******the number of edges in each network is %s", numEdges),sep="\n")    


countInteractions=numeric(length=dim(Opt$stringsTol)[2])
for(i in 1:dim(Opt$stringsTol)[2]){
  countInteractions[i]=sum(Opt$stringsTol[,i])
}
averageNw=countInteractions/dim(Opt$stringsTol)[1]*100

# *************************** calculate RMSE for each signal across all experiments

#midas@signals contains the value at 0 and 5, hence we take the second
#then we take the second again, since rows are stimuli and columns signals
numSignals=dim(midas@signals[[2]])[2]
numStimuli=dim(midas@signals[[2]])[1]
rmseAll=numeric(length=numSignals) 

for(i in 1:numSignals){
  rmseAll[i]=sqrt((sum((midas@signals[[2]][,i] - prova$simResults[[1]][,i])^2)/numStimuli))
}

Opt_results_list=list("averageNw"=averageNw,
                      "rmseAll"=rmseAll,
                      "countInteractions"=countInteractions,
                      "countNws"=dim(Opt$stringsTol)[1],
                      "OptT1"=Opt,
                      "OptT2"=OptT2,
                      "OptT3"=OptT3)   


saveRDS(Opt_results_list,
        file = paste0('03_run_CellNOptR/data/results_', cell.line, '.RDS'))