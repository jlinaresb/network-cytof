# ===
# Create MIDAS file
# ===
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(CellNOptR)

# Arguments
# ===
plotting = F
setwd('~/git/network-cytof/')

# Load data
# ===
inputPath = 'extdata/cytof/DATOS medianas nuevos SIDT1.xlsx'
data = readxl::read_excel(inputPath,
                          col_names = T,
                          skip = 1)


# Remove pTyrosine and p38 (not included in the network)
# ===
data = data[, - grep('pTyrosi', names(data))]
data = data[, - grep('p38', names(data))]

# Change values names
# ===
names(data)[which(names(data) == 'time of stimulation')] = 'time'
data$stimulation[which(data$stimulation == 'Imiquimob (TLR7 agonist)')] = 'TLR7'
data$stimulation[which(data$stimulation == 'ODN (TLR9 agonist)')] = 'TLR9'
data = data[which(data$stimulation != 'NO STIM/NO Anticuerpos'),]         # removing technical controls
data$stimulation[which(data$stimulation == 'NO STIM/Anticuerpos')] = 'medium'


# Define variables
# ===
celline_data = data
markers = names(celline_data)[grep('^p', names(celline_data))]
timepoints = unique(data$time)
cell_lines = unique(data$cell_line)
stimulations = unique(data$stimulation)


# Melt data
# ===
celline_data_melted = melt(celline_data,
                           measure.vars = markers,
                           variable.name = "markers",
                           value.name = "signal")


# Plotting
# ===
# Plotting time course 
if(plotting == T){
  ggplot(ext, aes(as.numeric(time), scaled_signal, col=stimulation)) +
    geom_point() +
    geom_line() +
    facet_wrap(~markers) + 
    ggtitle('Cell line') +
    theme_bw() +
    xlab("time")
  
  # Signal distributions
  ggplot(ext) + 
    geom_violin(aes(x = stimulation, y = signal, fill = stimulation), alpha = .3) +
    facet_wrap(~markers, scales = 'free_y') +
    ggtitle('Signal distributions') +
    theme_bw() +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
}


# Create CNO object
# ===
cno_data = celline_data_melted

# Cues
cno = list()
cno$namesCues = stimulations[2:3]
cno$namesStimuli = cno$namesCues
cno$namesInhibitors = c()
cno$namesSignals = markers
cno$timeSignals = timepoints

cno$valueCues = diag(3)
cno$valueCues = cno$valueCues[,-1]
colnames(cno$valueCues) = cno$namesCues

cno$valueInhibitors = matrix(nrow = 3, ncol = 0)
cno$valueStimuli = cno$valueCues[, cno$namesStimuli, drop=FALSE]


# Signals
cno_data = cno_data[order(cno_data$markers), 2:6]
cno_data = cno_data[,-1]

valueSignals = vector('list', length = length(timepoints))
names(valueSignals) = timepoints

# control 
cno_data_t = dplyr::filter(cno_data, time==0)
control_Signal = dcast(cno_data_t, 
                   time ~ markers, 
                   mean,
                   value.var = "signal")
control_Signal = tibble::add_column(control_Signal, stimulation = 'medium', .after = 1)

# timepoints
cno_data_no_ctrls = cno_data[-which(cno_data$time == 0), ]
for(time_ind in 2:length(timepoints)){
  
  act_time = timepoints[[time_ind]]
  cno_data_filter = dplyr::filter(cno_data_no_ctrls, time==act_time)
  formated_S = reshape2::dcast(cno_data_filter, 
                     time + stimulation ~ markers, 
                     value.var = "signal",
                     fill = NA)
  formated_S = rbind(control_Signal, formated_S)
  # adds the time 0 to each timepoint
  # -- end MIDAS v3
  
  valueSignals[[time_ind]] = as.matrix(formated_S[,markers])
}
cno$valueSignals = valueSignals
cno$valueSignals$`0` = as.matrix(rbind(control_Signal[,3:12], control_Signal[,3:12], control_Signal[,3:12]))



# Save MIDAS file
# ===
outDir = '02_preprocess_cytoff/data/'
setwd(outDir)
file = paste0(cell_lines, '.csv')
cnolist = CNOlist(cno)

writeMIDAS(cnolist, filename = file, overwrite = T)


