# ===
# Create MIDAS file
# ===
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(CellNOptR)


# Load data
# ===
inputPath = '~/projects/networks-cytof/data/cytof/cytof_jlb_v2.xlsx'
data = readxl::read_excel(inputPath,
                          col_names = T,
                          skip = 1)


# Remove pTyrosine (not included in the network)
# ===
data = data[, - grep('pTyrosi', names(data))]

# Change values names
# ===
names(data)[which(names(data) == 'time of stimulation')] = 'time'
data$stimulation[which(data$stimulation == 'Imiquimob (TLR7 agonist)')] = 'TLR7'
data$stimulation[which(data$stimulation == 'ODN (TLR9 agonist)')] = 'TLR9'
data = data[which(data$stimulation != 'control'),]         # removing technical controls


# Define variables
# ===
celline_data = data
markers = names(celline_data)[grep('^p', names(celline_data))]
timepoints = unique(data$time)
cell_lines = unique(data$cell_line)


# Normalize (??)
# ===
celline_data[, markers] = apply(celline_data[, markers], 2, function(x) log2(x + 1))


# Melt data
# ===
celline_data_melted = melt(celline_data,
                           measure.vars = markers,
                           variable.name = "markers",
                           value.name = "signal")


# Scale data
# ===
ext = ddply(celline_data_melted, .(markers), function(df){
  
  q99 = quantile(df$signal, probs = c(0.005, 0.995))
  df$scaled_signal = (df$signal - q99[[1]]) / (q99[[2]] - q99[[1]])
  df$scaled_signal[df$scaled_signal > 1] = 1
  df$scaled_signal[df$scaled_signal < 0] = 0
  return(df)

})

# save data
# saveRDS(ext, "~/projects/networks-cytof/data/cytof/cellline_scaled.rds")

# Plotting
# ===
# Plotting time course 
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


# Create CNO object
# ===
cell_lines = unique(celline_data_melted$cell_line)
stimulations = unique(celline_data_melted$stimulation)

cno_data = ext

cno = list()
cno$namesCues = c('TLR7', 'TLR9')
cno$namesStimuli = c('TLR7', 'TLR9')
cno$namesInhibitors = c()
cno$namesSignals = markers
cno$timeSignals = timepoints

cno$valueCues = diag(3)
cno$valueCues = cno$valueCues[,-1] 
colnames(cno$valueCues) = cno$namesCues

cno$valueInhibitors = matrix(nrow = 3, ncol = 0)
cno$valueStimuli = cno$valueCues[, cno$namesStimuli, drop=FALSE]


# Values by time
# ===
cno_data = cno_data[order(cno_data$markers), 2:6]

valueSignals = vector('list', length = length(timepoints))
names(valueSignals) = timepoints

for(time_ind in 1:length(timepoints)){
  
  act_time = timepoints[[time_ind]]
  
  cno_data_t = filter(cno_data, time==act_time)
  formated_S = dcast(cno_data_t, 
                     time + stimulations ~ markers, 
                     value.var = "signal",
                     fill = NA)
  
  # MIDAS v3 ---
  if(time_ind == 1){
    control_Signal = formated_S[1,]
    # save time point 0 and use this values  along each timepoint.
  }
  formated_S = rbind(control_Signal, formated_S)
  # adds the time 0 to each timepoint
  # -- end MIDAS v3
  
  valueSignals[[time_ind]] = as.matrix(formated_S[,markers])
}
cno$valueSignals = valueSignals


# Save MIDAS file
# ===
outDir = '~/projects/networks-cytof/data/MIDAS/'
setwd(outDir)
file = paste0(cell_lines, '.csv')
writeMIDAS(CNOlist(cno), filename = file, overwrite = T)
