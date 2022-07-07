# Annotation
# ===

setwd('~/git/network-cytof/')

annot = read.csv('extdata/cytof/annotation_markers.csv', header = T, skip = 1)
annot = annot[,-1]

# Add p38 protein as MAPK14

annot = rbind.data.frame(annot,
                 c('pp38', 'MAPK14', 1432, 'Q16539', 'measured'))

saveRDS(annot, 'extdata/cytof/annotation_markers.rds')
