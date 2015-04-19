library(limma)

# Files paths
organism <- 1 # 1: human, 2: mouse

setwd('/Users/emanuel/Projects/projects/pymist/resources/fh_cells')

# Read rma normalised transcriptomics files
data_df <- read.table(c(
  'human_transcriptomics_annot_rma.tab',
  'mouse_transcriptomics_annot_rma.tab'
)[organism], sep='\t', header=T, row.names=1)

design_df <- list(
  cbind(WT=rep(c(1, 0), each=3), KO=rep(c(0, 1), each=3)),
  cbind(KO=rep(c(1, 0), each=2), WT=rep(c(0, 1), each=2))
)[[organism]]

