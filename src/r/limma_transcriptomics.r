library(limma)
library(piano)

# Files paths
organism <- 2 # 1: human, 2: mouse

setwd('/Users/emanuel/Projects/projects/pymist/resources/fh_cells/')

# Read rma normalised transcriptomics files
data_df <- read.table(c(
  'files/human_transcriptomics_annot_rma.tab',
  'files/mouse_transcriptomics_annot_rma.tab'
)[organism], sep='\t', header=T, row.names=1)

result_file <- c(
  'files/human_transcriptomics.tab',
  'files/mouse_transcriptomics.tab'
)[organism]

design_df <- list(
  cbind(WT=rep(c(1, 0), each=3), KO=rep(c(0, 1), each=3)),
  cbind(KO=rep(c(1, 0), each=2), WT=rep(c(0, 1), each=2))
)[[organism]]

# Perform limma differential analysis
fit_1 <- lmFit(data_df, design_df)
cont_matrix <- makeContrasts(KOvsWT=KO-WT, levels=design_df)
fit_2 <- contrasts.fit(fit_1, cont_matrix) 
fit_2 <- eBayes(fit_2)

result <- as.data.frame(topTable(fit_2, adjust.method='fdr', n=Inf))

# Store reults
write.table(result, result_file, sep='\t', quote=F, row.names=T)

# Gene-set enrichment analysis
emt <- read.table(c(
  'files/human_byers_emt.tab',
  'files/mouse_byers_emt.tab'
)[organism], sep='\t', header=T, check.names=F, stringsAsFactors=F)[c('gene', 'signature')]

piano_results <- runGSA(result['logFC'], geneSetStat='gsea', gsc=loadGSC(emt), nPerm=100)
GSAsummaryTable(piano_results)
geneSetSummary(piano_results, 'EMT')

index <- list(
  emt_up = rownames(result['logFC']) %in% emt[which(emt$signature == 'EMT'), ]$gene,
#   emt_down = rownames(result['logFC']) %in% emt[which(emt$signature == 'EMT-DOWN'), ]$gene
)

camera_results <- camera(data_df, index, design_df, nrot=10000)

pdf(paste('plots/emt_', c('human', 'mouse')[organism], '.pdf', sep=''))
barcodeplot(
  as.numeric(unlist(result['logFC'])),
  rownames(result['logFC']) %in% emt[which(emt$signature == 'EMT'), ]$gene, 
#   index$emt_down, 
#   labels=c(
#     paste('up-regulated (p-value = ', signif(camera_results['emt_up',]$FDR, 3), ')', sep=''),
#     paste('down-regulated (p-value = ', signif(camera_results['emt_down',]$FDR, 3), ')', sep='')
#   ), 
  main='EMT - KO vs WT'
)
dev.off()