library(DESeq2)
library(piano)

setwd('/Users/emanuel/Projects/projects/tcga_rna_seq/')
# setwd('/Users/emanuel/Projects/projects/pymist/resources/fh_cells/')

# Data-sets
tumours <- c('KIRC', 'KIRP', 'KICH')

# Preprocess and run gsea
lapply(tumours, function (type) {
  # Datasets file paths
  tumour <- paste('data_preprocessed/', type, '_rsem_raw_counts_tumour.tab', sep='')
  normal <- paste('data_preprocessed/', type, '_rsem_raw_counts_normal.tab', sep='')
  
  # Import data
  tumour <- read.table(tumour, sep='\t', header=T, row.names=1, check.names=F, stringsAsFactors=F)
  normal <- read.table(normal, sep='\t', header=T, row.names=1, check.names=F, stringsAsFactors=F)
  
  # Assemble test data-set
  dataset <- cbind(tumour, normal)
  
  # Remove rows with zero or very low counts
  dataset <- dataset[rowSums(dataset) > 20, ]
  
  colData <- data.frame(row.names=colnames(dataset), condition=c(rep('tumour', dim(tumour)[2]), rep('normal', dim(normal)[2])))
  
  dds <- DESeqDataSetFromMatrix(countData=round(dataset), colData=colData, design=~condition)
  
  dds <- DESeq(dds)
  
  result <- results(dds, pAdjustMethod='fdr')
  
  # Gene-set enrichment analysis
  emt <- read.table('/Users/emanuel/Projects/projects/pymist/resources/fh_cells/files/human_byers_emt.tab', sep='\t', header=T, check.names=F, stringsAsFactors=F)[c('gene', 'signature')]
  
  piano_results <- runGSA(as.data.frame(result['stat']), geneSetStat='gsea', gsc=loadGSC(emt), nPerm=100)
  
  GSAsummaryTable(piano_results)
  
  pdf(paste('/Users/emanuel/Projects/projects/pymist/resources/fh_cells/plots/emt_', type, '.pdf', sep=''))
  barcodeplot(
    result[,'stat'],
    rownames(result['stat']) %in% emt[which(emt$signature == 'EMT'), ]$gene, 
    #   labels=c(
    #     paste('up-regulated (p-value = ', signif(camera_results['emt_up',]$FDR, 3), ')', sep=''),
    #     paste('down-regulated (p-value = ', signif(camera_results['emt_down',]$FDR, 3), ')', sep='')
    #   ), 
    main='EMT - Tumour vs Normal'
  )
  dev.off()
})