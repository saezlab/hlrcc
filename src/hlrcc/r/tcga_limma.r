# Configure workspace
library(limma)

# Data-sets
tumours <- c('KIRC', 'KIRP', 'KICH')

setwd('/Users/emanuel/Projects/projects/tcga_rna_seq/')

# Preprocess and run gsea
lapply(tumours, function (type) {
  # Datasets file paths
  normal <- paste('data_preprocessed/', type, '_rsem_raw_counts_normal.tab', sep='')
  tumour <- paste('data_preprocessed/', type, '_rsem_raw_counts_tumour.tab', sep='')
  
  # Import data
  normal <- read.table(normal, sep='\t', header=T, row.names=1, check.names=F, stringsAsFactors=F)[, 1:3]
  tumour <- read.table(tumour, sep='\t', header=T, row.names=1, check.names=F, stringsAsFactors=F)[, 1:3]
  
  # Assemble test data-set
  dataset <- cbind(tumour, normal)
  
  # Remove rows with zero or very low counts
  dataset <- dataset[rowSums(dataset) > 20, ]
  
  ### voom + limma  
  # Design
  design <- cbind(tumour=c(rep(1, length(tumour)), rep(0, length(normal))), normal=c(rep(0, length(tumour)), rep(1, length(normal))))
  
  # VOOM normalisation
  v <- voom(dataset, design, plot=F, normalize='quantile')
  
  # Differential analysis
  fit_1 <- lmFit(v, design)
  contrast_matrix <- makeContrasts(tumour_normal=tumour-normal, levels=design)
  fit_2 <- contrasts.fit(fit_1, contrast_matrix) 
  fit_2 <- eBayes(fit_2)
  
  # Get result
  result <- as.data.frame(topTable(fit_2, adjust.method='bonferroni', n=Inf))
  
  # Store reults
  result_file <- paste('/Users/emanuel/Projects/projects/pymist/resources/fh_cells/files/', 'tcga_', type, '_limma_rsem_raw_counts.tab', sep='')
  write.table(result, result_file, sep='\t', quote=F)
})