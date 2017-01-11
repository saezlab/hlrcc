# Copyright (C) 2017 Emanuel Goncalves

library(limma)

wd <- '~/Projects/projects/hlrcc/'

tp_file <- paste(wd, 'data/uok262_proteomics_processed.txt', sep='')
tp_file_res <- paste(wd, 'data/uok262_proteomics_logfc.txt', sep='')

# Import data-sets
tp <- read.table(tp_file, sep='\t', header=T, stringsAsFactors=F, check.names=F, row.names=1)

# Run differential analysis
tp_design <- cbind(KO=rep(c(1,0), each=3), WT=rep(c(0,1), each=3))

tp_fit <- lmFit(tp, tp_design)

tp_cont_matrix <- makeContrasts(KOvsWT=KO-WT, levels=tp_design)

tp_fit_2 <- contrasts.fit(tp_fit, tp_cont_matrix) 
tp_fit_2 <- eBayes(tp_fit_2)

tp_result <- as.data.frame(topTable(tp_fit_2, adjust.method='fdr', n=Inf))

tp_result$p.value.log10 <- -log10(tp_result$adj.P.Val)

# Store reults
write.table(tp_result, tp_file_res, sep='\t', quote=F)
message('[INFO] Differential analysis done')
