# Copyright (C) 2017 Emanuel Goncalves
library(limma)

setwd('~/Projects/hlrcc/')

# Import proteomics
proteomics <- read.table('data/proteomics_qe.csv', sep=',', header=T, stringsAsFactors=F, check.names=F, row.names=1)

# Define design table
conditions <- factor(unlist(lapply(colnames(proteomics), function (x) {substr(x, 1, nchar(x) - 3)})))
replicates <- factor(unlist(lapply(colnames(proteomics), function (x) {substr(x, 5, nchar(x) - 1)})))

design <- model.matrix(~0 + conditions + replicates)
rownames(design) <- colnames(proteomics)

# Fit
fit <- lmFit(proteomics, design)

# Constrasts fit
m_contrasts <- makeContrasts(conditions262vsconditionspfH=conditions262-conditionspfH, levels=design)

fit2 <- contrasts.fit(fit, m_contrasts)
fit2 <- eBayes(fit2)

# Export results
res <- as.data.frame(topTable(fit2, adjust.method='fdr', n=Inf))
write.table(res, 'data/proteomics_qe_fc.csv', quote=F, sep=',')
