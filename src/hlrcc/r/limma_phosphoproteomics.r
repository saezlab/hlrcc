library(limma)

pp_file <- '/Users/emanuel/Projects/projects/hlrcc/data/b1368p100_phospho_human_processed.tsv'
pp_file_res <- '/Users/emanuel/Projects/projects/hlrcc/data/b1368p100_phospho_human_limma.tsv'

# Import data-sets
pp <- read.table(pp_file, sep='\t', header=T, stringsAsFactors=F, check.names=F, row.names=1)

# Run differential analysis
pp_design <- cbind(KO=rep(c(1,0), each=6), WT=rep(c(0,1), each=6))

pp_fit <- lmFit(pp, pp_design)

pp_cont_matrix <- makeContrasts(KOvsWT=KO-WT, levels=pp_design)

pp_fit_2 <- contrasts.fit(pp_fit, pp_cont_matrix) 
pp_fit_2 <- eBayes(pp_fit_2)

pp_result <- as.data.frame(topTable(pp_fit_2, adjust.method='fdr', n=Inf))

pp_result$p.value.log10 <- -log10(pp_result$adj.P.Val)

# Store reults
write.table(pp_result, pp_file_res, sep='\t', quote=F)