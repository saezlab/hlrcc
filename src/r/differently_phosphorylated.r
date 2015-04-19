library(limma)

# Files paths
organism = 1 # 1: human, 2: mouse

pp_human_file = c(
  '~/Projects/data/fh_cells/human_phosphoproteomics/b1368p100_phosho_human_processed.tab', 
  '~/Projects/data/fh_cells/mouse_phosphoproteomics/b1368p100_phosho_mouse_processed.tab'
)[organism]

tp_human_file = c(
  '~/Projects/data/fh_cells/human_proteomics/b1368p100_protein_human_processed.tab', 
  '~/Projects/data/fh_cells/mouse_proteomics/b1368p100_protein_mouse_processed.tab'
)[organism]

pp_result_file = c(
  '~/Projects/projects/pymist/resources/fh_cells/files/b1368p100_phospho_human_limma.tsv',
  '~/Projects/projects/pymist/resources/fh_cells/files/b1368p100_phospho_mouse_limma.tsv'
)[organism]

tp_result_file = c(
  '~/Projects/projects/pymist/resources/fh_cells/files/b1368p100_protein_human_limma.tsv',
  '~/Projects/projects/pymist/resources/fh_cells/files/b1368p100_protein_mouse_limma.tsv'
)[organism]

# Import data-sets
pp <- read.table(pp_human_file, sep='\t', header=T, stringsAsFactors=F, check.names=F)
tp <- read.table(tp_human_file, sep='\t', header=T, stringsAsFactors=F, check.names=F)

# ---- Phosphoproteomics
# Preprocess data-set
pp_ss <- pp[,1:3]
pp <- pp[,-(1:3)]

# Run differential analysis
pp_design <- cbind(KO=rep(c(1,0), each=6), WT=rep(c(0,1), each=6))
pp_fit <- lmFit(pp, pp_design)
pp_cont_matrix <- makeContrasts(KOvsWT=KO-WT, levels=pp_design)
pp_fit_2 <- contrasts.fit(pp_fit, pp_cont_matrix) 
pp_fit_2 <- eBayes(pp_fit_2)

pp_result <- as.data.frame(topTable(pp_fit_2, adjust.method='fdr', n=Inf))
pp_result <- cbind(pp_ss[rownames(pp_result),], pp_result)
pp_result$p.value.log10 <- -log10(pp_result$adj.P.Val)

# Store reults
write.table(pp_result, pp_result_file, sep='\t', quote=F, row.names=F)

# ---- Proteomics
# Preprocess data-set
tp_ss <- data.frame(tp[,1:3])
tp <- tp[,-(1:3)]

# Run differential analysis
tp_design <- cbind(KO=rep(c(1,0), each=3), WT=rep(c(0,1), each=3))
tp_fit <- lmFit(tp, tp_design)
tp_cont_matrix <- makeContrasts(KOvsWT=KO-WT, levels=tp_design)
tp_fit_2 <- contrasts.fit(tp_fit, tp_cont_matrix)
tp_fit_2 <- eBayes(tp_fit_2)

tp_result <- as.data.frame(topTable(tp_fit_2, adjust.method='fdr', n=Inf))
tp_result <- cbind(tp_ss[rownames(tp_result), ], tp_result)
tp_result$p.value.log10 <- -log10(tp_result$adj.P.Val)

# Store reults
write.table(tp_result, tp_result_file, sep='\t', quote=F, row.names=F)