library('pheatmap')

# Human heatmap
human_gl = read.table('/Users/emanuel/Projects/projects/pymist/resources/fh_cells/etc_proteins_human.tab', sep='\t', header=T, check.names=F, stringsAsFactors=F)
human_tp = read.table('/Users/emanuel/Projects/data/fh_cells/petros_tmt/human_proteomics.tab', sep='\t', header=T, row.names=1, check.names=F, stringsAsFactors=F)

df = t(human_tp[human_gl[,'accession'],])
colnames(df) <- human_gl[,'gene_symbol']
annot = data.frame(row.names=human_gl[,'gene_symbol'], group=human_gl[,'signature'], stringsAsFactors=F)
pheatmap(df, cluster_rows=F, cluster_cols=T, annotation=annot)

# Mouse heatmap
mouse_gl = read.table('/Users/emanuel/Projects/projects/pymist/resources/fh_cells/etc_proteins_mouse_v2.tab', sep='\t', header=T, check.names=F, stringsAsFactors=F)
mouse_tp = read.table('/Users/emanuel/Projects/data/fh_cells/petros_tmt/mouse_proteomics_v2.tab', sep='\t', header=T, row.names=1, check.names=F, stringsAsFactors=F)

df <- t(mouse_tp[mouse_gl[,'uniprot'],])
df <- df[,complete.cases(t(df))]
colnames(df) <- mouse_gl[,'gene_symbol']
annot = data.frame(row.names=mouse_gl[,'gene_symbol'], group=mouse_gl[,'signatue'], stringsAsFactors=F)
pheatmap(df, cluster_rows=F, cluster_cols=F, annotation=annot)
