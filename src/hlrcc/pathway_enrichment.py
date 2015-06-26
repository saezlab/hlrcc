import re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from hlrcc import wd
from scipy.stats.stats import pearsonr
from statsmodels.stats.multitest import multipletests
from pymist.enrichment.gsea import gsea
from pandas import DataFrame, Series, read_csv
from bioservices import KEGG, KEGGParser, UniProt

sns.set_style('white')

# ---- Import Kinase enrichment
kinases_es = read_csv(wd + '/files/kinase_activity.tab', sep='\t', index_col=0, header=None).to_dict()[1]

# ---- Perform GSEA on the kinase activities
# Set-up KEGG bioservice
kegg, kegg_parser = KEGG(cache=True), KEGGParser()
kegg.organism = 'hsa'
print '[INFO] KEGG service configured'

kegg_pathways = {p: kegg.parse_kgml_pathway(p) for p in kegg.pathwayIds}
print '[INFO] KEGG pathways extracted: ', len(kegg_pathways)

# Convert KEGG pathways Gene Name to UniProt
kegg_to_uniprot = kegg.conv('uniprot', 'hsa')
kegg_pathways = {p: {kegg_to_uniprot[x].split(':')[1] for i in kegg_pathways[p]['entries'] if i['type'] == 'gene' for x in i['name'].split(' ') if x in kegg_to_uniprot} for p in kegg_pathways}
kegg_pathways = {p: kegg_pathways[p] for p in kegg_pathways if len(kegg_pathways[p].intersection(kinases_es)) > 2}
kegg_pathways_names = {p: re.findall('NAME\s*(.*) - Homo sapiens\n?', kegg.get(p))[0] for p in kegg_pathways}
print '[INFO] KEGG pathways ids converted to UniProt: ', len(kegg_pathways)

# ---- Pathway enrichment analysis
pathways_es = {p: gsea(kinases_es, kegg_pathways[p], 10000) for p in kegg_pathways}
pathways_es = DataFrame(pathways_es, index=['es', 'pvalue']).T
pathways_es['adj.pvalue'] = multipletests(pathways_es['pvalue'], method='fdr_bh')[1]
pathways_es['name'] = [kegg_pathways_names[i] for i in pathways_es.index]
pathways_es['intersection'] = [len(kegg_pathways[i].intersection(kinases_es)) for i in pathways_es.index]
pathways_es = pathways_es.sort('adj.pvalue', ascending=False)
print '[INFO] Pathways enrichment: ', len(pathways_es)

# ---- Plot enrichment
plot_df = pathways_es[pathways_es['pvalue'] < .05]
colours, y_pos = sns.color_palette('Paired', 2), [x + 1.5 for x in range(len(plot_df['name']))]

plt.barh(y_pos, -np.log10(plot_df['pvalue']), lw=0, align='center', height=.5, color=colours[0], label='p-value')
plt.barh(y_pos, -np.log10(plot_df['adj.pvalue']), lw=0, align='center', height=.5, color=colours[1], label='FDR')
plt.yticks(y_pos, plot_df['name'])

plt.axvline(-np.log10(0.05), ls='--', lw=0.4, c='gray')
plt.axvline(-np.log10(0.01), ls='--', lw=0.4, c='gray')

plt.text(-np.log10(0.05) * 1.01, .5, '5%', ha='left', color='gray', fontsize=9)
plt.text(-np.log10(0.01) * 1.01, .5, '1%', ha='left', color='gray', fontsize=9)

sns.despine()
plt.xlabel('-log10')
plt.title('KEGG pathways enrichment')
plt.legend(loc=4)
plt.savefig('%s/reports/kegg_pathway_kinase_activity_enrichment.pdf' % wd, bbox_inches='tight')
plt.close()
print '[INFO] Pathways enrichment plotted'