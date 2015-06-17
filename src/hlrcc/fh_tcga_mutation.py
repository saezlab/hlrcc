import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import DataFrame, Series, read_csv
from pymist.enrichment.mutual_exclusivity import mutual_exclusivity

types = ['KICH', 'KIRC', 'KIRP']

wd, dd = '/Users/emanuel/Projects/projects/pymist/resources/fh_cells/', '/Users/emanuel/Projects/data/tcga/mutation/'

genes_to_highlight = {'FH', 'VHL', 'KDM6A', 'KDM6B'}

for tumour in types:
    print '[INFO] Tumour type: ', tumour

    tumour_mutations = read_csv(dd + tumour + '.maf', sep='\t')
    print {gene: tumour_mutations[tumour_mutations['Hugo_Symbol'] == gene].shape[0] for gene in genes_to_highlight}


tumour_mutations = read_csv(dd + 'KIRP.maf', sep='\t')

barcodes = set(tumour_mutations['Tumor_Sample_Barcode'])

array_a = {barcode: tumour_mutations.loc[tumour_mutations['Tumor_Sample_Barcode'] == barcode, 'Hugo_Symbol'].isin(['KDM6A']).any() for barcode in barcodes}
array_b = {barcode: tumour_mutations.loc[tumour_mutations['Tumor_Sample_Barcode'] == barcode, 'Hugo_Symbol'].isin(['FH']).any() for barcode in barcodes}

mutual_exclusivity(np.array(array_a.values()), np.array(array_b.values()))