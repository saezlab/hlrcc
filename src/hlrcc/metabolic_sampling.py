__author__ = 'emanuel'

import sys
sys.path.extend(['/Users/emanuel/Projects/projects/pymist', '/Users/emanuel/Projects/projects/pymist/pymist'])

import time
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
from reader.sbml_reader import read_sbml_model
from sampler import sample, fix_futile_cycles
from reduce import reduce_model
from pandas import read_csv, DataFrame, Series, melt
from simulation import min_differences, pFBA, FVA
from balance import elements_atoms
from utils.plot_utils import save_plot

# General variables
data_dir = '/Users/emanuel/Projects/data/fh_cells/'
conditions = ['UOK262', 'UOK262pFH']

# Import metabolic model
model = read_sbml_model('/Users/emanuel/Projects/resources/metabolic_models/recon1.xml')
model.remove_b_metabolites()

# Import medium
medium_elements = read_csv('/Users/emanuel/Projects/resources/mediums/DMEM_41966.txt', sep='\t').dropna()

# ---- Recon1 small fixes
# Reactions not used
for reaction in ['R_EX_retpalm', 'R_EX_retpalm_e_']:
    model.set_constraint(reaction, 0)
####

# Import medium
medium = read_csv(data_dir + 'human_metabolomics/human_metabolomics.txt', sep='\t', index_col=0)

sampling, reduced_models, fitted_mediums, fitted_solutions = {}, {}, {}, {}
# Fit metabolomics data and perform model reduction
for condition in conditions:
    # condition = conditions[0]
    print '[INFO] Condition: ', condition

    # Deepcopy metabolic model
    reduced_model = model.deepcopy()

    # Add FH KO to tumour cell lines
    if condition == 'UOK262':
        reduced_model.set_constraint('R_FUM', 0, 0)
        reduced_model.set_constraint('R_FUMm', 0, 0)
        print '[INFO] R_FUM blocked'

    # Minimise differences between simulations and measurements
    for reaction in reduced_model.get_exchanges():
        if reaction != 'R_EX_o2_e_':
            if reaction in medium.index:
                reduced_model.set_constraint(reaction, medium.loc[reaction, condition + '_p25'], medium.loc[reaction, condition + '_p75'])

            elif not reaction in medium_elements['exchange'].values:
                reduced_model.set_constraint(reaction, lower_bound=0)

    num_atoms = {r: np.sum(elements_atoms(reduced_model.formulas[reduced_model.get_reactants(r)[0]]).values()) for r in reduced_model.get_exchanges()}
    fitted_solution = min_differences(reduced_model, {r: 0 for r in reduced_model.get_exchanges() if not r in medium.index or r == 'R_EX_o2_e_'}, num_atoms, '')
    fitted_media = fitted_solution.get_net_conversions(reduced_model)

    # Constraint metabolic model lower bound
    for reaction in fitted_media:
        if not reaction in medium.index:
            reduced_model.set_constraint(reaction, lower_bound=fitted_media[reaction])

    # # Reduce metabolic model
    # print model
    # reduced_model = reduce_model(reduced_model)
    # print reduced_model
    #
    # # Store reduced model and fitted medium
    # reduced_models[condition] = reduced_model.deepcopy()
    fitted_solutions[condition] = fitted_solution
    fitted_mediums[condition] = fitted_media

print '[INFO] Metabolomics fitting and model reduction DONE!'

# Plot fitted medium
mediums_df = DataFrame(fitted_mediums)
mediums_df_order = [reduced_model.metabolites[reduced_model.get_reactants(m)[0]] for m in np.abs(mediums_df['UOK262'] - mediums_df['UOK262pFH']).sort(inplace=False, ascending=False).index]
mediums_df['metabolite'] = [reduced_model.metabolites[reduced_model.get_reactants(e)[0]] for e in mediums_df.index]
mediums_df['type'] = ['Measured' if e in medium.index else 'Predicted' for e in mediums_df.index]
mediums_df = melt(mediums_df, id_vars=['metabolite', 'type'])
mediums_df['combination'] = mediums_df['variable'] + ' ' + mediums_df['type']

pallete = sns.color_palette('Paired')
pallete = [pallete[0], pallete[2], pallete[1], pallete[3]]

sns.set(style='white')
sns.set_context('paper')
g = sns.factorplot('metabolite', 'value', data=mediums_df, hue='combination', x_order=mediums_df_order, legend=False, aspect=2, palette=pallete)
sns.despine()
plt.axhline(0, c='#95a5a6', lw=.5, alpha=.15)
plt.ylabel('mol / gDW / h')
plt.legend(loc='lower right')
g.set_xticklabels(rotation=90)
ax = plt.gca()
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
plt.savefig('/Users/emanuel/Projects/projects/pymist/reports/fh_cells/' + 'fitted_media' + '.pdf', bbox_inches='tight')
plt.close('all')

# Sample reduced models
for condition in conditions:
    # Sample
    samples = sample(reduced_models[condition], 1000, 350, verbose=1)
    samples = DataFrame(samples, columns=reduced_model.reactions.keys())

    # Store sampling
    sampling[condition] = samples
    print '[INFO] Sampling finished!'

    # Store reduced model
    reduced_models[condition] = reduced_model

print '[INFO] Sampling DONE!'

# Plot reactions by pathway
reactions = model.get_reactions_by_pathways(['Citric Acid Cycle'])
reactions.extend(['R_PDHm', 'R_ACONT'])
reactions = ['R_PDHm', 'R_ACONT', 'R_SUCD1m', 'R_ICDHy', 'R_GLNS', 'R_GND']

reactions = list(set(reactions).intersection(set(sampling[conditions[0]].columns).intersection(set(sampling[conditions[1]].columns))))

ko_samples = sampling['UOK262'][reactions]
ko_samples['condition'] = 'UOK262'

wt_samples = sampling['UOK262pFH'][reactions]
wt_samples['condition'] = 'UOK262pFH'

df = ko_samples.append(wt_samples)

sns.set(style='white')
sns.pairplot(df, hue='condition', palette=sns.color_palette('Paired'), diag_kind='kde')
plt.show()
plt.savefig('/Users/emanuel/Projects/projects/pymist/reports/fh_cells/' + conditions[0] + '_vs_' +  conditions[1] + '_sampling_histogram' + '.png', bbox_inches='tight')

df = melt(df, id_vars='condition')

sns.set(style='white')
g = sns.factorplot('variable', 'value', 'condition', data=df, kind='box', palette=sns.color_palette('Paired'), legend=False)
g.set_xticklabels(rotation=90)
plt.axhline(0, c='#95a5a6', lw=.5, alpha=.15)
plt.legend(loc='upper right')
plt.ylabel('mol / gDW / h')
plt.savefig('/Users/emanuel/Projects/projects/pymist/reports/fh_cells/' + conditions[0] + '_vs_' +  conditions[1] + '_sampling' + '.pdf', bbox_inches='tight')