import time
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
from hlrcc import wd
from pymist.reader.sbml_reader import read_sbml_model
from pymist.sampler import sample, fix_futile_cycles
from pymist.reduce import reduce_model
from pandas import read_csv, DataFrame, Series, melt
from pymist.simulation import min_differences, pFBA, FVA
from pymist.balance import elements_atoms
from pymist.utils.plot_utils import save_plot
from pandas.stats.misc import zscore


# -- Imports
# Metabolic model
model = read_sbml_model('/Users/emanuel/Projects/resources/metabolic_models/recon1.xml')
model.remove_b_metabolites()

# CORE metabolomics (mmol/gDW/h)
core = read_csv('%s/data/uok262_metabolomics_core_processed.txt' % wd, sep='\t', index_col=0)
core.columns = [c.split('_')[0] for c in core]
core = DataFrame({c: core[c].median(1) for c in set(core)})
core = core[core.abs() > .1]
core = core.dropna(how='all')

# Medium
medium = read_csv('%s/files/DMEM_41966_medium.txt' % wd, sep='\t').dropna()

# O2 consumption (mmol/gDW/h)
o2_exch = 'R_EX_o2_e_'
core_o2 = read_csv('%s/data/uok262_metabolomics_core_o2_processed.txt' % wd, sep='\t', index_col=0).to_dict()[o2_exch]


# -- Fit medium
# Set condition
conditions = ['UOK262', 'UOK262pFH']

sampling_results = {}
for condition in conditions:
    # condition = conditions[0]
    print '[INFO] Condition: ', condition

    # -- Set model constraint variables
    exchange_rates = core[condition].dropna()

    # metabolite = exchange_rates.index[0]
    # Reactions not used
    for reaction in ['R_EX_retpalm', 'R_EX_retpalm_e_']:
        model.set_constraint(reaction, 0)

    # Deepcopy metabolic model
    reduced_model = model.deepcopy()

    # Add FH KO to tumour cell lines
    if condition == 'UOK262':
        reduced_model.set_constraint('R_FUM', 0, 0)
        reduced_model.set_constraint('R_FUMm', 0, 0)
        print '[INFO] R_FUM blocked'

    # -- Fit CORE and medium
    # Block exhange fluxes not present in the medium or CORE experiment
    for reaction in reduced_model.get_exchanges(check_matrix=True):
        if (reaction not in exchange_rates.index) and (reaction != o2_exch) and (reaction not in medium['exchange'].values):
            reduced_model.set_constraint(reaction, lower_bound=0)
    print '[INFO] Exchange reactions constrained'

    # Minimise differences between simulations and measurements
    constrains = {}
    for r in reduced_model.get_exchanges(check_matrix=True):
        if r in exchange_rates.index:
            constrains[r] = exchange_rates.ix[r, condition]
        elif r == 'R_EX_o2_e_':
            constrains['R_EX_o2_e_'] = core_o2[condition]

    fitted_medium = min_differences(reduced_model, constrains).get_net_conversions(reduced_model, check_matrix=True)

    # Constraint metabolic model lower bound
    for reaction in fitted_medium:
        if reaction in exchange_rates.index:
            reduced_model.set_constraint(reaction, lower_bound=fitted_medium[reaction], upper_bound=fitted_medium[reaction])

    # -- Metabolic sampling
    samples = sample(reduced_model, 2000, 350, verbose=1)
    samples = DataFrame(samples, columns=reduced_model.reactions.keys())

    # Store sampling
    sampling_results[condition] = samples
    print '[INFO] Sampling finished: ', condition

print '[INFO] Sampling DONE!'


# -- Plot
reactions = ['R_PDHm', 'R_ACONT', 'R_SUCD1m', 'R_ICDHy', 'R_GLNS', 'R_GND']

ko_samples = sampling_results['UOK262'][reactions]
ko_samples['condition'] = 'UOK262'

wt_samples = sampling_results['UOK262pFH'][reactions]
wt_samples['condition'] = 'UOK262pFH'

df = ko_samples.append(wt_samples)

sns.set(style='ticks')
sns.pairplot(df, hue='condition', palette=sns.color_palette('Paired'), diag_kind='kde')
plt.savefig('%s/reports/sampling_pairplot.pdf' % wd, bbox_inches='tight')
plt.close('all')
print '[INFO] Plot done'
