import time
import seaborn as sns
import matplotlib.pylab as plt
import numpy as np
from hlrcc import wd
from pymist.reader.sbml_reader import read_sbml_model
from pymist.sampler import sample
from pandas import read_csv, DataFrame
from pymist.simulation import min_differences


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
conditions = ['UOK262', 'UOK262pFH']

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
    samples = sample(reduced_model, 10000, 2500, verbose=1)
    samples = DataFrame(samples, columns=reduced_model.reactions.keys())

    # -- Export sampling
    samples.to_csv('%s/data/%s_sampling.txt' % (wd, condition), sep='\t', index=False)
    print '[INFO] Sampling finished: ', condition
<<<<<<< HEAD
=======

print '[INFO] Sampling DONE!'

# -- Export sampling
[sampling_results[c].to_csv('%s/data/%s_sampling.txt' % (wd, c), sep='\t') for c in conditions]
print '[INFO] Sampling matrices exported'

# # -- Plot
# reactions = ['R_PDHm', 'R_ACONT', 'R_SUCD1m', 'R_ICDHy', 'R_GLNS', 'R_GND']
#
# ko_samples = sampling_results['UOK262'][reactions]
# ko_samples['condition'] = 'UOK262'
#
# wt_samples = sampling_results['UOK262pFH'][reactions]
# wt_samples['condition'] = 'UOK262pFH'
#
# plot_df = ko_samples.append(wt_samples)
#
# sns.set(style='ticks')
# sns.pairplot(plot_df, hue='condition', palette=sns.light_palette('#34495e', 3)[1:], diag_kind='kde')
# plt.savefig('%s/reports/sampling_pairplot.pdf' % wd, bbox_inches='tight')
# plt.close('all')
# print '[INFO] Plot done'
>>>>>>> ec222437a72c155a7c909e47f04a42389ec57662
