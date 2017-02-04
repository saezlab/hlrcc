#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import seaborn as sns
import matplotlib.pyplot as plt
from framed import load_cbmodel, simplify
from hlrcc.metabolomics.sampler import sample, fix_futile_cycles

# Import model
model = load_cbmodel('./files/ecoli_core_model.xml', flavor='cobra')
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(model.metabolites), len(model.reactions), len(model.genes))

# Simplify
simplify(model, inplace=True)
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(model.metabolites), len(model.reactions), len(model.genes))

# Environmental conditions
env = {'R_EX_glc_e': (-5, 0)}

# Sampler
samples = sample(model, n_samples=1000, n_steps=250, n_steps_projection=50, constraints=env, verbose=1)
print samples

# Fix futile cycles
samples = fix_futile_cycles(model, samples)

# Plotting
reactions = ['R_EX_glc_e', 'R_PGI', 'R_FBA', 'R_FBP', 'R_FUM', 'R_SUCDi', 'R_FRD7']

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.pairplot(samples[reactions], diag_kws={'bins': 100}, plot_kws={'s': 5})
plt.gcf().set_size_inches(8, 8)
plt.savefig('./reports/sampler_test_ecoli_core.png', bbox_inches='tight', dpi=200)
plt.close('all')
print '[INFO] Plot done'
