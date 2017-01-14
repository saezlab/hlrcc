#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import seaborn as sns
import matplotlib.pyplot as plt
from framed import load_cbmodel, simplify
from hlrcc.metabolomics.sampler import sample

# Import model
model = load_cbmodel('./files/ecoli_core_model.xml', flavor='cobra')
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(model.metabolites), len(model.reactions), len(model.genes))

# Simplify
simplify(model, inplace=True)
print 'Metabolites: %d, Reactions: %d, Genes: %d' % (len(model.metabolites), len(model.reactions), len(model.genes))

# Sampler
solution = sample(model, n_samples=1000, n_steps=200, verbose=1)
print solution

# Plotting
reactions = ['R_ACALD', 'R_MDH', 'R_FUM', 'R_FBA', 'R_FBP', 'R_ADK1']

sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})
sns.pairplot(solution[reactions], diag_kws={'bins': 100}, plot_kws={'s': 5})
plt.gcf().set_size_inches(8, 8)
plt.savefig('./reports/sampler_test_ecoli_core.png', bbox_inches='tight', dpi=200)
plt.close('all')
print '[INFO] Plot done'
