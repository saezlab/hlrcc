#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

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
solution = sample(model, n_samples=100)
print solution
