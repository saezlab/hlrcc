#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import matplotlib.pyplot as plt
from pandas import read_csv, DataFrame, Series


# -- Import
transcriptomics = Series.from_csv('./data/UOK262_rnaseq_preprocessed.csv')

proteomics = Series.from_csv('./data/uok262_proteomics_tmt_preprocessed.csv')

phosphoproteomics = Series.from_csv('./data/uok262_phosphoproteomics_tmt_preprocessed.csv')

fluxomics = read_csv('./data/pfba_atp.csv', index_col=0)

