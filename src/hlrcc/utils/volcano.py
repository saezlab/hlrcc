#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def signif(value):
    if value < 0.01:
        return '*'
    elif value < 0.05:
        return '**'
    else:
        return '-'


def volcano(f, d, x_label, y_label, adj_p_value_label, title='', to_highlight=None, ids=None):
    # Define params
    sns.set(style='ticks', context='paper', font_scale=.75, rc={'axes.linewidth': .3, 'xtick.major.width': .3, 'ytick.major.width': .3})

    # Descritise significance
    d['signif'] = [signif(v) for v in d[adj_p_value_label]]

    # Define pallete
    colour_pallete = sns.light_palette('#34495e', 4, reverse=True)[:-1]
    g = sns.lmplot(x=x_label, y=y_label, data=d, hue='signif', fit_reg=False, palette=colour_pallete, legend=False)
    g.axes[0, 0].set_ylim(0,)

    # Add FDR threshold lines
    plt.text(plt.xlim()[0]*.98, -np.log10(0.01) * 1.01, 'FDR 1%', ha='left', color=colour_pallete[0], alpha=0.65, fontsize=5)
    plt.axhline(-np.log10(0.01), c=colour_pallete[0], ls='--', lw=.7, alpha=.7)

    plt.text(plt.xlim()[0]*.98, -np.log10(0.05) * 1.01, 'FDR 5%', ha='left', color=colour_pallete[1], alpha=0.65, fontsize=5)
    plt.axhline(-np.log10(0.05), c=colour_pallete[1], ls='--', lw=.7, alpha=.7)

    # Add axis lines
    plt.axvline(0, c='#95a5a6', lw=.3, ls='-', alpha=.3)

    # Add axis labels and title
    plt.title(title, fontsize=10, fontname='Arial')
    plt.xlabel('Fold-change (log2)', fontsize=8, fontname='Arial')
    plt.ylabel('FDR (-log10)', fontsize=8, fontname='Arial')

    # Add text to highlighted genes
    if to_highlight is not None:
        for i, r in d.iterrows():
            if r[ids] in to_highlight:
                plt.text((r[x_label] * 1.01), (r[y_label] * 1.01), r[ids], ha='left', alpha=0.75, fontsize=8)

    # Adjust axis lines thickness
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.xaxis.set_tick_params(width=0.5)
    ax.yaxis.set_tick_params(width=0.5)

    # Save plot
    plt.gcf().set_size_inches(3., 5., forward=True)
    plt.savefig(f, bbox_inches='tight')

    print '[INFO] Volcano generated: ' + f
