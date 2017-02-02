#!/usr/bin/env python
# Copyright (C) 2017 Emanuel Goncalves

import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pandas import read_csv

uniprot_fasta = './files/uniprot_sprot.fasta'


# --
def read_gmt(file_path):
    with open(file_path) as f:
        signatures = {l.split('\t')[0]: set(l.strip().split('\t')[2:]) for l in f.readlines()}
    return signatures


def get_complexes_dict(organism='Human', corum_file='./files/allComplexesCore.csv'):
    # Import
    complexes = read_csv(corum_file, sep=';').dropna(subset=['Complex name'])

    # Select by organism
    complexes = complexes[complexes['organism'] == organism]

    # Build dict
    complexes = {c: {x.replace('(', '').replace(')', '') for p in complexes.loc[complexes['Complex id'] == c, 'subunits (UniProt IDs)'] for x in p.split(',')} for c in complexes['Complex id']}

    # Map uniprot to genename
    with open('./files/uniprot_to_genename.pickle', 'rb') as handle:
        gmap = pickle.load(handle)

    complexes = {c: {gmap[g][0] for g in complexes[c] if g in gmap and len(gmap[g]) == 1} for c in complexes}

    return complexes


def get_complexes_name(organism='Human', corum_file='./files/allComplexesCore.csv'):
    complexes = read_csv(corum_file, sep=';').dropna(subset=['Complex name'])

    complexes = complexes[complexes['organism'] == organism]

    complexes_name = complexes.groupby('Complex id')['Complex name'].first().to_dict()

    return complexes_name


def get_ktargets(mapfile='./files/Kinase_Substrate_Dataset.txt', organism='human', in_vivo_only=False):
    # Import
    t = read_csv(mapfile, sep='\t')

    # Select organism
    t = t[(t['KIN_ORGANISM'] == organism) & (t['SUB_ORGANISM'] == organism)]

    # P-site id
    t['PSITE'] = t['SUB_GENE'] + '_' + t['SUB_MOD_RSD']

    # Select only in vitro associations
    if in_vivo_only:
        t = t[t['IN_VIVO_RXN'] == 'X']

    # Build dict
    t = t.groupby('GENE')['PSITE'].agg(lambda x: set(x)).to_dict()

    return t


def get_ktargets_omnipath(mapfile='./files/omnipath_ptms_all.txt', ref=['Signor', 'PhosphoSite']):
    # ref = ['Li2012', 'MIMP', 'PhosphoNetworks', 'dbPTM', 'PhosphoSite', 'DEPOD', 'phosphoELM', 'Signor', 'HPRD']

    # Import
    t = read_csv(mapfile, sep='\t')

    # Select data-bases
    t = t[[len(set(i.split(';')).intersection(ref)) != 0 for i in t['sources']]]

    # Map uniprot to genename
    with open('./files/uniprot_to_genename.pickle', 'rb') as handle:
        gmap = pickle.load(handle)

    t = t[[(k in gmap) and (s in gmap) for k, s in t[['enzyme', 'substrate']].values]]
    t = t[[(len(gmap[k]) == 1) and (len(gmap[s]) == 1) for k, s in t[['enzyme', 'substrate']].values]]
    t['enzyme'] = [gmap[i][0] for i in t['enzyme']]
    t['substrate'] = [gmap[i][0] for i in t['substrate']]

    # P-site id
    t['psite'] = t['substrate'] + '_' + t['residue_type'] + t['residue_offset'].astype(str)

    # Build dict
    t = t.groupby('enzyme')['psite'].agg(lambda x: set(x)).to_dict()

    return t


# -- Plotting
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


# -- Sequence and mapping
def read_fasta(fasta_file=None, os='Homo sapiens'):
    fasta_file = uniprot_fasta if fasta_file is None else fasta_file

    sequences = {}

    with open(fasta_file) as f:
        lines = f.readlines()

        for i in range(len(lines)):
            if lines[i].startswith('>sp') and (('OS='+os) in lines[i]):
                uniprot = lines[i].split('|')[1].strip()
                sequence = ''

                i += 1
                while (i < len(lines)) and (not lines[i].startswith('>sp')):
                    sequence += lines[i].strip()
                    i += 1

                sequences[uniprot] = sequence

    return sequences


def match_sequence(sequences, sequence):
    return [k for k, v in sequences.items() if sequence in v]


def read_uniprot_genename(fasta_file=None, os='Homo sapiens'):
    fasta_file = uniprot_fasta if fasta_file is None else fasta_file

    uniprot2genename = {}

    with open(fasta_file) as f:
        lines = f.readlines()

        for i in range(len(lines)):
            if lines[i].startswith('>sp') and (('OS='+os) in lines[i]) and ('GN=' in lines[i]):
                uniprot = lines[i].split('|')[1].strip()
                genename = lines[i].split(' GN=')[1].split(' ')[0]
                accname = lines[i].split('|')[2].strip().split(' ')[0].strip()
                uniprot2genename[uniprot] = (genename, accname)

    return uniprot2genename


# -- ID conversion maps
def build_uniprot_genesymbol_map(mapfile='./files/HUMAN_9606_idmapping.txt'):
    # Build map
    gmap = read_csv(mapfile, sep='\t')
    gmap = gmap[gmap['database'] == 'Gene_Name']
    gmap = gmap.groupby('uniprot')['value'].agg(lambda x: list(x)).to_dict()

    # Store as a pickle object
    with open('./files/uniprot_to_genename.pickle', 'wb') as handle:
        pickle.dump(gmap, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return gmap
