from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO as SeqIO
import pickle as pkl
import itertools as itls
from collections import Counter 
import seaborn as sns
import matplotlib.pyplot as plt
import Levenshtein as LV
from sklearn.manifold import MDS
from tqdm import tqdm
import re
import warnings

from skbio.stats import distance

def result_overview(result, fig_path):
    """
    Examine Kraken2 classification of each read across taxons.
        rank: U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies

    Compares species-level vs. family/genus-level classification counts: if there are more 
    higher-rank (G, F) assignments than species-level ones (S\\d), the reads don't divide well on species 
    Saves a bar plot of rank proportions to fig_path.
    """
    fig_path.mkdir(exist_ok= True)
    
    species_counts = (
        result
        [result['rank'].str.match(r"^S\d*")]
        .groupby("name")
        .size()
        .sort_values(ascending = False)
    )
    high_rank_counts = (
        result
        [result['rank'].str.match(r"G|F")]
        .groupby('name')
        .size()
        .sort_values(ascending = False)
    )

    species_to_highrank_count_ratio = pd.DataFrame(species_counts).sum().div(high_rank_counts.sum()).mean()

    print(f"There are {species_to_highrank_count_ratio:.0f}:1 species:family/genus kraken classifications.")

    if species_to_highrank_count_ratio < 1:
        warnings.warn(
            "Species:family/genus ratio < 1 — no LCA resolved below family/genus rank. "
            "Family/genus may explain read diversity better than species.",
            UserWarning
        )
        
    #\\\
    # plot rank abundance
    #\\\
        
    # order ranks
    ranks = ["R", "K", "P", "C", "O", "F", "G", "S"]
    granular_ranks = [r+str(i) for r in ranks for i in [""] + list(range(1,5))]
    observed_ranks = [r for r in granular_ranks if r in result['rank'].unique()]
    idx_key = {r: i for i,r in enumerate(observed_ranks)}
    rank_prop = (
        result
        .groupby("rank")
        .size().div(len(result))
        .sort_index(key = lambda idx:
            idx.map(idx_key)
            )
    )

    label_key = {"R":"Root", "K":"Kingdom", "P":"Phylum","C":"Class", "O":"Order", "F":"Family", "G":"Genus", "S":"Species"}
    labels = [label_key[R[0]] + R[1] 
            if len(R)>1 else label_key[R[0]]
            for R in observed_ranks 
            ]

    rank_prop.index = labels

    fig,ax = plt.subplots(figsize = (8,8))
    sns.barplot(y = rank_prop.index, x =rank_prop, ax = ax, color = "#5a0d87")
    ax.set_xscale("log")
    ax.set_ylabel("")
    ax.set_xlabel("Log-Abundance")
    ax.set_title(f"Kraken2 Rank Classification Proportions")
    plt.tight_layout()
    plt.savefig(fig_path / "genus_abundance_bar.png", dpi = 150)
    
    
def result_abundances(result,fig_path, percentile = .90):
    """
    Calculate and plot species/genus/family abundance proportions and Shannon diversity.

    rank: U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies

    Saves bar plots of abundances to fig_path and returns the top species
    by proportion (above the given percentile) along with that percentile.
    """
    
    fig_path.mkdir(exist_ok= True)
    
    def get_prop(res: pd.DataFrame, rank = 'S1'):
        if isinstance(rank, re.Pattern):
            filt = result['name'][result['rank'].str.contains(rank)]
        else:
            filt = result['name'][result['rank'] == rank]
        cnts = filt.value_counts()
        return (cnts / len(filt)).sort_values(ascending = False)

    species_prop = get_prop(result, rank = re.compile(r"^S\d"))
    genus_prop = get_prop(result, rank = re.compile(r"^G\d"))
    family_prop = get_prop(result, rank = re.compile(r"^F\d"))

    # –––– Diversity 
        
    # calculate shannon diversity
    def shannon_diversity(P):
        return - np.sum(P * np.log(P))
    
    P = percentile
    percentile = np.exp(np.quantile(np.log(species_prop),P)) # 2.2e-5
    top_species_by_prop = species_prop[(species_prop > percentile )]

    fig,ax = plt.subplots(figsize = (8,8))
    sns.barplot(y = top_species_by_prop.index, x = top_species_by_prop, ax = ax, color = "#5a0d87")
    ax.set_xscale("log")
    ax.set_ylabel("")
    ax.set_xlabel("Log-Abundance")
    ax.set_title(f"Top {P*100:.0f}% classifications\nShannon Diversity (all) = {shannon_diversity(species_prop):.2f}")
    plt.tight_layout()
    plt.savefig(fig_path / "species_abundance_bar.png", dpi = 150)

    fig,ax = plt.subplots(figsize = (8,8))
    sns.barplot(y = genus_prop.index, x =genus_prop, ax = ax, color = "#5a0d87")
    ax.set_xscale("log")
    ax.set_ylabel("")
    ax.set_xlabel("Log-Abundance")
    ax.set_title(f"Genus Classifications\nShannon Diversity = {shannon_diversity(genus_prop):.2f}")
    plt.tight_layout()
    plt.savefig(fig_path / "genus_abundance_bar.png", dpi = 150)

    fig,ax = plt.subplots(figsize = (8,8))
    sns.barplot(y = family_prop.index, x =family_prop, ax = ax, color = "#5a0d87")
    ax.set_xscale("log")
    ax.set_ylabel("")
    ax.set_xlabel("Log-Abundance")
    ax.set_title(f"Family Classifications\nShannon Diversity = {shannon_diversity(family_prop):.2f}")
    plt.tight_layout()
    plt.savefig(fig_path / "family_abundance_bar.png", dpi = 150)
    
    return top_species_by_prop

def read_str_distance_analysis(result, top_species_by_prop, out_dir, fig_path, read_sample_N  = 10_000):
    """
    Analyze pairwise Levenshtein distances between reads of the top-abundance taxa.

    Stratified-samples up to read_sample_N reads (by taxon proportion), and computes pairwise string distances and   runs a PERMANOVA to test for
    separation between taxa on the stratified read sample string distance matrix.
    Then plots the MDS embedding of the distance matrix to fig_path and saves PERMANOVA
    result to out_dir / "permanova_res.csv" 
    
    :param read_sample_N: number of reads to sample stratified on species and calculate string distance between.
    :return: MDS embedding (X).
    """
    reads = result['reads'][result['name'].isin(top_species_by_prop.index)]
    N = read_sample_N # 2min w/ 10_000 for strdist
    
    # # uniform sample
    # samp_inds = pd.Series(reads.index).sample(N, replace = False, random_state = 139)
    
    # stratified sample 
    samp_inds = (
        result.loc[reads.index]
        .groupby('name', group_keys=False)
        .apply(lambda g: g.sample(n=min(len(g), max(1, int(N * len(g) / len(result) ) )), random_state=139, replace = False), include_groups = False)
        .index
    )

    samp_reads = reads.loc[samp_inds]
    n = len(samp_inds)
    str_dist = np.zeros((n, n))

    for i in tqdm(range(len(samp_reads)), desc = "Calculating pairwise distance"):
        a = samp_reads.iloc[i]
        for j, b in enumerate(samp_reads):
            str_dist[i,j] = LV.distance(a, b)

    # \\\\
    # –––– Permanova
    # \\\\

    str_DIST = distance.DistanceMatrix(str_dist)
    perm = distance.permanova(str_DIST, result.loc[samp_inds,"name"].values)
    print(f"*** Permanova ***")      
    for idx in perm.index[1:]:
        print(f"\t{idx} = {perm[idx]}")
    perm_path = out_dir / "permanova_res.csv"
    print(f"Saving permanova result to {str(perm_path)}")
    perm.to_csv(perm_path, header=False)

    # \\\\
    # –––– MDS
    # \\\\
    mds = MDS(dissimilarity = "precomputed",random_state = 139, verbose = 1,n_init=2, max_iter=100, n_jobs=2)

    print("Computing MDS...")
    X = mds.fit_transform(str_dist)
    print("Done.")

    if len(samp_inds) != len(X):
        raise ValueError("Shapes do not match")

    # plot mds
    names = result["name"].loc[samp_inds]
    cats = pd.Categorical(names)
    unique_cats = cats.categories

    plt.figure(figsize = (13,10))
    scatter = plt.scatter(X[:,0], X[:,1], alpha = .5, c = cats.codes, cmap = 'tab10')
    elements = scatter.legend_elements()
    elements = (elements[0],[unique_cats[i] for el in elements[1] for i in [int(re.findall(r"\d+",el)[0])]])
    plt.legend(*elements, title = "Taxon", bbox_to_anchor=(1.0, 1.0) )
    plt.title("MDS plot of read string distances")
    plt.tight_layout()
    plt.savefig(fig_path / f"top_perc_read_str_dist_mds.png",dpi = 150)
    
    return X    

