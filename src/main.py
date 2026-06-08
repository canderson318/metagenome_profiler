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
from sklearn.cluster import KMeans
from tqdm import tqdm
import re
import warnings
from skbio.stats import distance

from src.lib.qc.filter_reads import filter_reads
from src.lib.assembly.DBG import DBG
from src.lib.classification.classify_reads import classify_reads
from src.lib.qc.kmers import  get_num_reads

in_dir = Path("in/")
out_dir = Path("out/")

fig_path = out_dir / 'figs'
fig_path.mkdir(exist_ok= True)

# Data files
ref_file = in_dir / "bat_genome/GCF_004115265.2/GCF_004115265.2_mRhiFer1_v1.p_genomic.fna"
samp_file = in_dir / "SRR12464727.fastq"

# \\\\
# \\\\
# filter reads
# \\\\
# \\\\

filtered_read_path = filter_reads(samp_file, ref_file, out_dir = out_dir, K = 31,
                                  out_file = "filtered.fasta", 
                                  aln_thresh= .5,
                                  force_phred = False,
                                  force_aln = False,
                                #   which_host_scaffolds = [0]
                                  which_host_scaffolds = "all"
                                  )
# ^^out/host_filtered.fasta

# \\\\
# \\\\
# assemble contigs
# \\\\
# \\\\

# G = DBG(SeqIO.parse(filtered_read_path, 'fasta'), K = 85)
# dbg_dir = out_dir / "DBG"
# dbg_dir.mkdir(exist_ok=True)

# with open(dbg_dir / "dbg.pkl", 'wb') as fout: 
#     pkl.dump(G,fout)
    
# contigs = G.make_contigs(True, join_sub_contigs=True)
# with open(dbg_dir / "dbg_contigs.pkl", 'wb') as fout: 
#     pkl.dump(G,fout)

# with open(out_dir / "contigs.fasta" , 'wt') as fout:
#     for i,seq in enumerate(contigs):
#         fout.write(f">contig_{i}\n{seq}\n")

# covs = np.array([ sum(c)/len(c) for s,n,c in contigs])
# lens = np.array([len(s) for s,n,c in contigs])
# import numpy as np
# import seaborn as sns
# sns.histplot((covs/covs.max()) * 100, kde = True)
# sns.histplot((lens/lens.max()) * 100, kde = True)
# p = 70
# len_perc = np.percentile(lens, p)
# cov_perc = np.percentile(covs, p)
# top_contigs = [s for s,n,c in contigs if (sum(c)/len(c) > cov_perc) and (len(s) > len_perc)]


# \\\\
# \\\\
# reference viral genome
# \\\\
# \\\\

# filtered_read_path = out_dir / "filtered.fasta"
viral_db_path = in_dir / "kraken_2_viral_db"
result_dir = out_dir / "kraken"

# load kraken query results
results  = classify_reads(filtered_read_path, viral_db_path,result_dir, 
                          force = False,
                          ncores = 5)
# results.set_index("readid",inplace = True, drop =  False)

cnt_reads = len(results)

# filter out unclassifieds
results = results[results['classified']== 'C'].copy()
print(f"{len(results):,} reads classified of {cnt_reads:,} total reads")

# add reads to each 
reads = {}
target_ids = set(results['readid'])
with open(filtered_read_path) as fin:
    while True:
        h, read = fin.readline().lstrip(">").split(" ")[0], fin.readline().strip()
        if not h:
            break
        if h in target_ids:
            reads.update({h:read})
results['reads'] = results['readid'].map(reads)
del reads

# \\\\
# \\\\
# –––– Analyze results
# \\\\
# \\\\

# \\\\
# –––– Look at ranks
# \\\\
    
# look at how kraken classified
### rank: U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies
## if there are more higher rank (G,F) assignments than species level ones (S\d), 
## then the reads don't divide well on species, 
## otherwise species explains diversity well

species_counts = (
    results
    [results['rank'].str.match(r"^S\d*")]
    .groupby("name")
    .size()
    .sort_values(ascending = False)
)
high_rank_counts = (
    results
    [results['rank'].str.match(r"G|F")]
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

# order ranks
ranks = ["R", "K", "P", "C", "O", "F", "G", "S"]
granular_ranks = [r+str(i) for r in ranks for i in [""] + list(range(1,5))]
observed_ranks = [r for r in granular_ranks if r in results['rank'].unique()]
idx_key = {r: i for i,r in enumerate(observed_ranks)}
rank_prop = (
    results
    .groupby("rank")
    .size().div(len(results))
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
sns.barplot(y = rank_prop.index, x =rank_prop, ax = ax)
ax.set_xscale("log")
ax.set_ylabel("")
ax.set_xlabel("Abundance")
ax.set_title(f"Kraken2 Rank Classification Proportions")
plt.tight_layout()
plt.savefig(fig_path / "genus_abundance_bar.png", dpi = 150)

# \\\\
# –––– Look at abundances
# \\\\

# calculate abundance
### rank: U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies
def get_prop(res: pd.DataFrame, rank = 'S1'):
    if isinstance(rank, re.Pattern):
        filt = results['name'][results['rank'].str.contains(rank)]
    else:
        filt = results['name'][results['rank'] == rank]
    cnts = filt.value_counts()
    return (cnts / len(filt)).sort_values(ascending = False)

species_prop = get_prop(results, rank = re.compile(r"^S\d"))
genus_prop = get_prop(results, rank = re.compile(r"^G\d"))
family_prop = get_prop(results, rank = re.compile(r"^F\d"))

# –––– Diversity 
    
# calculate shannon diversity
def shannon_diversity(P):
    return - np.sum(P * np.log(P))

species_shann = shannon_diversity(species_prop)

P = .90
percentile = np.exp(np.quantile(np.log(species_prop),P)) # 2.2e-5
top_species_prop = species_prop[(species_prop > percentile )]

fig,ax = plt.subplots(figsize = (8,8))
sns.barplot(y = top_species_prop.index, x = top_species_prop, ax = ax)
ax.set_xscale("log")
ax.set_ylabel("")
ax.set_xlabel("Abundance")
ax.set_title(f"Top {P*100:.0f}% classifications")
plt.tight_layout()
plt.savefig(fig_path / "species_abundance_bar.png", dpi = 150)

fig,ax = plt.subplots(figsize = (8,8))
sns.barplot(y = genus_prop.index, x =genus_prop, ax = ax)
ax.set_xscale("log")
ax.set_ylabel("")
ax.set_xlabel("Abundance")
ax.set_title(f"Genus Classifications")
plt.tight_layout()
plt.savefig(fig_path / "genus_abundance_bar.png", dpi = 150)

fig,ax = plt.subplots(figsize = (8,8))
sns.barplot(y = family_prop.index, x =family_prop, ax = ax)
ax.set_xscale("log")
ax.set_ylabel("")
ax.set_xlabel("Abundance")
ax.set_title(f"Family Classifications")
plt.tight_layout()
plt.savefig(fig_path / "family_abundance_bar.png", dpi = 150)


# \\\\
# \\\\
# ––– String distance
# \\\\
# \\\\
def H(string):
    prob = [ float(string.count(c)) / len(string) for c in list(string) ]
    # calculate the entropy
    return - sum([ p * np.log(p) / np.log(2.0) for p in prob ])

def euclid(x,y):
    return np.sqrt((x-y)**2)

reads = results['reads'][results['name'].isin(top_species_prop.index)]
np.random.seed(1034)
N = 10_000 # 2min w/ 10_000 for strdist
# N = 1000

# # uniform sample
# samp_inds = pd.Series(reads.index).sample(N, replace = False, random_state = 139)
# stratified sample 
samp_inds = (
      results.loc[reads.index]
      .groupby('name', group_keys=False)
      .apply(lambda g: g.sample(n=min(len(g), max(1, int(N * len(g) / len(results) ) )), random_state=139, replace = False), include_groups = False)
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

# samp_inds.shape
# str_dist.shape

str_DIST = distance.DistanceMatrix(str_dist)
perm = distance.permanova(str_DIST, results.loc[samp_inds,"name"].values)
perm.to_csv(out_dir / "permanova_res.csv", header=False)

print(f"*** Permanova ***")      
for idx in perm.index[1:]:
    print(f"{idx} = {perm[idx]}")
dir(perm)


# \\\\
# –––– MDS
# \\\\

mds = MDS(dissimilarity = "precomputed",
          random_state = 139, 
          verbose = 2,
          n_init=2, 
          max_iter=100, 
          n_jobs=2)

print("Computing MDS...")
X = mds.fit_transform(str_dist)
print("Done.")

if len(samp_inds) != len(X):
    raise ValueError("Shapes do not match")

# plot mds
names = results["name"].loc[samp_inds]
cats = pd.Categorical(names)
unique_cats = cats.categories

plt.figure(figsize = (10,10))
scatter = plt.scatter(X[:,0], X[:,1], alpha = .5, c = cats.codes, cmap = 'tab10')
elements = scatter.legend_elements()
elements = (elements[0],
                 [unique_cats[i] for el in elements[1] for i in [int(re.findall(r"\d+",el)[0])]]
)
plt.legend(*elements, title = "Taxon")
plt.title("MDS plot of read string distances")
plt.tight_layout()
plt.savefig(fig_path / f"top{P*100}perc_read_str_dist_mds.png",dpi = 150)


