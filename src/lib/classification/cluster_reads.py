from pathlib import Path
from Bio import SeqIO
from scipy.spatial.distance import pdist
import Levenshtein as LV
from rapidfuzz import process
import numpy as np
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
from src.lib.qc.kmers import get_num_reads
from sklearn.cluster import KMeans

in_dir = Path("in/")
out_dir = Path("out/")


filtered_read_path = out_dir / "host_filtered.fasta"

rng = np.random.default_rng(123)
inds = rng.choice(get_num_reads(filtered_read_path,'fasta'), 500, replace = False)
all_reads = [str(rec.seq)
         for i,rec in tqdm(enumerate(SeqIO.parse(filtered_read_path, 'fasta'))) 
         if i in inds
         ]
reads = all_reads[:100]


def H(string):
    prob = [ float(string.count(c)) / len(string) for c in list(string) ]
    # calculate the entropy
    return - sum([ p * np.log(p) / np.log(2.0) for p in prob ])

def euclid(x,y):
    return np.sqrt((x-y)**2)


n = len(reads)
str_dist = np.zeros((n, n))
entropy_dist = str_dist.copy()

for i, a in enumerate(reads):
    for j, b in enumerate(reads):
        str_dist[i,j] = LV.distance(a, b)
        entropy_dist[i,j] = euclid(H(a),H(b))

cumm_dist = entropy_dist  + str_dist 
mds = MDS(dissimilarity = "precomputed",random_state = 139)
X = mds.fit_transform(cumm_dist)
# X = mds.fit_transform(str_dist)
plt.scatter(X[:,0], X[:,1], alpha = .5)


KM = KMeans()
KM.fit(X)
clust = KM.predict(X)
scatter = plt.scatter(X[:,0], X[:,1], alpha = .5, c = clust, cmap = "tab10", label = 'cluster')
plt.legend(*scatter.legend_elements(), title="Cluster")

{c:len(clust[clust == c]) for c in np.unique(clust)}

[str(reads[i]) for i in np.where(clust == 7)[0]]