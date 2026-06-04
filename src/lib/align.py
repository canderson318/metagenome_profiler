from difflib import SequenceMatcher
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import os
import numpy as np 
import pandas as pd
from itertools import islice
from  multiprocessing import Pool
from functools import partial

#\\\\
#\\\\
# Load Data
#\\\\
#\\\\
in_dir = Path("in/")
ref_dat_dir = [x for x in os.listdir(in_dir / "data") if "GCF" in x][0]
ref_file = in_dir / "data/" / ref_dat_dir 
ref_file = ref_file / os.listdir(ref_file)[0]

samp_file = in_dir / "SRR12464727.fastq"


#\\\\
#\\\\
# 
#\\\\
#\\\\

count = 0
lens = []
for rec in SeqIO.parse(ref_file, "fasta"):
    count +=1
    L = len(rec.seq)
    lens.append(L)
    # print(rec.id, L) 
print(f"{count} scaffolds.")
# 134
lens_d = pd.DataFrame({"len" : lens})
descr  = lens_d.describe()
outliers = (lens_d > descr.loc['75%']  + ((descr.loc['75%'] - descr.loc['50%'] )* 1.5) )
whr_outlier = np.where(outliers)[0]

#\\\\
#\\\\
# build index from bat genome
#\\\\
#\\\\

# load host genome into memory where large outlier
SEQS = [str(S.seq).upper() 
        for i, S in enumerate(SeqIO.parse(ref_file, "fasta")) 
        if i in whr_outlier
        ]

from src.lib.kmers import extract_kmers
K = 31

# which_SEQS = list(range(len(SEQS)))
which_SEQS = [0]
if __name__ == "__main__":
      with Pool(processes=8) as pool:
          res = pool.map(
              partial(extract_kmers,k = int(31)), 
              [s for i, s in enumerate(SEQS) if i in which_SEQS]
              )
      genome_kmers = set().union(*res)

#\\\\
#\\\\
# classify reads
#\\\\
#\\\\

SEQS[0].__len__()
    
reads = islice(SeqIO.parse(samp_file, 'fastq'), 500)

def is_host(read, threshold=0.2):
    seq = str(read.seq).upper()
    kmers = [seq[i:i+k] for i in range(len(seq) - k)]
    hits = sum(1 for km in kmers if km in genome_kmers)
    return hits / len(kmers) > threshold




# chunk_size = 500

# for i in range(len(reads)-chunk_size+1):
#     for 