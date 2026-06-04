from Bio import SeqIO
import os 
from pathlib import Path
import numpy as np
import seaborn as sns

wd = Path.home() / "Documents/school/prelims/day3/"

# in_dir = wd / "in"
# paths = sorted([in_dir / x for x in os.listdir(in_dir)])

reads = []
with SeqIO.parse(Path('/Users/canderson/Documents/school/prelims/day3/in/SRR12464727.fastq'), "fastq") as iterator:
    for i in range(5_000):
       reads.append(next(iterator)) 
       
means = np.array([np.mean(x.letter_annotations['phred_quality']) for x in reads])

sns.kdeplot(means)

