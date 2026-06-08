from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import numpy as np
from src.lib.qc.kmers import  get_num_reads
from numpy.typing import NDArray
from tqdm import tqdm

def iterate_phred_lines(samp_file):
    phreds = []
    with open(samp_file) as f:
        while True:
            f.readline() # h1
            f.readline() # seq 
            f.readline() # h2
            phred = f.readline()
            if not phred:
                break
            yield phred # phred

# Setting up the ASCII Phred score conversion table
phred_table = {}
for i in range(33, 127):
    phred_table[chr(i)] = i - 33
    
def get_avg_phred(samp_file)-> NDArray[float]:
    phreds = []
    ITER = iterate_phred_lines(samp_file)
    for _ in tqdm(range(get_num_reads(samp_file)), desc = "Calculating average Phred per read", unit = 'line'):
        line = next(ITER)
        phreds.append(np.mean([ phred_table[p] for p in line.strip() ]))
    return np.array( phreds, dtype = np.float32)

def phred_filter(samp_file, out_dir, thresh = 20, force = False) -> np.array:
    
    hit_file = out_dir / "phred_hits.txt"
    phred_file = out_dir / "phreds.txt"
    
    if not hit_file.exists() or force:
        phreds = get_avg_phred(samp_file)
        print(f"Done.")
        np.savetxt(phred_file,phreds, delimiter="\n", fmt="%d")
        print(f"Phreds saved to {phred_file}")
        
        np.savetxt(hit_file, hit_inds, delimiter="\n", fmt="%d")
        print(f"Phred hit indices (p < {thresh}) saved to {str(hit_file)}")
    else:
        print("Phreds already calculated so using those.\nUse `force` to run again.")
        phreds = np.loadtxt(phred_file)
        hit_inds = np.where([x < thresh for x in phreds])[0]
    
    return hit_inds

