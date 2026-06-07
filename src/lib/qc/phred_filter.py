from src.lib.par_utils import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import numpy as np

def get_avg_phred(rec:SeqRecord) -> float:
    if not isinstance(rec, SeqRecord):
        raise TypeError("Record should be of type SeqRecord")
    return np.mean([p for p in rec.letter_annotations['phred_quality']])

def phred_filter(samp_file, out_dir, thresh = 90, force = False):
    
    hit_file = out_dir / "phred_hits.txt"
    
    print("Calculating average per read Phred scores...")
    if (not hit_file.exists()) or force:
        phreds = [get_avg_phred(rec) for rec in SeqIO.parse(samp_file, "fastq")]
        print(f"Done.")
        hit_inds = np.where([x > thresh for x in phreds])
        print(f"Phred hit indices saved to {str(hit_file)}")
        np.savetxt(hit_file, hit_inds, delimiter="\n", fmt="%d")
    else:
        print("Phreds already calculated so using those.\nUse `force` to run again.")
        hit_inds = np.loadtxt(hit_file)
    
    return hit_inds

