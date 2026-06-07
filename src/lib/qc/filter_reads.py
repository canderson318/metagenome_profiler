from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from src.lib.qc.align_reads import align
from src.lib.qc.phred_filter import phred_filter

def filter_reads(samp_file, ref_file,out_dir, K = 31,  out_file = "host_filtered.fasta", test = False):
    out_dir = Path(out_dir)
    
    print("Aligning reads to host genome...")
    
    if test:
        aln_hits = set(np.loadtxt(out_dir / "dummy_reads_from_host.inds", dtype=int)) 
    else:
        aln_hits = align(samp_file = samp_file, ref_file=ref_file, K = K, out_dir = out_dir, ncores = 8, test = test)
    
    print("Removing low quality reads...")
    qual_hits = phred_filter(samp_file = samp_file, ref_file=ref_file, K = K, out_dir = out_dir, ncores = 8, test = test)
    
    # \\\\
    # filter out host reads
    # \\\\
    hits_set = set(aln_hits) + set(qual_hits)

    out_file_path = out_dir / out_file
    print(f"Writing non-host reads to {out_file_path}")
    with open(samp_file) as fin, open(out_file_path, 'wt') as fout:
        # for every four lines, check if ith line is a hit, save the header and sequence of that line
        for i, (h1, seq, h2, qual) in enumerate(zip(*[iter(fin)]*4)): 
            if i not in aln_hits_set:
                fout.write(f">{h1[1:].strip()}\n{seq}")  # @ → >, skip qual
    print(f"Done.")
    
    return out_file_path