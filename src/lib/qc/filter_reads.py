from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from src.lib.qc.align_reads import align
from src.lib.qc.phred_filter import phred_filter

def filter_reads(samp_file, ref_file,out_dir, K = 31,  out_file = "host_filtered.fasta", phred_thresh = 20, aln_thresh = .2,force_index = False,force_aln = False,force_phred = False, which_host_scaffolds = [0] ):
    out_dir = Path(out_dir)
    
    # \\\\
    # Align reads to host index: indexes host to K-length kmers and compares these to read K-mers
    # \\\\
    aln_hits = align(samp_file = samp_file, ref_file=ref_file, K = K,thresh=aln_thresh, out_dir = out_dir, ncores = 8, force_index = force_index,force_aln = force_aln,which_scaffolds=which_host_scaffolds)
    
    # \\\\
    # Calculate Phred scores and record where lower than threshold (default 20)
    # \\\\
    print("Removing low quality reads...")
    qual_hits = phred_filter(samp_file, out_dir, force = force_phred)
    
    # \\\\
    # filter out host reads
    # \\\\
    hits_set = set(aln_hits) | set(qual_hits)
    type(qual_hits)
    out_file_path = out_dir / out_file
    print(f"Writing non-host reads to {out_file_path}")
    with open(samp_file) as fin , open(out_file_path, 'wt') as fout:
        # for every fourth line, check if hit, save the header and sequence of that line
        i = 0
        while True:
        # for i, (h1, seq, h2, qual) in enumerate(zip(*[iter(fin)]*4)): 
            h1 = fin.readline() # h1
            read = fin.readline() # read 
            fin.readline(); fin.readline() # h2, phred
            if not h1: 
                break
            if i not in hits_set:
                fout.write(f">{h1[1:]}{read}")  # @ → >, skip qual
            i+=1
            
    print(f"Done.")
    
    return out_file_path