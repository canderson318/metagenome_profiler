from Bio import SeqIO
from pathlib import Path
import numpy as np
import pickle as pkl
from multiprocessing import Pool
from functools import partial
from src.lib.qc.kmers import extract_kmers

def outlier_scaffold_indices(ref_file):
    lens = [len(rec) for rec in SeqIO.parse(ref_file, "fasta")]
    print(f"{len(lens)} scaffolds.")
    lens = np.array(lens)
    q1, q3 = np.percentile(lens, [25, 75])
    threshold = q3 + 1.5 * (q3 - q1)
    indices = set(np.where(lens > threshold)[0])
    print(f"{len(indices)} main scaffolds (chromosomes).")
    return indices


def index_host(K, ref_file, which_scaffolds, force = False):
    """
    compare K-length read mers to K-length reference genome mers. 
    :param which_scaffolds: fasta headers of scaffolds to align against as set
    :return: path to index and the reference genome kmers
    """

    out_file = Path("out/") / f"host_index_k{K}.pkl"

    if (not out_file.exists()) or force:
        print("Indexing host genome...")
        seqs = [
            str(rec.seq).upper()
            for rec in SeqIO.parse(ref_file, "fasta")
            if rec.id in which_scaffolds
        ]
        
        seqs_len = len(seqs)
        if  seqs_len == 0:
            raise ValueError("Scaffold sequences not found.")
        
        if seqs_len > 1:
            with Pool(processes=8) as pool:
                res = pool.map(partial(extract_kmers, k=K), seqs)
            genome_kmers = set().union(*res)
        else:
            genome_kmers = extract_kmers(seqs[0], k=K)

        with open(out_file, "wb") as f:
            pkl.dump(genome_kmers, f)

        print("Done.")
        print(f"Host index saved to {str(out_file)}")

    else:
        with open(out_file, "rb") as f:
            genome_kmers = pkl.load(f)

    return str(out_file), genome_kmers