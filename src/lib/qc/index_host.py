from Bio import SeqIO
from pathlib import Path
import numpy as np
import pickle as pkl
from multiprocessing import Pool
from functools import partial
from src.lib.qc.kmers import extract_kmers


def _outlier_scaffold_indices(ref_file):
    lens = [len(rec) for rec in SeqIO.parse(ref_file, "fasta")]
    print(f"{len(lens)} scaffolds.")
    lens = np.array(lens)
    q1, q3 = np.percentile(lens, [25, 75])
    threshold = q3 + 1.5 * (q3 - q1)
    indices = set(np.where(lens > threshold)[0])
    print(f"{len(indices)} main scaffolds (chromosomes).")
    return indices


def index_host(K, ref_file, which_scaffolds="all", filter_scaffolds=True, force = False):

    if which_scaffolds == "all":
        out_file = Path("out/") / f"host_index_k{K}_all.pkl"
    else:
        out_file = Path("out/") / f"host_index_k{K}.pkl"

    if (not out_file.exists()) or force:

        if not filter_scaffolds:
            raise RuntimeWarning("Running on all scaffolds will take much too long.")

        outlier_indices = _outlier_scaffold_indices(ref_file)

        seqs = [
            str(rec.seq).upper()
            for i, rec in enumerate(SeqIO.parse(ref_file, "fasta"))
            if i in outlier_indices
        ]

        if which_scaffolds != "all":
            seqs = [seqs[i] for i in which_scaffolds]

        if len(seqs) > 1:
            with Pool(processes=8) as pool:
                res = pool.map(partial(extract_kmers, k=K), seqs)
            genome_kmers = set().union(*res)
        else:
            genome_kmers = extract_kmers(seqs[0], k=K)

        with open(out_file, "wb") as f:
            pkl.dump(genome_kmers, f)

        print(f"Host index saved to {str(out_file)}")

    else:
        with open(out_file, "rb") as f:
            genome_kmers = pkl.load(f)

    return str(out_file), genome_kmers