from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import os
import numpy as np 
from itertools import islice
import pickle as pkl
import multiprocessing
multiprocessing.set_start_method('fork', force=True)
from  multiprocessing import Pool

from src.lib.qc.kmers import extract_kmers, get_num_reads
from src.lib.qc.index_host import index_host

def host_match(read:SeqIO.SeqRecord,ref_kmers: set, K = None) -> float:
    """
    Returns average kmer hit against genome reference kmers. 
    :param k: size of kmers
    :param read: read to align to reference
    :return: proportion of a read's kmers that match reference kmers
    """
    
    ref_k = len(next(iter(ref_kmers)))

    if K is None:
        K = ref_k
    elif ref_k != K:
        raise ValueError("Reference Kmer size does not match specified K.")
    seq = str(read.seq).upper()
    kmers = extract_kmers(seq, K)
    if not kmers:
        return 0.0
    hits = len(kmers & ref_kmers)
    return hits / len(kmers)

# host_match(SeqIO.SeqRecord(Seq("ALKDJFLSKDJFSLDKALKDJFLSKDJFSLDKFJSALKDJFLSKDJFSLDKFJSALKDJFLSKDJFSLDKFJSALKDJFLSKDJFSLDKFJSALKDJFLSKDJFSLDKFJSALKDJFLSKDJFSLDKFJSFJS")), ref_kmers)
# host_match(SeqIO.SeqRecord(Seq("TCATCTATTTGGGTCTTTAAACTGTATTTAA"*10)), ref_kmers)

def find_chunk_offsets(filepath, n_chunks, N_reads):
    chunk_size = N_reads // n_chunks
    # byte offset for every chunk's first record
    target_lines = [i * chunk_size * 4 for i in range(n_chunks)]

    offsets = []
    with open(filepath, 'rb') as f:
        line_num = 0
        target_idx = 0
        while target_idx < len(target_lines):
            if line_num == target_lines[target_idx]:
                offsets.append(f.tell())
                target_idx += 1
            f.readline()
            line_num += 1
    return offsets, chunk_size

def init_worker(kmers):
    global _ref_kmers
    _ref_kmers = kmers

def process_chunk(args):
    samp_file, chunk_size, offset, K = args
    res = []
    with open(samp_file) as f:
        f.seek(offset)
        for read in islice(SeqIO.parse(f, "fastq"), chunk_size):
            res.append(host_match(read, _ref_kmers, K=K))
    return res


def align_reads(samp_file, ref_index_file, K,thresh, n_chunks = 16,ncores = 8):
    """
    Align each read to reference by comparing each read's kmers to reference kmers
    """
    ref_kmers = pkl.load(open(ref_index_file, 'rb'))

    # make chunks of indices to parallelize over
    N_reads = get_num_reads(samp_file)
    print(f"\tFinding {n_chunks} chunk offsets...")
    offsets, chunk_size = find_chunk_offsets(samp_file, n_chunks, N_reads)
    chunk_args = [(samp_file, chunk_size, offset, K) for offset in offsets]
    print("Done.")

    print(f"\tProcessing {len(chunk_args)} chunks...")
    with Pool(processes=ncores, initializer=init_worker, initargs=(ref_kmers,)) as pool:
        res = pool.map(process_chunk, chunk_args)

    print(f"Done.")
    alignment = [x for chunk in res for x in chunk]
    hit_inds = np.where([x > thresh for x in alignment ])[0]
    
    if not hit_inds.size:
        return np.array([])
    else:
        return hit_inds
    
def align(samp_file, ref_file, K, out_dir,thresh = .2,ncores = 8, force_aln = False,force_index = False, which_scaffolds = None):
    """
    Align a samples reads against a reference genome using kmer comparison
    
    :param thresh: kmer string similarity threshold. reads more similar than thresh are considered hits against the reference
    :param K: reads/reference genome are cut up into K-length strings
    :param force_aln: should alignment run or should it grab a previously saved one
    :param force_index: should the reference be indexed run or should it grab a previously saved one
    :param which_scaffolds: a set or list of fasta type headers for which scaffold of the reference genome is used in indexing->alignment
    
    :return: array of read indices where alignment with genome was greather than thresh
    """
    if which_scaffolds is None:
        raise ValueError("Please provide reference genome scaffolds to `which_scaffolds` in the format of a set of fasta headers (e.g. '>XXX_1234.1').")
    
    ref_index_file_path, _ = index_host(K=K, ref_file=ref_file, which_scaffolds=which_scaffolds , force=force_index)

    index_file = out_dir / "reads_from_host.inds"
    
    if (not index_file.exists()) or force_aln:
        print("Aligning reads to host genome...")
        hit_inds = align_reads(samp_file, ref_index_file_path, K=K, thresh = thresh, n_chunks = ncores * 2, ncores = ncores) 
        print(f"Hit indices (similarity > {thresh}) saved to {str(index_file)}")
        np.savetxt(index_file, hit_inds, delimiter="\n", fmt="%d")
    else:
        print("Reads already aligned to host so using those.\nUse `force` to align again.")
        hit_inds = np.loadtxt(index_file)
    
    if hit_inds.size == 0:
        print("No host reads found.")
        return np.array([])
    return hit_inds

        