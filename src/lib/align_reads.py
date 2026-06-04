from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import os
import numpy as np 
from itertools import islice
from  multiprocessing import Pool
import pickle as pkl
from src.lib.kmers import extract_kmers

def host_match(read:SeqIO.SeqRecord,ref_kmers: set, K = None):
    """
    Returns average kmer hit against genome reference kmers. 
    :param k: size of kmers
    :param read: read to align to reference
    """
    ref_k = len(next(iter(ref_kmers)))

    if K is None:
        K = ref_k
    elif ref_k != K:
        raise ValueError("Reference Kmer size does not match specified K.")
    seq = str(read.seq).upper()
    kmers = extract_kmers(seq, K)
    hits = len(set(kmers) & ref_kmers)
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
    results = []
    with open(samp_file) as f:
        f.seek(offset)
        for read in islice(SeqIO.parse(f, "fastq"), chunk_size):
            results.append(host_match(read, _ref_kmers, K=K))
    return results


def align_reads(samp_file, ref_index_file, K, test = False):
    """
    Align each read to reference by comparing each read's kmers to reference kmers
    """
    ref_kmers = pkl.load(open(ref_index_file, 'rb'))

    # make chunks of indices to parallelize over
    print("Counting reads...")
    # N_reads = sum(1 for _ in SeqIO.parse(samp_file, 'fastq')) # 39_257_492
    N_reads = 39_257_492
    n_cores = 8
    n_chunks = n_cores * 2
    offsets, chunk_size = find_chunk_offsets(samp_file, n_cores * 2, N_reads)
    chunk_args = [(samp_file, offset, chunk_size, K) for offset in offsets]

    print(f"Processing {len(CHUNKS)} chunks...")
    if test:
        xchunk_args = chunk_args[:3]
        Xref_kmers = set(list(ref_kmers)[:100])
        xn_cores = 3

        with Pool(processes=xn_cores, initializer=init_worker, initargs=(Xref_kmers,)) as pool:
            results = pool.map(process_chunk, xchunk_args)
    else:
        chunk_args = [(samp_file, start, stop, K) for start, stop in CHUNKS]
        with Pool(processes=n_cores, initializer=init_worker, initargs=(ref_kmers,)) as pool:
            results = pool.map(process_chunk, chunk_args)

    print(f"Done.")
    return [x for chunk in results for x in chunk]
    