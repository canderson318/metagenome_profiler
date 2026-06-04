from pathlib import Path
from Bio import SeqIO

def get_num_reads(samp_file):
    samp_file = Path(samp_file)
    size_file = samp_file.parent / ("."+samp_file.stem + "_size")
    if not size_file.exists():
        size = sum(1 for _ in SeqIO.parse(samp_file, 'fastq'))
        with open(size_file, 'wt') as f:
            f.write(str(size))
    else:
        with open(size_file, 'rt') as f:
            size = int(f.read())
    return size

complement = str.maketrans('ACGT', 'TGCA')
def reverse_complement(seq):
    return seq.translate(complement)[::-1]

def extract_kmers(seq,k):
    kmers = set()
    for i in range(len(seq) - k):
        kmer = seq[i:i+k]
        # take smaller of two strand kmers (identical but read from opposite strands)
        rc = reverse_complement(kmer)
        kmers.add(min(kmer, rc))  
    return kmers

