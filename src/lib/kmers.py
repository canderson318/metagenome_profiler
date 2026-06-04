from Bio.Seq import Seq

def extract_kmers(seq,k):
    kmers = set()
    for i in range(len(seq) - k):
        kmer = seq[i:i+k]
        # take smaller of two strand kmers (identical but read from opposite strands)
        rc = str(Seq(kmer).reverse_complement())
        kmers.add(min(kmer, rc))  
    return kmers

