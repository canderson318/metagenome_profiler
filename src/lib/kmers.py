

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

