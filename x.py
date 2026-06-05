# np.random.default_rng(123)
# reads = sim_reads("".join(np.random.choice(["A", "T", "C", "G"], 500)), 1000, 31)
reads = sim_reads(load_hobbit(236), 1000, 15)
k = 5
G = DBG(reads, k = k)
G.calculate_degrees()
G.prune_low_coverage_nodes(coverage_thresh=2)
G.prune_low_coverage_connections(seq_divergence_ceiling=.2)

import matplotlib.pyplot as plt
_,ax = plt.subplots(figsize = (15,15))
# G.plotG(weighted = False, ax = ax)

# \\\\
# –––– Contig Creation
# \\\\
    
G.make_contigs()

contigs = [seq for seq, kmers, covs in G.contigs]
lens = [len(seq) for seq in contigs ]
top_contigs = [seq for seq in contigs if  len(seq) > np.percentile(lens,85)]
