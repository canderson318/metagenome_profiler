from typing import List, Dict, Tuple
import numpy as np
from  tqdm import tqdm
import matplotlib.pyplot as plt
from src.lib.assembly.DBG import *

# \\\\
# \\\\
# –––– Debrujin graph construction
# \\\\
# \\\\
    
# np.random.default_rng(123)
# reads = sim_reads("".join(np.random.choice(["A", "T", "C", "G"], 500)), 1000, 31)
reads = sim_reads(load_hobbit(236), 1000, 15)
k = 5
G = DBG(reads, k = k)
G.calculate_incoming_weights()
G.calculate_outgoing_weights()
G.prune_low_coverage_nodes(coverage_thresh=2)
G.prune_low_coverage_connections(pruned_seq_size_ceiling=3, seq_divergence_ceiling=.8)
G.plotG(weighted = False)

# \\\\
# \\\\
# –––– Contig Creation
# \\\\
# \\\\

def make_contigs(graph: DBG):
    """
    """
    contigs = []
    visited = set()

    in_weights = graph.incoming_weights
    G = graph.G
    
    for u in tqdm(G, desc="Finding Contigs"):
        # get nodes connected to u
        out_edges = G[u]
        # check if there are outgoing edges 
        start_branches = (len(out_edges) != 1) or (in_weights.get(u, 0) != 1)
        if not start_branches:
            continue
        # for each vertex in outgoing edges
        for v in out_edges:
            # if this pair has been visited, this is a loop, pass on 
            if (u, v) in visited:
                continue
            path_nodes = [u, v]
            # edge coverages outgoing from u
            coverages = [G[u][v]]
            visited.add((u, v))
            curr = v
            # while no alternate routes to take (no branches) and only one previous edge
            while len(G.get(curr, {})) == 1 and in_weights.get(curr, 0) == 1:
                nxt = next(iter(G[curr]))
                path_nodes.append(nxt)
                coverages.append(G[curr][nxt])
                visited.add((curr, nxt))
                curr = nxt
            # construct contiguous sequence: first node + last char of all next
            seq = path_nodes[0] + ''.join(n[-1] for n in path_nodes[1:])
            contigs.append((seq, path_nodes, coverages))
        
        # traverse again, this time exploring loops
        for u in G:
            for v in G[u]:
                if (u, v) in visited:
                    continue
                path_nodes = [u, v]
                coverages = [G[u][v]]
                visited.add((u, v))
                curr = v
                # while not at begining of loop again
                while curr != u and len(G[curr]):
                    nxt = next(iter(G[curr]))
                    path_nodes.append(nxt)
                    coverages.append(G[curr][nxt])
                    visited.add((curr, nxt))
                    curr = nxt
                seq = path_nodes[0] + ''.join(n[-1] for n in path_nodes[1:])
                contigs.append((seq, path_nodes, coverages))
                    
    return contigs

    # selected_contigs = []
    # for seq, nodes, covs in tqdm(contigs, desc="Selecting long contigs", unit="unitig"):
    #     if len(seq) >= length_thresh and (sum(covs)/len(covs)) >= cov_thresh:
    #         selected_contigs.append((seq, nodes, covs))

    # return selected_contigs

# def select_linear_contigs(contigs: List[Tuple[str, List[str], List[int]]], graph: DBG, length_thresh: int = 1000, cov_thresh: int = 5):
#     """
#     """
#     in_weights = graph.in_degrees
#     G = graph.G
    
#     contigs = []
#     for seq, nodes, covs in tqdm(contigs, desc="Selecting linear contigs", unit="contig"):
#         start, end = nodes[0], nodes[-1]
#         if in_weights.get(start, 0) == 0 and out_deg.get(end, 0) == 0:
#             if len(seq) >= length_thresh and (sum(covs)/len(covs)) >= cov_thresh:
#                 contigs.append((seq, sum(covs)/len(covs)))

#     return contigs

contigs = make_contigs(G)


