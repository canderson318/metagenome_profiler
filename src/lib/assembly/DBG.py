from collections import defaultdict
from typing import List, Dict
from difflib import SequenceMatcher
import numpy as np
from  tqdm import tqdm
import networkx as nx
import matplotlib.pyplot as plt
import re

class DBG:
    def __init__(self, reads: List[str], k: int):
        self.k = k
        self.__graph = None
        self._construct_graph(reads, k)
        self.__incoming_weight = None
        self.__outgoing_weight = None

    def _construct_graph(self, reads: List[str], k: int):
        """
        Generate debrujin graph from reads. 
        Works by ticking nested dictionary at index of [kmer[-end], kmer[-start]].
        Each top-level key is preceding, items are successor nodes
        
        :param reads: list of gene reads
        :param k: Kmer size
        """
        # make dict where new keys make new element so indexing at nonexistent key doesnt keyerror
        self.__graph = defaultdict(lambda: defaultdict(int))
        for read in tqdm(reads, desc = "Constructing graph"):
            for i in range(len(read) - k + 1):
                kmer = read[i:i+k]
                # add index to dict at u, v k-1mer
                u, v = kmer[:-1], kmer[1:]
                self.__graph[u][v] += 1
        
    def calculate_incoming_weights(self):
        """
        Compute incoming edge weights between each node. 
        """
        self.__incoming_weight = defaultdict(int)
        for u, outs in tqdm(self.__graph.items(), desc = 'Calculating incoming edge weights'):
            for v,n in outs.items():
                self.__incoming_weight[v] += n

    def calculate_outgoing_weights(self):
        self.__outgoing_weight = defaultdict(int)
        for u, outs in tqdm(self.__graph.items(), desc = 'Calculating outgoing edge weights'):
            self.__outgoing_weight[u] = len(self.__graph.get(u, {}))
        
    def prune_low_coverage_nodes(self, coverage_thresh: int = 2):
        """
        For each node, cut edges (and nodes) where coverage < threshold.
        """
        for u in tqdm(list(self.__graph), desc = "Pruning nodes"):
            for v, w in list(self.__graph[u].items()): 
                if w < coverage_thresh: 
                    del self.__graph[u][v] # remove edges to infrequent nodes
            if not self.__graph[u]:
                del self.__graph[u] # remove nodes with no edges
                
    @staticmethod
    def _calc_divergence(seq1, seq2):
        matcher = SequenceMatcher(None, seq1, seq2)
        ed = 0
        for tag, i1, i2, j1, j2 in matcher.get_opcodes():
            if tag != 'equal':
                ed += max(i2 - i1, j2 - j1)
        return ed / min(len(seq1), len(seq2))

    def prune_low_coverage_connections(self,pruned_seq_size_ceiling: int = 1000, seq_divergence_ceiling: float = 0.2):
        """
        For each node, compare forward connected nodes and pick larger of pair, swallowing up hte smaller's coverage 
        where sequences are small enough (< pruned_seq_size_ceiling) and similar enough (< seq_divergence_ceiling) 
        
        :param pruned_seq_size_celing: sequence pairs smaller than this will be compared
        :param seq_divergence_ceiling: sequences more divergent than this (where increasing to 1 is more different) will not be collapsed into each other
        """
        for u, outs in tqdm(list(self.__graph.items()), desc = 'Pruning connections'):
            vs = list(outs.keys())
            # if no bulge
            if len(vs) < 2:
                continue
            # for each pair of v_i & v_i+1
            for i in range(len(vs)):
                for j in range(i+1, len(vs)):
                    # get node sequence
                    v1, v2 = vs[i], vs[j] 
                    # reconstruct original sequence for both: prev + last char of node
                    seq1 = u + v1[-1] 
                    seq2 = u + v2[-1]
                    # if sequences small enough 
                    if len(seq1) <= pruned_seq_size_ceiling and len(seq2) <= pruned_seq_size_ceiling:
                        # calculate string difference between both sequences
                        divergence = self._calc_divergence(seq1,seq2) # closer to 1 is more divergent
                        # if they are divergent enough
                        if divergence < seq_divergence_ceiling: 
                            # pick more covered of two, adding its coverage and removing it
                            if self.__graph[u][v1] >= self.__graph[u][v2]: # if v2 seen more, swallow v2 coverage
                                self.__graph[u][v1] += self.__graph[u].pop(v2)
                            else:
                                self.__graph[u][v2] += self.__graph[u].pop(v1)  # if v1 seen more, etc.    

    @property
    def incoming_weights(self):
        if self.__incoming_weight is None:
            self.calculate_incoming_weights()
        return self.__incoming_weight
    
    @property
    def outgoing_weights(self):
        if self.__outgoing_weight is None:
            self.calculate_outgoing_weights()
        return self.__outgoing_weight
    
    @property
    def G(self):
        return self.__graph

    def plotG(self, weighted = True, ax = None):
        # make graph object
        if weighted:
            graph_list = [(u,v,w) for u,vs in self.__graph.items() for v,w in vs.items()]
            G = nx.Graph()
            G.add_weighted_edges_from(graph_list)
        else:
            graph_list = [(u,v) for u,vs in self.__graph.items() for v,_ in vs.items()]
            G = nx.Graph()
            G.add_edges_from(graph_list)
        # plot graph
        _,ax = plt.subplots(figsize = (10,10)) if ax is None else None, ax
        pos = nx.spring_layout(G, seed = 2026)
        nx.draw_networkx_nodes(G, ax = ax, pos = pos,node_size=800, node_color = "#9c5ac2ff")
        nx.draw_networkx_labels(G, ax = ax, pos = pos, font_size=9,font_color = "black")
        nx.draw_networkx_edges(G, edge_color="black", pos = pos)

def load_hobbit(up_to = None):
    """
    Load First paragraph of hobbit and return cleaned sequence.
    :param up_to: Size of returned sequence from start
    """
    with open("in/hobbit.txt") as fin:
        seq = fin.read()
    seq = re.sub(r"[.,:;\-\"'?)(),]", "", seq)
    seq = seq.lower()
    seq = re.sub(r"\s", "", seq)
    seq = re.sub(r"_", "", seq)
    if up_to is None:
        return seq
    return seq[:up_to]

def sim_reads(seq, N, read_len):
    '''
    Simulate reads from a sequence
    
    :param seq: string of characters
    :param N: number of reads to generate
    :param read_len: size of each read
    '''
    np.random.default_rng(123)
    return [seq[strt:strt+read_len]  
            for _ in range(N) 
                for strt in [np.random.choice(range(len(seq) - read_len+1))]
                ]
    
