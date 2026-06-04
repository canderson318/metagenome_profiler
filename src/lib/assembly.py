from collections import defaultdict
from typing import List, Dict
from difflib import SequenceMatcher

class DBG:
    
    def __init__(self, reads: List[str], k: int):
        self.k = k
        self.__graph = None
        self.__in_deg = None
        self._build(reads, k)

    def _build(self, reads: List[str], k: int):
        self.__graph = defaultdict(lambda: defaultdict(int))
        for read in reads:
            for i in range(len(read) - k + 1):
                kmer = read[i:i+k]
                u, v = kmer[:-1], kmer[1:]
                self.__graph[u][v] += 1

    def compute_in_degrees(self):
        self.__in_deg = defaultdict(int)
        for u, outs in self.__graph.items():
            for v in outs:
                self.__in_deg[v] += outs[v]

    def prune_low_coverage_edges(self, min_coverage: int = 2):
        for u in list(self.graph):
            for v, w in list(self.graph[u].items()):
                if w < min_coverage:
                    del self.graph[u][v]
            if not self.graph[u]:
                del self.graph[u]

    def collapse_bulge(self,max_bulge: int = 1000, max_div: float = 0.2):
        for u, outs in list(self.__graph.items()):
            vs = list(outs.keys())
            if len(vs) < 2:
                continue
            for i in range(len(vs)):
                for j in range(i+1, len(vs)):
                    v1, v2 = vs[i], vs[j]
                    seq1 = u + v1[-1]
                    seq2 = u + v2[-1]
                    if len(seq1) <= max_bulge and len(seq2) <= max_bulge:
                        matcher = SequenceMatcher(None, seq1, seq2)
                        ed = 0
                        for tag, i1, i2, j1, j2 in matcher.get_opcodes():
                            if tag != 'equal':
                                ed += max(i2 - i1, j2 - j1)
                        divergence = ed / min(len(seq1), len(seq2))
                        if divergence < max_div:
                            if self.__graph[u][v1] >= self.__graph[u][v2]:
                                self.__graph[u][v1] += self.__graph[u].pop(v2)
                            else:
                                self.__graph[u][v2] += self.__graph[u].pop(v1)    

    @property
    def in_degrees(self) -> Dict[str, int]:
        if self.__in_deg is None:
            self._compute_in_degrees()
        return self.__in_deg