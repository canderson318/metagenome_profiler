import pytest
from src.lib.assembly.DBG import *

# \\\\
# –––– Debrujin graph construction
# \\\\
    
# np.random.default_rng(123)
# reads = sim_reads("".join(np.random.choice(["A", "T", "C", "G"], 500)), 1000, 31)
reads = sim_reads(load_hobbit(236), 1000, 15)
k = 5
G = DBG(reads, k = k)
G.calculate_degs()
G.prune_low_coverage_nodes(coverage_thresh=2)
G.prune_low_coverage_connections(pruned_seq_size_ceiling=10, seq_divergence_ceiling=.2)

import matplotlib.pyplot as plt
_,ax = plt.subplots(figsize = (15,15))
G.plotG(weighted = False, ax = ax)

# \\\\
# –––– Contig Creation
# \\\\
    
G.make_contigs()

contigs = [seq for seq, kmers, covs in G.contigs]
lens = [len(seq) for seq in contigs ]
top_contigs = [seq for seq in contigs if  len(seq) > np.percentile(lens,85)]


# \\\
# –––– Fixtures
# \\\

@pytest.fixture
def linear_reads():
    # simple linear sequence: ab→bc→cd→de
    return ["abcde"] * 5

@pytest.fixture
def linear_graph(linear_reads):
    g = DBG(linear_reads, k=3)
    g.calculate_incoming_weights()
    return g

@pytest.fixture
def branching_reads():
    # two paths diverging at 'ab': abcde and abfde
    return ["abcde"] * 5 + ["abfde"] * 5

@pytest.fixture
def branching_graph(branching_reads):
    g = DBG(branching_reads, k=3)
    g.calculate_incoming_weights()
    return g

@pytest.fixture
def cycle_reads():
    # "abcab" produces cycle: ab→bc→ca→ab
    return ["abcab"] * 5

@pytest.fixture
def cycle_graph(cycle_reads):
    g = DBG(cycle_reads, k=3)
    g.calculate_incoming_weights()
    return g


# \\\
# –––– DBG construction
# \\\

def test_graph_contains_expected_edges(linear_graph):
    G = linear_graph.G
    assert "bc" in G["ab"]
    assert "cd" in G["bc"]
    assert "de" in G["cd"]

def test_graph_edge_coverage_matches_read_count(linear_reads, linear_graph):
    n = len(linear_reads)
    G = linear_graph.G
    assert G["ab"]["bc"] == n
    assert G["bc"]["cd"] == n
    assert G["cd"]["de"] == n

def test_incoming_weights_source_node_is_zero(linear_graph):
    # "ab" is the source — nothing points to it
    assert linear_graph.incoming_weights.get("ab", 0) == 0

def test_incoming_weights_interior_node_correct(linear_graph):
    assert linear_graph.incoming_weights["bc"] == 5
    assert linear_graph.incoming_weights["cd"] == 5


# \\\
# –––– prune_low_coverage_nodes
# \\\

def test_prune_removes_edges_below_threshold():
    g = DBG(["abcde"] * 1, k=3)  # coverage 1
    g.prune_low_coverage_nodes(coverage_thresh=2)
    assert "bc" not in g.G.get("ab", {})

def test_prune_keeps_edges_above_threshold(linear_graph):
    linear_graph.prune_low_coverage_nodes(coverage_thresh=2)
    G = linear_graph.G
    assert "bc" in G.get("ab", {})
    assert "cd" in G.get("bc", {})
    assert "de" in G.get("cd", {})

def test_prune_removes_empty_nodes():
    g = DBG(["abcde"] * 1, k=3)
    g.prune_low_coverage_nodes(coverage_thresh=2)
    assert "ab" not in g.G


# \\\
# –––– make_contigs
# \\\

def test_make_contigs_returns_list(linear_graph):
    contigs = linear_graph.make_contigs()
    assert isinstance(contigs, list)

def test_make_contigs_stored_on_graph(linear_graph):
    linear_graph.make_contigs()
    assert linear_graph.contigs is not None

def test_make_contigs_tuples_have_three_elements(linear_graph):
    contigs = linear_graph.make_contigs()
    for entry in contigs:
        assert len(entry) == 3

def test_make_contigs_linear_recovers_sequence(linear_graph):
    contigs = linear_graph.make_contigs()
    seqs = [seq for seq, _, _ in contigs]
    assert "abcde" in seqs

def test_make_contigs_linear_returns_one_contig(linear_graph):
    contigs = linear_graph.make_contigs()
    assert len(contigs) == 1

def test_make_contigs_coverages_nonempty(linear_graph):
    contigs = linear_graph.make_contigs()
    for _, _, covs in contigs:
        assert len(covs) > 0

def test_make_contigs_path_nodes_match_seq_length(linear_graph):
    contigs = linear_graph.make_contigs()
    for seq, nodes, _ in contigs:
        # seq = first node (k-1 chars) + 1 char per subsequent node
        k_minus_1 = len(nodes[0])
        expected_len = k_minus_1 + (len(nodes) - 1)
        assert len(seq) == expected_len

def test_make_contigs_branching_returns_two_contigs(branching_graph):
    contigs = branching_graph.make_contigs()
    seqs = [seq for seq, _, _ in contigs]
    assert len(seqs) == 2
    assert "abcde" in seqs
    assert "abfde" in seqs

def test_make_contigs_cycle_recovered(cycle_graph):
    # pure cycle ab→bc→ca→ab, no branching nodes — requires cycle pass to run
    contigs = cycle_graph.make_contigs()
    seqs = [seq for seq, _, _ in contigs]
    assert len(contigs) >= 1
    assert any(len(s) >= 4 for s in seqs)
