import pytest
from src.lib.assembly.DBG import *

    
# \\\
# –––– Fixtures
# \\\

@pytest.fixture
def linear_reads():
    # simple linear sequence: ab→bc→cd→de
    return ["abcde"] * 5
# linear_reads = linear_reads()

@pytest.fixture
def linear_graph(linear_reads):
    g = DBG(linear_reads, k=3)
    g.calculate_degrees()
    return g
# linear_graph = linear_graph(linear_reads)

@pytest.fixture
def branching_reads():
    # two paths diverging at 'ab': abcde and abfde
    return ["axbcde"] * 5 + ["abfdeghij"] * 5
# branching_reads = branching_reads()

@pytest.fixture
def branching_graph(branching_reads):
    g = DBG(branching_reads, k=3)
    g.calculate_degrees()
    return g
# branching_graph = branching_graph(branching_reads)
# branching_graph.plotG()

@pytest.fixture
def cycle_reads():
    # "abcab" produces cycle: ab→bc→ca→ab
    return ["abcab"] * 5
# cycle_reads = cycle_reads()

@pytest.fixture
def cycle_graph(cycle_reads):
    g = DBG(cycle_reads, k=3)
    g.calculate_degrees()
    return g
# cycle_graph = cycle_graph(cycle_reads)
# cycle_graph.plotG()

# \\\
# –––– DBG construction
# \\\

def test_graph_contains_expected_edges(linear_graph):
    G = linear_graph.G
    assert "bc" in G["ab"]
    assert "cd" in G["bc"]
    assert "de" in G["cd"]
# G = linear_graph.G
# "bc" in G["ab"]
# "cd" in G["bc"]
# "de" in G["cd"]

def test_graph_edge_coverage_matches_read_count(linear_reads, linear_graph):
    n = len(linear_reads)
    G = linear_graph.G
    assert G["ab"]["bc"] == n
    assert G["bc"]["cd"] == n
    assert G["cd"]["de"] == n
# n = len(linear_reads)
# G = linear_graph.G
# G["ab"]["bc"] == n
# G["bc"]["cd"] == n
# G["cd"]["de"] == n


def test_incoming_degree_source_node_is_zero(linear_graph):
    # "ab" is the source — nothing points to it
    assert linear_graph._DBG__in_deg.get("ab", 0) == 0
# linear_graph._DBG__in_deg.get("ab", 0) == 0

def test_incoming_degree_interior_node_correct(linear_graph):
    assert linear_graph._DBG__in_deg["bc"] == 1
    assert linear_graph._DBG__in_deg["cd"] == 1
# linear_graph._DBG__in_deg["bc"] == 1
# linear_graph._DBG__in_deg["cd"] == 1

# \\\
# –––– prune_low_coverage_nodes
# \\\

def test_prune_removes_edges_below_threshold():
    g = DBG(["abcde"] * 1, k=3)  # coverage 1
    g.prune_low_coverage_nodes(coverage_thresh=2)
    assert "bc" not in g.G.get("ab", {})
# g = DBG(["abcde"] , k=3)  # coverage 1
# g.prune_low_coverage_nodes(coverage_thresh=2)
# "bc" not in g.G.get("ab", {})

def test_prune_keeps_edges_above_threshold(linear_graph):
    linear_graph.prune_low_coverage_nodes(coverage_thresh=2)
    G = linear_graph.G
    assert "bc" in G.get("ab", {})
    assert "cd" in G.get("bc", {})
    assert "de" in G.get("cd", {})
# linear_graph.prune_low_coverage_nodes(coverage_thresh=2)
# G = linear_graph.G
# "bc" in G.get("ab", {})
# "cd" in G.get("bc", {})
# "de" in G.get("cd", {})

def test_prune_removes_empty_nodes():
    g = DBG(["abcde"] * 1, k=3)
    g.prune_low_coverage_nodes(coverage_thresh=2)
    assert "ab" not in g.G
# g = DBG(["abcde"] * 1, k=3)
# g.prune_low_coverage_nodes(coverage_thresh=2)
# "ab" not in g.G


# \\\
# –––– make_contigs
# \\\

def test_make_contigs_tuples_have_three_elements(linear_graph):
    linear_graph.make_contigs()
    contigs = linear_graph.contigs
    for entry in contigs:
        assert len(entry) == 3
# linear_graph.make_contigs()
# contigs = linear_graph.contigs
# for entry in contigs:
#     print(len(entry) == 3)

def test_make_contigs_linear_recovers_sequence(linear_graph):
    linear_graph.make_contigs()
    contigs = linear_graph.contigs
    seqs = [seq for seq, _, _ in contigs]
    assert "abcde" in seqs
# seqs = [seq for seq, _, _ in contigs]
# "abcde" in seqs

def test_make_contigs_linear_returns_one_contig(linear_graph):
    linear_graph.make_contigs()
    contigs = linear_graph.contigs
    assert len(contigs) == 1
# len(contigs) == 1

def test_make_contigs_coverages_nonempty(linear_graph):
    linear_graph.make_contigs()
    contigs = linear_graph.contigs
    for _, _, covs in contigs:
        assert len(covs) > 0
# for _, _, covs in contigs:
#     print(len(covs) > 0)

def test_make_contigs_path_nodes_match_seq_length(linear_graph):
    linear_graph.make_contigs()
    contigs = linear_graph.contigs
    for seq, nodes, _ in contigs:
        # seq = first node (k-1 chars) + 1 char for each subsequent node
        k_minus_1 = len(nodes[0])
        expected_len = k_minus_1 + (len(nodes) - 1)
        assert len(seq) == expected_len
# for seq, nodes, _ in contigs:
#     # seq = first node (k-1 chars) + 1 char for each subsequent node
#     k_minus_1 = len(nodes[0])
#     expected_len = k_minus_1 + (len(nodes) - 1)
#     print(len(seq) == expected_len)

def test_make_contigs_branching_returns_two_contigs(branching_graph):
    branching_graph.make_contigs()
    contigs = branching_graph.contigs
    seqs = [seq for seq, _, _ in contigs]
    assert len(seqs) == 5
    for true_contig in ["axbcde","abfde","axbcdeghij","deghij","abfdeghij"]:
        assert true_contig in seqs
# branching_graph.make_contigs()
# contigs = branching_graph.contigs
# seqs = [seq for seq, _, _ in contigs]
# assert len(seqs) == 5
# for true_contig in ["axbcde","abfde","axbcdeghij","deghij","abfdeghij"]:
#     true_contig in seqs
    
def test_make_contigs_cycle_recovered(cycle_graph):
    # pure cycle ab→bc→ca→ab, no branching nodes — requires cycle pass to run
    cycle_graph.make_contigs()
    contigs = cycle_graph.contigs
    seqs = [seq for seq, _, _ in contigs]
    assert len(contigs) >= 1
    assert any(len(s) >= 4 for s in seqs)
# cycle_graph.make_contigs()
# contigs = cycle_graph.contigs
# seqs = [seq for seq, _, _ in contigs]
# len(contigs) >= 1
# any(len(s) >= 4 for s in seqs)
