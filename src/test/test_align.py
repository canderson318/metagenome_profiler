import pytest
import numpy as np
from pathlib import Path
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# funcs to test
from src.lib.align_reads import host_match, find_chunk_offsets, align_reads
from src.lib.kmers import extract_kmers


# \\\
# –––– fixtures 
# \\\

@pytest.fixture
def ref_seq():
    return "ACGTACGTACGTACGTACGTACGTACGTACGT" * 10  # ~320bp known sequence

@pytest.fixture
def ref_kmers(ref_seq):
    return extract_kmers(ref_seq, k=31)

@pytest.fixture
def matching_read(ref_seq):
    return SeqRecord(Seq(ref_seq[:100]), id="read_host", description="")

@pytest.fixture
def nonmatching_read():
    rng = np.random.default_rng(42)
    seq = "".join(rng.choice(["A", "T", "C", "G"], 100))
    return SeqRecord(Seq(seq), id="read_nonhost", description="")

@pytest.fixture
def tmp_fastq(tmp_path):
    s = "ACGTACGTACGTACGTACGTACGTACGTACGT"*3
    records = [
        SeqRecord(Seq(s), id=f"r{i}",
                  description="", letter_annotations={"phred_quality": [40] * len(s) })
        for i in range(30)
    ]
    path = tmp_path / "test.fastq"
    with open(path, "w") as f:
        SeqIO.write(records, f, "fastq")
    return path

# \\\
# –––– host_match 
# \\\

def test_host_match_high_for_matching_read(matching_read, ref_kmers):
    score = host_match(matching_read, ref_kmers, K=31)
    assert score > 0.5
# host_match(matching_read(ref_seq()), ref_kmers(ref_seq()), K=31)

def test_host_match_low_for_nonmatching_read(nonmatching_read, ref_kmers):
    score = host_match(nonmatching_read, ref_kmers, K=31)
    assert score < 0.2
# host_match(nonmatching_read(), ref_kmers(ref_seq()), K=31)

def test_host_match_short_read_returns_zero(ref_kmers):
    short = SeqRecord(Seq("ACGT"), id="short", description="")
    score = host_match(short, ref_kmers, K=31)
    assert score == 0.0
# short = SeqRecord(Seq("ACGT"), id="short", description="")
# host_match(short, ref_kmers(ref_seq()), K=31)

def test_host_match_k_mismatch_raises(matching_read, ref_kmers):
    with pytest.raises(ValueError):
        host_match(matching_read, ref_kmers, K=15)
# host_match(matching_read(ref_seq()), ref_kmers(ref_seq()), K=15)

# \\\
# –––– find_chunk_offsets 
# \\\
def test_find_chunk_byte_offsets_returns_correct_count(tmp_fastq):
    offsets, chunk_size = find_chunk_offsets(tmp_fastq, n_chunks=3, N_reads=30)
    assert len(offsets) == 3
# tmp_path = Path("/tmp/test_align")
# tmp_path.mkdir(exist_ok=True)
# offsets, chunk_size = find_chunk_offsets(tmp_fastq(tmp_path), n_chunks=3, N_reads=30)
# len(offsets)

def test_find_chunk_byte_offsets_first_offset_is_0_and_last_4010(tmp_fastq):
    offsets, _ = find_chunk_offsets(tmp_fastq, n_chunks=3, N_reads=30)
    assert offsets[0] == 0
    assert offsets[2] == 4010
# offsets, _ = find_chunk_offsets(tmp_fastq(tmp_path), n_chunks=3, N_reads=30)
# (offsets[0] == 0) & (offsets[2] == 4010)

def test_byte_offset_line_content(tmp_fastq):
    offsets, _ = find_chunk_offsets(tmp_fastq, n_chunks=3, N_reads=30)
    offset = offsets[1] # 2000
    with open(tmp_fastq) as fin:
        fin.seek(offset)
        seq = fin.readline()
    assert seq.strip() == "@r10"
# offsets, _ = find_chunk_offsets(tmp_fastq(tmp_path), n_chunks=3, N_reads=30)
# offset = offsets[1] # 2000
# with open(tmp_fastq(tmp_path)) as fin:
#     fin.seek(offset)
#     seq = fin.readline()
# seq.strip() == "@r10"

# \\\
# –––– align_reads 
# \\\
def test_align_reads_returns_two_arrays(tmp_fastq, ref_kmers, tmp_path):
    import pickle
    index_file = tmp_path / "index.pkl"
    with open(index_file, "wb") as f:
        pickle.dump(ref_kmers, f)

    hit_inds, non_hit_inds = align_reads(tmp_fastq, index_file, K=31, test=False, ncores = 3)
    assert hit_inds is not None
    assert non_hit_inds is not None
# import pickle
# index_file = tmp_path / "index.pkl"
# with open(index_file, "wb") as f:
#     pickle.dump(ref_kmers(ref_seq()), f)
# hit_inds, non_hit_inds = align_reads(tmp_fastq(tmp_path), index_file, K=31, test=False, ncores = 3)
# hit_inds is not None
# non_hit_inds is not None

def test_align_reads_hit_and_nonhit_partition_reads(tmp_fastq, ref_kmers, tmp_path):
    import pickle
    index_file = tmp_path / "index.pkl"
    with open(index_file, "wb") as f:
        pickle.dump(ref_kmers, f)

    hit_inds, non_hit_inds = align_reads(tmp_fastq, index_file, K=31, test=False, ncores = 3)
    if hit_inds is not None:
        combined = np.sort(np.concatenate([hit_inds, non_hit_inds]))
        assert np.array_equal(combined, np.arange(len(combined)))
# import pickle
# index_file = tmp_path / "index.pkl"
# with open(index_file, "wb") as f:
#     pickle.dump(ref_kmers(ref_seq()), f)
# hit_inds, non_hit_inds = align_reads(tmp_fastq(tmp_path), index_file, K=31, test=False)
# if hit_inds is not None:
#     combined = np.sort(np.concatenate([hit_inds, non_hit_inds]))
#     np.array_equal(combined, np.arange(len(combined)))

