import pytest
import numpy as np
import tempfile
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from src.lib.align_reads import host_match, find_chunk_offsets, align_reads
from src.lib.kmers import extract_kmers


# ---- fixtures ----

@pytest.fixture
def ref_seq():
    return "ACGTACGTACGTACGTACGTACGTACGTACGT" * 10  # 320bp known sequence

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
    records = [
        SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGTACGT" * 3), id=f"r{i}",
                  description="", letter_annotations={"phred_quality": [40] * 96})
        for i in range(30)
    ]
    path = tmp_path / "test.fastq"
    with open(path, "w") as f:
        SeqIO.write(records, f, "fastq")
    return path

