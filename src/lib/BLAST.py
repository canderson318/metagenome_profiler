pip install biopython matplotlib pandas numpy requests
# For local BLAST (optional but recommended for large datasets):
# conda install -c bioconda blast

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Align
import pandas as pd
import numpy as np
import requests
import time


# An unknown sequence — let's see what BLAST says
unknown_seq = """
ATGAGAGCCTTGTCCTTTATAGCACAAATGGCAACCAAATGCATCAAAATGAGAGAATCAGTCTCAATGG
ATGGAAATGAATAAGCCTGTCAAGCTGGATGGAAACAACATTACAACAGAAGATCAAATGGGGAAGATGG
AGAAAGAAGTGGTCATCAGAAGCAACAGCTTATCATCAGCAGCAGCCATCAGCAGGATAGAAGCGGCAGC
""".replace("\n", "").strip()

Entrez.email = "your@email.com"

print(f"Query length: {len(unknown_seq)} bp")
print("Running BLAST against NT database...")

result_handle = NCBIWWW.qblast(
    "blastn",
    "nt",
    unknown_seq,
    hitlist_size=10,
    format_type="XML"
)

# Save result to avoid re-running
with open("blast_result.xml", "w") as f:
    f.write(result_handle.read())
result_handle.close()
print("BLAST complete.")

    
    
def parse_blast_results(xml_file: str) -> pd.DataFrame:
    with open(xml_file) as f:
        blast_record = next(NCBIXML.parse(f))

    rows = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps[:1]:  # best HSP per hit
            rows.append({
                "title": alignment.title[:80],
                "accession": alignment.accession,
                "score": hsp.score,
                "e_value": hsp.expect,
                "identity_pct": 100 * hsp.identities / hsp.align_length,
                "alignment_length": hsp.align_length,
                "query_coverage_pct": 100 * (hsp.query_end - hsp.query_start + 1) / blast_record.query_length,
            })

    return pd.DataFrame(rows)

df = parse_blast_results("blast_result.xml")
print(df[["title", "identity_pct", "e_value", "query_coverage_pct"]].to_string(index=False))



def classify_blast_hit(identity: float, coverage: float) -> str:
    """Classify the BLAST hit quality."""
    if identity >= 99 and coverage >= 95:
        return "Near-identical match (same strain or minor variant)"
    elif identity >= 95 and coverage >= 90:
        return "High confidence match (closely related strain)"
    elif identity >= 80 and coverage >= 70:
        return "Moderate match (same species, different strain)"
    elif identity >= 70:
        return "Distant homolog (same genus or family)"
    else:
        return "Weak match — possibly novel or contamination"

for _, row in df.head(5).iterrows():
    classification = classify_blast_hit(row["identity_pct"], row["query_coverage_pct"])
    print(f"\n{row['title'][:60]}")
    print(f"  Identity: {row['identity_pct']:.1f}%, Coverage: {row['query_coverage_pct']:.1f}%")
    print(f"  {classification}")
    
      from Bio.Blast import NCBIWWW, NCBIXML

  def blast_contig(contig_seq):
      result_handle = NCBIWWW.qblast(
          program="blastn",
          database="nt",          # or "refseq_rna"
          sequence=contig_seq,
          entrez_query="Viruses[Organism]",  # restrict to viruses only
          hitlist_size=5
      )
      return NCBIXML.read(result_handle)
  

############
############
############
############

from Bio import Entrez

Entrez.email = "your@email.com"

def get_metadata(accession):
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="gb",
        retmode="text"
    )
    record = SeqIO.read(handle, "genbank")
    return {
        "accession": accession,
        "organism": record.annotations.get("organism"),
        "taxonomy": record.annotations.get("taxonomy"),  # full lineage list
        "description": record.description
    }

results = []
for contig in contigs:
    blast_rec = blast_contig(str(contig.seq))
    if not blast_rec.alignments:
        results.append({"contig": contig.id, "hit": None})
        continue

    top_hit = blast_rec.alignments[0]
    top_hsp = top_hit.hsps[0]
    accession = top_hit.accession

    meta = get_metadata(accession)
    results.append({
        "contig": contig.id,
        "hit": top_hit.title,
        "accession": accession,
        "pct_identity": top_hsp.identities / top_hsp.align_length * 100,
        "evalue": top_hsp.expect,
        "organism": meta["organism"],
        "taxonomy": meta["taxonomy"]
    })
