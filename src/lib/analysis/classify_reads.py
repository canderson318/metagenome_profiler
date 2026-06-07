import subprocess as sp
from pathlib import Path
import pandas as pd

def classify_reads(fasta_path: str or Path, force = False):
    """
    Queries fasta reads against viral genome database using kraken2. 
    Returns a pandas DataFrame of classified reads including their fasta accession number, read index, 
    and the location of alignment with the viral taxon.
    
    :param fasta_path: Path to fasta file.
    :return: pandas DataFrame `result`
    """
    fasta_path = Path(fasta_path)
    
    if not fasta_path.exists():
        raise FileNotFoundError

    viral_db_path = in_dir / "kraken_2_viral_db"
    if not viral_db_path.exists():
        raise FileNotFoundError

    result_dir = out_dir / "kraken"
    result_dir.mkdir(exist_ok=True)

    if not (result_dir / "report.txt").exists() or force:
        print("Querying viral reference database...")
        sp.run([
            "kraken2", "--db", str(viral_db_path),
            "--output", str(result_dir / "results.txt"),
            "--report", str(result_dir / "report.txt"),
            str(fasta_path)
        ], stderr=sp.DEVNULL, stdout = sp.DEVNULL)  
        print("Done.")
    else:
        print(f"Viral taxon classification already exists at {str(result_dir)} so using those.\nUse `force` to query again.")
    
    print("Loading results...")
    # columns: pct_reads, reads_covered, reads_direct, rank, taxid, name
    ### rank: U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies
    report = pd.read_csv(result_dir / "report.txt", sep="\t", header=None, names=["pct","reads_cov","reads_direct","rank","taxid","name"])
    result = pd.read_csv(result_dir / "results.txt", sep = '\t', names = ["classified", "readid", "taxid", "length", "loc"],header = None)
    print("Done.")
    
    # clean ws
    report = report.apply(lambda col: col.str.strip() if col.dtype == object else col)
    result = result.apply(lambda col: col.str.strip() if col.dtype == object else col)
    # filter for classified reads
    result = result[result['classified'] == "C"]
    # make accession column
    result[['acc','idx']] = [x.split('.') for x in result['readid']]
    # assign taxon name for each read
    result["name"] = result['taxid'].map(
        {taxid: report.loc[report['taxid'] == taxid]["name"].iloc[0] 
         for taxid in report['taxid'].unique()}
    )
    
    return result    
