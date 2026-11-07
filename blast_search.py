#!/usr/bin/env python3
"""Fetch a public RefSeq mRNA and run a remote NCBI BLASTN search."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML


def fetch_fasta(accession: str, output: Path):
    """Fetch a nucleotide accession from NCBI and save it as FASTA."""
    email = os.environ.get("NCBI_EMAIL")
    if not email:
        raise ValueError("Set NCBI_EMAIL to your email address before querying NCBI.")
    Entrez.email = email
    with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
        record = SeqIO.read(handle, "fasta")
    output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(record, output, "fasta")
    return record


def run_blast(record, xml_output: Path, database: str = "nt", hitlist_size: int = 10):
    """Submit BLASTN to NCBI and save the returned XML report."""
    xml_output.parent.mkdir(parents=True, exist_ok=True)
    result_handle = NCBIWWW.qblast(
        program="blastn",
        database=database,
        sequence=record.format("fasta"),
        expect=10.0,
        hitlist_size=hitlist_size,
        megablast=True,
        format_type="XML",
    )
    try:
        xml_output.write_text(result_handle.read(), encoding="utf-8")
    finally:
        result_handle.close()


def parse_blast(xml_path: Path, csv_path: Path, query_length: int):
    """Convert BLAST alignments into an annotation table."""
    rows = []
    with xml_path.open(encoding="utf-8") as handle:
        for result in NCBIXML.parse(handle):
            for alignment in result.alignments:
                for hsp in alignment.hsps:
                    identity = (hsp.identities / hsp.align_length * 100) if hsp.align_length else 0
                    coverage = (hsp.align_length / query_length * 100) if query_length else 0
                    rows.append({
                        "subject_accession": alignment.accession,
                        "subject_title": alignment.title,
                        "alignment_length": hsp.align_length,
                        "percent_identity": round(identity, 2),
                        "query_coverage_approx_percent": round(min(coverage, 100), 2),
                        "e_value": hsp.expect,
                        "bit_score": round(hsp.bits, 2),
                    })
    rows.sort(key=lambda row: (row["e_value"], -row["bit_score"]))
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "subject_accession", "subject_title", "alignment_length",
        "percent_identity", "query_coverage_approx_percent", "e_value", "bit_score",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", default="NM_001101.5", help="NCBI nucleotide accession")
    parser.add_argument("--database", default="nt", help="NCBI BLAST database (default: nt)")
    parser.add_argument("--hits", type=int, default=10, help="Maximum number of hits to retain")
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    args = parser.parse_args()

    fasta_path = Path("data/query_sequence.fasta")
    xml_path = args.output_dir / "blast_results.xml"
    csv_path = args.output_dir / "annotation_table.csv"

    print(f"Fetching {args.accession} from NCBI Nucleotide...")
    record = fetch_fasta(args.accession, fasta_path)
    print(f"Fetched {record.id}: {len(record.seq)} nucleotides.")
    print("Submitting BLASTN search to NCBI. This may take several minutes...")
    run_blast(record, xml_path, args.database, args.hits)
    rows = parse_blast(xml_path, csv_path, len(record.seq))
    print(f"Saved {len(rows)} alignment row(s) to {csv_path}")
    print("Interpret the accession, identity, coverage, and E-value using the README.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
