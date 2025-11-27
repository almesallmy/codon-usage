"""
codon_usage
===========

Public API for the codon-usage package.

For v1, the main entry point is `analyze_directory`, which:

- reads all .fasta and .cds.fasta files from a directory,
- parses headers into species_name and record_id,
- builds species-level composite CDS sequences,
- computes codon counts and frequencies,
- returns a tidy pandas DataFrame ready to write to CSV.
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

from .parsing import load_fasta_directory
from .analysis import compute_codon_usage


__all__ = ["analyze_directory"]

def analyze_directory(input_dir: str | Path) -> pd.DataFrame:
    """
    Analyze all FASTA/CDS FASTA files in the given directory.

    Parameters
    ----------
    input_dir:
        Path to a directory containing .fasta and/or .cds.fasta files.

    Returns
    -------
    pandas.DataFrame
        Tidy table with one row per (species, codon) pair and columns:
        - species
        - codon
        - count
        - total_codons
        - freq_fraction
        - freq_percent
    """
    input_path = Path(input_dir).expanduser().resolve()
    records = load_fasta_directory(input_path)
    df = compute_codon_usage(records)
    return df