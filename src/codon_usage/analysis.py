"""
analysis.py

Codon counting and usage statistics.

Given a list of CDSRecord objects, this module will:
- build species-level composite CDS sequences,
- split sequences into codons,
- count codons per species,
- compute frequencies and percentages,
- return a tidy pandas DataFrame.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from typing import Dict, List

import pandas as pd

from .parsing import CDSRecord


# Canonical set of 64 RNA codons (U/C/A/G)
_NUCLEOTIDES = "UCAG"
ALL_CODONS: List[str] = [
    a + b + c for a in _NUCLEOTIDES for b in _NUCLEOTIDES for c in _NUCLEOTIDES
]


def _normalize_sequence_to_rna(seq: str) -> str:
    """
    Convert a DNA-like sequence to uppercase RNA (T->U).

    Any character outside A/C/G/T (or a/c/g/t) is left as-is for now;
    codons containing non-A/C/G/U characters will be ignored downstream.
    """
    seq = seq.upper()
    return seq.replace("T", "U")


def _count_codons_for_species(records: Iterable[CDSRecord]) -> Dict[str, Dict[str, int]]:
    """
    Aggregate codon counts per species from a collection of CDSRecord objects.

    Returns
    -------
    dict[str, dict[str, int]]
        Mapping: species -> { codon -> count }.
    """
    species_counts: Dict[str, Dict[str, int]] = {
        species: {codon: 0 for codon in ALL_CODONS} for species in set(r.species for r in records)
    }

    # We iterate again because we needed the set of species above.
    for rec in records:
        rna_seq = _normalize_sequence_to_rna(rec.sequence)
        # Ignore trailing nucleotides that don't form a full codon.
        limit = len(rna_seq) - (len(rna_seq) % 3)

        for i in range(0, limit, 3):
            codon = rna_seq[i:i+3]
            if len(codon) != 3:
                continue

            # Skip codons that contain gaps or ambiguous bases.
            if any(base not in _NUCLEOTIDES for base in codon):
                continue

            species_counts[rec.species][codon] += 1

    return species_counts


def compute_codon_usage(records: Iterable[CDSRecord]) -> pd.DataFrame:
    """
    Compute codon-usage statistics from a collection of CDSRecord objects.

    For v1, the result has one row per (species, codon) pair with columns:
    - species
    - codon
    - count
    - total_codons
    - freq_fraction
    - freq_percent
    """
    # Materialize records into a list so we can iterate multiple times.
    record_list = list(records)

    if not record_list:
        # Empty input: return an empty DataFrame with the expected columns.
        return pd.DataFrame(
            columns=["species", "codon", "count", "total_codons", "freq_fraction", "freq_percent"]
        )

    species_counts = _count_codons_for_species(record_list)

    rows: List[dict] = []
    for species, codon_counts in species_counts.items():
        total_codons = sum(codon_counts.values())
        if total_codons == 0:
            # No valid codons for this species; still emit rows with zero counts.
            for codon in ALL_CODONS:
                rows.append(
                    {
                        "species": species,
                        "codon": codon,
                        "count": 0,
                        "total_codons": 0,
                        "freq_fraction": 0.0,
                        "freq_percent": 0.0,
                    }
                )
            continue

        for codon in ALL_CODONS:
            count = codon_counts.get(codon, 0)
            freq_fraction = count / total_codons if total_codons > 0 else 0.0
            freq_percent = round(freq_fraction * 100.0, 3)
            rows.append(
                {
                    "species": species,
                    "codon": codon,
                    "count": count,
                    "total_codons": total_codons,
                    "freq_fraction": freq_fraction,
                    "freq_percent": freq_percent,
                }
            )

    df = pd.DataFrame(rows, columns=["species", "codon", "count", "total_codons", "freq_fraction", "freq_percent"])
    return df