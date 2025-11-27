from __future__ import annotations

from pathlib import Path

import pandas as pd

from codon_usage import analyze_directory


def test_basic_pipeline(tmp_path: Path) -> None:
    """
    Minimal end-to-end test:

    - Create a temporary directory.
    - Write one small FASTA file with two species.
    - Run analyze_directory on that directory.
    - Check that codon counts and total_codons behave as expected.
    """

    fasta_content = """>gam2@rec1
ATGCCGACTT
>gam5@rec2
ATGCCGACTT
"""

    fasta_path = tmp_path / "tiny_test.cds.fasta"
    fasta_path.write_text(fasta_content, encoding="utf-8")

    # Run the main analysis on the temp directory.
    df = analyze_directory(tmp_path)

    # We expect two species: gam2 and gam5.
    species_set = set(df["species"])
    assert species_set == {"gam2", "gam5"}

    # Each sequence is 10 nt -> 3 codons (ATG CCG ACT), 1 base leftover ignored.
    # So per species: 3 codons total.
    per_species_totals = df.groupby("species")["count"].sum().to_dict()
    assert per_species_totals["gam2"] == 3
    assert per_species_totals["gam5"] == 3

    # Check that total_codons column matches 3 for each species-row.
    for species in species_set:
        sub = df[df["species"] == species]
        # All rows for that species should have the same total_codons value.
        assert sub["total_codons"].nunique() == 1
        assert int(sub["total_codons"].iloc[0]) == 3

    # Check that freq_fraction sums to ~1 per species.
    for species in species_set:
        sub = df[df["species"] == species]
        freq_sum = sub["freq_fraction"].sum()
        assert abs(freq_sum - 1.0) < 1e-6