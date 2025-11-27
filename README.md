# codon-usage

A lightweight tool for batch codon-usage analysis of CDS FASTA files. It reads all `.fasta` and `.cds.fasta` files in an input directory, groups sequences by species, builds per-species composite CDS sequences, counts codons, and outputs a single tidy CSV ready for downstream analysis.

---

## Overview

Each FASTA header is interpreted in the form:
```
>species_name@record_id
```

- Everything before the separator is treated as the species name.
- Everything after the separator is treated as the record identifier.
- Future versions will allow customizing the separator or disabling it entirely.

All sequences belonging to the same species are concatenated to form that species’ composite CDS. Codons are counted on this composite sequence, and usage frequencies are calculated from those counts.

The tool supports DNA FASTA files (A/C/G/T). Internally all T characters are converted to U before codon parsing so codons are counted in the RNA alphabet.

---

## Input directory

Point the tool at any directory containing `.fasta` / `.cds.fasta` files.

This repository includes an example_data/ folder with sample files. These are for demonstration only; place your own files (or copies of the examples) into any directory you wish to analyze.

Example structure:
```
example_data/
    cog04046.cds.fasta
    cog04047.cds.fasta
```

---

## Output format

Version 1 produces a single CSV containing species-level composite codon usage. Each row corresponds to one (species, codon) pair.

Columns:

- species        — species name
- codon          — codon in RNA alphabet (A/C/G/U)
- count          — number of occurrences of that codon in the species composite CDS
- total_codons   — total number of codons in the species composite CDS
- freq_fraction  — count / total_codons
- freq_percent   — (count / total_codons) * 100 (typically rounded to three decimals)

This long-format table is immediately usable with pandas, R/tidyverse and most statistical workflows.

---

## Installation (development)

Create and activate a development environment, then install in editable mode.

Example commands:

conda create -n ankhbio-env python=3.11
conda activate ankhbio-env
pip install -e .

---

## Basic usage

The package provides a simple API. A typical usage pattern is:

```python

from codon_usage import analyze_directory

df = analyze_directory("path/to/input_dir")
df.to_csv("codon_usage.csv", index=False)
```

analyze_directory:

- reads all `.fasta` / `.cds.fasta` files from the given input directory,
- parses headers into species_name and record_id,
- builds species-level composite CDS sequences,
- computes codon counts and frequencies,
- returns a DataFrame in the CSV format described above.

---

## Project layout
```
codon-usage/
    src/
        codon_usage/
            __init__.py
            parsing.py
            analysis.py
    example_data/
    notebooks/
        01_explore_codon_usage.ipynb
    README.md
    pyproject.toml
```
---

## Roadmap and TODO

Initial release (v1.0) includes:

- Species-level composite CDS per species.
- Codon counts plus fractional and percent usage.
- Single tidy CSV output.
- T to U conversion for codon parsing.

Planned for v1.x:

- Configurable header separator and header modes (species + record id, species-only, raw header).
- Optional per-record codon-usage CSV with the same structure as the species-level CSV plus record_id.
- Ability to toggle composite and record-level outputs independently.

Planned for v2:

- GC content and GC3 per species.
- Relative synonymous codon usage (RSCU).
- Per–amino-acid codon usage tables.

Planned for v3:

- Effective number of codons (ENC).
- Codon adaptation index (CAI), using a user-defined reference set of highly expressed genes.

---

## License

MIT License. Refer to the LICENSE file for details.
