"""
parsing.py

Functions for reading FASTA / CDS FASTA files and turning them into
a simple in-memory representation that downstream analysis can use.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

DEFAULT_HEADER_SEPARATOR = "@"  # v1: fixed; later we make this configurable.


@dataclass
class CDSRecord:
    """
    Representation of a single CDS record from a FASTA file.

    Attributes
    ----------
    species : str
        Species name parsed from the FASTA header (before the separator).
    record_id : str
        Record identifier parsed from the FASTA header (after the separator).
    sequence : str
        Nucleotide sequence as DNA (A/C/G/T and possibly other symbols).
        We will handle T->U conversion and codon parsing in analysis.py.
    """
    species: str
    record_id: str
    sequence: str


def parse_header(header_line: str, separator: str = DEFAULT_HEADER_SEPARATOR) -> tuple[str, str]:
    """
    Parse a FASTA header into (species, record_id).

    The input header_line should include the leading '>' character.
    For example: '>gam2@c52334_g1_i1'.

    If the separator is not found, the entire header (minus '>') is used
    as both species and record_id.
    """
    if not header_line.startswith(">"):
        raise ValueError(f"Header line does not start with '>': {header_line!r}")

    text = header_line[1:].strip()

    if separator in text:
        species, record_id = text.split(separator, 1)
    else:
        species = text
        record_id = text

    return species.strip(), record_id.strip()


def iter_fasta_file(path: Path, separator: str = DEFAULT_HEADER_SEPARATOR) -> Iterable[CDSRecord]:
    """
    Yield CDSRecord objects from a single FASTA file.
    """
    current_species: str | None = None
    current_record_id: str | None = None
    seq_chunks: list[str] = []

    with path.open("r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()

            if not line:
                continue  # skip empty lines

            if line.startswith(">"):
                # Flush previous record if any
                if current_species is not None and seq_chunks:
                    sequence = "".join(seq_chunks)
                    yield CDSRecord(
                        species=current_species,
                        record_id=current_record_id or "",
                        sequence=sequence,
                    )

                # Start new record
                current_species, current_record_id = parse_header(line, separator=separator)
                seq_chunks = []
            else:
                # Sequence line
                seq_chunks.append(line)

        # Flush last record at EOF
        if current_species is not None and seq_chunks:
            sequence = "".join(seq_chunks)
            yield CDSRecord(
                species=current_species,
                record_id=current_record_id or "",
                sequence=sequence,
            )


def load_fasta_directory(input_dir: str | Path) -> list[CDSRecord]:
    """
    Load all .fasta and .cds.fasta files from a directory.

    Parameters
    ----------
    input_dir:
        Directory containing .fasta / .cds.fasta files.

    Returns
    -------
    list[CDSRecord]
        One CDSRecord per FASTA entry across all files.
    """
    input_path = Path(input_dir).expanduser().resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_path}")
    if not input_path.is_dir():
        raise NotADirectoryError(f"Input path is not a directory: {input_path}")

    # Collect both *.fasta and *.cds.fasta files, avoiding duplicates.
    fasta_patterns = ["*.fasta", "*.cds.fasta"]
    files: set[Path] = set()
    for pattern in fasta_patterns:
        files.update(input_path.glob(pattern))

    records: list[CDSRecord] = []
    for path in sorted(files):
        records.extend(iter_fasta_file(path))

    return records