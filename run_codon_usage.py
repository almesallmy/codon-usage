from __future__ import annotations

import sys
import time
from pathlib import Path

from codon_usage import analyze_directory


def main(argv: list[str] | None = None) -> None:
    """
    Simple command-line entry point for codon-usage v1.

    Usage:
        python run_codon_usage.py <input_dir> [<output_csv>]

    - <input_dir>  : directory containing .fasta / .cds.fasta files
    - <output_csv> : optional output path (default: codon_usage.csv in CWD)
    """
    if argv is None:
        argv = sys.argv[1:]

    if not argv:
        print("Usage: python run_codon_usage.py <input_dir> [<output_csv>]")
        raise SystemExit(1)

    input_dir = Path(argv[0]).expanduser().resolve()
    if len(argv) >= 2:
        output_csv = Path(argv[1]).expanduser().resolve()
    else:
        output_csv = Path("codon_usage.csv").resolve()

    if not input_dir.exists() or not input_dir.is_dir():
        print(f"Error: input_dir does not exist or is not a directory: {input_dir}")
        raise SystemExit(1)

    print(f"[codon-usage] Reading FASTA files from: {input_dir}")

    start = time.perf_counter()
    df = analyze_directory(input_dir)
    elapsed_analysis = time.perf_counter() - start

    print(f"[codon-usage] Parsed {len(df)} rows of codon-usage data "
          f"in {elapsed_analysis:.2f} seconds.")

    start_write = time.perf_counter()
    df.to_csv(output_csv, index=False)
    elapsed_write = time.perf_counter() - start_write

    total_time = elapsed_analysis + elapsed_write

    print(f"[codon-usage] Wrote CSV to: {output_csv}")
    print(f"[codon-usage] Write time: {elapsed_write:.2f} seconds.")
    print(f"[codon-usage] Total time (analysis + write): {total_time:.2f} seconds.")


if __name__ == "__main__":
    main()