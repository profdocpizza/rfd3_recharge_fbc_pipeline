#!/usr/bin/env python3
import argparse
import csv
import shutil
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["full", "filter-only"], required=True)
    parser.add_argument("--input-complex", required=True)
    parser.add_argument("--pass1-scores", default=None)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    did = Path(args.input_complex).stem.replace("_recharged", "")

    if args.mode == "full":
        fc = outdir / "filtered_candidates"
        fc.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(args.input_complex, fc / f"{did}.pdb")

        with (outdir / "final_design_stats.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["design_id", "fbc_pass1_score"])
            w.writeheader()
            w.writerow({"design_id": did, "fbc_pass1_score": 0.82})

        (outdir / "designed_sequences.fasta").write_text(f">{did}\nMSEQUENCE\n")
        return

    with (outdir / "revalidation_scores.csv").open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["design_id", "fbc_revalidation_score", "fbc_revalidation_pass"],
        )
        w.writeheader()
        w.writerow(
            {
                "design_id": did,
                "fbc_revalidation_score": 0.79,
                "fbc_revalidation_pass": "true",
            }
        )


if __name__ == "__main__":
    main()

