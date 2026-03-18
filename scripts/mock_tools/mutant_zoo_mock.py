#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-annotations", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    ann_dir = outdir / "annotated_structures"
    ann_dir.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(open(args.input_annotations, newline="")))
    designs = []
    for row in rows:
        did = row.get("design_id")
        if not did:
            continue
        designs.append({"design_id": did, "mutation_count": 3})
        (ann_dir / f"{did}.pdb").write_text("END\n")

    (outdir / "mutation_sets.json").write_text(json.dumps({"designs": designs}, indent=2))


if __name__ == "__main__":
    main()

