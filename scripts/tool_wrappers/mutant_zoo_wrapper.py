#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

from mutant_zoo import MutantZoo


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-annotations", required=True)
    p.add_argument("--input-complex", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    ann = outdir / "annotated_structures"
    ann.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(open(args.input_annotations, newline="")))
    zoo = MutantZoo(seed=42)
    zoo.add_file("input_annotations", str(Path(args.input_annotations).resolve()))

    designs = []
    input_complex = Path(args.input_complex)
    has_pass_columns = any(("filter_pass" in r or "near_pass" in r) for r in rows)
    for row in rows:
        did = row.get("design_id") or row.get("name") or row.get("id")
        if not did:
            continue
        if has_pass_columns:
            keep = str(row.get("filter_pass", "")).lower() == "true" or str(row.get("near_pass", "")).lower() == "true"
        else:
            keep = True
        if not keep:
            continue
        designs.append(
            {
                "design_id": did,
                "mutation_count": 1,
                "filter_pass": row.get("filter_pass", "true" if not has_pass_columns else "false"),
                "near_pass": row.get("near_pass", "true" if not has_pass_columns else "false"),
            }
        )
        if input_complex.exists():
            (ann / f"{did}.pdb").write_text(input_complex.read_text())
        else:
            (ann / f"{did}.pdb").write_text("END\n")

    zoo_yaml = outdir / "zoo.yaml"
    zoo.save_class_state(str(zoo_yaml))
    (outdir / "mutation_sets.json").write_text(json.dumps({"designs": designs, "zoo_yaml": str(zoo_yaml)}, indent=2))


if __name__ == "__main__":
    main()
