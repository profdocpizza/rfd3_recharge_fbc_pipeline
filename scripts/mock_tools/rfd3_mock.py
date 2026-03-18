#!/usr/bin/env python3
import argparse
import json
from pathlib import Path


PDB_TEMPLATE = """\
ATOM      1  N   ALA A   1      11.104  13.207  12.011  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.104  13.207  12.011  1.00 20.00           C
ATOM      3  N   GLY A   2      13.104  13.207  12.011  1.00 20.00           N
ATOM      4  CA  GLY A   2      14.104  13.207  12.011  1.00 20.00           C
ATOM      5  N   SER B   1      21.104  23.207  22.011  1.00 20.00           N
ATOM      6  CA  SER B   1      22.104  23.207  22.011  1.00 20.00           C
ATOM      7  N   THR B   2      23.104  23.207  22.011  1.00 20.00           N
ATOM      8  CA  THR B   2      24.104  23.207  22.011  1.00 20.00           C
TER
END
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True)
    parser.add_argument("--hotspots", required=True)
    parser.add_argument("--config", default=None)
    parser.add_argument("--num-designs", type=int, default=3)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    designs = []

    for i in range(1, args.num_designs + 1):
        did = f"design_{i:03d}"
        pdb = outdir / f"{did}.pdb"
        pdb.write_text(PDB_TEMPLATE)
        designs.append({"design_id": did, "rfd3_score": round(0.5 + i * 0.1, 3)})

    (outdir / "score_metadata.json").write_text(json.dumps({"designs": designs}, indent=2))


if __name__ == "__main__":
    main()

