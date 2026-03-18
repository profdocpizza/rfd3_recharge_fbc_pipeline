#!/usr/bin/env python3
import argparse
import csv
import json
import shutil
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-complex", required=True)
    parser.add_argument("--config", default=None)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    src = Path(args.input_complex)
    did = src.stem.replace("_standardized", "")
    dst = outdir / f"{did}_recharged.pdb"
    shutil.copyfile(src, dst)

    report = {"designs": [{"design_id": did, "charge_net": 1.5}]}
    (outdir / "charge_report.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

