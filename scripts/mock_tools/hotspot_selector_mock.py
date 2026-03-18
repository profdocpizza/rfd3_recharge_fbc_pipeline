#!/usr/bin/env python3
import argparse
import json
import shutil
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_pdb")
    parser.add_argument("--sasa-threshold")
    parser.add_argument("--support-dist")
    parser.add_argument("--anchor")
    parser.add_argument("--surface-radius")
    parser.add_argument("--graph-step")
    parser.add_argument("--design-spec", default=None)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    trimmed_pdb = outdir / "trimmed_target_001.pdb"
    shutil.copyfile(args.input_pdb, trimmed_pdb)

    hotspot = {
        "target_id": Path(args.input_pdb).stem,
        "anchors": args.anchor.split() if args.anchor else [],
        "params": {
            "sasa_threshold": args.sasa_threshold,
            "support_dist": args.support_dist,
            "surface_radius": args.surface_radius,
            "graph_step": args.graph_step,
        },
    }
    (outdir / "hotspot_metadata.json").write_text(json.dumps(hotspot, indent=2))

    residue_map = {"residue_map": {"B": {str(i): str(100 + i) for i in range(1, 6)}}}
    (outdir / "target_numbering_map.json").write_text(json.dumps(residue_map, indent=2))


if __name__ == "__main__":
    main()

