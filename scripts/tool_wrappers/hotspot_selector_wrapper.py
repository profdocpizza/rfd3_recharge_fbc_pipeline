#!/usr/bin/env python3
import argparse
import json
import subprocess
from pathlib import Path


def run(cmd):
    subprocess.run(cmd, check=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input_pdb")
    p.add_argument("--sasa-threshold", required=True)
    p.add_argument("--support-dist", required=True)
    p.add_argument("--hotspot", default="")
    p.add_argument("--surface-radius", required=True)
    p.add_argument("--graph-step", required=True)
    p.add_argument("--design-spec", default=None)
    p.add_argument("--outdir", required=True)
    p.add_argument(
        "--hotspot-script",
        default="/home/tadas/code/hotspot_selector/hotspot_selector.py",
    )
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python",
        args.hotspot_script,
        args.input_pdb,
        "--sasa-threshold",
        str(args.sasa_threshold),
        "--support-dist",
        str(args.support_dist),
        "--surface-radius",
        str(args.surface_radius),
        "--graph-step",
        str(args.graph_step),
        "--output-dir",
        str(outdir),
    ]
    if args.hotspot:
        cmd += ["--anchor", *args.hotspot.split()]
    run(cmd)

    input_stem = Path(args.input_pdb).stem
    hotspot_pdb = outdir / f"{input_stem}_hotspot.pdb"
    if not hotspot_pdb.exists():
        raise FileNotFoundError(f"Expected hotspot output not found: {hotspot_pdb}")

    hotspot_index_txt = outdir / f"{input_stem}_hotspot_residue_index.txt"
    target_contig = None
    if hotspot_index_txt.exists():
        target_contig = hotspot_index_txt.read_text().strip()

    standardized = outdir / "trimmed_target_001.pdb"
    standardized.write_text(hotspot_pdb.read_text())

    meta = {
        "source": str(hotspot_pdb),
        "anchors": args.hotspot.split() if args.hotspot else [],
        "hotspot_residue_index_file": str(hotspot_index_txt) if hotspot_index_txt.exists() else None,
        "target_contig": target_contig,
    }
    (outdir / "hotspot_metadata.json").write_text(json.dumps(meta, indent=2))

    numbering = {"residue_map": {"B": {}}}
    (outdir / "target_numbering_map.json").write_text(json.dumps(numbering, indent=2))


if __name__ == "__main__":
    main()
