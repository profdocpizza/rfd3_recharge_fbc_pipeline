#!/usr/bin/env python3
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path


def run(cmd):
    subprocess.run(cmd, check=True)


def parse_optional_int(value):
    if value is None:
        return None
    if isinstance(value, int):
        return value
    text = str(value).strip()
    if text == "" or text.lower() in {"null", "none"}:
        return None
    try:
        return int(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid --max_residues value '{value}'; expected integer or null/None."
        ) from exc


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input_pdb")
    p.add_argument("--sasa-threshold", required=True)
    p.add_argument("--support-dist", required=True)
    p.add_argument("--hotspot", default="")
    p.add_argument("--max_residues", dest="max_residues", default=None)
    p.add_argument("--surface-radius", required=True)
    p.add_argument("--graph-step", required=True)
    p.add_argument("--design-spec", default=None)
    p.add_argument("--outdir", required=True)
    p.add_argument(
        "--hotspot-script",
        default="/home/tadas/code/hotspot_selector/hotspot_selector.py",
    )
    args = p.parse_args()
    args.max_residues = parse_optional_int(args.max_residues)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    debug_dir = outdir.parent / "debug"
    debug_dir.mkdir(parents=True, exist_ok=True)

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
    if args.max_residues is not None:
        cmd.extend(["--max_residues", str(args.max_residues)])
    attempted_cmd = None
    used_anchor = False
    fallback_no_anchor = False
    anchors = args.hotspot.split() if args.hotspot else []
    if anchors:
        cmd_with_anchor = [*cmd, "--anchor", *anchors]
        attempted_cmd = cmd_with_anchor
        used_anchor = True
        try:
            run(cmd_with_anchor)
        except subprocess.CalledProcessError as exc:
            print(
                "hotspot_selector failed with anchors; retrying without --anchor. "
                f"anchors={anchors} returncode={exc.returncode}",
                file=sys.stderr,
            )
            fallback_no_anchor = True
            attempted_cmd = cmd
            run(cmd)
    else:
        attempted_cmd = cmd
        run(cmd)

    input_stem = Path(args.input_pdb).stem
    hotspot_pdb = outdir / f"{input_stem}_hotspot.pdb"
    if not hotspot_pdb.exists():
        raise FileNotFoundError(f"Expected hotspot output not found: {hotspot_pdb}")

    hotspot_index_txt = outdir / f"{input_stem}_hotspot_residue_index.txt"
    annotated_pdb = outdir / f"{input_stem}_annotated.pdb"
    target_contig = None
    if hotspot_index_txt.exists():
        target_contig = hotspot_index_txt.read_text().strip()

    standardized = outdir / "trimmed_target_001.pdb"
    standardized.write_text(hotspot_pdb.read_text())

    meta = {
        "source": str(hotspot_pdb),
        "anchors": anchors,
        "hotspot_residue_index_file": str(hotspot_index_txt) if hotspot_index_txt.exists() else None,
        "target_contig": target_contig,
    }
    (outdir / "hotspot_metadata.json").write_text(json.dumps(meta, indent=2))

    numbering = {"residue_map": {"B": {}}}
    (outdir / "target_numbering_map.json").write_text(json.dumps(numbering, indent=2))

    # Persist hotspot run configuration/args for debugging and reproducibility.
    run_args = {
        "input_pdb": args.input_pdb,
        "sasa_threshold": args.sasa_threshold,
        "support_dist": args.support_dist,
        "hotspot": args.hotspot,
        "max_residues": args.max_residues,
        "surface_radius": args.surface_radius,
        "graph_step": args.graph_step,
        "design_spec": args.design_spec,
        "outdir": str(outdir),
        "hotspot_script": args.hotspot_script,
        "anchors": anchors,
        "used_anchor": used_anchor,
        "fallback_no_anchor": fallback_no_anchor,
        "executed_command": attempted_cmd,
    }
    (debug_dir / "hotspot_run_args.json").write_text(json.dumps(run_args, indent=2))

    if annotated_pdb.exists():
        shutil.copy2(annotated_pdb, debug_dir / annotated_pdb.name)


if __name__ == "__main__":
    main()
