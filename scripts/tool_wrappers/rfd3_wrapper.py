#!/usr/bin/env python3
import argparse
import json
import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path

import yaml
from Bio.PDB import MMCIFParser, PDBIO


def run(cmd):
    subprocess.run(cmd, check=True)


def collect_residues(pdb_path: Path):
    residues = {}
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21]
        try:
            resi = int(line[22:26])
        except ValueError:
            continue
        residues.setdefault(chain, set()).add(resi)
    return {k: sorted(v) for k, v in residues.items()}


def contiguous_segments(vals):
    if not vals:
        return []
    segs = []
    start = prev = vals[0]
    for v in vals[1:]:
        if v == prev + 1:
            prev = v
            continue
        segs.append((start, prev))
        start = prev = v
    segs.append((start, prev))
    return segs


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--target", required=True)
    p.add_argument("--hotspots", required=True)
    p.add_argument("--config", default=None)
    p.add_argument("--outdir", required=True)
    p.add_argument("--num-designs", type=int, default=3)
    p.add_argument("--binder-length", required=True)
    p.add_argument("--hotspot", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    run_dir = outdir / "rfd3_raw"
    run_dir.mkdir(parents=True, exist_ok=True)

    if shutil.which("rfd3") is None:
        raise RuntimeError(
            "rfd3 CLI not found in PATH. Install rc-foundry[rfd3] and checkpoints."
        )

    if not args.config:
        raise ValueError("RFD3 wrapper requires --config YAML path.")

    config_in = yaml.safe_load(Path(args.config).read_text())
    if not isinstance(config_in, dict) or not config_in:
        raise ValueError("RFD3 config must contain at least one design entry.")
    first_key = next(iter(config_in.keys()))
    entry = config_in[first_key] or {}
    entry["input"] = str(Path(args.target).resolve())

    hotspot_meta = json.loads(Path(args.hotspots).read_text())
    target_contig = (hotspot_meta.get("target_contig") or "").strip()
    if not target_contig:
        index_file = hotspot_meta.get("hotspot_residue_index_file")
        if index_file and Path(index_file).exists():
            target_contig = Path(index_file).read_text().strip()

    if target_contig:
        # Use hotspot selector's exact compact range format for the target contig.
        entry["contig"] = f"{args.binder_length},/0,{target_contig}"
        anchors = args.hotspot.split() if args.hotspot else hotspot_meta.get("anchors") or []
        if anchors:
            entry["select_hotspots"] = {
                a.replace(":", ""): "CA,CB" for a in anchors if ":" in a
            }
    else:
        residues = collect_residues(Path(args.target))
        best_chain = None
        best_seg = None
        for chain, vals in residues.items():
            for seg in contiguous_segments(vals):
                if seg[1] - seg[0] + 1 >= 20:
                    best_chain = chain
                    best_seg = seg
                    break
            if best_chain:
                break
        if not best_chain or not best_seg:
            raise RuntimeError("Could not find a contiguous 20-residue target segment for RFD3.")

        seg_start = best_seg[0]
        seg_end = min(best_seg[1], seg_start + 19)
        hs1 = seg_start + (seg_end - seg_start) // 3
        hs2 = seg_start + 2 * (seg_end - seg_start) // 3
        entry["contig"] = f"{args.binder_length},/0,{best_chain}{seg_start}-{seg_end}"
        entry["select_hotspots"] = {
            f"{best_chain}{hs1}": "CA,CB",
            f"{best_chain}{hs2}": "CA,CB",
        }

    normalized = {first_key: entry}
    run_cfg = outdir / "rfd3_runtime.yaml"
    run_cfg.write_text(yaml.safe_dump(normalized, sort_keys=False))

    cmd = [
        "rfd3",
        "design",
        f"out_dir={run_dir}",
        f"inputs={run_cfg}",
        "skip_existing=False",
        "prevalidate_inputs=True",
        "dump_trajectories=False",
    ]
    run(cmd)

    designs = []
    produced = sorted(run_dir.rglob("*.pdb"))
    if not produced:
        produced = sorted(run_dir.rglob("*.cif")) + sorted(run_dir.rglob("*.cif.gz"))
    if not produced:
        raise RuntimeError("RFD3 completed but no structure outputs found.")

    for idx, f in enumerate(produced[: args.num_designs], start=1):
        did = f"design_{idx:03d}"
        dst = outdir / f"{did}.pdb"
        if f.suffix == ".pdb":
            shutil.copyfile(f, dst)
        elif f.suffix == ".cif":
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure(did, str(f))
            io = PDBIO()
            io.set_structure(structure)
            io.save(str(dst))
        elif f.suffix == ".gz" and f.name.endswith(".cif.gz"):
            with tempfile.NamedTemporaryFile(suffix=".cif", delete=False) as tmp:
                with gzip.open(f, "rt") as fin:
                    tmp.write(fin.read().encode())
                tmp_path = Path(tmp.name)
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure(did, str(tmp_path))
            io = PDBIO()
            io.set_structure(structure)
            io.save(str(dst))
            tmp_path.unlink(missing_ok=True)
        else:
            continue
        designs.append({"design_id": did, "path": str(dst), "rfd3_score": None})

    (outdir / "score_metadata.json").write_text(json.dumps({"designs": designs}, indent=2))


if __name__ == "__main__":
    main()
