#!/usr/bin/env python3
import argparse
import json
import gzip
import shutil
import subprocess
import tempfile
import time
import sys
from pathlib import Path
from collections import defaultdict

import yaml
from Bio.PDB import MMCIFParser, PDBIO


# Residue type → important atom names (for hotspot selection constraints)
RESIDUE_ATOM_MAP = {
    "ARG": "NH1,NH2",  # Arginine - guanidinium group
    "ASP": "OD1,OD2",  # Aspartate - carboxyl oxygens
    "ASN": "OD1,ND2",  # Asparagine - amide group
    "GLU": "OE1,OE2",  # Glutamate - carboxyl oxygens
    "GLN": "OE1,NE2",  # Glutamine - amide group
    "HIS": "ND1,NE2",  # Histidine - imidazole nitrogens
    "LYS": "NZ",       # Lysine - amino group
    "SER": "OG",       # Serine - hydroxyl oxygen
    "THR": "OG1",      # Threonine - hydroxyl oxygen
    "TYR": "OH",       # Tyrosine - phenol oxygen
    "CYS": "SG",       # Cysteine - thiol sulfur
    "TRP": "NE1",      # Tryptophan - indole nitrogen
}


def parse_bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def run(cmd, max_retries=0, retry_delay_seconds=30, log_path: Path | None = None):
    attempts = max_retries + 1
    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, attempts + 1):
        attempt_log = None
        if log_path:
            attempt_log = log_path.parent / f"rfd3_design_attempt_{attempt}.tmp.log"
            attempt_log.write_text("")
        try:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                bufsize=0,
            )
            assert proc.stdout is not None
            with open(log_path, "ab", buffering=0) if log_path else open("/dev/null", "wb") as live_f:
                with open(attempt_log, "ab", buffering=0) if attempt_log else open("/dev/null", "wb") as attempt_f:
                    while True:
                        chunk = proc.stdout.read(4096)
                        if not chunk:
                            break
                        live_f.write(chunk)
                        live_f.flush()
                        attempt_f.write(chunk)
                        attempt_f.flush()
                        sys.stdout.buffer.write(chunk)
                        sys.stdout.buffer.flush()
            rc = proc.wait()
            if rc != 0:
                raise subprocess.CalledProcessError(rc, cmd)
            if attempt_log and attempt_log.exists():
                attempt_log.unlink(missing_ok=True)
            return
        except subprocess.CalledProcessError as e:
            if log_path and attempt_log and attempt_log.exists():
                crash_copy = log_path.parent / f"rfd3_design_{attempt}.log"
                shutil.copyfile(attempt_log, crash_copy)
                attempt_log.unlink(missing_ok=True)
            if attempt >= attempts:
                raise e
            msg = (
                f"[rfd3_wrapper] attempt {attempt}/{attempts} failed; "
                f"retrying in {retry_delay_seconds}s\n"
            )
            print(msg, flush=True)
            if log_path:
                with open(log_path, "a") as f:
                    f.write(msg)
            time.sleep(retry_delay_seconds)

def write_json(path: Path, payload: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def collect_residues(pdb_path: Path):
    """Collect residues and their types from PDB."""
    residues = {}
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21]
        try:
            resi = int(line[22:26])
        except ValueError:
            continue
        resname = line[17:20].strip()
        key = (chain, resi)
        if key not in residues:
            residues[key] = {"chain": chain, "resi": resi, "resname": resname}
    
    # Group by chain for backwards compatibility
    by_chain = defaultdict(set)
    for (chain, resi), info in residues.items():
        by_chain[chain].add(resi)
    
    return {k: sorted(v) for k, v in by_chain.items()}, residues


def get_atoms_for_residue(resname: str) -> str:
    """Get important atoms for a residue type, default to CA if not found."""
    return RESIDUE_ATOM_MAP.get(resname, "CA")


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
    p.add_argument("--num-designs", type=int, default=2)
    p.add_argument("--rfd3-env", default="rfd3")
    p.add_argument("--binder-length", required=True)
    p.add_argument("--hotspot", required=True)
    p.add_argument("--max-retries", type=int, default=20)
    p.add_argument("--retry-delay-seconds", type=int, default=30)
    p.add_argument("--skip-existing", default="true")
    args = p.parse_args()
    if args.num_designs <= 0:
        raise ValueError("--num-designs must be a positive integer.")
    if args.max_retries < 0:
        raise ValueError("--max-retries must be >= 0.")
    if args.retry_delay_seconds < 0:
        raise ValueError("--retry-delay-seconds must be >= 0.")
    diffusion_batch_size = 1
    n_batches = args.num_designs
    skip_existing = parse_bool(args.skip_existing)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    run_dir = outdir / "rfd3_raw"
    run_dir.mkdir(parents=True, exist_ok=True)
    process_debug_dir = outdir.parent / "debug"
    process_debug_dir.mkdir(parents=True, exist_ok=True)

    if not args.config:
        raise ValueError("RFD3 wrapper requires --config YAML path.")

    cfg_path = Path(args.config)
    cfg_exists = cfg_path.exists()
    cfg_size = cfg_path.stat().st_size if cfg_exists else None

    write_json(process_debug_dir / "rfd3_wrapper_invocation.json", {
        "target": str(Path(args.target).resolve()),
        "hotspots": str(Path(args.hotspots).resolve()),
        "config": str(cfg_path.resolve()),
        "config_exists": cfg_exists,
        "config_size_bytes": cfg_size,
        "outdir": str(outdir.resolve()),
        "num_designs": args.num_designs,
        "n_batches": n_batches,
        "diffusion_batch_size": diffusion_batch_size,
        "rfd3_env": args.rfd3_env,
        "binder_length": args.binder_length,
        "hotspot": args.hotspot,
        "max_retries": args.max_retries,
        "retry_delay_seconds": args.retry_delay_seconds,
        "skip_existing": skip_existing,
    })

    config_in = yaml.safe_load(cfg_path.read_text())
    if not isinstance(config_in, dict) or not config_in:
        raise ValueError("RFD3 config must contain at least one design entry.")
    first_key = next(iter(config_in.keys()))
    entry = config_in[first_key] or {}
    entry["input"] = str(Path(args.target).resolve())
    entry.setdefault("is_non_loopy", True)
    entry.setdefault("infer_ori_strategy", "hotspots")
    entry.setdefault("dialect", 2)

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
            # Load residue info to map residue types to important atoms
            _, res_info = collect_residues(Path(args.target))
            select_hotspots = {}
            for anchor in anchors:
                if ":" in anchor:
                    chain, resi_str = anchor.split(":")
                    try:
                        resi = int(resi_str)
                        key = (chain, resi)
                        if key in res_info:
                            resname = res_info[key]["resname"]
                            atoms = get_atoms_for_residue(resname)
                            select_hotspots[anchor.replace(":", "")] = atoms
                        else:
                            # Fallback if residue not found
                            select_hotspots[anchor.replace(":", "")] = "CA"
                    except ValueError:
                        pass
            if select_hotspots:
                entry["select_hotspots"] = select_hotspots
    else:
        by_chain, res_info = collect_residues(Path(args.target))
        best_chain = None
        best_seg = None
        for chain, vals in by_chain.items():
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
        
        # Get residue types for automatic atom selection
        hs1_key = (best_chain, hs1)
        hs2_key = (best_chain, hs2)
        hs1_atoms = get_atoms_for_residue(res_info[hs1_key]["resname"]) if hs1_key in res_info else "CA"
        hs2_atoms = get_atoms_for_residue(res_info[hs2_key]["resname"]) if hs2_key in res_info else "CA"
        
        entry["select_hotspots"] = {
            f"{best_chain}{hs1}": hs1_atoms,
            f"{best_chain}{hs2}": hs2_atoms,
        }

    normalized = {first_key: entry}
    run_cfg = outdir / "rfd3_runtime.yaml"
    run_cfg.write_text(yaml.safe_dump(normalized, sort_keys=False))
    write_json(process_debug_dir / "rfd3_wrapper_pre_run.json", {
        "runtime_config_path": str(run_cfg.resolve()),
        "resolved_design_key": first_key,
        "resolved_design_entry": entry,
        "hotspot_meta_path": str(Path(args.hotspots).resolve()),
        "hotspot_meta": hotspot_meta,
        "target_contig": target_contig,
    })

    cmd = [
        "conda",
        "run",
        "--no-capture-output",
        "-n",
        args.rfd3_env,
        "env",
        "PYTHONUNBUFFERED=1",
        "rfd3",
        "design",
        f"out_dir={run_dir}",
        f"inputs={run_cfg}",
        f"n_batches={n_batches}",
        f"diffusion_batch_size={diffusion_batch_size}",
        f"skip_existing={'True' if skip_existing else 'False'}",
        "prevalidate_inputs=True",
        "dump_trajectories=False",
    ]
    write_json(process_debug_dir / "rfd3_design_command.json", {
        "command": cmd
    })
    run(
        cmd,
        max_retries=args.max_retries,
        retry_delay_seconds=args.retry_delay_seconds,
        log_path=process_debug_dir / "rfd3_retry.log",
    )

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
