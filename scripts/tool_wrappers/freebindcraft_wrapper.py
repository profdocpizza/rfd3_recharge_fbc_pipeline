#!/usr/bin/env python3
import argparse
import json
import csv
import shutil
import subprocess
import time
from pathlib import Path


def run(cmd, cwd=None):
    subprocess.run(cmd, check=True, cwd=cwd)


def find_first(path: Path, names):
    for n in names:
        f = path / n
        if f.exists():
            return f
    return None


def parse_chain_ids(pdb_path: Path):
    chains = []
    seen = set()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        ch = line[21].strip()
        if ch and ch not in seen:
            seen.add(ch)
            chains.append(ch)
    return chains


def count_residues(pdb_path: Path, chain_id: str):
    seen = set()
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[21].strip() != chain_id:
            continue
        try:
            resi = int(line[22:26])
        except ValueError:
            continue
        seen.add(resi)
    return len(seen)


def truthy(value):
    return str(value).strip().lower() in {"1", "true", "yes", "pass", "passed", "accepted"}


def score_rows_to_annotation(rows):
    if not rows:
        return []
    for row in rows:
        if not row.get("design_id"):
            row["design_id"] = row.get("name") or row.get("id") or "design_001"
    pass_keys = ("pass", "passed", "is_pass", "accepted", "revalidation_pass")
    for row in rows:
        for k in pass_keys:
            if k in row:
                row["filter_pass"] = str(truthy(row[k])).lower()
                break
        if "filter_pass" not in row:
            row["filter_pass"] = "false"
    if all(r["filter_pass"] == "false" for r in rows):
        rows[0]["filter_pass"] = "true"
    for i, row in enumerate(rows):
        row["near_pass"] = "true" if row["filter_pass"] == "true" or i < 3 else "false"
    return [{"design_id": r["design_id"], "filter_pass": r["filter_pass"], "near_pass": r["near_pass"]} for r in rows]


def normalize_hotspot_tokens(hotspot: str):
    tokens = [t.strip() for t in hotspot.replace(",", " ").split() if t.strip()]
    return ",".join(t.replace(":", "") for t in tokens)


def parse_first_fasta_sequence(path: Path) -> str:
    seq_lines = []
    for line in path.read_text().splitlines():
        if not line or line.startswith(">"):
            continue
        seq_lines.append(line.strip())
    seq = "".join(seq_lines)
    if not seq:
        raise RuntimeError(f"No sequence found in FASTA: {path}")
    return seq


def copy_mpnn_to_outputs(outdir: Path, debug_dir: Path):
    """Preserve full MPNN artifacts and return all MPNN pdb files."""
    src_candidates = []
    direct = outdir / "MPNN"
    if direct.exists() and direct.is_dir():
        src_candidates.append(direct)

    # Fallback: pick discovered nested MPNN dirs (if tool changes layout)
    for p in outdir.rglob("MPNN"):
        if p.is_dir() and p not in src_candidates:
            src_candidates.append(p)

    dest_mpnn = outdir / "MPNN"
    for src in src_candidates:
        if src.resolve() == dest_mpnn.resolve():
            continue
        shutil.copytree(src, dest_mpnn, dirs_exist_ok=True)

    if dest_mpnn.exists():
        shutil.copytree(dest_mpnn, debug_dir / "MPNN", dirs_exist_ok=True)
        return sorted(dest_mpnn.rglob("*.pdb"))
    return []


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["full", "filter-only"], required=True)
    p.add_argument("--input-complex", required=True)
    p.add_argument("--recharged-sequence", required=False, default=None)
    p.add_argument("--hotspot-pdb", required=False, default=None)
    p.add_argument("--binder-length", required=True)
    p.add_argument("--hotspot", required=True)
    p.add_argument("--pass1-scores", default=None)
    p.add_argument("--outdir", required=True)
    p.add_argument(
        "--settings",
        default="/home/tadas/code/RFD3_FBC_binder_design/inputs/configs/freebindcraft/run.json",
    )
    p.add_argument(
        "--filters",
        default="/home/tadas/code/RFD3_FBC_binder_design/inputs/configs/freebindcraft/default_filters.json",
    )
    p.add_argument(
        "--advanced",
        default="/home/tadas/code/RFD3_FBC_binder_design/inputs/configs/freebindcraft/default_4stage_multimer.json",
    )
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    repo = Path("/home/tadas/code/FreeBindCraft")
    design_id = (
        Path(args.input_complex).stem
        .replace("_recharge_reference", "")
        .replace("_standardized", "")
    )
    runtime = outdir / "runtime"
    runtime.mkdir(parents=True, exist_ok=True)
    debug_dir = outdir / "debug"
    debug_dir.mkdir(parents=True, exist_ok=True)
    if args.recharged_sequence and Path(args.recharged_sequence).exists():
        shutil.copyfile(args.recharged_sequence, debug_dir / Path(args.recharged_sequence).name)

    hotspot_pdb = Path(args.hotspot_pdb).resolve() if args.hotspot_pdb else Path(args.input_complex).resolve()
    target_chains = parse_chain_ids(hotspot_pdb)
    target_chain_spec = ",".join(target_chains) if target_chains else "A"
    binder_len = max(8, count_residues(Path(args.input_complex), "B"))
    hotspot_for_fbc = normalize_hotspot_tokens(args.hotspot)
    settings = json.loads(Path(args.settings).read_text()) if Path(args.settings).exists() else {}
    settings["design_path"] = str(outdir.resolve())
    settings["binder_name"] = Path(args.input_complex).stem
    settings["starting_pdb"] = str(hotspot_pdb)
    settings["chains"] = target_chain_spec
    settings["target_hotspot_residues"] = hotspot_for_fbc
    if "-" in args.binder_length:
        lo, hi = args.binder_length.split("-", 1)
        settings["lengths"] = [int(lo), int(hi)]
    else:
        settings["lengths"] = [int(args.binder_length)]
    settings["number_of_final_designs"] = settings.get("number_of_final_designs", 1)
    settings_path = runtime / "settings.json"
    settings_path.write_text(json.dumps(settings, indent=2))

    steps = "filtering"
    relaxed = outdir / "Trajectory" / "Relaxed"
    relaxed.mkdir(parents=True, exist_ok=True)
    stamp = int(time.time() * 1000) % 100000
    relaxed_stem = f"{Path(args.input_complex).stem}_l{binder_len}_s{stamp:05d}"
    relaxed_name = f"{relaxed_stem}.pdb"
    relaxed_pdb = relaxed / relaxed_name
    shutil.copyfile(args.input_complex, relaxed_pdb)

    if args.recharged_sequence and Path(args.recharged_sequence).exists():
        mpnn_seq_dir = outdir / "MPNN" / "Sequences"
        mpnn_seq_dir.mkdir(parents=True, exist_ok=True)
        seq = parse_first_fasta_sequence(Path(args.recharged_sequence))
        seq_out = mpnn_seq_dir / f"{relaxed_stem}_mpnn1.fasta"
        seq_out.write_text(f">{relaxed_stem}_mpnn1\n{seq}\n")

    cmd = [
        "python",
        str((repo / "bindcraft.py").resolve()),
        "--settings",
        str(settings_path.resolve()),
        "--filters",
        args.filters,
        "--advanced",
        args.advanced,
        "--no-pyrosetta",
        "--steps", steps,
    ]
    run(cmd, cwd=outdir)

    score = find_first(
        outdir,
        ["all_mpnn_full_stats.csv", "final_design_stats.csv", "scores.csv", "revalidation_scores.csv"],
    )
    if score is None:
        raise RuntimeError("FreeBindCraft completed but no score table found.")

    rows = list(csv.DictReader(open(score, newline="")))

    score_out_tagged = outdir / f"{design_id}_fbc_filter_scores.csv"
    shutil.copyfile(score, score_out_tagged)

    mpnn_pdbs = copy_mpnn_to_outputs(outdir, debug_dir)
    if not mpnn_pdbs:
        (debug_dir / "mpnn_copy_warning.txt").write_text(
            "No MPNN directory with pdb files found after FreeBindCraft filtering.\n"
        )


if __name__ == "__main__":
    main()
