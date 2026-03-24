#!/usr/bin/env python3
import argparse
import json
import csv
import shutil
import subprocess
import time
import sys
import os
from pathlib import Path


def run(cmd, cwd=None, log_path: Path | None = None, attempt_log_path: Path | None = None):
    """
    Run command and persist combined stdout/stderr to a log file while streaming
    to the wrapper stdout so Nextflow-level logs still receive the same output.
    """
    env = dict(**os.environ)
    env["PYTHONUNBUFFERED"] = "1"

    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)
    if attempt_log_path:
        attempt_log_path.parent.mkdir(parents=True, exist_ok=True)
    with (open(log_path, "ab", buffering=0) if log_path else open(os.devnull, "wb")) as log_f:
        proc = subprocess.Popen(
            cmd,
            cwd=cwd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=0,
            env=env,
        )
        assert proc.stdout is not None
        with (open(attempt_log_path, "ab", buffering=0) if attempt_log_path else open(os.devnull, "wb")) as attempt_f:
            while True:
                chunk = proc.stdout.read(4096)
                if not chunk:
                    break
                if log_path:
                    log_f.write(chunk)
                    log_f.flush()
                if attempt_log_path:
                    attempt_f.write(chunk)
                    attempt_f.flush()
                sys.stdout.buffer.write(chunk)
                sys.stdout.buffer.flush()

        rc = proc.wait()
        if rc != 0:
            msg = f"\n[freebindcraft_wrapper] Command failed with exit code {rc}: {' '.join(cmd)}\n"
            if log_path:
                log_f.write(msg.encode("utf-8", errors="replace"))
                log_f.flush()
            sys.stdout.write(msg)
            sys.stdout.flush()
            raise subprocess.CalledProcessError(rc, cmd)

def parse_bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def run_with_retries(
    cmd,
    cwd=None,
    log_path: Path | None = None,
    max_retries: int = 20,
    retry_delay_seconds: int = 30,
):
    attempts = max_retries + 1
    for attempt in range(1, attempts + 1):
        attempt_log = None
        if log_path:
            attempt_log = log_path.parent / f"freebindcraft_filtering_attempt_{attempt}.tmp.log"
            attempt_log.write_text("")
        try:
            run(cmd, cwd=cwd, log_path=log_path, attempt_log_path=attempt_log)
            if attempt_log and attempt_log.exists():
                attempt_log.unlink(missing_ok=True)
            return
        except subprocess.CalledProcessError:
            if log_path and attempt_log and attempt_log.exists():
                crash_copy = log_path.parent / f"freebindcraft_filtering_{attempt}.log"
                shutil.copyfile(attempt_log, crash_copy)
                attempt_log.unlink(missing_ok=True)
            if attempt >= attempts:
                raise
            msg = (
                f"[freebindcraft_wrapper] attempt {attempt}/{attempts} failed; "
                f"retrying in {retry_delay_seconds}s\n"
            )
            if log_path:
                with open(log_path, "ab", buffering=0) as log_f:
                    log_f.write(msg.encode("utf-8", errors="replace"))
                    log_f.flush()
            sys.stdout.write(msg)
            sys.stdout.flush()
            time.sleep(retry_delay_seconds)


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


def design_id_from_complex(path: Path) -> str:
    return path.stem.replace("_recharge_reference", "").replace("_standardized", "")


def design_id_from_sequence(path: Path) -> str:
    return path.stem.replace("_recharged", "")


def collect_design_pairs(args):
    if args.input_complex_dir:
        complex_files = sorted(Path(args.input_complex_dir).glob("*_standardized.pdb"))
    elif args.input_complex:
        complex_files = [Path(args.input_complex)]
    else:
        complex_files = []

    if args.recharged_sequence_dir:
        seq_files = sorted(Path(args.recharged_sequence_dir).glob("*_recharged.fasta"))
    elif args.recharged_sequence:
        seq_files = [Path(args.recharged_sequence)]
    else:
        seq_files = []

    seq_by_id = {p.stem.replace("_recharged", ""): p for p in seq_files}
    pairs = []
    for complex_path in complex_files:
        design_id = design_id_from_complex(complex_path)
        seq_path = seq_by_id.get(design_id)
        if seq_path:
            pairs.append((complex_path.resolve(), seq_path.resolve(), design_id))
    if not pairs:
        raise RuntimeError("No matched *_standardized.pdb / *_recharged.fasta pairs found for filtering.")
    return pairs


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
    p.add_argument("--input-complex", required=False, default=None)
    p.add_argument("--recharged-sequence", required=False, default=None)
    p.add_argument("--input-complex-dir", required=False, default=None)
    p.add_argument("--recharged-sequence-dir", required=False, default=None)
    p.add_argument("--hotspot-pdb", required=False, default=None)
    p.add_argument("--binder-length", required=True)
    p.add_argument("--hotspot", required=False, default="")
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
    p.add_argument("--max-retries", type=int, default=20)
    p.add_argument("--retry-delay-seconds", type=int, default=30)
    p.add_argument("--reuse", default="true")
    args = p.parse_args()
    if args.max_retries < 0:
        raise ValueError("--max-retries must be >= 0.")
    if args.retry_delay_seconds < 0:
        raise ValueError("--retry-delay-seconds must be >= 0.")
    reuse = parse_bool(args.reuse)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    repo = Path("/home/tadas/code/FreeBindCraft")
    design_id = "batch"
    runtime = outdir / "runtime"
    runtime.mkdir(parents=True, exist_ok=True)
    debug_dir = outdir / "debug"
    debug_dir.mkdir(parents=True, exist_ok=True)

    pairs = collect_design_pairs(args)
    first_complex = pairs[0][0]
    hotspot_pdb = Path(args.hotspot_pdb).resolve() if args.hotspot_pdb else first_complex

    target_chains = parse_chain_ids(hotspot_pdb)
    target_chain_spec = ",".join(target_chains) if target_chains else "A"
    binder_len = max(8, count_residues(first_complex, "B"))
    hotspot_for_fbc = "" if args.mode == "filter-only" else normalize_hotspot_tokens(args.hotspot)
    settings = json.loads(Path(args.settings).read_text()) if Path(args.settings).exists() else {}
    settings["design_path"] = str(outdir.resolve())
    settings["binder_name"] = design_id
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
    (debug_dir / "effective_settings.json").write_text(json.dumps(settings, indent=2))

    steps = "filtering"
    relaxed = outdir / "Trajectory" / "Relaxed"
    relaxed.mkdir(parents=True, exist_ok=True)
    mpnn_seq_dir = outdir / "MPNN" / "Sequences"
    mpnn_seq_dir.mkdir(parents=True, exist_ok=True)

    prepared = []
    for i, (complex_path, seq_path, pair_id) in enumerate(pairs):
        shutil.copyfile(seq_path, debug_dir / seq_path.name)
        stamp = (int(time.time() * 1000) + i) % 100000
        pair_len = max(8, count_residues(complex_path, "B"))
        suffix = f"_l{pair_len}_s{stamp:05d}"
        design_stem = f"{pair_id}{suffix}"
        relaxed_pdb = relaxed / f"{design_stem}.pdb"
        shutil.copyfile(complex_path, relaxed_pdb)
        seq = parse_first_fasta_sequence(seq_path)
        # FreeBindCraft filtering expects <trajectory_stem>.fasta for each Relaxed pdb.
        seq_out = mpnn_seq_dir / f"{design_stem}.fasta"
        seq_out.write_text(f">{design_stem}\n{seq}\n")
        # Keep legacy alias for compatibility with old readers.
        seq_out_mpnn1 = mpnn_seq_dir / f"{design_stem}_mpnn1.fasta"
        seq_out_mpnn1.write_text(f">{design_stem}_mpnn1\n{seq}\n")
        prepared.append(
            {
                "design_id": pair_id,
                "complex": str(complex_path),
                "sequence": str(seq_path),
                "design_stem": design_stem,
                "fbc_sequence_file": str(seq_out),
                "stamp": f"{stamp:05d}",
                "binder_length": pair_len,
            }
        )

    (debug_dir / "prepared_pairs.json").write_text(json.dumps(prepared, indent=2))

    cmd = [
        "python",
        "-u",
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
    if reuse:
        cmd.append("--reuse")
    (debug_dir / "freebindcraft_runtime_args.json").write_text(
        json.dumps(
            {
                "mode": args.mode,
                "work_dir": str(outdir.resolve()),
                "bindcraft_script": str((repo / "bindcraft.py").resolve()),
                "starting_pdb": settings.get("starting_pdb"),
                "chains": settings.get("chains"),
                "binder_name": settings.get("binder_name"),
                "target_hotspot_residues": settings.get("target_hotspot_residues"),
                "lengths": settings.get("lengths"),
                "settings_file": str(settings_path.resolve()),
                "filters_file": str(Path(args.filters).resolve()),
                "advanced_file": str(Path(args.advanced).resolve()),
                "max_retries": args.max_retries,
                "retry_delay_seconds": args.retry_delay_seconds,
                "reuse": reuse,
                "command": cmd,
            },
            indent=2,
        )
    )
    pipeline_live_log = outdir.parent / "debug" / "freebindcraft_filtering.log"
    with open(pipeline_live_log, "ab", buffering=0) as log_f:
        log_f.write(
            b"[freebindcraft_wrapper] Starting FreeBindCraft program execution.\n"
        )
        log_f.flush()

    run_with_retries(
        cmd,
        cwd=outdir,
        log_path=pipeline_live_log,
        max_retries=args.max_retries,
        retry_delay_seconds=args.retry_delay_seconds,
    )

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
