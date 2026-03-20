#!/usr/bin/env python3
import argparse
import json
import shutil
import subprocess
from pathlib import Path

import yaml


def run(cmd, cwd=None):
    subprocess.run(cmd, check=True, cwd=cwd)


def can_import_pandas_in_env(env_name: str) -> bool:
    try:
        subprocess.run(
            ["conda", "run", "-n", env_name, "python", "-c", "import pandas"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return True
    except subprocess.CalledProcessError:
        return False


def pick_recharge_env(preferred_env: str) -> str | None:
    candidates = []
    if preferred_env:
        candidates.append(preferred_env)
    for fallback in ["recharge_env", "tadas311", "ligandmpnn_env"]:
        if fallback not in candidates:
            candidates.append(fallback)
    for env_name in candidates:
        if can_import_pandas_in_env(env_name):
            return env_name
    return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-complex", required=True)
    p.add_argument("--config", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    did = Path(args.input_complex).stem.replace("_standardized", "")

    cfg = yaml.safe_load(Path(args.config).read_text()) if Path(args.config).exists() else {}
    cfg["inputs"] = [{"input_pdb_path": str(Path(args.input_complex).resolve())}]
    cfg["out_folder"] = str(outdir.resolve())
    cfg.setdefault("redesign_chain", "B")
    cfg.setdefault("ligandmpnn_run_py", "/home/tadas/code/LigandMPNN/run.py")
    cfg.setdefault("ligandmpnn_conda_env", "ligandmpnn_env")
    cfg.setdefault("protein_recharge_conda_env", "recharge_env")
    desired_charge = cfg.get("desired_chain_charge", None)

    generated_cfg = outdir / f"{did}_recharge_runtime.yaml"
    generated_cfg.write_text(yaml.safe_dump(cfg, sort_keys=False))

    repo = Path("/home/tadas/code/ProteinRecharge")
    preferred_env = str(cfg.get("protein_recharge_conda_env") or "").strip()
    recharge_env = pick_recharge_env(preferred_env)
    if not recharge_env:
        raise RuntimeError(
            "No conda environment with pandas found for ProteinRecharge. "
            f"Tried: {preferred_env or '(none)'}, recharge_env, tadas311, ligandmpnn_env."
        )

    run(
        ["conda", "run", "-n", recharge_env, "python", "recharge.py", str(generated_cfg.resolve())],
        cwd=repo,
    )

    selection = outdir / "selection"
    fasta_files = sorted(selection.glob(f"{did}*_recharged.fasta")) or sorted(selection.glob("*_recharged.fasta"))
    if not fasta_files:
        raise RuntimeError("ProteinRecharge finished but no recharged FASTA found in selection/.")
    fasta_src = max(fasta_files, key=lambda p: p.stat().st_mtime)
    fasta_dst = outdir / f"{did}_recharged.fasta"
    shutil.copyfile(fasta_src, fasta_dst)

    check_files = sorted(selection.glob("*_check_redesigned.cif"))
    debug_dir = outdir / "debug"
    debug_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(generated_cfg, debug_dir / generated_cfg.name)
    (debug_dir / "protein_recharge_env.json").write_text(
        json.dumps(
            {
                "preferred_env": preferred_env,
                "selected_env": recharge_env,
            },
            indent=2,
        )
    )
    for f in check_files:
        shutil.copyfile(f, debug_dir / f.name)

    charge_val = None
    header = fasta_src.read_text().splitlines()[0] if fasta_src.exists() else ""
    if "_charge_" in header:
        try:
            charge_val = float(header.split("_charge_")[-1].strip())
        except ValueError:
            charge_val = None

    charge = {
        "designs": [
            {
                "design_id": did,
                "charge_net": charge_val,
                "desired_charge": desired_charge,
                "sequence_fasta": str(fasta_dst),
                "input_complex": str(Path(args.input_complex).resolve()),
                "runtime_config": str(generated_cfg),
            }
        ]
    }
    (outdir / "charge_report.json").write_text(json.dumps(charge, indent=2))
    (outdir / f"{did}_charge_report.json").write_text(json.dumps(charge, indent=2))


if __name__ == "__main__":
    main()
