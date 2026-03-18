#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path
from typing import Dict, List


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def file_entries(paths: List[Path]) -> List[Dict[str, str]]:
    return [{"path": str(p), "sha256": sha256(p)} for p in paths if p.is_file()]


def main() -> None:
    parser = argparse.ArgumentParser(description="Create run-level manifest from stage outputs.")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    stage_dirs = sorted([p for p in args.results_dir.iterdir() if p.is_dir()])
    manifests = []
    for stage in stage_dirs:
      stage_files = [p for p in stage.rglob("*") if p.is_file()]
      manifests.append({
          "stage": stage.name,
          "files": file_entries(stage_files)
      })

    run_manifest = {
        "run_id": args.run_id,
        "results_dir": str(args.results_dir),
        "stages": manifests
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(run_manifest, indent=2))


if __name__ == "__main__":
    main()

