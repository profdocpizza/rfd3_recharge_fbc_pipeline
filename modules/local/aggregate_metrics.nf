process AGGREGATE_METRICS {
  label 'small'
  container params.container_aggregate
  publishDir "${params.outdir}/08_aggregate", mode: 'copy'
  publishDir "${params.outdir}/debug/08_aggregate", mode: 'copy', pattern: "debug/*"

  input:
    path rfd3_scores
    path fbc_filter_scores
    path charge_report
    path mutation_sets

  output:
    path "aggregate/candidates.csv", emit: aggregate_csv
    path "aggregate/candidates.json", emit: aggregate_json
    path "debug/*", emit: debug_files

  script:
  """
  mkdir -p aggregate debug
  python - <<'PY'
import csv
import json
from pathlib import Path

def load_json(path):
    return json.loads(Path(path).read_text())

def load_csv(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))

rfd3 = load_json("${rfd3_scores}")
fbc = load_csv("${fbc_filter_scores}")
charge = load_json("${charge_report}")
mut = load_json("${mutation_sets}")

rows = []
charge_by_id = {r.get("design_id"): r for r in charge.get("designs", [])}
mut_by_id = {r.get("design_id"): r for r in mut.get("designs", [])}
rfd3_by_id = {r.get("design_id"): r for r in rfd3.get("designs", [])}

for row in fbc:
    did = row.get("design_id")
    merged = {}
    merged.update(rfd3_by_id.get(did, {}))
    merged.update(row)
    merged.update({f"charge_{k}": v for k, v in charge_by_id.get(did, {}).items()})
    merged.update({f"mut_{k}": v for k, v in mut_by_id.get(did, {}).items()})
    merged["design_id"] = did
    rows.append(merged)

if rows:
    fields = sorted({k for r in rows for k in r.keys()})
    with open("aggregate/candidates.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
else:
    with open("aggregate/candidates.csv", "w", newline="") as f:
        f.write("design_id\\n")

Path("aggregate/candidates.json").write_text(json.dumps({"candidates": rows}, indent=2))
PY
  echo "aggregate complete" > debug/aggregate_metrics.log
  """
}
