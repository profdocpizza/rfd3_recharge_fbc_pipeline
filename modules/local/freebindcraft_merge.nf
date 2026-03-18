process FREEBINDCRAFT_MERGE {
  label 'small'
  container params.container_aggregate
  publishDir "${params.outdir}/05_fbc_filtering", mode: 'copy'

  input:
    path tagged_scores

  output:
    path "merged/fbc_filter_scores.csv", emit: fbc_filter_scores

  script:
  """
  mkdir -p merged
  python - <<'PY'
import csv
from pathlib import Path

def merge_csvs(files, out_path):
    rows = []
    field_order = []
    seen = set()
    for f in files:
        with open(f, newline="") as fh:
            reader = csv.DictReader(fh)
            for fn in (reader.fieldnames or []):
                if fn not in seen:
                    seen.add(fn)
                    field_order.append(fn)
            rows.extend(list(reader))
    if not field_order:
        field_order = ["design_id"]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=field_order)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in field_order})

score_files = sorted(Path(".").glob("*_fbc_filter_scores.csv"))
merge_csvs(score_files, Path("merged/fbc_filter_scores.csv"))
PY
  """
}
