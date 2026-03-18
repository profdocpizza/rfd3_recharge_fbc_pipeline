process FREEBINDCRAFT_REVALIDATE {
  label 'medium'
  container params.container_fbc
  publishDir "${params.outdir}/06_fbc_revalidate", mode: 'copy'

  input:
    path recharged_complex
    path pass1_scores

  output:
    path "fbc_revalidate/revalidation_scores.csv", emit: revalidation_scores
    path "fbc_revalidate/pass_fail_annotation.csv", emit: revalidated_annotations

  script:
  """
  mkdir -p fbc_revalidate
  ${params.freebindcraft_cmd} \\
    --mode filter-only \\
    --input-complex ${recharged_complex} \\
    --pass1-scores ${pass1_scores} \\
    --outdir fbc_revalidate

  SCORE_FILE=""
  for CANDIDATE in revalidation_scores.csv final_design_stats.csv all_mpnn_full_stats.csv scores.csv; do
    if [[ -f "fbc_revalidate/\${CANDIDATE}" ]]; then
      SCORE_FILE="\${CANDIDATE}"
      break
    fi
  done

  if [[ -z "\${SCORE_FILE}" ]]; then
    echo "No revalidation score table found in fbc_revalidate output." >&2
    exit 1
  fi

  if [[ "\${SCORE_FILE}" != "revalidation_scores.csv" ]]; then
    cp "fbc_revalidate/\${SCORE_FILE}" "fbc_revalidate/revalidation_scores.csv"
  fi

  if [[ ! -f "fbc_revalidate/pass_fail_annotation.csv" ]]; then
    python - <<'PY'
import csv
from pathlib import Path

pass1 = list(csv.DictReader(open("${pass1_scores}", newline="")))
reval = list(csv.DictReader(open("fbc_revalidate/revalidation_scores.csv", newline="")))

def ids(rows):
    keys = ("design_id", "name", "id")
    out = set()
    for row in rows:
        for k in keys:
            if row.get(k):
                out.add(row[k])
                break
    return out

pass1_ids = ids(pass1)
reval_ids = ids(reval)
rows = [{"design_id": did, "revalidation_pass": str(did in reval_ids).lower()} for did in sorted(pass1_ids)]

with open("fbc_revalidate/pass_fail_annotation.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["design_id", "revalidation_pass"])
    w.writeheader()
    w.writerows(rows)
PY
  fi
  """
}
