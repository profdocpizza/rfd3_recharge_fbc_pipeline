process DASHBOARD_DATA {
  label 'small'
  container params.container_aggregate
  publishDir "${params.outdir}/09_dashboard", mode: 'copy'
  publishDir "${params.outdir}/debug/09_dashboard", mode: 'copy', pattern: "debug/*"

  input:
    path aggregate_csv
    path aggregate_json

  output:
    path "dashboard/candidates.csv", emit: dashboard_csv
    path "dashboard/candidates.json", emit: dashboard_json
    path "debug/*", emit: debug_files

  script:
  """
  mkdir -p dashboard debug
  cat > debug/run_args.json <<EOF
  {
    "step": "DASHBOARD_DATA",
    "work_dir": "\${PWD}",
    "aggregate_csv": "${aggregate_csv}",
    "aggregate_json": "${aggregate_json}",
    "outputs": ["dashboard/candidates.csv", "dashboard/candidates.json"]
  }
EOF
  cp ${aggregate_csv} dashboard/candidates.csv
  cp ${aggregate_json} dashboard/candidates.json
  echo "dashboard export complete" > debug/dashboard_data.log
  """
}
