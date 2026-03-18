process RUN_MANIFEST {
  label 'small'
  container params.container_aggregate
  publishDir "${params.outdir}/manifests", mode: 'copy'
  publishDir "${params.outdir}/debug/manifests", mode: 'copy', pattern: "debug/*"

  input:
    path dashboard_csv
    path dashboard_json

  output:
    path "run_manifest.json", emit: run_manifest
    path "debug/*", emit: debug_files

  script:
  """
  mkdir -p debug
  python ${params.manifest_script} \\
    --run-id ${params.run_id} \\
    --results-dir ${projectDir}/${params.outdir} \\
    --output run_manifest.json > debug/run_manifest.log 2>&1
  """
}
