process FREEBINDCRAFT_FILTERING {
  label 'large'
  container params.container_fbc
  publishDir "${params.outdir}/05_fbc_filtering", mode: 'copy'
  publishDir "${params.outdir}/debug/05_fbc_filtering", mode: 'copy', pattern: "debug/*"

  input:
    tuple path(input_complexes), path(recharged_sequences), path(hotspot_pdb)

  output:
    path "fbc_filtering/*_fbc_filter_scores.csv", emit: fbc_filter_scores_tagged, optional: true
    path "fbc_filtering/MPNN/**", emit: mpnn_tree, optional: true
    path "debug/MPNN/**", emit: debug_mpnn_tree, optional: true
    path "debug/*", emit: debug_files

  script:
  """
  mkdir -p fbc_filtering debug
  cat > debug/run_args.json <<EOF
  {
    "step": "FREEBINDCRAFT_FILTERING",
    "work_dir": "\${PWD}",
    "input_complexes": "${input_complexes}",
    "recharged_sequences": "${recharged_sequences}",
    "hotspot_pdb": "${hotspot_pdb}",
    "binder_length": "${params.binder_length}",
    "hotspot": "${params.hotspot}",
    "settings": "${params.fbc_settings}",
    "filters": "${params.fbc_filters}",
    "advanced": "${params.fbc_advanced}",
    "freebindcraft_cmd": "${params.freebindcraft_cmd}",
    "outdir": "fbc_filtering"
  }
EOF
  ${params.freebindcraft_cmd} \\
    --mode filter-only \\
    --input-complex-dir . \\
    --recharged-sequence-dir . \\
    --hotspot-pdb ${hotspot_pdb} \\
    --binder-length "${params.binder_length}" \\
    --hotspot "" \\
    --max-retries ${params.fbc_max_retries} \\
    --retry-delay-seconds ${params.fbc_retry_delay_seconds} \\
    --reuse "${params.fbc_reuse}" \\
    --settings "${params.fbc_settings}" \\
    --filters "${params.fbc_filters}" \\
    --advanced "${params.fbc_advanced}" \\
    --outdir fbc_filtering

  if [[ -d fbc_filtering/debug ]]; then
    cp -r fbc_filtering/debug/* debug/ || true
  fi
  """
}
