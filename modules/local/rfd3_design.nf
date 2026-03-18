process RFD3_DESIGN {
  label 'large'
  container params.container_rfd3
  publishDir "${params.outdir}/02_rfd3", mode: 'copy'
  publishDir "${params.outdir}/debug/02_rfd3", mode: 'copy', pattern: "debug/*"
  maxForks params.max_parallel_designs

  input:
    path trimmed_pdb
    path hotspot_metadata

  output:
    path "rfd3/*.pdb", emit: designed_complexes
    path "rfd3/score_metadata.json", emit: rfd3_scores
    path "debug/*", emit: debug_files

  script:
  """
  mkdir -p rfd3 debug
  RFD3_CONFIG_ARG=""
  if [[ -f "${params.rfd3_config}" && -s "${params.rfd3_config}" ]]; then
    RFD3_CONFIG_ARG="--config ${params.rfd3_config}"
  fi

  ${params.rfd3_cmd} \\
    --target ${trimmed_pdb} \\
    --hotspots ${hotspot_metadata} \\
    --binder-length "${params.binder_length}" \\
    --hotspot "${params.hotspot}" \\
    \${RFD3_CONFIG_ARG} \\
    --outdir rfd3 > debug/rfd3_design.log 2>&1
  """
}
