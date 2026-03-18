process HOTSPOT_TRIM {
  label 'medium'
  container params.container_hotspot
  publishDir "${params.outdir}/01_hotspot_trim", mode: 'copy'
  publishDir "${params.outdir}/debug/01_hotspot_trim", mode: 'copy', pattern: "debug/*"

  input:
    path target_pdb

  output:
    path "trimmed/trimmed_target_001.pdb", emit: trimmed_pdbs
    path "trimmed/hotspot_metadata.json", emit: hotspot_metadata
    path "trimmed/*_hotspot_residue_index.txt", emit: hotspot_contig_ranges, optional: true
    path "trimmed/target_numbering_map.json", emit: target_numbering_map
    path "debug/*", emit: debug_files

  script:
  """
  mkdir -p trimmed debug
  DESIGN_SPEC_ARG=""
  if [[ -n "${params.design_spec}" && -f "${params.design_spec}" ]]; then
    DESIGN_SPEC_ARG="--design-spec ${params.design_spec}"
  fi

  ${params.hotspot_selector_cmd} \\
    ${target_pdb} \\
    --sasa-threshold ${params.hotspot_sasa_threshold} \\
    --support-dist ${params.hotspot_support_dist} \\
    --hotspot "${params.hotspot}" \\
    --surface-radius ${params.hotspot_surface_radius} \\
    --graph-step ${params.hotspot_graph_step} \\
    \${DESIGN_SPEC_ARG} \\
    --outdir trimmed > debug/hotspot_trim.log 2>&1
  """
}
