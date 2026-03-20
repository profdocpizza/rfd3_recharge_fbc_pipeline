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
  EFFECTIVE_TRIM_HOTSPOT="${params.hotspot}"
  if [[ -n "${params.hotspot_trimming_residues}" ]]; then
    EFFECTIVE_TRIM_HOTSPOT="${params.hotspot_trimming_residues}"
  fi
  if [[ -n "${params.design_spec}" && -f "${params.design_spec}" ]]; then
    DESIGN_SPEC_ARG="--design-spec ${params.design_spec}"
  fi
  MAX_RESIDUES_ARG=""
  if [[ -n "${params.hotspot_trimming_max_residues}" && "${params.hotspot_trimming_max_residues}" != "null" && "${params.hotspot_trimming_max_residues}" != "None" ]]; then
    MAX_RESIDUES_ARG="--max_residues ${params.hotspot_trimming_max_residues}"
  fi
  cat > debug/run_args.json <<EOF
  {
    "step": "HOTSPOT_TRIM",
    "work_dir": "\${PWD}",
    "target_pdb": "${target_pdb}",
    "sasa_threshold": "${params.hotspot_sasa_threshold}",
    "support_dist": "${params.hotspot_support_dist}",
    "hotspot": "${params.hotspot}",
    "hotspot_trimming_residues": "${params.hotspot_trimming_residues}",
    "effective_trim_hotspot": "\${EFFECTIVE_TRIM_HOTSPOT}",
    "surface_radius": "${params.hotspot_surface_radius}",
    "graph_step": "${params.hotspot_graph_step}",
    "hotspot_trimming_max_residues": "${params.hotspot_trimming_max_residues}",
    "max_residues_arg": "\${MAX_RESIDUES_ARG}",
    "design_spec_arg": "\${DESIGN_SPEC_ARG}",
    "hotspot_selector_cmd": "${params.hotspot_selector_cmd}",
    "outdir": "trimmed"
  }
EOF

  ${params.hotspot_selector_cmd} \\
    ${target_pdb} \\
    --sasa-threshold ${params.hotspot_sasa_threshold} \\
    --support-dist ${params.hotspot_support_dist} \\
    --hotspot "\${EFFECTIVE_TRIM_HOTSPOT}" \\
    --surface-radius ${params.hotspot_surface_radius} \\
    --graph-step ${params.hotspot_graph_step} \\
    \${MAX_RESIDUES_ARG} \\
    \${DESIGN_SPEC_ARG} \\
    --outdir trimmed > debug/hotspot_trim.log 2>&1
  """
}
