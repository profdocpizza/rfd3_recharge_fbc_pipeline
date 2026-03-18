process FREEBINDCRAFT_FILTERING {
  label 'large'
  container params.container_fbc
  publishDir "${params.outdir}/05_fbc_filtering", mode: 'copy'
  publishDir "${params.outdir}/debug/05_fbc_filtering", mode: 'copy', pattern: "debug/*"

  input:
    tuple path(input_complex), path(recharged_sequence), path(hotspot_pdb)

  output:
    path "fbc_filtering/*_fbc_filter_scores.csv", emit: fbc_filter_scores_tagged, optional: true
    path "fbc_filtering/MPNN/**", emit: mpnn_tree, optional: true
    path "debug/MPNN/**", emit: debug_mpnn_tree, optional: true
    path "debug/*", emit: debug_files

  script:
  """
  mkdir -p fbc_filtering debug
  ${params.freebindcraft_cmd} \\
    --mode filter-only \\
    --input-complex ${input_complex} \\
    --recharged-sequence ${recharged_sequence} \\
    --hotspot-pdb ${hotspot_pdb} \\
    --binder-length "${params.binder_length}" \\
    --hotspot "${params.hotspot}" \\
    --settings "${params.fbc_settings}" \\
    --filters "${params.fbc_filters}" \\
    --advanced "${params.fbc_advanced}" \\
    --outdir fbc_filtering > debug/freebindcraft_filtering.log 2>&1
  """
}
