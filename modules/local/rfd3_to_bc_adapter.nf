process RFD3_TO_BC_ADAPTER {
  label 'small'
  container params.container_adapter
  publishDir "${params.outdir}/03_adapter", mode: 'copy'
  publishDir "${params.outdir}/debug/03_adapter", mode: 'copy', pattern: "debug/*"

  input:
    tuple path(rfd3_complex), path(target_numbering_map)

  output:
    path "adapter/*_standardized.pdb", emit: standardized_complexes
    path "adapter/*_chain_mapping.json", emit: chain_map_metadata
    path "debug/*", emit: debug_files

  script:
  """
  mkdir -p adapter debug
  base_name=\$(basename "${rfd3_complex}" .pdb)
  python ${params.adapter_script} \\
    --input-pdb ${rfd3_complex} \\
    --target-numbering-map ${target_numbering_map} \\
    --output-pdb adapter/\${base_name}_standardized.pdb \\
    --output-chain-map adapter/\${base_name}_chain_mapping.json > debug/rfd3_to_bc_adapter.log 2>&1
  """
}
