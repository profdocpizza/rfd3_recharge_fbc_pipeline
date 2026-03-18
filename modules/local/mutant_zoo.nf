process MUTANT_ZOO {
  label 'medium'
  container params.container_mutant_zoo
  publishDir "${params.outdir}/07_mutant_zoo", mode: 'copy'
  publishDir "${params.outdir}/debug/07_mutant_zoo", mode: 'copy', pattern: "debug/*"

  input:
    path pass_fail_annotations
    path recharged_complex

  output:
    path "mutant_zoo/mutation_sets.json", emit: mutation_sets
    path "mutant_zoo/annotated_structures/*.pdb", emit: annotated_structures
    path "debug/*", emit: debug_files

  script:
  """
  mkdir -p mutant_zoo/annotated_structures debug
  ${params.mutant_zoo_cmd} \\
    --input-annotations ${pass_fail_annotations} \\
    --input-complex ${recharged_complex} \\
    --outdir mutant_zoo > debug/mutant_zoo.log 2>&1
  """
}
