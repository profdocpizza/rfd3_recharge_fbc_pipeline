process PROTEIN_RECHARGE {
  label 'medium'
  container params.container_recharge
  publishDir "${params.outdir}/04_protein_recharge", mode: 'copy'
  publishDir "${params.outdir}/debug/04_protein_recharge", mode: 'copy', pattern: "debug/*"

  input:
    path standardized_complex

  output:
    path "recharge/*_recharged.fasta", emit: recharged_sequences
    path "recharge/charge_report.json", emit: charge_report
    path "recharge/*_charge_report.json", emit: charge_reports, optional: true
    path "debug/*", emit: debug_files

  script:
  """
  mkdir -p recharge debug
  RECHARGE_CONFIG_ARG=""
  if [[ -f "${params.recharge_config}" && -s "${params.recharge_config}" ]]; then
    RECHARGE_CONFIG_ARG="--config ${params.recharge_config}"
  fi

  ${params.protein_recharge_cmd} \\
    --input-complex ${standardized_complex} \\
    \${RECHARGE_CONFIG_ARG} \\
    --outdir recharge > debug/protein_recharge.log 2>&1

  if [[ -d recharge/debug ]]; then
    cp -r recharge/debug/* debug/ || true
  fi
  """
}
