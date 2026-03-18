process FREEBINDCRAFT_PASS1 {
  label 'large'
  container params.container_fbc
  publishDir "${params.outdir}/05_fbc_pass1", mode: 'copy'

  input:
    path recharged_complex

  output:
    path "fbc_pass1/filtered_candidates/*.pdb", emit: filtered_candidates
    path "fbc_pass1/pass1_scores.csv", emit: pass1_scores
    path "fbc_pass1/designed_sequences.fasta", emit: designed_sequences

  script:
  """
  mkdir -p fbc_pass1
  ${params.freebindcraft_cmd} \\
    --mode full \\
    --input-complex ${recharged_complex} \\
    --outdir fbc_pass1

  SCORE_FILE=""
  for CANDIDATE in final_design_stats.csv all_mpnn_full_stats.csv scores.csv; do
    if [[ -f "fbc_pass1/\${CANDIDATE}" ]]; then
      SCORE_FILE="\${CANDIDATE}"
      break
    fi
  done

  if [[ -z "\${SCORE_FILE}" ]]; then
    echo "No score table found in fbc_pass1 output." >&2
    exit 1
  fi

  cp "fbc_pass1/\${SCORE_FILE}" "fbc_pass1/pass1_scores.csv"
  """
}
