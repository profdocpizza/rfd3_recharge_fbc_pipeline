include { HOTSPOT_TRIM } from './modules/local/hotspot_trim'
include { RFD3_DESIGN } from './modules/local/rfd3_design'
include { RFD3_TO_BC_ADAPTER } from './modules/local/rfd3_to_bc_adapter'
include { PROTEIN_RECHARGE } from './modules/local/protein_recharge'
include { FREEBINDCRAFT_FILTERING } from './modules/local/freebindcraft_filtering'
include { FREEBINDCRAFT_MERGE } from './modules/local/freebindcraft_merge'
include { MUTANT_ZOO } from './modules/local/mutant_zoo'
include { AGGREGATE_METRICS } from './modules/local/aggregate_metrics'
include { DASHBOARD_DATA } from './modules/local/dashboard_data'
include { RUN_MANIFEST } from './modules/local/run_manifest'

workflow {
  if (!params.target_pdb || !params.binder_length || !params.hotspot) {
    error "Required params missing: --target_pdb --binder_length --hotspot"
  }

  target_pdb_ch = Channel.fromPath(params.target_pdb, checkIfExists: true)

  HOTSPOT_TRIM(target_pdb_ch)

  RFD3_DESIGN(
    HOTSPOT_TRIM.out.trimmed_pdbs.flatten(),
    HOTSPOT_TRIM.out.hotspot_metadata
  )

  adapter_input_ch = RFD3_DESIGN.out.designed_complexes
    .flatten()
    .combine(HOTSPOT_TRIM.out.target_numbering_map)
    .map { row -> tuple(row[0], row[1]) }

  RFD3_TO_BC_ADAPTER(adapter_input_ch)

  PROTEIN_RECHARGE(RFD3_TO_BC_ADAPTER.out.standardized_complexes.flatten())

  adapter_complex_by_id = RFD3_TO_BC_ADAPTER.out.standardized_complexes
    .flatten()
    .map { p -> tuple(p.baseName.replace('_standardized',''), p) }
  recharge_seq_by_id = PROTEIN_RECHARGE.out.recharged_sequences
    .flatten()
    .map { p -> tuple(p.baseName.replace('_recharged',''), p) }
  recharge_pairs = adapter_complex_by_id
    .join(recharge_seq_by_id)
    .map { id, complex, seq -> tuple(complex, seq) }
  fbc_input_ch = recharge_pairs
    .combine(HOTSPOT_TRIM.out.trimmed_pdbs)
    .map { row -> tuple(row[0], row[1], row[2]) }

  FREEBINDCRAFT_FILTERING(fbc_input_ch)
  FREEBINDCRAFT_MERGE(FREEBINDCRAFT_FILTERING.out.fbc_filter_scores_tagged.collect())

  MUTANT_ZOO(
    FREEBINDCRAFT_MERGE.out.fbc_filter_scores,
    RFD3_TO_BC_ADAPTER.out.standardized_complexes.flatten()
  )

  AGGREGATE_METRICS(
    RFD3_DESIGN.out.rfd3_scores,
    FREEBINDCRAFT_MERGE.out.fbc_filter_scores,
    PROTEIN_RECHARGE.out.charge_report,
    MUTANT_ZOO.out.mutation_sets
  )

  DASHBOARD_DATA(AGGREGATE_METRICS.out.aggregate_csv, AGGREGATE_METRICS.out.aggregate_json)

  RUN_MANIFEST(DASHBOARD_DATA.out.dashboard_csv, DASHBOARD_DATA.out.dashboard_json)
}

workflow FBC_ONLY {
  if (!params.fbc_input_complex || !params.fbc_recharged_sequence || !params.fbc_hotspot_pdb || !params.binder_length || !params.hotspot) {
    error "FBC_ONLY requires: --fbc_input_complex --fbc_recharged_sequence --fbc_hotspot_pdb --binder_length --hotspot"
  }

  recharged_complex_ch = Channel.fromPath(params.fbc_input_complex, checkIfExists: true)
  recharged_sequence_ch = Channel.fromPath(params.fbc_recharged_sequence, checkIfExists: true)
  hotspot_pdb_ch = Channel.fromPath(params.fbc_hotspot_pdb, checkIfExists: true)
  fbc_only_input_ch = recharged_complex_ch
    .combine(recharged_sequence_ch)
    .combine(hotspot_pdb_ch)
    .map { row -> tuple(row[0], row[1], row[2]) }

  FREEBINDCRAFT_FILTERING(fbc_only_input_ch)
  FREEBINDCRAFT_MERGE(FREEBINDCRAFT_FILTERING.out.fbc_filter_scores_tagged.collect())
}

workflow ADAPTER_RECHARGE_FBC {
  if (!params.rerun_rfd3_designs_glob || !params.rerun_target_numbering_map || !params.rerun_hotspot_pdb || !params.binder_length || !params.hotspot) {
    error "ADAPTER_RECHARGE_FBC requires: --rerun_rfd3_designs_glob --rerun_target_numbering_map --rerun_hotspot_pdb --binder_length --hotspot"
  }

  rfd3_designs_ch = Channel.fromPath(params.rerun_rfd3_designs_glob, checkIfExists: true).flatten()
  target_map_ch = Channel.value(file(params.rerun_target_numbering_map, checkIfExists: true))
  hotspot_pdb_ch = Channel.value(file(params.rerun_hotspot_pdb, checkIfExists: true))

  adapter_rerun_input_ch = rfd3_designs_ch
    .combine(target_map_ch)
    .map { row -> tuple(row[0], row[1]) }

  RFD3_TO_BC_ADAPTER(adapter_rerun_input_ch)

  PROTEIN_RECHARGE(RFD3_TO_BC_ADAPTER.out.standardized_complexes.flatten())

  rerun_adapter_by_id = RFD3_TO_BC_ADAPTER.out.standardized_complexes
    .flatten()
    .map { p -> tuple(p.baseName.replace('_standardized',''), p) }
  rerun_seq_by_id = PROTEIN_RECHARGE.out.recharged_sequences
    .flatten()
    .map { p -> tuple(p.baseName.replace('_recharged',''), p) }
  rerun_recharge_pairs = rerun_adapter_by_id
    .join(rerun_seq_by_id)
    .map { id, complex, seq -> tuple(complex, seq) }
  rerun_fbc_input_ch = rerun_recharge_pairs
    .combine(hotspot_pdb_ch)
    .map { row -> tuple(row[0], row[1], row[2]) }

  FREEBINDCRAFT_FILTERING(rerun_fbc_input_ch)
  FREEBINDCRAFT_MERGE(FREEBINDCRAFT_FILTERING.out.fbc_filter_scores_tagged.collect())
}
