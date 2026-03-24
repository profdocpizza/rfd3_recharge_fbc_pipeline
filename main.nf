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

def parseSteps(raw) {
  if (!raw) return []
  return raw.toString()
    .split(',')
    .collect { it.trim().toLowerCase() }
    .findAll { it }
    .unique()
}

def stepEnabled(List steps, String key) {
  return steps.isEmpty() || steps.contains(key)
}

workflow {
  def steps = parseSteps(params.steps)
  def validSteps = [
    'hotspot_trim',
    'rfd3_design',
    'rfd3_to_bc_adapter',
    'protein_recharge',
    'freebindcraft_filtering',
    'freebindcraft_merge',
    'mutant_zoo',
    'aggregate_metrics',
    'dashboard_data',
    'run_manifest'
  ]
  def unknownSteps = steps.findAll { !validSteps.contains(it) }
  if (unknownSteps) {
    error "Unknown --steps value(s): ${unknownSteps.join(',')}. Valid: ${validSteps.join(',')}"
  }
  def stepIndex = validSteps.withIndex().collectEntries { k, i -> [(k): i] }
  def maxRequestedIndex = steps ? steps.collect { stepIndex[it] }.max() : (validSteps.size() - 1)

  outdirPath = file("${projectDir}/${params.outdir}")
  def stepPath = { String p -> Channel.fromPath("${outdirPath}/${p}", checkIfExists: true) }

  if (!params.target_pdb || !params.binder_length || !params.hotspot || !params.num_designs) {
    error "Required params missing: --target_pdb --binder_length --hotspot --num_designs"
  }

  target_pdb_ch = Channel.fromPath(params.target_pdb, checkIfExists: true)

  def trimmedPdbCh
  def hotspotMetaCh
  def targetMapCh
  if (stepEnabled(steps, 'hotspot_trim')) {
    HOTSPOT_TRIM(target_pdb_ch)
    trimmedPdbCh = HOTSPOT_TRIM.out.trimmed_pdbs.flatten()
    hotspotMetaCh = HOTSPOT_TRIM.out.hotspot_metadata
    targetMapCh = HOTSPOT_TRIM.out.target_numbering_map
  } else {
    trimmedPdbCh = stepPath("01_hotspot_trim/trimmed/trimmed_target_001.pdb")
    hotspotMetaCh = stepPath("01_hotspot_trim/trimmed/hotspot_metadata.json")
    targetMapCh = stepPath("01_hotspot_trim/trimmed/target_numbering_map.json")
  }

  def rfd3DesignsCh
  def rfd3ScoresCh
  if (maxRequestedIndex >= 1) {
    if (stepEnabled(steps, 'rfd3_design')) {
      RFD3_DESIGN(trimmedPdbCh, hotspotMetaCh)
      rfd3DesignsCh = RFD3_DESIGN.out.designed_complexes.flatten()
      rfd3ScoresCh = RFD3_DESIGN.out.rfd3_scores
    } else {
      rfd3DesignsCh = stepPath("02_rfd3/rfd3/*.pdb").flatten()
      rfd3ScoresCh = stepPath("02_rfd3/rfd3/score_metadata.json")
    }
  }

  def standardizedComplexesCh = null
  if (maxRequestedIndex >= 2) {
    adapter_input_ch = rfd3DesignsCh
      .combine(targetMapCh)
      .map { row -> tuple(row[0], row[1]) }

    if (stepEnabled(steps, 'rfd3_to_bc_adapter')) {
      RFD3_TO_BC_ADAPTER(adapter_input_ch)
      standardizedComplexesCh = RFD3_TO_BC_ADAPTER.out.standardized_complexes.flatten()
    } else {
      standardizedComplexesCh = stepPath("03_adapter/adapter/*_standardized.pdb").flatten()
    }
  }

  def rechargedSeqCh = null
  def chargeReportCh = null
  if (maxRequestedIndex >= 3) {
    if (stepEnabled(steps, 'protein_recharge')) {
      PROTEIN_RECHARGE(standardizedComplexesCh)
      rechargedSeqCh = PROTEIN_RECHARGE.out.recharged_sequences.flatten()
      chargeReportCh = PROTEIN_RECHARGE.out.charge_report
    } else {
      rechargedSeqCh = stepPath("04_protein_recharge/recharge/*_recharged.fasta").flatten()
      chargeReportCh = stepPath("04_protein_recharge/recharge/charge_report.json")
    }
  }

  def fbcTaggedScoresCh = null
  if (maxRequestedIndex >= 4) {
    fbc_complexes_batch_ch = standardizedComplexesCh.collect().map { tuple('batch', it) }
    fbc_sequences_batch_ch = rechargedSeqCh.collect().map { tuple('batch', it) }
    fbc_hotspot_batch_ch = trimmedPdbCh.first().map { tuple('batch', it) }
    fbc_input_ch = fbc_complexes_batch_ch
      .join(fbc_sequences_batch_ch)
      .join(fbc_hotspot_batch_ch)
      .map { id, complexes, sequences, hotspot -> tuple(complexes, sequences, hotspot) }

    if (stepEnabled(steps, 'freebindcraft_filtering')) {
      FREEBINDCRAFT_FILTERING(fbc_input_ch)
      fbcTaggedScoresCh = FREEBINDCRAFT_FILTERING.out.fbc_filter_scores_tagged.collect()
    } else {
      fbcTaggedScoresCh = stepPath("05_fbc_filtering/fbc_filtering/*_fbc_filter_scores.csv").collect()
    }
  }

  def mergedFbcScoresCh = null
  if (maxRequestedIndex >= 5) {
    if (stepEnabled(steps, 'freebindcraft_merge')) {
      FREEBINDCRAFT_MERGE(fbcTaggedScoresCh)
      mergedFbcScoresCh = FREEBINDCRAFT_MERGE.out.fbc_filter_scores
    } else {
      mergedFbcScoresCh = stepPath("05_fbc_filtering/merged/fbc_filter_scores.csv")
    }
  }

  if (maxRequestedIndex >= 6) {
    def mutationSetsCh
    if (stepEnabled(steps, 'mutant_zoo')) {
      MUTANT_ZOO(mergedFbcScoresCh, standardizedComplexesCh)
      mutationSetsCh = MUTANT_ZOO.out.mutation_sets
    } else {
      mutationSetsCh = stepPath("07_mutant_zoo/mutant_zoo/mutation_sets.json")
    }

    if (maxRequestedIndex >= 7) {
      def aggregateCsvCh
      def aggregateJsonCh
      if (stepEnabled(steps, 'aggregate_metrics')) {
        AGGREGATE_METRICS(rfd3ScoresCh, mergedFbcScoresCh, chargeReportCh, mutationSetsCh)
        aggregateCsvCh = AGGREGATE_METRICS.out.aggregate_csv
        aggregateJsonCh = AGGREGATE_METRICS.out.aggregate_json
      } else {
        aggregateCsvCh = stepPath("08_aggregate/aggregate/candidates.csv")
        aggregateJsonCh = stepPath("08_aggregate/aggregate/candidates.json")
      }

      if (maxRequestedIndex >= 8) {
        def dashboardCsvCh
        def dashboardJsonCh
        if (stepEnabled(steps, 'dashboard_data')) {
          DASHBOARD_DATA(aggregateCsvCh, aggregateJsonCh)
          dashboardCsvCh = DASHBOARD_DATA.out.dashboard_csv
          dashboardJsonCh = DASHBOARD_DATA.out.dashboard_json
        } else {
          dashboardCsvCh = stepPath("09_dashboard/dashboard/candidates.csv")
          dashboardJsonCh = stepPath("09_dashboard/dashboard/candidates.json")
        }

        if (maxRequestedIndex >= 9 && stepEnabled(steps, 'run_manifest')) {
          RUN_MANIFEST(dashboardCsvCh, dashboardJsonCh)
        }
      }
    }
  }
}

workflow FBC_ONLY {
  if (!params.fbc_input_complex || !params.fbc_recharged_sequence || !params.fbc_hotspot_pdb || !params.binder_length || !params.hotspot) {
    error "FBC_ONLY requires: --fbc_input_complex --fbc_recharged_sequence --fbc_hotspot_pdb --binder_length --hotspot"
  }

  recharged_complex_ch = Channel.fromPath(params.fbc_input_complex, checkIfExists: true).flatten().collect().map { tuple('batch', it) }
  recharged_sequence_ch = Channel.fromPath(params.fbc_recharged_sequence, checkIfExists: true).flatten().collect().map { tuple('batch', it) }
  hotspot_pdb_ch = Channel.fromPath(params.fbc_hotspot_pdb, checkIfExists: true).first().map { tuple('batch', it) }
  fbc_only_input_ch = recharged_complex_ch
    .join(recharged_sequence_ch)
    .join(hotspot_pdb_ch)
    .map { id, complexes, sequences, hotspot -> tuple(complexes, sequences, hotspot) }

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

  rerun_complexes_batch_ch = RFD3_TO_BC_ADAPTER.out.standardized_complexes.flatten().collect().map { tuple('batch', it) }
  rerun_sequences_batch_ch = PROTEIN_RECHARGE.out.recharged_sequences.flatten().collect().map { tuple('batch', it) }
  rerun_hotspot_batch_ch = hotspot_pdb_ch.first().map { tuple('batch', it) }
  rerun_fbc_input_ch = rerun_complexes_batch_ch
    .join(rerun_sequences_batch_ch)
    .join(rerun_hotspot_batch_ch)
    .map { id, complexes, sequences, hotspot -> tuple(complexes, sequences, hotspot) }

  FREEBINDCRAFT_FILTERING(rerun_fbc_input_ch)
  FREEBINDCRAFT_MERGE(FREEBINDCRAFT_FILTERING.out.fbc_filter_scores_tagged.collect())
}
