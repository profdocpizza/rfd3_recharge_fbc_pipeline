   cd /home/tadas/code/rfd3_recharge_fbc_pipeline
   # --steps value 'hotspot_selection'. Valid: hotspot_trim,rfd3_design,rfd3_to_bc_adapter,protein_recharge,freebindcraft_filtering,freebindcraft_merge,mutant_zoo,aggregate_metrics,dashboard_data,run_manifest
# (chain A and resi 159+108) or (chain D and resi 339+309)
   ./nextflow run main.nf -profile local -c conf/smoke.config \
     --run_id "A_hotspot_run100" \
     --target_pdb "/home/tadas/code/startup_lanternfish_binders/inputs/TUBB4B_TUBA1A/hexamer/TUBB4B_TUBA1A_hexamer.pdb" \
     --binder_length "80-150" \
     --num_designs "100" \
     --hotspot "A:159 D:339 A:156 D:309" \
     --hotspot_trimming_max_residues "470" \
     --hotspot_trimming_residues "A:159 A:108 D:339 D:309" \
     --hotspot_selector_cmd "conda run -n hotspot python /home/tadas/code/rfd3_recharge_fbc_pipeline/scripts/tool_wrappers/hotspot_selector_wrapper.py" \
     --rfd3_cmd "conda run -n ligandmpnn_env python /home/tadas/code/rfd3_recharge_fbc_pipeline/scripts/tool_wrappers/rfd3_wrapper.py --rfd3-env rfd3" \
     --protein_recharge_cmd "conda run -n ligandmpnn_env python /home/tadas/code/rfd3_recharge_fbc_pipeline/scripts/tool_wrappers/protein_recharge_wrapper.py" \
     --freebindcraft_cmd "conda run -n FreeBindCraft python /home/tadas/code/rfd3_recharge_fbc_pipeline/scripts/tool_wrappers/freebindcraft_wrapper.py" \
     --mutant_zoo_cmd "conda run -n mutantzoo_env python /home/tadas/code/rfd3_recharge_fbc_pipeline/scripts/tool_wrappers/mutant_zoo_wrapper.py"
