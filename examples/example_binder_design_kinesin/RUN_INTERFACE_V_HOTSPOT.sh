   cd /home/tadas/code/rfd3_recharge_fbc_pipeline
   # --steps value 'hotspot_selection'. Valid: hotspot_trim,rfd3_design,rfd3_to_bc_adapter,protein_recharge,freebindcraft_filtering,freebindcraft_merge,mutant_zoo,aggregate_metrics,dashboard_data,run_manifest

   ./nextflow run main.nf -profile local \
     --steps "protein_recharge,freebindcraft_filtering" \
     --run_id "V_hotspot_run100" \
     --target_pdb "/home/tadas/code/startup_lanternfish_binders/inputs/TUBB4B_TUBA1A/hexamer/TUBB4B_TUBA1A_hexamer.pdb" \
     --binder_length "80-150" \
     --num_designs "100" \
     --hotspot "A:159 F:397" \
     --hotspot_trimming_max_residues "410" \
     --hotspot_trimming_residues "F:106 F:111 A:159 A:263" \
     --hotspot_selector_cmd "conda run -n hotspot python $(pwd)/scripts/tool_wrappers/hotspot_selector_wrapper.py" \
     --rfd3_cmd "conda run -n ligandmpnn_env python $(pwd)/scripts/tool_wrappers/rfd3_wrapper.py --rfd3-env rfd3" \
     --protein_recharge_cmd "conda run -n ligandmpnn_env python $(pwd)/scripts/tool_wrappers/protein_recharge_wrapper.py" \
     --freebindcraft_cmd "conda run -n FreeBindCraft python $(pwd)/scripts/tool_wrappers/freebindcraft_wrapper.py" \
     --mutant_zoo_cmd "conda run -n mutantzoo_env python $(pwd)/scripts/tool_wrappers/mutant_zoo_wrapper.py"
