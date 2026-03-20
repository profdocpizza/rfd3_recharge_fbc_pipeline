   cd /home/tadas/code/rfd3_recharge_fbc_pipeline
   # --steps Valid: hotspot_trim,rfd3_design,rfd3_to_bc_adapter,protein_recharge,freebindcraft_filtering,freebindcraft_merge,mutant_zoo,aggregate_metrics,dashboard_data,run_manifest

   ./nextflow run main.nf -profile local -c conf/smoke.config \
     --steps mutant_zoo \
     --run_id "kinesin_hotspot_run10" \
     --target_pdb "/home/tadas/code/startup_lanternfish_binders/inputs/TUBB4B_TUBA1A/hexamer/TUBB4B_TUBA1A_hexamer.pdb" \
     --binder_length "80-150" \
     --num_designs "10" \
     --hotspot "B:160 A:409 B:260" \
     --hotspot_trimming_max_residues "470" \
     --hotspot_trimming_residues "A:409 A:113 A:108 B:160 B:157" \
     --hotspot_selector_cmd "conda run -n hotspot python $(pwd)/scripts/tool_wrappers/hotspot_selector_wrapper.py" \
     --rfd3_cmd "conda run -n ligandmpnn_env python $(pwd)/scripts/tool_wrappers/rfd3_wrapper.py --rfd3-env rfd3" \
     --protein_recharge_cmd "conda run -n ligandmpnn_env python $(pwd)/scripts/tool_wrappers/protein_recharge_wrapper.py" \
     --freebindcraft_cmd "conda run -n FreeBindCraft python $(pwd)/scripts/tool_wrappers/freebindcraft_wrapper.py" \
     --mutant_zoo_cmd "conda run -n mutantzoo_env python $(pwd)/scripts/tool_wrappers/mutant_zoo_wrapper.py"

# A:260 i think did not work bcs was burried