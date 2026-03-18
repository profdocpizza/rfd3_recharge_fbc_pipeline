# RFD3 → Recharge → FBC Binder Pipeline

Nextflow DSL2 pipeline (active dev mode), with reusable step wrappers and debug artifacts.

## Pipeline stages

1. `HOTSPOT_TRIM`
2. `RFD3_DESIGN`
3. `RFD3_TO_BC_ADAPTER`
4. `PROTEIN_RECHARGE`
5. `FREEBINDCRAFT_FILTERING` (`--steps filtering`)
6. `FREEBINDCRAFT_MERGE` (merges per-design FBC score CSVs)
7. `MUTANT_ZOO`
8. `AGGREGATE_METRICS`
9. `DASHBOARD_DATA`
10. `RUN_MANIFEST`

## Install on a new system (from git)

### 1) Prerequisites

- Linux
- Conda (or Mamba)
- Java 17+
- Nextflow 25+
- Git

### 2) Clone pipeline repo

```bash
mkdir -p ~/code
cd ~/code
git clone <YOUR_PIPELINE_REPO_URL> RFD3_FBC_binder_design
cd RFD3_FBC_binder_design
```

### 3) Clone tool repos used by wrappers

```bash
cd ~/code
git clone https://github.com/profdocpizza/hotspot_selector
git clone https://github.com/profdocpizza/ProteinRecharge
git clone https://github.com/profdocpizza/FreeBindCraft
git clone https://github.com/profdocpizza/MUTANT_ZOO
```

### 4) Create per-tool conda envs (example names)

```bash
conda create -n hotspot_env python=3.10 -y
conda create -n rfd3_env python=3.12 -y
conda create -n recharge_env python=3.10 -y
conda create -n FreeBindCraft python=3.10 -y
conda create -n mutantzoo_env python=3.10 -y
```

Install each tool in its env per upstream docs.

### 5) Put your inputs/configs in place

- Target PDB: `inputs/targets/...`
- RFD3 config: `inputs/configs/rfd3/...`
- Recharge config: `inputs/configs/recharge/...`
- FBC configs:
  - `inputs/configs/freebindcraft/run.json`
  - `inputs/configs/freebindcraft/default_filters.json`
  - `inputs/configs/freebindcraft/default_4stage_multimer.json`

## Runtime model (important)

- Pipeline uses **existing local repos and envs**.
- It does **not** auto-clone tools and does **not** auto-create conda envs.
- You control tool commands via params in `nextflow.config` or CLI:
  - `--hotspot_selector_cmd`
  - `--rfd3_cmd`
  - `--protein_recharge_cmd`
  - `--freebindcraft_cmd`
  - `--mutant_zoo_cmd`

## Required runtime params

- `--target_pdb`
- `--binder_length` (example: `"80-150"`)
- `--hotspot` (example: `"A:159 F:106 E:307 D:337 B:405"`)

## Full run command (dev mode)

```bash
cd ~/code/RFD3_FBC_binder_design

./nextflow run main.nf -profile local -c conf/smoke.config \
  --run_id debug_run_full \
  --target_pdb "$(pwd)/inputs/targets/TUBB4B_TUBA1A_hexamer.pdb" \
  --binder_length "80-150" \
  --hotspot "A:159 F:106 E:307 D:337 B:405" \
  --hotspot_selector_cmd "conda run -n hotspot_env python $(pwd)/scripts/tool_wrappers/hotspot_selector_wrapper.py" \
  --rfd3_cmd "conda run -n rfd3_env python $(pwd)/scripts/tool_wrappers/rfd3_wrapper.py --num-designs 2" \
  --protein_recharge_cmd "conda run -n recharge_env python $(pwd)/scripts/tool_wrappers/protein_recharge_wrapper.py" \
  --freebindcraft_cmd "conda run -n FreeBindCraft python $(pwd)/scripts/tool_wrappers/freebindcraft_wrapper.py" \
  --mutant_zoo_cmd "conda run -n mutantzoo_env python $(pwd)/scripts/tool_wrappers/mutant_zoo_wrapper.py"
```

## Rerun only adapter → recharge → FBC (overwrite)

Use helper script:

```bash
bash /home/tadas/code/random/run_RFD3_FBC_pipeline.sh debug_run_full
```

This entrypoint uses existing `02_rfd3` outputs and reruns:
- `RFD3_TO_BC_ADAPTER`
- `PROTEIN_RECHARGE`
- `FREEBINDCRAFT_FILTERING`
- `FREEBINDCRAFT_MERGE`

## Where to inspect each step

- Published outputs: `results/<run_id>/0X_*`
- Step debug logs/artifacts: `results/<run_id>/debug/0X_*`
- Exact task command/env: `work/<hash>/.command.sh`
- Task stdout/stderr: `work/<hash>/.command.log` / `.command.err`
- Workflow trace/report: `results/<run_id>/manifests/`

## Current FBC filtering contract

Per design, wrapper prepares:

- `Trajectory/Relaxed/<design>_l{chainBlen}_s{5digits}.pdb`
- `MPNN/Sequences/<same_stem>_mpnn1.fasta` (from recharge sequence)

Then runs:

```bash
python /home/tadas/code/FreeBindCraft/bindcraft.py \
  --settings <runtime settings.json> \
  --filters <default_filters.json> \
  --advanced <default_4stage_multimer.json> \
  --steps filtering
```

Pipeline output for FBC is a single merged score file:

- `results/<run_id>/05_fbc_filtering/merged/fbc_filter_scores.csv`
