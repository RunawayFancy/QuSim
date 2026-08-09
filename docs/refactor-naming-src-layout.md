# Refactor run record — naming & structure pass (Phase 0)

Date: 2026-08-09 · Branch: `refactor/naming-src-layout` · Status: **complete, verified statically**

This is Phase 0 of the qutip-4 → qutip-5 upgrade: a **pure naming & structure pass**.
No qutip API usage or numerical behavior was changed. Written so the next-phase
(qutip-5) agent does not have to re-derive the environment or the layout.

## Environment (verified, do not re-discover)
- Conda env `qusim_qutip5`: **Python 3.10.16, qutip 5.1.1, numpy 2.2.1, scipy 1.15.1**.
- Interpreter: `E:\conda\envs\qusim_qutip5\python.exe`.
- The package is **not** pip-installed in the env. Import it for testing with
  `PYTHONPATH=src` (Bash) / `$env:PYTHONPATH=(Join-Path (pwd) 'src')` (PowerShell).
  Avoid `pip install -e .` without `--no-deps` — `setup.py` still pins `qutip==4.7.5`
  and would downgrade the env.

## What changed
### Layout: flat package → src layout
`qusim/` → `src/qusim/` (still imported as `import qusim`). `setup.py` now uses
`package_dir={'': 'src'}` + `find_packages(where='src')`.

Sub-packages renamed to lowercase snake_case (git-tracked as renames):
`System→system`, `PulseGen→pulse_gen`, `DataPlot→data_plot`, `DataView→data_view`,
`Instruments→instruments`, `Scans→scans`, `Utils→utils`.

Module files: `scan1D.py→scan_1d.py`, `scan2D.py→scan_2d.py`, `scanND.py→scan_nd.py`.

### Renames applied (strict PEP-8 violations + obvious typos only)
- **Classes:** `scan1D→Scan1D`, `qsviewer→QSViewer`, `INSTR→Instr`.
- **`system/transmon_system.py` (TransmonSys):** `get_Hq_Ha→get_hq_ha`,
  `get_H_inter→get_h_inter`, `H_XY_drive→h_xy_drive`, `H_Z_bias→h_z_bias`.
- **`system/arb_qubit_system.py` (ArbQubitSys):** `get_H_0→get_h_0`,
  `get_H_extra→get_h_extra`, `get_H_inter→get_h_inter`, `get_V→get_v`,
  `get_H_XY_drive→get_h_xy_drive`, `get_H_Z_bias→get_h_z_bias`,
  `get_H_int_bias→get_h_int_bias`, `get_H_d→get_h_d`, `get_Hd_channel→get_hd_channel`,
  `H_XY_drive→h_xy_drive`, `H_Z_bias→h_z_bias`, `H_int_bias→h_int_bias`
  (co-located local vars renamed for consistency).
- **`pulse_gen/pulse_config.py`:** method `DRAG→drag`, param `DRAG_config→drag_config`.
  Kept: attribute `DRAG_config_list`, class `DRAGConfig` (out of scope this pass).
- **`instruments/avoid_crossing.py`:** `get_E_diff→get_e_diff`.
- **`utils/noise_trafofn_tdbasefn.py` (ChrgNoiseExchangeQD):** `Jexchange→j_exchange`,
  `J2Vbarrier→j2v_barrier`, `dJdV→djdv`, and typo `tranfofn_charge_noise→trafofn_charge_noise`.
- **`instruments/tools.py`:** `get_XY_element→get_xy_element`, `get_Z_element→get_z_element`.
- **`data_plot/plot_lib.py`:** `plot_Elevel_dynamics→plot_elevel_dynamics`.
- **`data_view/tracer_device.py`:** `RunOne→run_one`.
- **`pulse_gen/edges.py`:** typo `rasing_t→rising_t`.
- **Cross-cutting typo:** `NoiseTimeConfig.tranfofn→trafofn` field
  (`pulse_gen/noise_config.py` + `pulse_gen/noise_gen.py`), matching the module spelling.

All internal call sites and absolute imports (`qusim.PulseGen→qusim.pulse_gen`, etc.)
were updated across `src/`, plus the top-level `src/qusim/__init__.py`.

### Cleanup
- Deleted dead `noise_gen_old.py`.
- De-duplicated `SwapAvoidCrossing`: canonical copy in `instruments/avoid_crossing.py`;
  `instruments/tools.py` now re-exports it (`from .avoid_crossing import SwapAvoidCrossing`),
  so both import paths still work (verified `tools.SwapAvoidCrossing is avoid_crossing.SwapAvoidCrossing`).
- Removed stale on-disk `build/` and `dist/` artifacts (regenerable; already gitignored).

### Notebooks
Updated import paths + `plot_elevel_dynamics` in the 8 **active** notebooks
(`Tutorial/transmon_tutorial/{DRAG,Energy_level,Setup_system,System_dynamics,ZZ_coupling,data_saving,iswap_like}.ipynb`
and `Tutorial/arb_qubit_tutorial/arb_qubit_demo.ipynb`). `Tutorial/**/archive/` left untouched.

## Verification (static — full runtime is deferred to the qutip-5 phase)
- `python -m compileall src/qusim` → exit 0.
- Grep sweep: **zero** old tokens remain in `src/` and active notebooks.
- `import qusim` succeeds under qutip 5.1; renamed-symbol assertions all pass
  (methods/classes/functions resolve; `SwapAvoidCrossing` is a single class object).
- 8 notebooks re-validated as JSON.

## Known qutip-5 breakages already observed (head-start for next phase — NOT naming bugs)
1. **`Options` deprecation** — every solver default arg `option = Options(rtol=1e-8)`
   (`system/transmon_system.py`, `system/arb_qubit_system.py`) triggers a qutip
   `FutureWarning`: options should be passed as a plain `dict`. `mesolve`/`propagator`
   signatures also changed.
2. **`qusim.data_plot.plot_tomo`** fails to import: `No module named 'qutip.qobj'`
   (qutip moved `Qobj`; update to `from qutip import Qobj`). Note `plot_tomo.py` is
   gitignored via a broad `plot_tomo.py` rule in `.gitignore`.

## Deferred to the qutip-5 phase (out of scope here)
qutip-5 API migration (`Options`→dict/`SolverOptions`, `mesolve`/`propagator` args,
`enr_destroy`, `Qobj`/`qutip.qobj`, numpy-2 compat), dependency pins + version in
`setup.py`, semantic/clarity renames (`TransmonSys→TransmonSystem`, `co_list`, …), and
attribute-name cleanup (`DRAG_config_list`, `H_q`/`H_inter`, etc.).
