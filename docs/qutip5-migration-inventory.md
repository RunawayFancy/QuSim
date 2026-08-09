# qutip 4 → 5 breakage inventory (Phase 1)

Date: 2026-08-09 · Builds on `docs/refactor-naming-src-layout.md` (Phase 0, src layout done).
This is the "review" deliverable: every place in `src/` that must change for qutip 5.1 /
numpy 2. Fixes are applied by subagents; results recorded in
`docs/qutip5-migration-core.md` and `docs/qutip5-migration-plot-tomo.md`.

## Environment (verified)
- `qusim_qutip5`: Python 3.10.16, **qutip 5.1.1**, numpy 2.2.1, scipy 1.15.1,
  **lmfit 1.3.4**, tqdm 4.67.1. Interpreter `E:\conda\envs\qusim_qutip5\python.exe`.
- Run tests with `PYTHONPATH=src`. Do NOT `pip install -e .` without `--no-deps`
  (setup.py still pins qutip 4.7.5).

## Verified qutip-5 signatures / facts
- `mesolve(H, rho0, tlist, c_ops=None, *, e_ops=None, args=None, options=None)` —
  `e_ops`/`args`/`options` are **keyword-only**; positional slots are deprecated shims.
- `propagator(H, t, c_ops=None, args=None, options=None, **kwargs)` — **no** `parallel`
  or `progress_bar` kwargs. Parallelism/progress are controlled via the `options` dict.
- `options` must be a plain **dict** (e.g. `{'rtol': 1e-8}`). `qutip.Options` still exists
  as a deprecated shim (emits FutureWarning) — replace it.
- Still top-level in qutip 5: `enr_destroy`, `isket`, `Qobj`, `vector_to_operator`, `settings`.
- Gone in qutip 5: modules `qutip.qobj`, `qutip.superoperator` (path), `qutip.superop_reps`.

## Breakages

### A. Core solver API — HIGH (exercised by active notebooks)
**`system/transmon_system.py`**
- L225, L248, L254: default arg `option = Options(rtol=1e-8)` → `option = {'rtol': 1e-8}`
  (define a module constant, e.g. `DEFAULT_OPTIONS = {'rtol': 1e-8}`, and use `None`-guard).
- L249 `mesolve(..., c_ops=self.co_list(), options=option)` — OK once `option` is a dict.
- L261 `propagator(H_d, tlist, self.co_list(), {}, option)` — positional c_ops/args/options
  OK once `option` is a dict.

**`system/arb_qubit_system.py`**
- L499, L524, L535: `Options(rtol=1e-8)` defaults → dict.
- L525 `mesolve(H_d, initial_state, tlist, self.co_list, [], options=option)`:
  (i) `[]` is the deprecated positional `e_ops` slot → pass `e_ops=[]` or drop;
  (ii) **`self.co_list` is passed without `()`** — `co_list` is a method (def at L387),
  so this passes a bound method as `c_ops`. transmon calls `self.co_list()`. Confirm intended
  behavior and fix to `self.co_list()` (or make `co_list` a `@property`). Same at L546.
- L546 `propagator(H_d, tlist, self.co_list, {}, option, parallel=do_parallel,
  progress_bar=do_progress_bar)`: **`parallel=` removed**; **`progress_bar=` removed** →
  move progress into `options` (`{'progress_bar': 'tqdm'}` when requested). Decide how to
  honor `do_parallel`/`do_progress_bar` params (map to options or drop with a note).

### B. Result / Qobj indexing — MEDIUM (verify at runtime, don't assume)
- transmon L302 `np.abs((states[ii].dag()*state))**2`; arb L631
  `((states[ii]).dag()*state)[0][0][0]` — Qobj `*` → 1×1 Qobj; qutip-5 scalar extraction
  differs (`.full()[0,0]` / complex()). Verify against a running notebook and fix minimally.
- transmon L313 / arb L644: density-matrix `.tr()` expressions — likely OK, verify.
- `instruments/tunablec.py` L30-31 `state.dag()*result.states[-1].data` — `.data` is now a
  qutip-5 `Data` object, not a numpy/scipy matrix. Verify/adjust.
- `instruments/angle.py` uses `isket` (still valid) + `res.states[-1]` — verify.

### C. Module-path breakages — LOW for active harness
**`data_plot/plot_tomo.py`** (gitignored; imported only by *archive* notebooks, none active):
- L24 `from qutip.qobj import Qobj` → `from qutip import Qobj`.
- L25 `from qutip.superoperator import vector_to_operator` → `from qutip import vector_to_operator`.
- L26 `from qutip.superop_reps import _super_to_superpauli, _isqubitdims` → **no qutip-5
  equivalent** (private, module removed). Options: reimplement the two helpers locally, or
  delegate `hinton` to qutip-5's built-in `qutip.hinton`. Preserve current call signature.
- L28 `from qutip import settings` — still valid.

### D. Packaging / deps
- `src/qusim/processing/` has **no `__init__.py`** → add one.
- `setup.py`: bump `install_requires` (`qutip>=5.0`, `numpy>=2.0`, add `lmfit>=1.3`), and
  bump the `Programming Language :: Python` classifier as needed. (Phase 0 deferred this.)

### E. Clean (no action)
- numpy-2 scalar aliases (`np.float`/`np.int`/`np.complex`/`np.bool`): none present.
- `utils/floquet.py`: already qutip-5 native (`qt.solver.Propagator`, `parallel_map`,
  `.data_as('ndarray')`, dict options).
- `utils/units.py`, `utils/general_utils.py`, `processing/fitting.py`: no qutip; fine.

## New modules — status
- `utils/floquet.py`: qutip-5 native; imports `qusim.processing.fitting` (needs `__init__.py`).
- `processing/fitting.py`: lmfit-based, no qutip; needs `processing/__init__.py`; `lmfit` present.
- `utils/units.py`, `utils/general_utils.py`: pure Python, OK.
- None of the four are imported by active notebooks yet — verification = clean import.

## Acceptance for Phase 1
1. `import qusim` clean, no FutureWarning from `Options`.
2. `import qusim.utils.floquet`, `qusim.processing.fitting`, `qusim.data_plot.plot_tomo` all OK.
3. Active transmon + arb notebooks execute end-to-end in `qusim_qutip5` (real acceptance test).
4. Grep: no `Options(`, no `qutip.qobj`/`qutip.superop_reps`, no `parallel=` in propagator.
