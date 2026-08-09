# qutip 5 migration — final summary (Phase 1)

Date: 2026-08-09 · Branch: `refactor/naming-src-layout` (continues Phase 0 src-layout work).
Integrates: `qutip5-migration-inventory.md` (review), `qutip5-migration-core.md`
(core subagent), `qutip5-migration-plot-tomo.md` (plot_tomo subagent), + host packaging.
Env: `qusim_qutip5` — Python 3.10.16, qutip 5.1.1, numpy 2.2.1, scipy 1.15.1, lmfit 1.3.4.

## Done & verified
### Core solvers (`system/transmon_system.py`, `system/arb_qubit_system.py`)
- `Options(rtol=1e-8)` defaults → `DEFAULT_SOLVER_OPTIONS = {'rtol': 1e-8}` dict + `None`-guard.
- `mesolve`: `options`/`e_ops` passed as keywords; dropped deprecated positional `[]` e_ops (arb).
- `propagator`: `options` is a dict; **added mandatory `tlist=` keyword** (required for the
  `[H0,[H1,array],...]` array-coefficient Hamiltonian in qutip 5); removed `parallel=`/
  `progress_bar=` kwargs. `do_progress_bar`→`options['progress_bar']='tqdm'`; `do_parallel`
  kept for API compatibility but documented as a no-op (qutip-5 `propagator` has no parallel knob).
- numpy-2 fix: `np.complex_` → `np.complex128` (was blocking *all* ArbQubitSys construction).
- qutip-5 Qobj-indexing fixes in `get_state_index` / `get_data_list` (`[0][0][0]`→`[0]` etc.).
- Finding: the flagged `co_list` "bug" is NOT a bug — arb's `co_list` is a `@property`.

### Other files
- `instruments/tunablec.py` `fsim_zhen`: dropped unsupported `Qobj * Data` (`.states[-1].data`) + trailing `[0][0]`.
- `data_plot/plot_tomo.py`: qutip-4 internal imports → qutip-5 (`from qutip import Qobj,
  vector_to_operator, settings`; `qutip.core.superop_reps._to_superpauli` /public `isqubitdims`).
  `hinton` verified on 4×4, 2-qubit dm, superoperator, operket.
- `instruments/angle.py`, `pulse_gen/simulation_option.py`, `pulse_gen/pulse_buffer.py`:
  no qutip-5 change required.

### Packaging (host)
- Added `src/qusim/processing/__init__.py` (was missing; broke package discovery).
- `setup.py`: `qutip>=5.0`, `numpy>=1.26`, `scipy>=1.12`, added `lmfit>=1.2`, relaxed the rest.

### Verification evidence
- `import qusim` under `-W error::FutureWarning` → clean (no `Options` warning).
- Grep: zero real `Options(...)` constructors, no removed propagator kwargs, no `np.complex_`,
  no `qutip.qobj`/old `qutip.superop_reps` path.
- Subagent `verify_core.py`: 13/13 — real `mesolve` AND `propagator` on 2-qubit TransmonSys
  and ArbQubitSys (± collapse ops), `get_data_list`, `get_angle`, `tunablec.fsim_zhen`.
- New modules import clean: `processing.fitting`, `utils.floquet`, `utils.units`, `general_utils`.

### Active-notebook verification matrix (executed end-to-end in `qusim_qutip5`)
| Notebook | Status | Notes |
|---|---|---|
| transmon/DRAG | PASS | real end-to-end mesolve (host re-ran, 13 cells) |
| transmon/Setup_system | PASS | subagent |
| transmon/ZZ_coupling | PASS | subagent |
| transmon/Energy_level | PASS | subagent |
| transmon/System_dynamics | PASS | 26 cells; unblocked by `merge_pulse_chan` + `__eq__` fixes |
| transmon/iswap_like | PASS | 19 cells; iSWAP fidelity 99.98%; needed notebook fixes (below) |
| transmon/data_saving | NOT RUN | pure `qsave`/pickle, no qutip surface; writes into `Data/` so not executed per AGENT.md |
| arb_qubit/arb_qubit_demo | FAIL (pre-existing) | stale ArbQubitSys API; not a qutip-5 issue — see Open items |

### Follow-up fixes applied after the migration (to unblock notebooks)
- `pulse_gen/pulse_buffer.py::merge_pulse_chan` — rewrote the crashing
  `int(np.intersect1d(...))` as a clean type+qindex match loop (pre-existing bug; blocked
  System_dynamics + iswap_like).
- `pulse_gen/pulse_config.py::PulseConfig.__eq__`/`__lt__` — return `NotImplemented` for
  non-PulseConfig operands (a notebook helper compared a PulseConfig to a Qobj; standard
  Python-correctness fix, general library improvement).
- `Tutorial/transmon_tutorial/iswap_like.ipynb` — two pre-existing NOTEBOOK bugs:
  (a) scan loop never did `result.append(y)` → empty `imshow`; (b) qutip-4 idiom
  `uQobj.data[i,j]` → `uQobj.full()[i,j]` (qutip-5 `Qobj.data` is a `Data` object).
- `system/arb_qubit_system.py::ArbQubitSys.__init__` — **argument-order fix**. `extra_list`
  was the 4th positional param (`freq, inter, r, extra, gamma, driving, bias`), but every
  positional caller (`arb_debug.py`, all `arb_qubit_demo.ipynb` calls) uses
  `(freq, inter, r, gamma, driving, bias)`. The shift silently put the Z-bias dicts into
  `self.driving_list`, so the XY-drive operator built to all-zeros → **no drive / wrong
  transitions**. Moved `extra_list` to the end (keyword-only in practice; the only user,
  archive `arb_Ramsy.ipynb`, passes it by keyword). Verified: `arb_debug.py` block 1
  reproduces the 000→200 transfer (P₂₀₀ 0.378) and block 2 the 000↔100 Rabi (ρ₁₁ 0.905).
  NB: the signature was byte-identical to the qutip-4 original, so this predates (is not
  caused by) the migration; it also unblocks the constructor calls in `arb_qubit_demo.ipynb`.

## Open items — PRE-EXISTING bugs, NOT qutip-5 (need a decision)
1. **`pulse_gen/pulse_buffer.py::merge_pulse_chan`** — `int(np.intersect1d(index_type,
   index_qi))` crashes on multi-channel pulse sequences (empty/multi intersection; also
   `np.where` returns a truthy tuple so the guard misfires). Version-independent; fails
   *before* any solver call. **Blocks `System_dynamics.ipynb` and `iswap_like.ipynb`.**
   Intended behavior: find the buffer row matching BOTH pulse_type and qindex, accumulate
   `Hd_i[1]` there, else append a new row (a simple `zip` loop over the three parallel lists).
2. **`Tutorial/arb_qubit_tutorial/arb_qubit_demo.ipynb`** — targets a stale ArbQubitSys API
   (dict configs, old constructor order); fails at cell ~9 before any solver call. Needs a
   notebook update, not a library change.
3. **`instruments/angle.py::cal_angle`** — dead code that would break under qutip 5 if wired up.

## Suggested next steps
- Decide whether to fix #1 (unblocks 2 notebooks; clear, low-risk logic fix) now or later.
- #2 requires rewriting the arb demo notebook against the current API.
- Consider a follow-up numerical-equivalence check vs a qutip-4 baseline for a couple of gates.
