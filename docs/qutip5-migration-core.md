# qutip 4 → 5 core-solver migration (Phase 1, sections A & B)

Date: 2026-08-09 · Env: `qusim_qutip5` (Python 3.10.16, **qutip 5.1.1**, **numpy 2.2.1**,
scipy 1.15.1). Interpreter `E:\conda\envs\qusim_qutip5\python.exe`, run with `PYTHONPATH=src`.

Scope: `system/transmon_system.py`, `system/arb_qubit_system.py`, `instruments/angle.py`,
`instruments/tunablec.py`, `pulse_gen/{simulation_option,pulse_buffer}.py`. No other files touched.

All changes were **verified by execution** under qutip 5.1.1 — a real `mesolve` AND a real
`propagator` run on small (2-qubit) TransmonSys and ArbQubitSys systems, plus four notebooks
executed end-to-end. Evidence at the bottom.

---

## 1. Changes per file

### `system/transmon_system.py`
- Added module constant `DEFAULT_SOLVER_OPTIONS = {'rtol': 1e-8}`.
- Three solver methods: default arg `option = Options(rtol=1e-8)` → `option = None`.
  Each call site resolves `options = option if option is not None else DEFAULT_SOLVER_OPTIONS`.
- `master_eq_solver`: `mesolve(..., c_ops=self.co_list(), options=options)` (options now a dict).
- `system_dynamics_propagator`: **before**
  `propagator(H_d, sim_opts.tlist, self.co_list(), {}, option)` → **after**
  `propagator(H_d, sim_opts.tlist, self.co_list(), {}, options, tlist=sim_opts.tlist)`.
  The `tlist=` keyword is **required** in qutip 5 for the array-coefficient Hamiltonian
  (see §3). No section-B change needed here: `get_data_list`'s `np.abs(bra*ket)**2` works
  because in qutip 5 `bra*ket` returns a Python complex scalar (verified).

### `system/arb_qubit_system.py`
- Added `DEFAULT_SOLVER_OPTIONS = {'rtol': 1e-8}`; three `Options(rtol=1e-8)` defaults → `None`
  guard; removed the `option: Options` type annotation.
- `master_eq_solver`: **before** `mesolve(H_d, initial_state, tlist, self.co_list, [], options=option)`
  → **after** `mesolve(H_d, initial_state, tlist, self.co_list, options=options)`.
  The deprecated positional `[]` e_ops slot is dropped (empty e_ops is the qutip-5 default and
  e_ops is keyword-only). `self.co_list` is left **without `()`** — it is a `@property` (see §2).
- `system_dynamics_propagator`: **before**
  `propagator(H_d, sim_opts.tlist, self.co_list, {}, option, parallel=do_parallel, progress_bar=do_progress_bar)`
  → **after**
  ```python
  options = option if option is not None else DEFAULT_SOLVER_OPTIONS
  if do_progress_bar:
      options = dict(options)
      options['progress_bar'] = do_progress_bar if isinstance(do_progress_bar, str) else 'tqdm'
  result = propagator(H_d, sim_opts.tlist, self.co_list, {}, options, tlist=sim_opts.tlist)
  ```
  (`parallel=`/`progress_bar=` removed; see §4 for the do_parallel/do_progress_bar decision.)
- **numpy-2 fix (blocker):** `hermitian_conjugate` used `np.array(matrix, dtype=np.complex_)`.
  `np.complex_` was **removed in numpy 2.0** → changed to `np.complex128`. This blocked
  *construction* of every ArbQubitSys (called from `get_h_inter`/`get_h_extra`/`get_v`/
  `get_h_xy_drive`). (Inventory §E missed this — it checked `np.complex` but not `np.complex_`.)
- **Qobj indexing (section B):**
  - `get_state_index` L302 **before** `np.abs(arr[state_index][0][0])` → **after**
    `np.abs(arr[state_index][0])`. In qutip 5 Qobj indexing drops one nesting level
    (`arr[i]` → 1-elem array, `arr[i][0]` → scalar), so `[0][0]` over-indexed. Preserves the
    original scalar intent. (transmon's analogue at L118 uses `arr[state_index]` and worked.)
  - `get_state_index` degeneracy branch **before** `np.abs(row[0][0])` → **after** `np.abs(row[0])`
    (iterating a ket yields 1-elem arrays; `row[0]` is the scalar).
  - `get_data_list` L631 **before**
    `np.abs(((result_list.states[ii]).dag()*state)[0][0][0])**2` → **after**
    `np.abs((result_list.states[ii]).dag()*state)**2`. `bra*ket` is a complex scalar in qutip 5,
    so `[0][0][0]` raised `TypeError: 'complex' object is not subscriptable`.

### `instruments/tunablec.py`
- `fsim_zhen` **before**:
  ```python
  p0 = state_001.dag()*result_list[0].states[-1].data
  p1 = state_100.dag()*result_list[0].states[-1].data
  ...
  theta = np.arctan(p1/p0)[0][0]
  ```
  **after** (drop `.data`, drop trailing `[0][0]`):
  ```python
  p0 = state_001.dag()*result_list[0].states[-1]
  p1 = state_100.dag()*result_list[0].states[-1]
  ...
  theta = np.arctan(p1/p0)
  ```
  In qutip 5 `.data` is a `qutip.core.data.Data` object, and `Qobj * Data` is unsupported
  (`TypeError`). `bra * ket` now yields a complex scalar, so `np.abs(...)` gives a float and the
  `[0][0]` indexing is invalid. Numerical intent (θ = arctan(|⟨100|ψ⟩| / |⟨001|ψ⟩|)) preserved.
  Exercised by `iswap_like.ipynb`; verified via a real mesolve result (θ finite, see evidence).

### `instruments/angle.py` — **no change**
- `get_angle` (the used function) is fine: `np.angle(tstate.dag()*spt)` returns a float scalar in
  qutip 5, and `isket(...)` is still valid. Verified: returns floats for ket results, `None` for
  density-matrix results.
- `cal_angle` is **dead code** (defined, never called anywhere in `src/`). It *would* raise under
  qutip 5 (`np.angle(bra*ket)[0][0]` indexes a float scalar), but per "change only if it errors in
  practice / keep changes minimal" it is left untouched and flagged in §5.

### `pulse_gen/simulation_option.py`, `pulse_gen/pulse_buffer.py` — **no change**
- No qutip-5 breakage in either (they only `from qutip import Qobj` / `import *`). See §5 for a
  *pre-existing, out-of-scope* logic bug found in `pulse_buffer.merge_pulse_chan`.

---

## 2. The `co_list` investigation (inventory A said "latent bug")

**Finding: there is NO bug in the current code — it was already resolved by a `@property`.**

- `arb_qubit_system.py`: `co_list` is decorated `@property` (L386-387). Therefore `self.co_list`
  (used **without** `()` in mesolve/propagator) correctly evaluates to the collapse-operator list.
  Verified at runtime: `type(ArbQubitSys.__dict__['co_list']) is property`, and `sysA.co_list`
  returns `list[Qobj]` (or `[]` when `gamma_list is None`).
- `transmon_system.py`: `co_list` is a plain method (no decorator), so `self.co_list()` **with**
  `()` is correct. Verified: `sysT.co_list()` returns `list[Qobj]`.

The two classes use different mechanisms but are each internally consistent and both pass the
actual collapse-operator list as `c_ops`. The inventory's note (co_list "def at L387… passes a
bound method") predates the `@property`. **No change made**; I left `self.co_list` (property) in
arb and `self.co_list()` (method) in transmon, and proved both feed real collapse operators to
`mesolve` (density-matrix output with `tr == 1`) — see evidence.

---

## 3. propagator + array-coefficient Hamiltonian (`tlist=` is mandatory)

The time-dependent Hamiltonian is `[H0, [H1, coeff_array], ...]` sampled on `tlist`. `mesolve`
receives `tlist` positionally, so it interpolates fine. **`propagator` does not** — it builds
`QobjEvo(H, args=args, **kwargs)` and, without the sample grid, raises
`ValueError: tlist must be the same len as the array to interpolate`.

qutip 5's `propagator` docstring: *"the output times in `t` are not used for array time dependent
system. `tlist` must be passed as a keyword argument."* Fix = pass `tlist=sim_opts.tlist` as a
keyword (positional `t` still selects output times; here they are the same grid). Verified safe
both with empty `c_ops` (returns list of unitary Qobj) and with collapse operators (returns list
of superoperator Qobj).

---

## 4. do_parallel / do_progress_bar decision (arb propagator)

qutip 5 `propagator(H, t, c_ops=None, args=None, options=None, **kwargs)` has **no** `parallel`
and **no** `progress_bar` parameter.
- **`do_progress_bar`** → mapped into the options dict: when truthy, set
  `options['progress_bar'] = 'tqdm'` (or the caller's string if one is given). The options dict is
  copied first so a shared/default dict is never mutated. Verified: `do_progress_bar=True` prints a
  tqdm bar and runs with no kwarg error.
- **`do_parallel`** → **kept in the signature for backward-compatible call sites but is a documented
  no-op.** qutip 5's public `propagator` exposes no parallelism knob (parallel maps belong to the
  lower-level `Propagator`/`parallel_map` API, which this method does not use). A comment records
  this so callers passing `do_parallel=...` don't silently expect speedups.

---

## 5. Verification (acceptance test)

**nbconvert is NOT installed** in `qusim_qutip5` (and `pip install` is disallowed), so notebooks
were executed with a headless runner that strips IPython magics, forces the Agg matplotlib
backend, and `exec`s code cells. Notebook outputs were never written back over the originals.

1. **Import / FutureWarning gate (PASS):**
   `PYTHONPATH=src python -c "import warnings; warnings.simplefilter('error', FutureWarning); import qusim; print('ok')"`
   → `ok`. No `Options` FutureWarning. `py_compile` of all scoped modules OK.

2. **Core solver script `verify_core.py` — 13/13 assertions PASS** (real qutip 5.1 runs):
   - TransmonSys (2q, dim 3): `co_list()` → `list[Qobj]`; `mesolve` **with collapse ops** →
     density matrices (`isoper`, `tr=1.0000`); `mesolve` pure-state (default options + a
     user-supplied `{'rtol':1e-9,'atol':1e-9}` dict); `get_data_list` finite; **`propagator`**
     → list of 200 (9×9) Qobj; `tunablec.fsim_zhen` → θ=1.4663 (finite).
   - ArbQubitSys (2q, dim 2): `co_list` **property** → `list[Qobj]`; `mesolve` with collapse ops
     (c_ops from the property) → density matrices; `mesolve` pure-state; `get_data_list` finite
     (after `[0][0][0]` removal); **`propagator`** with `do_parallel=True` no-op → list of (4×4)
     Qobj; `propagator` with `do_progress_bar=True` → tqdm bar, no kwarg error.

3. **Notebooks executed end-to-end (headless):**
   - `Setup_system.ipynb` — **PASS** (6 cells).
   - `ZZ_coupling.ipynb` — **PASS** (6 cells).
   - `Energy_level.ipynb` — **PASS** (7 cells, 800-pt eigen scan, ~8 s).
   - `DRAG.ipynb` — **PASS** (13 cells; real `system_dynamics_mesolve` end-to-end, ~5 s).
   - `System_dynamics.ipynb` — **FAILS at cell 15**, but **not in solver code**: it raises in
     `pulse_gen/pulse_buffer.merge_pulse_chan` (`int(np.intersect1d(...))`) *before* any `mesolve`
     call. Pre-existing, out-of-scope bug — see below.
   - `iswap_like.ipynb` — **not run in full**: it writes to the real `../../Data/tunablec/` dir and
     runs 200-point mesolve scans. Its migrated code paths (`system_dynamics_mesolve`,
     `system_dynamics_propagator`, `tunablec.fsim_zhen`, `.tr()`) are all covered by `verify_core.py`.
   - `arb_qubit_demo.ipynb` — **FAILS at cell 9**, *before* any solver call: it passes a plain
     `dict` where `SimulationOption` is expected (`'dict' object has no attribute 'tlist'`) and
     builds `ArbQubitSys(...)` with a stale positional argument order. This is **API drift unrelated
     to qutip 4→5** (owned by the notebook/host agent). The arb solver itself is proven working via
     `verify_core.py` using the current `PulseConfig`/`SimulationOption` API.

4. **Grep gates (PASS):** zero `Options(` in any scoped file; the only `do_parallel`/
   `do_progress_bar` matches are the kept parameter *definitions*; no `parallel=`/`progress_bar=`
   kwargs on any `propagator` call.

---

## 6. Remaining risks / notes for the host

1. **`pulse_buffer.merge_pulse_chan` multi-channel bug (pre-existing, blocks full `System_dynamics`
   and `iswap_like`).** `if index_type and index_qi:` is *always* true when both a matching
   pulse-type slot and a matching qindex slot exist anywhere in the buffer (they are non-empty
   `np.where` tuples). When those are *different* slots (e.g. adding Z-on-q2 while the buffer holds
   XY-on-q2 and Z-on-q0), `np.intersect1d` is empty and `int(empty)` raises
   `TypeError: only length-1 arrays can be converted to Python scalars`. This is independent of
   qutip/numpy version and out of my scope (not a qutip-5 breakage). Suggested fix for the owning
   agent: `idx = np.intersect1d(index_type[0], index_qi[0]); merge only if idx.size == 1, else append`.
2. **`arb_qubit_demo.ipynb` targets a stale ArbQubitSys API** (dict configs, old constructor order).
   Needs updating to `PulseConfig`/`SimulationOption` objects and the current
   `system_dynamics_mesolve(pseq, sim_opts, ...)` signature. Notebook-owner task.
3. **`angle.cal_angle` is dead but qutip-5-broken** (`np.angle(bra*ket)[0][0]` on a float). Harmless
   while unused; fix if it is ever wired up.
4. **`get_data_list` (both classes) assumes pure-state (ket) results.** With collapse operators
   `mesolve` returns density matrices and `np.abs(bra·ket)**2` no longer computes a population
   (would need `expect(state*state.dag(), rho)`). Pre-existing design limitation, not a qutip-5
   regression, and not combined with `get_data_list` in any active notebook.
