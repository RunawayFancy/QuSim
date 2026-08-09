# qutip 4 → 5 migration: `data_plot/plot_tomo.py`

Date: 2026-08-09 · Scope: **only** `src/qusim/data_plot/plot_tomo.py`.
Addresses section C of `docs/qutip5-migration-inventory.md`.

## Problem
Under qutip 5.1 the module failed to import with `No module named 'qutip.qobj'`.
Three qutip-4 *internal* module paths used by the imports were removed in qutip 5:
`qutip.qobj`, `qutip.superoperator`, `qutip.superop_reps`.

## What changed (imports only — the `hinton(...)` body is untouched)
Replaced the four old import lines (L24–L28) with:

```python
from qutip import Qobj, vector_to_operator, settings
from qutip.core.superop_reps import (
    _to_superpauli as _super_to_superpauli,
    isqubitdims as _isqubitdims,
)
```

Nothing else in the file was modified. The public `hinton(...)` signature and
all internal call sites (`_super_to_superpauli(rho)`, `_isqubitdims(rho.dims)`,
`vector_to_operator(...)`, `settings.colorblind_safe`, `Qobj`) are unchanged
because the imports are re-aliased to the original names.

## Decision: reuse maintained qutip-5 symbols, do NOT reimplement or delegate
I introspected qutip 5.1 first (as instructed) and found drop-in equivalents,
so reimplementation was unnecessary and reuse is the smallest correct change:

- `Qobj`, `vector_to_operator`, `settings` are still **top-level** in qutip 5
  (verified: `hasattr(qutip, ...)` all True; `settings.colorblind_safe` present).
- The two superoperator helpers survive under **`qutip.core.superop_reps`** (only
  the old top-level `qutip.superop_reps` *path* is gone):
  - `_isqubitdims` → **`isqubitdims`** (now a *public* function; same
    "all dims are powers of two" qubit check).
  - `_super_to_superpauli` → renamed **`_to_superpauli`**; source confirmed to be
    the same algorithm (`to_super` → `isqubitdims` guard → project into the
    normalized super-Pauli basis `B.dag() @ sqobj @ B`). Output Qobj has identical
    dims/shape/real entries, so `sqobj.full().T` in `hinton` behaves as before.

I deliberately did **not** delegate to the built-in `qutip.hinton`. It has an
incompatible signature/output for this module — it uses `x_basis`/`y_basis`
instead of `xlabels`/`ylabels`, has no `title` or `xtk_rotate`, and manages the
colorbar differently — so delegating would not preserve this module's public API.

`_to_superpauli` remains private in qutip 5; if a future qutip release renames it
again, only this single import line needs revisiting. That risk is preferred over
duplicating qutip's Pauli-basis math here.

## Verification (env `qusim_qutip5`, qutip 5.1.1, numpy 2.2.1, `PYTHONPATH=src`)

Command (Bash tool):

```
PYTHONPATH=src "/e/conda/envs/qusim_qutip5/python.exe" -c "
import matplotlib; matplotlib.use('Agg')
import qutip, matplotlib.pyplot as plt
import qusim.data_plot.plot_tomo as pt
from qusim.data_plot.plot_tomo import hinton
rho  = qutip.rand_dm(4)                                   # single-index dm
rho2 = qutip.tensor(qutip.rand_dm(2), qutip.rand_dm(2))   # 2-qubit dm
S    = qutip.to_super(qutip.sigmax())                     # qubit superoperator (Pauli path)
opk  = qutip.operator_to_vector(qutip.rand_dm(2))         # operket (vector_to_operator)
Sq   = qutip.to_super(qutip.qeye(3))                      # qutrit super -> must ValueError
for label, obj in [('4x4', rho), ('2-qubit', rho2), ('super', S), ('operket', opk)]:
    fig, ax = hinton(obj); print(label, type(fig).__name__, type(ax).__name__); plt.close('all')
try:
    hinton(Sq); print('FAIL')
except ValueError as e:
    print('qutrit rejected:', str(e)[:40])
"
```

Output:

```
1. import qusim.data_plot.plot_tomo OK
2. hinton 4x4 ok Figure Axes
2b. hinton 2-qubit dm ok Figure Axes
2c. hinton superoperator (Pauli path) ok Figure Axes
2d. hinton operket (vector_to_operator) ok Figure Axes
2e. qutrit superoperator correctly rejected: Hinton plots of superoperators are currently only
ALL CHECKS PASSED
```

`hinton` returns `(matplotlib.figure.Figure, matplotlib.axes.Axes)` for operators,
operkets and qubit superoperators, and still raises `ValueError` for non-qubit
superoperators. The only warning printed is qutip's own `Options` `FutureWarning`
emitted on `import qutip` (from `qutip/solver/options.py`); it is unrelated to this
module and out of scope here.
