# N-Body Simulation Prototype Notebook

This folder holds the Python prototype of the dissertation project *Comparative Reliability Performance of N-Body Simulation under Hamiltonian Framework Applying Barnes-Hut Trees*. The notebook `main.ipynb` implements the gravitational N-body problem with direct summation and with the Barnes-Hut octree, six integration schemes (three symplectic, three non-symplectic), the conserved quantities used to judge them, and a JSON-driven experiment pipeline. It is a prototype: every feature is validated here in Python before it is transferred to the C++ and R framework.

This document explains how the notebook is organised, how its components fit together, and how to describe and run an experiment with a JSON file.

---

## Contents

1. [Folder layout](#1-folder-layout)
2. [Requirements and running order](#2-requirements-and-running-order)
3. [Structure of the notebook](#3-structure-of-the-notebook)
4. [Conventions](#4-conventions)
5. [Core components](#5-core-components)
6. [The hand-written experiments](#6-the-hand-written-experiments)
7. [The automated pipeline](#7-the-automated-pipeline)
8. [Logger](#8-logger)
9. [Plotter](#9-plotter)
10. [JSON experiment format](#10-json-experiment-format)
11. [Runtime](#11-runtime)
12. [Example experiments](#12-example-experiments)
13. [Extending the notebook](#13-extending-the-notebook)
14. [Known limitations](#14-known-limitations)

---

## 1. Folder layout

```
notebook/
├── main.ipynb              the notebook (theory, implementation, experiments)
├── README.md               this document
└── experiments/            JSON experiment files read by run_experiment
    ├── figure_eight.json           Experiment 1 as a JSON file (VV4, closure check)
    ├── figure_eight_euler.json     the same orbit with Forward Euler, reduced plot selection
    ├── plummer_barnes_hut.json     one run of Experiment 2 (Plummer sphere, Barnes-Hut)
    └── hot_swap_random.json        the hot-swap schedule of Experiment 3
```

`notebook-markdown.css` and `futura-notebook-markdown/` are local VS Code font helpers and are excluded from version control by `.gitignore`.

---

## 2. Requirements and running order

The notebook requires Python 3 (developed on 3.12) with NumPy and Matplotlib. The standard-library modules `random`, `time`, `json`, `copy` and `pathlib` are also used.

The notebook is written to be run from top to bottom. The code cells fall into two groups:

- **Definition cells** contain functions only and run instantly: the imports and constants, the helper functions, the conservation functions, the six integrators, the Barnes-Hut octree, and the logger, plotter, loader and runtime cells.
- **Experiment cells** (Experiments 1 to 3) define a few experiment-specific functions and then run their simulations. Experiment 1 takes a few seconds; Experiments 2 and 3 take several minutes in pure Python.

The JSON runtime depends only on the definition cells, with one exception: the randomised initial-condition generators `plummer` and `random_sphere` reuse `init_plummer`/`virialise` (defined in the Experiment 2 cell) and `init_random_3` (defined in the Experiment 3 cell). An experiment that uses an explicit list of bodies needs none of the experiment cells. If a required generator is missing, the runtime stops with a message naming the cell to run.

Relative paths such as `"experiments/figure_eight.json"` are resolved from the working directory of the kernel, which in Jupyter and VS Code is the folder that contains the notebook.

---

## 3. Structure of the notebook

The notebook is kept in dated entry order. Each feature is introduced by a markdown cell (theory, equations and references), followed by the code cell that implements it.

| Section | Cell content | Main names defined |
|---|---|---|
| Introduction | Aims of the project and of the notebook, naming conventions | – |
| Environment | Imports, random seed (`seed = 256`), Yoshida constants | `yoshida_x0`, `yoshida_x1` |
| Helpers | Distances, initial velocities, direct-summation acceleration (Eq. 4) | `acc` |
| Conservation | Kinetic and potential energy, Hamiltonian, linear and angular momentum and their relative errors | `conservation_state`, `h_rel`, `p_rel`, `l_rel` |
| Velocity Verlet 4th order | Yoshida composition of three Velocity Verlet steps | `run_vv_4` |
| Experiment 1 | Three-body figure-eight choreography | `Body` |
| Barnes-Hut octree | Tree construction, mass moments, tree walk | `acc_bh`, `build_tree_bh` |
| Experiment 2 | Plummer cluster: direct summation against Barnes-Hut | `init_plummer`, `virialise`, `simulate` |
| Leapfrog | Kick-drift-kick leapfrog, 2nd order | `run_l_2` |
| Forest-Ruth | Forest-Ruth, 4th order | `run_fr_4` |
| Forward Euler | Explicit Euler, 1st order | `run_fe_1` |
| RK2 | Explicit midpoint method, 2nd order | `run_rk_2`, `acc_stage` |
| RK4 | Classical Runge-Kutta, 4th order | `run_rk_4` |
| Experiment 3 | Integrator hot-swap in a random 20-body system | `init_random_3`, `simulate_hot_swap` |
| Automation | Logger and plotter | `init_log`, `log_state`, `plot_log` |
| | JSON loader | `load_experiment`, `integrator_registry` |
| | Runtime | `run_experiment`, `run_experiments` |
| | Example run | – |
| References | Full reference list (Author, Year) | – |

---

## 4. Conventions

**Function names.** Integrators and their helper functions follow the pattern `operation_integrator_order`: for example `displ_vv_4` is the displacement update of Velocity Verlet 4th order, and `run_fr_4` the driver of Forest-Ruth 4th order. Barnes-Hut functions end in `_bh`. Mathematical quantities take their literal names (`numerator`, `denominator`, `product`).

**Bodies.** A body is any object with four attributes: mass `m` (scalar), displacement `r`, velocity `v` and acceleration `a` (each a NumPy array of three components). The class `Body` is the minimal implementation. It is declared in Experiment 1, again in Experiment 3 and again in the JSON loader cell, so that each of these cells can be run without the others; the three declarations are identical.

**Integrator signature.** Every integrator has the signature

```python
run_*(bod_obj, G, e, dt, n_steps, acc_func=acc)
```

and advances the list `bod_obj` *in place* by `n_steps` steps of size `dt`, evaluating forces through `acc_func`. Every integrator recomputes the accelerations from the positions on entry and leaves `r` and `v` at the same instant on exit. Consequently, consecutive calls to the same or to different integrators continue a run seamlessly, which is what makes the hot-swap of Experiment 3 exact.

**Force-routine signature.** A force routine has the signature `acc_func(G, bod_obj, e)` and returns a list of acceleration vectors, one per body. `acc_bh` takes the opening angle `theta` as an additional fourth argument, so it is wrapped in a lambda that fixes `theta`.

**Units.** All experiments use natural units with `G = 1`. Experiment 1 uses the dimensionless figure-eight units (unit masses); Experiments 2 and 3 use Hénon units (total mass `M = 1`, total energy `E = -1/4`, crossing time `2√2`).

**Time discretisation.** Every experiment is controlled by three numbers: the total simulated time `t_total`, the number of timesteps `n_steps` (which alone sets `dt = t_total / n_steps`), and the number of stored frames `n_frames` (which sets only the resolution of the record and of the plots). `n_steps` must be divisible by `n_frames`, so that the run ends exactly at `t_total`.

---

## 5. Core components

### 5.1 Forces

| Routine | Method | Cost | Notes |
|---|---|---|---|
| `acc(G, bod_obj, e)` | Direct summation over all pairs (Eq. 4) | O(N²) | Reference force, exact up to rounding |
| `acc_bh(G, bod_obj, e, theta)` | Barnes-Hut octree (Barnes & Hut, 1986) | O(N log N) | `theta = 0` reproduces direct summation |

Both routines use Plummer softening with length `e`. The softening must be the same in the force routine and in the potential energy `v_abs`, as otherwise `H` is not the quantity the integrator actually conserves.

### 5.2 Conserved quantities

The conservation cell provides the absolute quantities and their relative errors with respect to the initial state:

| Quantity | Absolute | Relative error |
|---|---|---|
| Hamiltonian | `H = T + V` (`t_abs`, `v_abs`, `h_abs`) | `h_rel = abs(H_n − H_0) / abs(H_0)` |
| Linear momentum | `P = Σ m_i v_i` (`p_abs`) | `p_rel = abs(P_n − P_0) / P_ref` |
| Angular momentum | `L = Σ m_i (r_i × v_i)` (`l_abs`) | `l_rel = abs(L_n − L_0) / L_ref` |

`P_ref` and `L_ref` equal `|P_0|` and `|L_0|` unless the initial momentum is (numerically) zero, as it is in any centre-of-mass frame. In that case the relative error would be undefined, and the characteristic sizes `p_scale = Σ m_i |v_i|` and `l_scale = Σ m_i |r_i × v_i|` are used instead. `conservation_state(bod_obj, G, e)` returns the triple `(H, P, L)` in one call. All conserved quantities are always evaluated with the exact pair sums, never with the tree, so Barnes-Hut runs are judged against the true Hamiltonian.

### 5.3 Integrators

| Key | Function | Scheme | Order | Symplectic | Force evaluations per step |
|---|---|---|---|---|---|
| `L2` | `run_l_2` | Leapfrog (kick-drift-kick) | 2 | yes | 1 |
| `FR4` | `run_fr_4` | Forest-Ruth (Forest & Ruth, 1990) | 4 | yes | 3 |
| `VV4` | `run_vv_4` | Velocity Verlet, Yoshida composition (Yoshida, 1990) | 4 | yes | 3 |
| `FE1` | `run_fe_1` | Forward Euler | 1 | no | 1 |
| `RK2` | `run_rk_2` | Explicit midpoint (Runge, 1895) | 2 | no | 2 |
| `RK4` | `run_rk_4` | Classical Runge-Kutta (Kutta, 1901) | 4 | no | 4 |

The keys in the first column are the names used in JSON experiment files.

---

## 6. The hand-written experiments

**Experiment 1: figure-eight (28/10/26).** Three equal masses on the periodic figure-eight orbit (Moore, 1993; Chenciner & Montgomery, 2000), integrated with VV4 for three periods. The orbit is known exactly, so the closure error `|r(T) − r(0)|` after a whole number of periods measures the integration error directly.

**Experiment 2: Barnes-Hut against direct summation (30/10/26).** A 64-body Plummer sphere in virial equilibrium, integrated with VV4 three times with identical initial conditions: once with direct summation and once each with Barnes-Hut at `theta = 0.5` and `theta = 1.0`. Any difference between the runs is caused by the force approximation alone.

**Experiment 3: integrator hot-swap (05/11/26).** A randomised 20-body system with unequal masses. The integrator is replaced at the end of every crossing time within one continuous run (`L2 → RK2 → FR4 → RK4 → VV4 → FE1`), with Barnes-Hut and with direct summation. An identity-swap test first confirms that the swap itself introduces no error.

These cells are kept as they are: they document the development of the project and contain diagnostics (tree accuracy, per-segment tables, run comparisons) that go beyond a single run. New experiments should use the JSON pipeline.

---

## 7. The automated pipeline

The three experiments repeat the same pattern: set the parameters, build the bodies, loop over frames while recording the state, and plot. The automation cells factor this pattern out:

```
experiments/my_experiment.json
        │
        ▼
load_experiment ── checks the file, fills in defaults, computes dt and the frame layout
        │
        ▼
build_bodies, build_force ── initial conditions and force routine
        │
        ▼
init_log ── creates the np.array table and logs the initial state (row 0)
        │
        ▼
loop over integrator segments and frames:
    run_*(bodies, G, e, dt, steps_per_frame, force)
    log_state ── one new row per frame
        │
        ▼
log_summary, closure report ── printed results
        │
        ▼
plot_log ── only the panels selected in the file
        │
        ▼
result dictionary (config, log, final bodies, initial bodies, runtime, figure)
```

Running Experiment 1 through this pipeline (`experiments/figure_eight.json`) reproduces the printout of the Experiment 1 cell to the last digit, and the reduced Plummer and hot-swap runs reproduce `simulate` and `simulate_hot_swap` bit for bit. The pipeline therefore changes the organisation of the code but not its results.

---

## 8. Logger

### 8.1 The log table

The log of one run is a dictionary. Its central entry `log["table"]` is a single two-dimensional NumPy array with **one row per stored frame** (`n_frames + 1` rows, including the initial state) and the following columns:

| Column(s) | Name(s) | Content |
|---|---|---|
| 0 | `t` | Simulated time of the frame (simulation units) |
| 1 | `H` | Absolute Hamiltonian `T + V` |
| 2 | `P_abs` | Absolute value `abs(P)` of the total linear momentum |
| 3 | `L_abs` | Absolute value `abs(L)` of the total angular momentum |
| 4 | `H_rel` | Relative Hamiltonian error (drift) with respect to `t = 0` |
| 5 | `P_rel` | Relative linear momentum error |
| 6 | `L_rel` | Relative angular momentum error |
| 7, 8, 9 | `x_1`, `y_1`, `z_1` | Displacement of body 1 |
| … | … | … |
| 3N+4 … 3N+6 | `x_N`, `y_N`, `z_N` | Displacement of body N |

The absolute and the relative values are computed in the same call and written into the same row, so they always refer to the same instant. The table is preallocated with `NaN` for the whole run, so logging never resizes the array, and any row that has not yet been written is easy to recognise.

A flat table of this shape maps one-to-one onto a CSV file, so saving it later only requires `np.savetxt` (see [Section 13](#13-extending-the-notebook)).

### 8.2 Other entries of the log

| Key | Content |
|---|---|
| `columns` | List of the column names (the header of the table) |
| `index` | Dictionary from column name to column number |
| `n_rows` | Number of rows written so far |
| `N`, `G`, `e` | Number of bodies, gravitational constant and softening used for every row |
| `H_0`, `P_0`, `L_0` | Reference values at `t = 0` |
| `P_scale`, `L_scale` | Fallback references for zero initial momenta |
| `segments` | List of `(first row, integrator key)` pairs, one per hot-swap segment |

### 8.3 Logger functions

| Function | Purpose |
|---|---|
| `init_log(G, e, bod_obj, n_frames, t=0.0)` | Create the log, store the references and log the initial state as row 0 |
| `log_state(log, bod_obj, t)` | Append the current state as the next row |
| `log_segment(log, label)` | Mark the start of an integrator segment (used by hot-swap schedules) |
| `log_rows(log)` | The rows written so far (a view of the table) |
| `log_column(log, name)` | One column by name, e.g. `log_column(log, "H_rel")` |
| `log_trajectories(log)` | Trajectories as an array of shape `(N, frames, 3)`, the layout used in Experiments 1 to 3 |
| `log_summary(log)` | Print `H_0`, `abs(P_0)`, `abs(L_0)` and the maximum and final relative errors |

The logger can also be used on its own, outside the runtime:

```python
log = init_log(G, e, bodies, n_frames)
for frame in range(1, n_frames + 1):
    run_vv_4(bodies, G, e, dt, steps_per_frame)
    log_state(log, bodies, frame * steps_per_frame * dt)
```

`init_log` raises an error if `H_0 = 0`, because the relative Hamiltonian error is undefined in that case.

---

## 9. Plotter

`plot_log(log, plots=None, title="", t_unit=1.0, t_unit_name="simulation time units")` draws a figure from a finished log. The figure has up to two rows:

- **Orbit row:** the 3D view and any of the 2D projections `xy`, `xz`, `yz`. The start of each orbit is marked `x` and the end `o`; all views share the same cubic limits, and the 2D views use an equal aspect ratio.
- **Conservation row:** the relative Hamiltonian drift `H_rel` and the relative momentum errors `P_rel` and `L_rel` against time, measured in the unit `t_unit` (for example orbital periods or crossing times). For hot-swap schedules each integrator segment is shaded and labelled.

Each row is a Matplotlib sub-figure with its own number of columns, so any combination of panels is laid out without gaps. A row with no selected panels is left out, and if nothing is selected `plot_log` returns `None` without drawing an empty figure.

### Plot options

| Key | Default | Meaning |
|---|---|---|
| `orbit_3d` | `true` | Orbits in 3D |
| `orbit_2d` | `["xy"]` | List of projection planes; `true` is shorthand for `["xy"]`, and `false` or `[]` switches the 2D views off |
| `h_rel` | `true` | Relative Hamiltonian drift |
| `p_rel` | `true` | Relative linear momentum error |
| `l_rel` | `true` | Relative angular momentum error |
| `max_bodies` | `10` | Largest number of orbits drawn as lines |
| `highlight` | `null` | Explicit list of body numbers to draw (counted from 1); overrides `max_bodies` |
| `log_scale` | `true` | Logarithmic y-axis on the conservation panels |

**Body selection.** When the system has more bodies than `max_bodies` and no `highlight` list is given, the `max_bodies` bodies that stay closest to the origin over the whole run are drawn, and the final positions of all bodies are shown as grey dots. Escapers would otherwise stretch the axes and squash the rest of the system. The choice is deterministic, so the same run always produces the same figure.

**Overlapping orbits.** With up to five drawn bodies, the line widths decrease from the first body to the last, so that bodies sharing one curve (as in the figure-eight) remain visible as nested bands.

**Logarithmic axes.** Exact zeros (the relative errors at `t = 0`) are drawn at a floor of `1e-18`, as in Experiments 1 to 3. The y-axis starts one decade below the smallest non-zero value, so the floor does not compress the actual error curve.

---

## 10. JSON experiment format

### 10.1 Overview

An experiment file is a JSON object with the following top-level blocks:

```json
{
  "name": "...",
  "description": "...",
  "seed": 256,
  "G": 1.0,
  "softening": 0.0,
  "integrator": "VV4",
  "force": {"method": "direct"},
  "time": {"duration": 3, "unit": 6.32591398, "unit_name": "orbital periods", "n_steps": 30000, "n_frames": 1500},
  "bodies": {"type": "explicit", "list": [...]},
  "plots": {"orbit_3d": true, "orbit_2d": ["xy"], "h_rel": true, "p_rel": true, "l_rel": true},
  "report": {"closure": true}
}
```

**Comments.** JSON has no comment syntax, so every key that starts with an underscore (for example `"_comment"`) is ignored by the loader in every block. These keys can be used to annotate a file.

**Strict keys.** Any other key that the loader does not recognise is rejected with an error that lists the allowed keys. A typo such as `"orbit3d"` therefore stops the run instead of silently falling back on the default value.

### 10.2 Top-level keys

| Key | Required | Default | Meaning |
|---|---|---|---|
| `name` | yes | – | Name of the experiment (printed and used in the plot title) |
| `description` | no | `""` | Free-text description, printed under the name |
| `seed` | no | `256` | Seed of the random generator; matters only for randomised initial conditions |
| `G` | no | `1.0` | Constant of gravitational attraction (> 0) |
| `softening` | no | `0.0` | Plummer softening length `e` (≥ 0), used for both the forces and `H` |
| `integrator` | yes | – | One integrator key (`"VV4"`), or a list of keys for a hot-swap schedule (`["L2", "RK2", "VV4"]`) |
| `force` | no | `{"method": "direct"}` | Force routine, see [10.3](#103-force) |
| `time` | yes | – | Time discretisation, see [10.4](#104-time) |
| `bodies` | yes | – | Initial conditions, see [10.5](#105-bodies) |
| `plots` | no | all defaults | Plot selection, see [Section 9](#plot-options) |
| `report` | no | `{}` | Extra printed checks, see [10.6](#106-report) |

**Hot-swap schedules.** With a list of integrators, the frames are divided into equal segments, one per entry and in the given order. All segments advance the same bodies in place, exactly as in Experiment 3. An integrator may appear more than once.

### 10.3 `force`

| `method` | Further keys | Meaning |
|---|---|---|
| `"direct"` | – | Direct summation `acc` |
| `"barnes_hut"` | `theta` (default `0.5`, ≥ 0) | Barnes-Hut octree `acc_bh` with opening angle `theta` |

### 10.4 `time`

| Key | Required | Default | Meaning |
|---|---|---|---|
| `duration` | yes | – | Length of the run in units of `unit` (> 0) |
| `unit` | no | `1.0` | Time unit in simulation units, e.g. the orbital period `6.32591398` or the crossing time `2.8284271247461903` |
| `unit_name` | no | `"simulation time units"` | Name of the unit, used in the printout and on the time axes |
| `n_steps` | yes | – | Total number of timesteps (positive integer) |
| `n_frames` | yes | – | Number of logged frames after the initial state (positive integer) |

The loader derives `t_total = duration × unit`, `dt = t_total / n_steps`, `steps_per_frame = n_steps / n_frames` and `frames_per_segment = n_frames / (number of integrators)`. Both divisions must be exact, and the loader rejects the file otherwise. This guarantees that the run ends exactly at `t_total` and that every integrator in a schedule covers the same simulated time with the same `dt`.

### 10.5 `bodies`

The key `type` selects one of three kinds of initial conditions.

**`"explicit"`: a list of bodies**

| Key | Required | Default | Meaning |
|---|---|---|---|
| `list` | yes | – | At least two bodies, each `{"m": mass, "r": [x, y, z], "v": [vx, vy, vz]}` |
| `centre_of_mass_frame` | no | `false` | Shift positions and velocities to the centre-of-mass frame, so that `P_0 = 0` and the system does not drift |

**`"plummer"`: a Plummer sphere** (requires the Experiment 2 cell)

| Key | Required | Default | Meaning |
|---|---|---|---|
| `N` | yes | – | Number of bodies (≥ 2) |
| `M` | no | `1.0` | Total mass |
| `a` | no | `null` | Plummer scale length; `null` uses `3πGM²/16`, which gives `E = −1/4` (Hénon units) |
| `r_cut` | no | `10.0` | Radius cut-off in units of `a` |
| `virialise` | no | `true` | Rescale the velocities to exact virial equilibrium `2T + V = 0` |

**`"random_sphere"`: a uniform random sphere with unequal masses** (requires the Experiment 3 cell)

| Key | Required | Default | Meaning |
|---|---|---|---|
| `N` | yes | – | Number of bodies (≥ 2) |
| `M` | no | `1.0` | Total mass |
| `R` | no | `1.2` | Radius of the sphere (`R = 6/5` gives `E ≈ −1/4` in equilibrium) |
| `mass_spread` | no | `0.5` | Masses uniform in `[1 − spread, 1 + spread]` before normalisation to `M` |
| `q_virial` | no | `0.5` | Virial ratio `T/abs(V)`; `0.5` is equilibrium, smaller values make the system collapse |

The random generator is re-seeded with `seed` before the bodies are built, so a randomised experiment produces the same bodies every time, whichever cells ran before.

### 10.6 `report`

| Key | Default | Meaning |
|---|---|---|
| `closure` | `false` | Print the closure error `abs(r(t_total) − r(0))` of every body (meaningful for periodic orbits run for a whole number of periods, as in Experiment 1) |

### 10.7 Validation

`load_experiment` checks the file before anything is simulated and stops with a descriptive error if:

- a required key is missing, or an unknown key is present (in any block);
- an integrator key, force method, bodies type or projection plane is not recognised;
- a physical parameter has the wrong sign or type (e.g. a negative `theta` or a non-integer `n_steps`);
- `n_steps` is not divisible by `n_frames`, or `n_frames` is not divisible by the number of integrators;
- a body vector does not have exactly three components.

`run_experiment` additionally checks that every `highlight` number refers to an existing body.

`load_experiment` also accepts a Python dictionary with the same structure as a file, which is convenient for quick tests directly in the notebook:

```python
result = run_experiment({
    "name": "Quick test",
    "integrator": "L2",
    "time": {"duration": 1.0, "n_steps": 1000, "n_frames": 100},
    "bodies": {"type": "explicit", "centre_of_mass_frame": True, "list": [
        {"m": 1.0, "r": [-0.5, 0.0, 0.0], "v": [0.0, -0.5, 0.0]},
        {"m": 1.0, "r": [0.5, 0.0, 0.0], "v": [0.0, 0.5, 0.0]}
    ]},
    "plots": {"orbit_3d": False}
})
```

---

## 11. Runtime

### `run_experiment(source, show=True)`

Runs one experiment from a JSON path or a dictionary:

1. loads and checks the experiment (`load_experiment`);
2. builds the initial conditions and the force routine, and keeps a copy of the initial state;
3. prints a header with everything needed to reproduce the run (N, G, softening, force routine, schedule, `t_total`, `dt`, frames);
4. creates the log and records the initial state;
5. for each integrator of the schedule, marks a new segment and advances the bodies frame by frame, logging one row per frame. The time of each row is computed from the frame counter rather than accumulated, so no rounding error builds up in the `t` column. Only the time spent inside the integrators is measured as the runtime, as in Experiments 2 and 3;
6. prints the runtime, the conservation summary and, if requested, the closure errors;
7. draws the selected plots (unless `show=False`).

It returns a dictionary:

| Key | Content |
|---|---|
| `config` | The complete, checked experiment description, with all defaults and derived quantities (`dt`, `t_total`, `steps_per_frame`, `schedule`) |
| `log` | The log of the run; `log["table"]` is the results table |
| `bodies` | The bodies in their final state |
| `initial` | An independent copy of the initial conditions |
| `runtime` | Wall-clock time spent inside the integrators, in seconds |
| `figure` | The Matplotlib figure, or `None` |

Assign the result to a variable (`result = run_experiment(...)`). If a bare `run_experiment(...)` is the last line of a cell, Jupyter prints the whole dictionary, including the table.

### `run_experiments(sources, show=True)`

Runs several experiments in turn and returns their results in a dictionary keyed by experiment name. To run every file in the folder:

```python
results = run_experiments(sorted(Path("experiments").glob("*.json")))
```

---

## 12. Example experiments

| File | Bodies | Integrator | Force | Plots | Approx. runtime |
|---|---|---|---|---|---|
| `figure_eight.json` | explicit, 3 | VV4 | direct | all, closure report | a few seconds |
| `figure_eight_euler.json` | explicit, 3 | FE1 | direct | x-y projection and `H_rel` only | about 1 s |
| `plummer_barnes_hut.json` | Plummer, 64 | VV4 | Barnes-Hut, θ = 0.5 | 3D, x-y and x-z projections, all errors; 6 orbits | about 30 s |
| `hot_swap_random.json` | random sphere, 20 | L2 → RK2 → FR4 → RK4 → VV4 → FE1 | Barnes-Hut, θ = 0.5 | 3D, x-y, all errors (segments shaded); 8 orbits | about 15 s |

`figure_eight_euler.json` shows how to exclude plots: every panel except the x-y projection and the Hamiltonian drift is switched off in its `plots` block.

### Writing a new experiment

1. Copy the example closest to what you need into `experiments/` under a new name.
2. Change `name` and `description`, and use `_comment` keys to record the purpose of the run.
3. Choose the `integrator` (or a schedule), the `force` and the `bodies`.
4. Set `time`: choose a convenient `unit`, then `duration`, `n_steps` (accuracy) and `n_frames` (resolution of the record). Keep `n_steps` divisible by `n_frames`, and `n_frames` divisible by the number of integrators.
5. Select the plots.
6. Run it with `result = run_experiment("experiments/<file>.json")`.

---

## 13. Extending the notebook

**A new integrator.** Implement `run_<name>_<order>(bod_obj, G, e, dt, n_steps, acc_func=acc)` following the conventions in [Section 4](#4-conventions): advance in place, recompute `a` on entry, and leave `r` and `v` at the same instant. Then add an entry to `integrator_registry` in the loader cell. The new key is then available in every JSON file and in hot-swap schedules.

**A new kind of initial conditions.** Add an entry `"<type>": ([required keys], {optional keys: defaults})` to `schema_bodies`, and a matching branch in `build_bodies` that returns a list of `Body` objects.

**A new force routine.** Add an entry to `schema_force` and a branch in `build_force` that returns a function with the signature `(G, bod_obj, e)`.

**A new logged quantity.** Append its name to `log_scalar_columns` and its value to the scalar row written in `log_state`. The displacement columns move along automatically, because they always start after the scalar columns.

**Saving the log to CSV.** The table is already in CSV shape:

```python
np.savetxt("results.csv", log_rows(log), delimiter=",", header=",".join(log["columns"]), comments="")
```

---

## 14. Known limitations

- **Pure Python.** Forces are evaluated with Python loops, so runs with more than about a hundred bodies, or with many steps, are slow. This is acceptable for a prototype; performance is the task of the C++ implementation.
- **Generators live in experiment cells.** The `plummer` and `random_sphere` initial conditions depend on functions defined in the Experiment 2 and 3 cells (see [Section 2](#2-requirements-and-running-order)).
- **One run per file.** A JSON file describes one run. Comparisons, such as direct summation against Barnes-Hut in Experiment 2, are made by running several files and comparing their logs.
- **In-memory log only.** The log is not saved to disk yet. Saving takes one call (see [Section 13](#13-extending-the-notebook)).
- **Memory.** The table holds `(n_frames + 1) × (7 + 3N)` floating-point numbers, which is about 7 MB for 1000 frames of 300 bodies. Keep `n_frames` moderate for large N.
- **Fixed timestep.** All integrators use a constant `dt`, so close encounters are only tamed by the softening length.
