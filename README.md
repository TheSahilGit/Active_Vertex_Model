# Vertex model with cell activity and Rho–ROCK–myosin feedback

Simulation code accompanying

> S. Islam, A. Gupta and M. S. Rizvi, *Length scale of cellular activity determines signatures of epithelial remodeling* (submitted to Biophysical Journal).

A two-dimensional vertex model of a confluent epithelium, written in Fortran 90. Cells are polygons that share vertices and edges. Vertices move by overdamped dynamics under the forces of the vertex-model energy and of the cell activity, and the tissue remodels through T1 (neighbour exchange) and T2 (cell extrusion) transitions.

This repository is meant to show how the model is implemented. The parameter file shipped here is a basic example, not the settings of the simulations in the paper; the parameter values used there are given in the paper.

## Model

Energy (λ, A₀ set the units of energy and length; time is in units of η/(λA₀)):

U = Σ_c [ λ (A_c − A₀)² + β_c P_c² + γ P_c ]

with A_c and P_c the area and perimeter of cell c, λ the area elasticity, β_c the contractility of cell c (`beta` for a passive cell) and γ the interfacial tension (`gamm`) (Eq. 1 of the paper; `Calculate_Energy` in `src/System_Info.f90`). Every vertex v moves as η dr_v/dt = −∇_v U + f_v^active, integrated with an explicit Euler step of size `dt`.

### Cell activity

All activity types are switched on in `para_Simulation.dat`. They are implemented in `src/Force.f90`; see `PARAMETERS.md` for every parameter.

| Activity (paper) | Switch | Strength / parameters | Implementation |
|---|---|---|---|
| Active motility (AM): independent random force on every vertex | `if_motility` | `etas_max` | `Motile_Force_Calculation` |
| Polar motility (PM): one random force per cell, applied to its vertices, redrawn every step | `if_polar_motility` | `polar_motility_strength` | `Polar_Motile_Force_Calculation` |
| Persistent polar motility (PPM, active Brownian): self-propulsion v₀(cos θ_c, sin θ_c) per cell with rotationally diffusing θ_c | `if_ABP` | `vo`, `Dr` (persistence time 1/`Dr`), `abp_force_mode` | `ABP_Force_Calculation` |
| Fluctuating contractility (FC): per-cell β_c following an Ornstein–Uhlenbeck process about β₀ = `beta` | `if_active_contractility` | `active_contr_strength`, `tau_contr`; optional `fc_clip_negative` keeps β_c ≥ 0 | `Solve_Myosin_Act_Contr` |
| Mechanochemical regulation (MR): area-driven Rho → ROCK → myosin kinetics per cell, β_c = (α·Myosin + coupling noise)/(λA₀) | `if_RhoROCK` (+ `if_coupling_noise`) | `nhill`, `K_hill`, `A_*`, `D_*`, `Myosin_Coupling_Strength` (α), `coupling_noise_strength` (σ_MR), `mr_noise_mode` | `Solve_RhoROCK_Euler` / `_RK4`, `compute_RHS`, `Force_Calculation` |
| MR control without feedback: β_c = (constant + coupling noise)/(λA₀) | `mr_feedback_off` | `mr_beta_mean` | `Force_Calculation` |

### Tissue, boundaries and protocols

| Feature | Switch | Implementation |
|---|---|---|
| T1 transitions (edge shorter than `min_d_T1`) | `if_Do_T1` | `src/T1_swap.f90` |
| T2 transitions (cell area below `min_area_T2`) | `if_Do_T2` | `src/T2_swap.f90` |
| Open tissue (free boundary) or periodic box | `if_PBC` (mesh: `if_periodic` in `para_MeshGen.dat`) | `src/Geometry.f90` |
| Step (sudden) or oscillatory shear; Lees–Edwards boundary conditions in the periodic box | `if_Shear_tissue`, `if_Sudden_Shearing`, `if_Oscillatory_Shearing` | `src/Stress.f90` |
| Periodic box whose size (and, optionally, aspect ratio) relaxes under the tissue stress | `if_variable_box`, `box_aniso` | `Variable_Box_Step` in `src/Stress.f90` |
| Stress tensor (virial), with each cell's own β_c | always written | `Calculate_StressTensor` in `src/Stress.f90` |
| Excluded volume: no overlapping cells in the open tissue | `if_excluded_volume` | `Excluded_Volume_Force` in `src/Force.f90`, `src/T1_swap.f90` |
| T1 energy barrier: one edge shrunk quasi-statically through a T1 and back (after Bi et al., *Nat. Phys.* 2015), optionally with every cell's β_c frozen | `if_shrink_edge`, `shrink_freeze_beta` | `src/ShrinkEdge.f90` |
| Tissue energy (with per-cell β_c) | always written | `src/System_Info.f90` |

The code also contains features that are not used in the paper (motility gradients and hotspots, cell division, tissue squeeze/perturbation and limb forces, differential free-edge tension, myosin-ODE noise). They are off by default.

## Implementation notes

The exact discrete forms used in the code (dt is the time step, η the vertex friction, u a uniform random number on (0, 1), drawn independently for every vertex, cell and component as stated):

- **Units.** β is divided by λA₀ and γ by λA₀^{3/2} once at start-up while λ(A − A₀)² is kept as given, so the parameters are in reduced units when λ = A₀ = 1 (as in the paper).
- **Time stepping.** Explicit Euler–Maruyama: deterministic forces and the PPM self-propulsion enter as dt·f/η, the AM and PM random forces as √dt·f/η.
- **AM.** Each vertex, each component, every step: f = √(2 M_AM η) (2u − 1). The variance per component is (2/3) M_AM η.
- **PM.** Each cell, every step: p_c with components M_PM (u − ½) (variance M_PM²/12), white in time. Every vertex receives the sum of p_c over the cells it belongs to.
- **PPM.** θ_c ← θ_c + √(2 D_r dt) n, n Gaussian with unit variance; persistence time τ_p = 1/D_r. Each vertex moves with the average of v₀(cos θ_c, sin θ_c) over its cells (`abp_force_mode = 0`) or with the sum (`abp_force_mode = 1`).
- **FC.** β_c λA₀ (stored as `Myosin`) ← β_c λA₀ − dt (β_c λA₀ − β₀)/τ + √dt σ_β (u − ½); stationary standard deviation σ_β √(τ/24). With `fc_clip_negative`, a negative value is set to 0 after the step.
- **MR.** The Rho, ROCK and myosin equations are integrated per cell (Euler or RK4) before the forces of each step. With `mr_noise_mode = 0`, β_c = (α M_c + σ_MR (u − ½))/(λA₀) with a fresh u every step and no dt scaling, so the effect of the coupling noise depends on dt. `mr_noise_mode = 1` replaces it by a per-cell Ornstein–Uhlenbeck process.
- **T1.** An edge shorter than `min_d_T1` is replaced by an edge of length 2·`min_d_T1` on its perpendicular bisector. With `if_excluded_volume`, a flip whose new edge would cross another edge is retried with the new edge half and a quarter as long, and skipped if all three cross.
- **T2.** A triangular cell with area below `min_area_T2` is removed and replaced by a single vertex.
- **Stress.** Per cell, σ_c = 2λ(A_c − A₀) I + (2β_c P_c + γ)/(2A_c) Σ_{e∈c} l_e ⊗ l_e / l_e, averaged over cells weighted by area. With only the cell's own tension in the second term, its contribution (and therefore the shear stress) is half the virial stress of the tissue; quantities normalised by their initial value are unaffected.
- **Variable box.** The box lengths relax at the friction the vertices feel, driven by the energy derivative with respect to the box size plus the active virial. With `box_aniso` and a nonzero Lees–Edwards strain, the box relaxes to zero normal stresses in the lab frame at fixed strain (the cross terms of order strain × shear stress are not included in the driving force).
- **Output names.** Snapshot file names carry the step number in 8 digits; split runs longer than 10⁸ steps into restarts (`nrun = 2`).

## Building

Requires `gfortran`.

```
make                    # bin/vertexmain (the simulation)
make initial_config     # bin/generate_initial_mesh (initial mesh)
make fbounds            # bin/vertexmain with array-bounds checking (debugging)
```

## Running

Both programs read and write through paths relative to the current directory, so run them from the repository root.

1. Initial mesh: set the lattice size and boundary type in `para_MeshGen.dat` (`if_periodic` must match `if_PBC` in `para_Simulation.dat`), then `./bin/generate_initial_mesh` (or `./run_build_initial.sh`). It writes `mesh/`. A 64 × 64 periodic mesh is included.
2. Simulation: set `para_Simulation.dat`, create `data/` (`mkdir -p data`), then `./bin/vertexmain` (or `./run.sh`, which runs it in the background).

`para_Simulation.dat` is read line by line in a fixed order (documented in `PARAMETERS.md`); only the value at the start of each line matters. The lines after `ev_k` are optional and keep their defaults when missing.

`nrun = 1` starts from the mesh in `mesh/`. `nrun = 2` restarts from snapshot `nrun2_initialTime` of an earlier run in the same `data/` directory and writes its output with the prefix `nrun2_`. In the paper, every run first relaxes the passive tissue (`nrun = 1`, no activity) and then switches the activity on with a restart (`nrun = 2`).

`rng_seed > 0` makes a run reproducible.

## Output (`data/`)

Every `it_dump` steps, a snapshot of the tissue: vertex positions (`v_<it>.dat`), cell–vertex connectivity (`inn_<it>.dat`, `num_<it>.dat`), per-cell Rho/ROCK/Myosin (`Myosin_<it>.dat`; for FC, Myosin = β_c λA₀), and class-specific state (e.g. ABP orientations). Time series: energy, shear stress and strain, T1 and T2 counts and events, and the box size for the variable box. Files are unformatted Fortran records (see `write_output` in `src/allocation.f90`).

## Source files

| File | Contents |
|---|---|
| `src/vertexmain.f90` | Main program: input, time loop |
| `src/allocation.f90` | Parameters and arrays, `read_input`, output |
| `src/Force.f90` | Mechanical and active forces, FC and MR dynamics, excluded volume |
| `src/Geometry.f90` | Areas, perimeters, centroids, boundary detection, periodic images |
| `src/T1_swap.f90`, `src/T2_swap.f90` | Topological transitions |
| `src/Stress.f90` | Stress tensor, shear protocols, Lees–Edwards, variable box |
| `src/System_Info.f90` | Energy, cell centres |
| `src/ShrinkEdge.f90` | T1 energy-barrier protocol |
| `src/Proliferation.f90` | Cell division (not used in the paper) |
| `src/array_info.f90` | Small array utilities |
| `src/Generate_Initial_Mesh.f90` | Initial mesh: jittered hexagonal lattice, open or periodic |

## Licence and citation

GPL-3.0 (see `LICENSE`). If you use this code, please cite the paper above (see also `CITATION.cff`).
