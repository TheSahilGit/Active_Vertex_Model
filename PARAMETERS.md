# Parameter reference

`para_MeshGen.dat` and `para_Simulation.dat` are both read **positionally** (plain
sequential `read()` calls — `Generate_Initial_Mesh.f90` and
`allocation.f90::read_input` respectively), so their inline comments are kept short
— just enough to identify each line, not to explain it. This file documents every
parameter in both, in the exact order they're read, with a little more explanation
where the one-line comment isn't the whole story. If a parameter's inline comment
already says everything there is to say, its entry here is just as short.

If you add a parameter complex enough to need more than a few words, give it the
same treatment: a short pointer in the `.dat` file (`-- see PARAMETERS.md`), the
real explanation here — and remember both files must stay positionally in sync
with their readers (`Generate_Initial_Mesh.f90`, `allocation.f90::read_input`).

## `para_MeshGen.dat` — read by `Generate_Initial_Mesh.f90`

| Parameter | Description |
|---|---|
| `Lx` | Lattice columns. Also doubles as the exact physical box width once the mesh is built (lattice spacing is 1) — used that way everywhere else in the codebase (`mot_Lc*Ly`, the stress-tensor core radius, `Lx_box`/`Ly_box` under PBC). |
| `Ly` | Lattice rows. Must be even (the hexagonal-offset lattice alternates row offset by parity). Same dual role as `Lx`. |
| `seed` | RNG seed for the lattice jitter. Fixed value → the exact same mesh every time; change it for a different random mesh. |
| `jitter_interior_x` | Max random x-offset applied to interior-cell seed points, as a fraction of lattice spacing. |
| `jitter_interior_y` | Same, y-offset. One random draw per point scales *both* its x and y jitter (correlated direction, not independent noise) — see the module header comment. |
| `jitter_boundary_x` | Max random x-offset for outer-ring seed points specifically. Normally smaller than the interior jitter to keep a free mesh's edge less ragged. Ignored (every cell uses the interior jitter instead) when `if_periodic=.true.` — see that entry. |
| `jitter_boundary_y` | Same, y-offset for the outer ring. |
| `cell_headroom` | Extra cell slots reserved beyond `Lx*Ly`, for cells created later by division. Sets `num_dim = Lx*Ly + cell_headroom` in the written `mesh/para_MeshDims.dat`. |
| `vertex_slot_headroom` | Extra per-cell vertex-list slots reserved beyond the largest polygon actually generated. Sets `inn_dim1`. |
| `vertex_pool_headroom` | Extra vertex-array rows reserved beyond the vertex count actually generated (new vertices come from T1 flips and division). Sets `v_dim2`; rule of thumb ~2×N. |
| `if_periodic` | `.false.` (default) — the ordinary free-boundary mesh every earlier feature assumes: ghost seed copies are used only to clip boundary-cell shapes against something, not to connect them to anything (see the module header — "a boundary-conditioning trick, not true periodic BCs"). `.true.` — a genuinely wrap-connected torus mesh: a post-pass identifies vertices independently computed by cells on opposite edges of the box (their raw positions differ by exactly one lattice period) and unifies their IDs via union-find, so opposite edges share real topology. Verified via Euler characteristic (V−E+F=0, the torus signature), uniform vertex degree 3, and exact area conservation. Required (and checked at simulation startup, `allocation.f90::read_data`) if you're going to run with `if_PBC=.true.` in `para_Simulation.dat`. |

## `para_Simulation.dat` — read by `allocation.f90::read_input`

### General / restart

| Parameter | Description |
|---|---|
| `nrun` | `1` = fresh start from `mesh/v_in.dat`/`mesh/inn_in.dat`/`mesh/num_in.dat`. `2` = restart from a mid-run snapshot at `nrun2_initialTime`; all restart output is written under an `nrun2_` filename prefix so it never overwrites the original run's data. |
| `nrun2_initialTime` | The snapshot `it` to resume from when `nrun=2`. Ignored when `nrun=1`. With `if_RhoROCK`, Rho/ROCK/Myosin are also restored from `data/Myosin_<it>.dat`. If that record is all zero (RhoROCK was off in the source run), they are initialised at random as in a fresh run. |
| `Ao` | Target (preferred) cell area in the area-elasticity energy term. |
| `Co` | Target (preferred) cell perimeter in the perimeter-elasticity energy term. |
| `lambda` | Area-elasticity modulus. |
| `beta` | Perimeter-elasticity modulus. Non-dimensionalized once at startup (divided by `lambda*Ao`) — the value here is the physical one, not the rescaled internal one. Silently overridden every step if `if_active_contractility` or `if_RhoROCK` is on (the two aren't mutually exclusive in code; if both are on, `if_active_contractility` wins). |
| `gamm` | Line-tension modulus (linear-in-perimeter energy term). Non-dimensionalized once at startup like `beta`. |
| `eta` | Damping/mobility coefficient — the friction constant relating force to velocity (`v += dt*F/eta`). Read into `eta_0`, the uniform baseline every vertex's `eta(:)` (a per-vertex array, `allocation.f90`) is initialized to; currently still spatially uniform in value (no flag varies it yet), but the array exists so a future flag/ODE can. |
| `totT` | Total number of timesteps. |
| `dt` | Timestep size. |

### T1 / T2 topological transitions

| Parameter | Description |
|---|---|
| `if_Do_T1` | Enable T1 neighbor-exchange transitions (short-edge flips). |
| `min_d_T1` | Edge-length threshold below which an edge is a T1 candidate. |
| `if_Do_T2` | Enable T2 cell-extrusion transitions (triangular-cell removal). |
| `min_area_T2` | Area threshold below which a 3-vertex (triangular) cell is extruded. Only 3-vertex cells can T2 — no support for extruding higher-order polygons. |

### Boundary conditions

| Parameter | Description |
|---|---|
| `if_Fixed_boundary` | Pin (zero the force on) every vertex on the tissue's free outer edge. |
| `if_bottom_borders_fixed` | Pin only the bottom-most border vertices. Independent of `if_Fixed_boundary` — can be used alone for a "one wall" setup. |
| `if_top_borders_fixed` | Pin only the top-most border vertices. |
| `if_PBC` | Periodic boundary conditions — periodic in both x and y (a torus). See `Generate_Initial_Mesh.f90`'s `if_periodic` above for the matching mesh requirement, and the note there on the safety check at startup. Mutually exclusive with all three flags above (no free edge exists to pin under full periodicity), with `if_squeeze_tissue`/`if_limb_force` (both assume a real edge/corners), with `if_motility_gradient`+`if_motility_Eulerian` together (that gradient formula isn't periodic in y), and with `if_free_edge_adhesion` (no free edge exists under full periodicity, so `gamm_free` would never be used). All checked at startup (`read_input` stops with an error if violated). |

### Output cadence

| Parameter | Description |
|---|---|
| `it_dump` | Write a full mesh snapshot (`v`/`inn`/`num`/`force`/`Myosin`/`cell_identity`) every this many steps (always also at `it=1,2`). Whole-run summary series (`Energy`/`ShearStress`/`T1_count`/`T2_count`) are rewritten periodically at a coarser, automatically-derived cadence (`summary_dump_interval`, capped to ~200 rewrites over the whole run) — no separate parameter needed. |
| `T1_time_interval` | How often (in units of `it`) the T1 attempt loop runs. |
| `T2_time_interval` | Declared but dead — `Do_T2` is always called from inside the `T1_time_interval` loop regardless of this value. |

### Motility

| Parameter | Description |
|---|---|
| `if_motility` | Enable apolar random motility (an independent random force per vertex each step, magnitude set by the per-vertex field `mot`). |
| `if_motility_gradient` | `.false.` — `mot` is uniform, equal to `etas_max` everywhere. `.true.` — `mot` varies with y (see `etas_max`/`mot_Lc`); *when* it's evaluated is controlled by `if_motility_Eulerian`. |
| `if_motility_Eulerian` | Only matters when `if_motility_gradient=.true.`. `.false.` (Lagrangian, default) — `mot` is assigned once at `it=1` from each vertex's position then, and carried by that vertex wherever it moves afterward (a persistent tag, not a standing field — measurably erodes over a long run, mainly from cell division creating new vertices whose motility must be estimated, not from passive drift). `.true.` (Eulerian) — re-derived from each vertex's *current* position every step, so the gradient stays spatially anchored regardless of tissue motion/mixing/division (~6% overhead). Not periodic in y — disallowed together with `if_PBC`. |
| `etas_max` | Apolar motility strength — the uniform value when `if_motility_gradient=.false.`, or the peak value (at y=0) when the gradient is on, or the peak value at a hotspot center when `if_motility_hotspot` is on. |
| `etas_min` | Read but never actually used in the gradient formula (only `etas_max`/`mot_Lc` matter) — vestigial. |
| `mot_Lc` | Motility-gradient decay length scale, as a fraction of `Ly`: `mot = etas_max*exp(-y/(mot_Lc*Ly))`. |
| `if_motility_decay` | Exponentially decay `mot` toward 0 over `motility_decay_timeScale`. Largely moot combined with `if_motility_Eulerian` — the next step's recompute discards the decayed value before it can accumulate. |
| `motility_decay_timeScale` | Decay time constant, in units of `dt`. |
| `if_motility_hotspot` | Enable localized Gaussian motility enhancement around 1, 2, or 4 fixed hotspot cells (hardcoded cell indices near mesh corners, not user-configurable positions), instead of the y-gradient. Reuses `etas_max` as the peak value. Assigned once at `it=1` (or at restart), like Lagrangian-mode gradient — not re-evaluated every step. |
| `number_of_hotspot` | Must be 1, 2, or 4 — anything else stops with an error. |
| `sigma_hotspot` | Gaussian width of each hotspot's motility enhancement. |

### Shear

| Parameter | Description |
|---|---|
| `if_Shear_tissue` | Master switch — must be on for either shear mode below (and for `Calculate_StressTensor` to run at all) to have any effect. With `if_PBC` the shear uses Lees–Edwards boundary conditions: the periodic image above the box is displaced in x by the applied strain times the box height (`LE_strain`, `LE_offset`; `Advance_LE_Strain` in `Stress.f90`), so the sheared tissue stays periodic. |
| `if_Sudden_Shearing` | One-time affine shear (`x += sudden_shearStrength*y`) applied at `it==sudden_shearWhen`. |
| `sudden_shearStrength` | Shear strain applied by the sudden-shear event. |
| `sudden_shearWhen` | Timestep the sudden shear is applied at. |
| `if_Oscillatory_Shearing` | Continuous oscillatory affine shear, applied every step once `it > Oscl_shearWhen`. |
| `Oscl_shearStrength` | Oscillatory shear strain amplitude. |
| `Oscl_shearWhen` | Timestep oscillatory shearing begins. |
| `Oscl_freq_wo` | Oscillatory shear angular frequency. |

### Tissue perturbation

| Parameter | Description |
|---|---|
| `if_perturb_tissue` | Master switch for the sinusoidal perturbation below. |
| `if_sin_perturb` | Apply a sinusoidal x-displacement, `dt*sin_perturb_strength*sin(2*pi*sin_perturb_waveNumber*y/Ly)`, to every vertex. **Note:** this currently does nothing at all unless `if_dirac_comb` is *also* `.true.` — the displacement is nested inside the comb gate (`Force.f90`), likely an unintended coupling between two conceptually independent flags. |
| `sin_perturb_when` | Timestep the perturbation begins (unit of `it`). |
| `sin_perturb_strength` | Perturbation displacement amplitude. |
| `sin_perturb_waveNumber` | Number of spatial oscillations across the tissue height. |
| `if_dirac_comb` | Gate the perturbation on/off in a square-wave pulse train (see `if_sin_perturb`'s note — currently required for any sinusoidal perturbation to apply at all, even a continuous one, via a long `comb_onPeriod`). |
| `comb_onPeriod` | Pulse "on" duration, unit of `it`. |
| `comb_offPeriod` | Pulse "off" duration, unit of `it`. |

### Limb force

| Parameter | Description |
|---|---|
| `if_limb_force` | Apply a fixed outward x-force at four hardcoded near-corner cell locations (pushes left-half cells left, right-half cells right) — assumes a bounded domain with real corners; incompatible with `if_PBC`. |
| `limb_force_strength` | Force magnitude for the above. |

### Cell division

| Parameter | Description |
|---|---|
| `if_cell_division` | Enable area-threshold-triggered division (split along the principal axis). Daughter cells correctly inherit split-vertex motility; `Rho`/`ROCK`/`Myosin`/`cell_identity` are **not** propagated to daughters — a known gap for any lineage-tracking or RhoROCK+division study. |
| `area_0` | Area threshold above which a cell is a division candidate. |

### Squeeze

| Parameter | Description |
|---|---|
| `if_squeeze_tissue` | One-time uniform compression of every vertex's y-coordinate toward the current bottom border, at `it==squeeze_when`. Assumes a real bottom edge to anchor to — incompatible with `if_PBC`. |
| `squeeze_when` | Timestep the squeeze is applied at. |
| `percent_squeeze` | Compression amount, percent. |

### Active contractility

| Parameter | Description |
|---|---|
| `if_active_contractility` | Ornstein-Uhlenbeck fluctuating contractility: a per-cell "Myosin" level relaxing to the basal `beta` with correlation time `tau_contr`, directly setting the perimeter-elasticity modulus `beta` per cell. Simpler and independent of the full RhoROCK system below (see `beta`'s entry for the interaction if both are on). |
| `tau_contr` | Correlation time for the OU process. |
| `active_contr_strength` | Noise strength for the OU process. |

### Active Brownian motility (ABP)

| Parameter | Description |
|---|---|
| `if_ABP` | Self-propelled cells: each **cell** has one heading θ that diffuses rotationally (dθ = √(2·Dr)·dW, one noise per cell). Each vertex moves with the **average** of `vo*(cos θ, sin θ)` over the cells it belongs to (Barton et al. 2017), so a cell whose neighbours share its heading moves at `vo`. No neighbour alignment — headings diffuse independently per cell. Headings move with the cells at T2, are inherited at division, and are written to `data/ABP_theta_<it>.dat`, which an `nrun=2` restart reads back. Péclet number Pe = `vo`/(`Dr`·√A₀). |
| `vo` | Self-propulsion speed. |
| `Dr` | Rotational diffusion constant for the orientation angle. |

### RhoROCK mechanochemical regulation

| Parameter | Description |
|---|---|
| `if_RhoROCK` | Per-cell Rho→ROCK→Myosin reaction cascade (area-dilation-driven Rho activation via a Hill function, linear production/decay down the chain), feeding back into the perimeter-elasticity modulus `beta` via `Myosin_Coupling_Strength`. **Despite the naming, `D_rho`/`D_ROCK`/`D_Myosin` are pure local decay rates, not diffusion** — there is no neighbor-coupling term anywhere, so Rho/ROCK/Myosin never actually spread between cells. See `beta`'s entry for the interaction with `if_active_contractility` if both are on. |
| `if_RK4` | Integrate the RhoROCK ODEs with 4th-order Runge-Kutta instead of Euler. |
| `nhill` | Hill-function exponent for area-driven Rho activation. |
| `K_hill` | Hill-function half-saturation constant. |
| `A_Rho` | Rho production rate constant. |
| `A_ROCK` | ROCK production rate constant (from Rho). |
| `A_Myosin` | Myosin production rate constant (from ROCK). |
| `D_rho` | Rho decay rate (not diffusion — see `if_RhoROCK`). |
| `D_ROCK` | ROCK decay rate (not diffusion). |
| `D_Myosin` | Myosin decay rate (not diffusion). |
| `Myosin_Coupling_Strength` | Scales Myosin's feedback into the perimeter-elasticity modulus `beta`. |
| `if_myosin_noise` | Add noise directly to the Myosin ODE. |
| `if_gaussian_noise` | Use Gaussian (vs. uniform) noise when `if_myosin_noise` is on. |
| `Myosin_noise_strength` | Myosin noise magnitude. |
| `if_coupling_noise` | Add noise to the Myosin→`beta` coupling itself, instead of/in addition to the Myosin ODE. |
| `coupling_noise_strength` | Coupling-noise magnitude. |

### Polar motility

| Parameter | Description |
|---|---|
| `if_polar_motility` | A random 2D kick per cell per step, applied to all its vertices. Despite the name, this has **no persistence/direction memory** (redrawn independently every step, unlike `if_ABP`'s diffusing orientation) — it's cell-scale white noise, not a polarized/directed active force. |
| `polar_motility_strength` | Kick magnitude. |

### Differential (free-edge) line tension

| Parameter | Description |
|---|---|
| `if_free_edge_adhesion` | Gives free edges (touched by exactly one cell — only possible with `if_PBC=.false.`) their own line-tension value, `gamm_free`, instead of the usual `gamm`. Every other edge (shared by two cells) keeps `gamm`, unchanged. Which edges are free is never stored — recomputed fresh from `inn`/`num` every time it's needed (`Get_Free_Edges`, `Geometry.f90`), so T1/T2/division topology changes need no bookkeeping. Mutually exclusive with `if_PBC` (checked at startup): a periodic mesh has no free edges, so `gamm_free` would never be used. |
| `gamm_free` | Line-tension modulus for free edges only, used only when `if_free_edge_adhesion` is on. Non-dimensionalized once at startup exactly like `gamm` (divided by `lambda*Ao^1.5`) — the value here is the physical one, not the rescaled internal one. |

### Edge shrinking (T1 energy path)

T1 path of one edge after Bi, Lopez, Schwarz & Manning, *Nat. Phys.* 11, 1074 (2015) — implemented in `src/ShrinkEdge.f90` (module header has the full description). It runs inside the ordinary overdamped time loop, so every other flag behaves as usual. The edge's length is held fixed: after every step its two vertices are put back exactly that far apart, about their freely moving midpoint and along their freely rotating direction. The target length is stepped through stages of `shrink_stage_steps` steps each:
- hold at the starting length l₀
- `shrink_nstages` stages down to `shrink_l_min`
- the T1 of that edge, done through the ordinary `Do_T1`
- `shrink_nstages` + 1 stages back up to l₀
- release: the constraint is switched off for the rest of the run

With activity off, each stage gradient-descends to the constrained minimum E(l), which is Bi's protocol. While the protocol runs, the automatic T1 skips this one edge, because it is driven below `min_d_T1` on purpose. Every other edge T1s normally. If the edge stops separating its expected two cells for any other reason (a neighbouring T1 or T2 involving its vertices), the protocol aborts. Output goes to `data/ShrinkEdge.dat` (`nrun2_` prefix for `nrun=2`): one row per stage end, plus a full snapshot at the same `it`. Its last column, `dE_dl_mean`, is the constraint force dE/dl averaged over the stage's second half; integrating it over the edge length gives the work W done on the edge. W equals the energy change when the forces derive from a fixed energy, but not with RhoROCK, whose β changes along the path. `data/ShrinkEdge_cellE.dat` holds every cell's energy averaged over the same window, two rows per stage end (each cell's own β, then passive β for all cells). Sum any patch of cells from it when single snapshots are too noisy. Make `totT` a stage end (`1 + k*shrink_stage_steps`) so the run's final snapshot and summary files are written.

| Parameter | Description |
|---|---|
| `if_shrink_edge` | Run the protocol above. |
| `shrink_v1`, `shrink_v2` | The edge's two vertex IDs. `0 0` auto-picks the first interior edge of the cell nearest the tissue centre (the box centre under PBC). Must be an interior edge: shared by two cells, both endpoints 3-fold. |
| `shrink_nstages` | Number of stages from l₀ down to `shrink_l_min`, and again back up. |
| `shrink_l_min` | Edge length at which the T1 is done. Bi et al. swap at exactly 0, but a zero-length edge has no direction. Its orientation also relaxes at a rate ~ force/l, and explicit Euler oscillates below l ~ 1e-2 at dt = 1e-3. Extrapolating E linearly from each side to l = 0 recovers the 4-fold-vertex energy. |
| `shrink_stage_steps` | Steps the edge is held at each length. |
| `shrink_relax_steps` | Steps of free relaxation before the protocol (phase −1): the edge is unconstrained, β frozen if `shrink_freeze_beta`, so the shrink starts from the local minimum of the (frozen) landscape. l₀ = the edge length at its end. Add these steps to `totT`. 0 = start at once. |

### Excluded volume

Optional trailing lines: a para file that ends at `shrink_relax_steps` reads as `if_excluded_volume = .false.`.

| Parameter | Meaning |
|---|---|
| `if_excluded_volume` | Stops edges from crossing (overlapping cells). (1) Each vertex is pushed away from every edge that does not touch it, U = (`ev_k`/2)(`ev_d0` − d)² for distance d < `ev_d0`; for an edge that starts at one of the vertex's neighbours only when the closest point is inside the edge (a corner folding shut). (2) A T1 is kept only if no edge at its two vertices crosses another edge; otherwise it is retried with the new edge ½ and ¼ as long, and skipped for that step if all three cross. Open boundaries only (stops with `if_PBC`). The repulsion is not in `Energy.dat`. |
| `ev_d0` | Range of the vertex-edge repulsion (default 0.05, half the T1 threshold). |
| `ev_k` | Stiffness of the vertex-edge repulsion (default 100). |

### Optional switches ()

Optional trailing lines after `ev_k`, in this order; reading stops at the first missing line and the rest keep their defaults, which reproduce the old behaviour.

| Parameter | Meaning |
|---|---|
| `series_stride` | The per-step series (`Energy.dat`, `ShearStress.dat`, `T1_count.dat`, `T2_count.dat`) get one entry per block of this many steps: Energy and shear stress at the block's last step, T1/T2 counts summed over the block. Their arrays used to have length `totT`, so a long run held ~1.2 GB whatever the tissue size; with 1000 a `totT` = 5e7 run needs ~5 MB. The energy is then only computed at block ends and snapshot steps. 1 = every step (default). Must be 1 with `if_shrink_edge`. Trajectories do not depend on it. |
| `rng_seed` | > 0: seeds the random numbers (same spreading as the mesh generator), so a run is reproducible. 0: unseeded, as before (gfortran picks a random seed). The generator state right after seeding is written to `data/rng_seed.dat` (`nrun2_rng_seed.dat` on a restart) in both cases. |
| `mr_noise_mode` | MR coupling noise (`if_coupling_noise`). 0: as before, β_c = (α·Myosin + σ_MR·U(−½,½))/(λA₀) with a fresh draw every step and no √dt, so its effect depends on dt and β_c < 0 in about half the steps. 1: β_c = (α·Myosin + ξ_c)/(λA₀) with a per-cell Ornstein–Uhlenbeck ξ_c (Gaussian, Euler–Maruyama), dξ = −ξ/τ dt + s√(2/τ) dW, which does not depend on dt. |
| `mr_noise_tau` | Mode 1: correlation time τ of ξ_c (default 1). |
| `mr_noise_sd` | Mode 1: stationary s.d. s of ξ_c. < 0 (default): s = σ_MR·√(10⁻⁴/(24τ)), which gives the same low-frequency noise intensity (2s²τ = σ_MR²·dt/12) as mode 0 at dt = 10⁻⁴ (0.0020 for σ_MR = 1, τ = 1). |
| `abp_force_mode` | `if_ABP`: 0 (default) each vertex moves with the average of its cells' self-propulsion v₀(cos θ, sin θ); 1 each vertex gets the sum over its cells, the same distribution as the PM force. |
| `if_variable_box` | PBC only. Rescale the periodic box (and every vertex, and the Lees-Edwards offset) isotropically every step, driven by the virial of all vertex forces (elastic and active), so the tissue's total area is free instead of pinned at Lx·Ly. Box length written to `data/BoxLength.dat`; an nrun=2 restart restores the box (and the Lees-Edwards strain, from `data/ShearStrain.dat`) of its start snapshot. See `Variable_Box_Step` (Stress.f90). |
| `box_mobility` | 0 (default): the box moves under the virial of every force the vertices feel (elastic, FC/RhoROCK via β, AM noise, PM, ABP), with the friction the vertices feel when every cell dilates about its centre, Γ = Σ_c Σ_{v∈c} η_v \|r_v − R_c\|²/n_v — no free parameter. > 0: fixed mobility, d ln s = μ·(generalized force)·dt/Nc (testing). See `Variable_Box_Step` (Stress.f90). |
| `box_pressure` | Imposed external pressure (0 = free tissue: the box settles where the tissue's elastic pressure vanishes). |
| `box_aniso` | With `if_variable_box`: `.true.` lets Lx_box and Ly_box relax separately, each driven by its own normal stress (dE/dln s_x = Σ_c [2λ(A−A₀)A + T_c Σ_e l_x²/l] plus the active virial, friction Σ w η (x−X_c)²), so the box reaches zero stress in every direction (a square 16×16 box cannot hold regular hexagons). The Lees–Edwards strain is held fixed while the box deforms. `.false.` (default): isotropic scaling. See `Variable_Box_Step` (Stress.f90) |
| `mr_feedback_off` | With `if_RhoROCK`: `.true.` replaces α·Myosin in β_c by the constant `mr_beta_mean`, keeping the coupling noise of the chosen `mr_noise_mode`: β_c = (`mr_beta_mean` + noise)/(λA₀) (noise-only control). The Rho–ROCK–myosin equations still run and are written out but do not act on the mechanics. |
| `mr_beta_mean` | The constant ⟨β⟩ of `mr_feedback_off` (same units as α·Myosin). |
| `fc_clip_negative` | With `if_active_contractility`: `.true.` keeps every cell's contractility non-negative — whenever an Ornstein–Uhlenbeck step makes Myosin (= β_c λA₀) negative, it is set to 0. `.false.` (default): unbounded OU process. |

| `shrink_freeze_beta` | With `if_RhoROCK` or `if_active_contractility`: during the protocol (hold, shrink, expand) Myosin is not updated (no Rho/ROCK/Myosin ODEs, no fluctuating-contractility OU step) and the RhoROCK coupling noise is off, so each cell keeps its β from the start snapshot. This is the energy landscape of the tissue at that instant, where W equals the energy change again. FC myosin is clipped at 0 when the freeze starts: a frozen negative β has no energy minimum, and the cell runs away. Everything resumes at release. No effect for other classes. |

