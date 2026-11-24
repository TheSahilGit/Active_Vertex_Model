module allocation

   implicit none
!   save

      integer :: Nc
      integer :: nrun, nrun2_initialTime
      real*8 :: nrun2_initialTime_r
      integer*4 :: i,j,k,it,jp,jm,n,nn,totT
      integer*4 :: ii,jj, iki, ik, jk, kke
      integer :: ic, jc
      real*8 :: totTr
      real*8 :: min_d_T1, min_area_T2
      integer :: count_T1, cellNoT1, verNoT1, verNoNextT1
      integer :: count_T2, cellNoT2, verNoT2
      real*8, allocatable, dimension(:) :: Total_T1_count, Total_T2_count

      integer :: next_idx, prev_idx, chosen_index
      real*8 :: len_d,len_d_sq, dx, dy, rr
      integer*4, dimension(:), allocatable :: cell_no, ver_no, ver_no_next
      real*8, dimension(:), allocatable :: d_val
      real*8, dimension(:), allocatable :: area_val


      logical :: if_Do_T1, if_Do_T2

      real*8 :: len_d_T1
      integer*4 :: Lx,Ly

      integer*4, dimension(:), allocatable :: num
      real*8, dimension(:,:), allocatable:: v
      integer*4, dimension(:,:), allocatable :: inn
      integer*4, dimension(:), allocatable :: borderver
      integer*4 :: borderdim
      integer*4, allocatable, dimension(:) :: coordNum, bound
      integer :: boundary_count
   
   
      real*8 :: lambda, beta, gamm, Ao, Co,  dt
      real*8 :: beta_0

      ! eta: per-vertex Langevin/friction coefficient . eta_0 is
      ! the single scalar read from para_Simulation.dat -- the uniform
      ! baseline every vertex starts at, and the fixed reference value any
      ! future spatially-varying eta (under its own flag, likely an ODE
      ! paralleling Rho/ROCK/Myosin) should read/relax back toward,
      ! mirroring how beta_0 already serves that exact role for Myosin's
      ! mean-reversion (Force.f90) and initial condition (vertexmain.f90).
      ! eta(:) is what every force/noise calculation actually uses.
      real*8 :: eta_0

      real*8, allocatable, dimension(:) :: coefficients
      real*8, dimension(:),  allocatable ::  fxx_temp, fyy_temp
      real*8, dimension(:), allocatable:: fxx, fyy
      real*8, dimension(:), allocatable:: fxx_ran, fyy_ran
      integer*4, dimension(:), allocatable :: workingzone
      real*8, dimension(:), allocatable :: vx,vy

      real*8, dimension(:), allocatable :: edgelength, edgelengthIn
      real*8 :: area, perimeter

      real*8 :: totalarea
      real*8 ::  pi, rann

      integer*4 :: num_dim, inn_dim1,inn_dim2, v_dim1, v_dim2

      integer*4, dimension(:), allocatable :: mainarea,leftP,rightP,topP,bottomP,corners,boundary, boundary_occ
      integer, dimension(:), allocatable :: inside1, inside2
      integer*4, dimension(:), allocatable :: GleftP,GrightP,GtopP,GbottomP,Gcorners,Gboundary

      real*8, allocatable, dimension(:) :: Energy, ShearStress
      ! Applied (cumulative) shear strain at each step, written next to
      ! ShearStress (data/ShearStrain.dat) so stress-strain Lissajous-
      ! Bowditch curves can be built from the exact strain the run applied
      ! rather than re-deriving it from the protocol parameters.
      real*8, allocatable, dimension(:) :: ShearStrain
      ! Per-cell perimeter modulus actually used by Force_Calculation this
      ! step (stress_relaxation_LE copy): Calculate_StressTensor read the
      ! module scalar beta, which the force loop overwrites cell by cell, so
      ! with if_RhoROCK / if_active_contractility every cell's stress used
      ! the LAST cell's beta. Equal to beta everywhere otherwise.
      real*8, allocatable, dimension(:) :: beta_cell
      character(100) :: fname_inn, fname_num, fname_v, fname_force
      character(100) :: fname_Myosin, fname_cell_identity
      character(100) :: fname_inn2, fname_num2, fname_v2
      integer*4 :: iunit_inn, iunit_num, iunit_v
      ! T1/T2 event log (spatial/cell-identity tracking). Opened once (Open_Event_Logs) before the main loop,
      ! written to incrementally (one line per actual event, from Do_T1/
      ! Do_T2), closed once (Close_Event_Logs) after it.
      integer, parameter :: iunit_T1events = 941, iunit_T2events = 942
      integer, parameter :: iunit_Divisionevents = 943

      integer*4 :: it_dump, T1_time_interval, T2_time_interval
      integer*4 :: summary_dump_interval

      real*8, allocatable, dimension(:) ::  mot, mot0
      real*8, allocatable, dimension(:) :: eta
      real*8 :: etas_max, etas_min, mot_Lc
      real*8, allocatable, dimension(:) :: msdt, cellCentInit
      logical :: if_Shear_tissue
      real*8 :: sudden_shearStrength
      integer :: sudden_shearWhen, Oscl_shearWhen
      real*8 :: strainRate, Oscl_freq_wo, Oscl_shearStrength

      integer :: sin_perturb_when
      real*8 :: sin_perturb_strength, sin_perturb_waveNumber
      logical :: if_Perturb_tissue, if_sin_perturb, if_dirac_comb
      integer :: comb_onPeriod, comb_offPeriod


      real*8 :: sigma(2,2)
      real*8, allocatable, dimension(:,:) :: TotalSigma

      logical :: if_motility

      logical :: if_motility_gradient, if_Fixed_boundary, if_Sudden_Shearing, if_Oscillatory_Shearing

      ! Eulerian (spatial-field) vs Lagrangian (per-vertex tag) motility
      ! gradient, per user request : if_motility_gradient assigns
      ! mot once, at it=1, from each vertex's position at that instant
      ! (Give_Motility_Gradient, called only before the main loop) -- once
      ! assigned, a vertex carries that value wherever it physically moves
      ! for the rest of the run (Lagrangian). Under active motility forces,
      ! the highest-motility vertices (originally at the gradient's low end)
      ! also move the most, so the visible gradient erodes over the run
      ! (confirmed directly: corr(cell mean mot, cell's current y) dropped
      ! from -0.66 at it=20000 to -0.31 at it=200000 in a real run). Setting
      ! if_motility_Eulerian=.true. instead re-derives mot from each
      ! vertex's CURRENT position every step (see the Give_Motility_Gradient
      ! call in vertexmain.f90's main loop), so the gradient stays anchored
      ! to physical location -- a true spatial field -- instead of drifting
      ! with the tissue. Defaults .false. (Lagrangian) to keep every
      ! existing para_Simulation.dat's behavior unchanged; has no effect
      ! unless if_motility_gradient is also .true.
      logical :: if_motility_Eulerian

      logical :: if_motility_decay
      real*8 :: motility_decay_timeScale

      logical :: if_bottom_borders_fixed, if_top_borders_fixed

      ! True runtime periodic boundary conditions .
      ! Placed here, right after the other boundary-condition flags (its
      ! read(112,*) counterpart is inserted right after if_top_borders_
      ! fixed in read_input, mirrored at the same position in
      ! para_Simulation.dat) since if_PBC is fundamentally a fourth
      ! boundary-condition choice, mutually exclusive with the three
      ! above -- there is no free edge left for if_Fixed_boundary/
      ! if_bottom_borders_fixed/if_top_borders_fixed to pin once the
      ! tissue is periodic in both x and y (enforced in read_input).
      ! Requires a mesh built with if_periodic=.true. in
      ! para_MeshGen.dat/Generate_Initial_Mesh.f90 (a genuinely
      ! wrap-connected torus, not the default free-boundary mesh) --
      ! nothing here checks that at runtime; using if_PBC with an
      ! ordinary free-boundary mesh would silently treat its real free
      ! edge as if it wrapped, which is wrong.
      logical :: if_PBC
      ! Exact periodic box lengths, set once from Lx,Ly (cast to real*8)
      ! at the end of read_data, once Lx,Ly are known -- deliberately a
      ! separate name from the integer lattice-count Lx,Ly (used all over
      ! the rest of the code for index-range/array-sizing logic) to avoid
      ! any accidental mixing of "cell count" and "physical length"
      ! meanings. The mesh generator's ghost-copy tiling places periodic
      ! images at exactly +-Lx,+-Ly, so this is exact, not approximate.
      real*8 :: Lx_box, Ly_box
      ! Lees-Edwards shear : under if_PBC + if_Shear_tissue, the
      ! periodic image of the box directly ABOVE the primary box is
      ! displaced in x by LE_offset (and the one below by -LE_offset), so
      ! an image translation is n_y*(LE_offset, Ly_box) + n_x*(Lx_box, 0).
      ! LE_strain is the cumulative applied strain (LE_offset = LE_strain*
      ! Ly_box, kept reduced modulo Lx_box since only that matters). Both
      ! stay exactly 0 when not shearing, so Unwrap_Relative/Wrap_Position
      ! reduce to the plain orthogonal-PBC result.
      real*8 :: LE_offset = 0.0d0, LE_strain = 0.0d0
      ! Variable periodic box (if_variable_box): under if_PBC the box
      ! (Lx_box, Ly_box, and LE_offset with it) is rescaled isotropically every
      ! step by an overdamped barostat, so the tissue's total area is free to
      ! change instead of being pinned at Lx*Ly -- see Variable_Box_Step
      ! (Stress.f90). box_dEdlns = dE/d(ln s) of the last step (s = isotropic
      ! scale factor); the tissue pressure is -box_dEdlns/(2*Lx_box*Ly_box).
      ! BoxLength(it_s) = Lx_box at the end of series block it_s (data/BoxLength.dat; Ly_box =
      ! BoxLength*Ly/Lx).
      ! Optional para lines (after abp_force_mode); defaults: box fixed.
      logical :: if_variable_box = .false.
      real*8 :: box_mobility = 0.0d0, box_pressure = 0.0d0
      ! box_aniso : Lx_box and Ly_box relax separately, each
      ! under its own normal stress; off = isotropic scaling .
      logical :: box_aniso = .false.
      ! mr_feedback_off, mr_beta_mean : MR noise-only control,
      ! beta_c = (mr_beta_mean + coupling noise)/(lambda Ao); the Rho-ROCK-
      ! myosin equations still run but no longer set beta_c.
      logical :: mr_feedback_off = .false.
      real*8 :: mr_beta_mean = 0.0d0
      ! fc_clip_negative: FC (if_active_contractility) contractility kept non-
      ! negative -- after each Ornstein-Uhlenbeck step a negative Myosin(ic)
      ! (= beta_c lambda Ao) is set to 0. Optional last para line; default off.
      logical :: fc_clip_negative = .false.
      real*8 :: box_dEdlns = 0.0d0
      real*8, allocatable, dimension(:) :: BoxLength, BoxLengthY

      real*8, allocatable, dimension(:,:) :: cellcen
      real*8 :: global_cellCenX, global_cellCenY

      integer, allocatable, dimension(:)  :: bottom_border, top_border, left_border, right_border, all_borders
      integer :: bottom_border_count, top_border_count, left_border_count, right_border_count
      integer :: all_border_count
      
      integer, allocatable, dimension(:) :: vertex_occurance_count

      logical :: if_limb_force
      real*8 :: limb_force_strength

      logical :: if_motility_hotspot
      integer :: number_of_hotspot
      real*8 :: sigma_hotspot

      logical :: if_cell_division
      real*8 :: area_0

      integer, allocatable, dimension(:) :: chosen_cell
      integer :: chosen_cell_count

      logical :: if_squeeze_tissue
      integer :: squeeze_when
      real*8 :: percent_squeeze


      logical :: if_active_contractility
      real*8 :: tau_contr, active_contr_strength


      logical :: if_ABP
      real*8 :: vo, Dr
      ! ABP : fxx_ABP/fyy_ABP are per vertex; theta_ABP/rot_noise
      ! are per CELL (each cell is one self-propelled particle).
      ! ABP_theta_restored: theta_ABP read back on an nrun=2 restart.
      real*8, dimension(:), allocatable :: fxx_ABP, fyy_ABP, theta_ABP, rot_noise
      logical :: ABP_theta_restored = .false.


      logical :: if_RhoROCK, if_RK4, if_myosin_noise, if_gaussian_noise
      real*8, dimension(:), allocatable :: Rho, ROCK, Myosin
      ! .true. once read_data has restored Rho/ROCK/Myosin from the nrun=2
      ! start snapshot (vertexmain.f90 then skips Initialize_RhoROCK, and
      ! FC's Myosin = beta_0 initialisation).
      logical :: myosin_restored = .false.
      real*8 :: A_Rho, nhill, K_hill, A_ROCK, A_Myosin, D_rho, D_ROCK, D_Myosin
      real*8 :: Myosin_Coupling_Strength, Myosin_noise_strength
      logical :: if_coupling_noise
      real*8 :: coupling_noise_strength

      character(500), dimension(:), allocatable :: cell_identity
      ! Division inheritance : monotonically-increasing counter for
      ! the NEXT never-before-used cell_identity number a new daughter cell
      ! (Proliferation.f90) gets -- must not just reuse Nc+1 (the new live
      ! array index), since a T2 extrusion can shrink Nc first, and reusing
      ! a small index would collide with an already-extruded cell's old
      ! identity. Initialized once in allocate_arrays.
      integer :: next_new_cell_id

      logical :: if_polar_motility
      real*8 :: polar_motility_strength
      real*8, dimension(:), allocatable :: fxx_Polar, fyy_Polar

      ! Differential line tension : gamm stays the "bulk"/shared-edge
      ! value; gamm_free is a second line-tension value used only for FREE
      ! edges (touched by exactly one cell -- only possible with if_PBC=
      ! .false.). Which edges are free is never stored -- Get_Free_Edges
      ! (Geometry.f90) recomputes it fresh from inn/num every time it's
      ! needed, so T1/T2/division topology changes need no bookkeeping here.
      logical :: if_free_edge_adhesion
      real*8 :: gamm_free

      ! Edge shrinking (ShrinkEdge.f90): T1 energy-landscape
      ! protocol after Bi et al., Nat. Phys. 11, 1074 (2015) -- one edge's
      ! length is held fixed and stepped down to shrink_l_min in
      ! shrink_nstages stages of shrink_stage_steps steps each, T1-swapped
      ! there, then stepped back up. shrink_v1/shrink_v2 = the edge's two
      ! vertex IDs (0 0 = auto-pick the interior edge nearest the tissue
      ! centre). shrink_va/shrink_vb/shrink_protect_edge/shrink_force_T1
      ! are the protocol's live state, kept here so find_T1 (T1_swap.f90)
      ! can see them.
      logical :: if_shrink_edge
      integer :: shrink_v1, shrink_v2
      integer :: shrink_nstages, shrink_stage_steps
      real*8 :: shrink_l_min
      ! shrink_freeze_beta: hold each cell's beta during the protocol (no
      ! RhoROCK ODEs or coupling noise, no FC OU update) -- the energy
      ! landscape of the tissue as it is at that instant. shrink_chem_frozen
      ! is the live switch (on in phases 0-2, off from release/abort).
      logical :: shrink_freeze_beta, shrink_chem_frozen = .false.
      ! shrink_relax_steps: steps of free relaxation (edge unconstrained,
      ! beta frozen if shrink_freeze_beta) before the protocol starts;
      ! l0 is the edge length at its end. 0 = start at once.
      integer :: shrink_relax_steps
      integer :: shrink_va = 0, shrink_vb = 0
      logical :: shrink_protect_edge = .false., shrink_force_T1 = .false.

      ! Excluded volume (branch no_overlap): short-range repulsion
      ! between every vertex and every edge not touching it or its neighbours,
      ! (ev_k/2)(ev_d0 - d)^2 for distance d < ev_d0, so edges cannot pass
      ! through each other (free boundary folding onto the tissue, cells
      ! overlapping). Off by default; older para files without these lines
      ! read as off.
      logical :: if_excluded_volume = .false.
      real*8 :: ev_d0 = 0.05d0, ev_k = 100.0d0

      ! Revision switches (; optional trailing para lines after
      ! ev_k, older files read as the defaults = old behaviour):
      ! series_stride: the per-step series (Energy, ShearStress, T1/T2 counts)
      !   are stored once per block of series_stride steps (Energy/ShearStress:
      !   the value at the block's last step; T1/T2 counts: summed over the
      !   block), entry it_s = (it-1)/series_stride + 1. They used to be
      !   allocated with length totT, so memory grew with run length (~1.2 GB
      !   at totT = 5e7). 1 = every step, as before.
      ! rng_seed: > 0 seeds the random numbers (random_seed put, spread as in
      !   Generate_Initial_Mesh); 0 = unseeded as before. The seed array used
      !   is written to data/[nrun2_]rng_seed.dat in both cases.
      ! mr_noise_mode: MR coupling noise. 0 = as before: a fresh uniform draw
      !   sigma_MR U(-1/2,1/2) added to beta every step, NOT scaled with dt, so
      !   its effect depends on dt and beta_c < 0 in ~half the steps.
      !   1 = per-cell Ornstein-Uhlenbeck noise xi_c with correlation time
      !   mr_noise_tau and stationary s.d. mr_noise_sd (Gaussian, Euler-
      !   Maruyama), dt-independent; mr_noise_sd < 0 picks the s.d. whose
      !   low-frequency intensity 2 sd^2 tau equals that of mode 0 at
      !   dt = 1e-4 (sigma_MR^2 1e-4/12): sd = sigma_MR sqrt(1e-4/(24 tau)).
      ! abp_force_mode: 0 = each vertex moves with the AVERAGE of its cells'
      !   self-propulsion (as before); 1 = each vertex gets the SUM over its
      !   cells, i.e. distributed exactly as the PM force.
      integer :: series_stride = 1, rng_seed = 0, mr_noise_mode = 0, abp_force_mode = 0
      integer :: it_s = 1
      real*8 :: mr_noise_tau = 1.0d0, mr_noise_sd = -1.0d0
      real*8, allocatable :: mr_xi(:)

      integer :: n_inside, n_outside
      real*8 :: radius_from_core
      real*8, allocatable, dimension(:,:) :: cell_centers
      integer, allocatable, dimension(:) :: inside_cells, outside_cells
      
        



   contains 

   subroutine read_input

     implicit none
     integer :: ios_ev
  
!     open(unit=121, file='para.in', status='old'); 
     open(112, file='para_Simulation.dat', status='old')
     open(unit=121, file='mesh/para_MeshDims.dat', status='old');
  
     read(121,*) Lx
     read(121,*) Ly
     read(121,*) num_dim
     read(121,*) v_dim1
     read(121,*) v_dim2
     read(121,*) inn_dim1
     read(121,*) inn_dim2

     ! Exact periodic box lengths -- see the Lx_box/Ly_box declaration
     ! comment. Computed here (Lx,Ly are already known) rather than at the
     ! end of read_input, so it's available regardless of where if_PBC
     ! itself is later validated.
     Lx_box = dble(Lx)
     Ly_box = dble(Ly)

     read(112,*) nrun
     read(112,*) nrun2_initialTime_r
     read(112,*) Ao
     read(112,*) Co
     read(112,*) lambda
     read(112,*) beta
     read(112,*) gamm
     read(112,*) eta_0
     read(112,*) totTr
     read(112,*) dt
     read(112,*) if_Do_T1
     read(112,*) min_d_T1
     read(112,*) if_Do_T2
     read(112,*) min_area_T2
     read(112,*) if_Fixed_boundary
     read(112,*) if_bottom_borders_fixed
     read(112,*) if_top_borders_fixed
     read(112,*) if_PBC

     ! Validation : if_PBC means the tissue is periodic in both x
     ! and y -- there is no free edge left, so a flag that means "pin/
     ! detect a free edge" is meaningless (and, worse, would silently
     ! mispin/misdetect on a periodic mesh's genuine interior vertices).
     ! Stop rather than run with an ambiguous/contradictory configuration.
     if (if_PBC) then
       if (if_Fixed_boundary .or. if_bottom_borders_fixed .or. if_top_borders_fixed) then
         write(*,*) 'read_input: if_PBC is incompatible with if_Fixed_boundary/', &
           'if_bottom_borders_fixed/if_top_borders_fixed (no free edge exists under PBC).'
         stop 1
       end if
     end if

     read(112,*) it_dump
     read(112,*) T1_time_interval
     read(112,*) T2_time_interval
     read(112,*) if_motility
     read(112,*) if_motility_gradient
     read(112,*) if_motility_Eulerian
     read(112,*) etas_max
     read(112,*) etas_min
     read(112,*) mot_Lc
     read(112,*) if_motility_decay
     read(112,*) motility_decay_timeScale
     read(112,*) if_motility_hotspot
     read(112,*) number_of_hotspot
     read(112,*) sigma_hotspot
     read(112,*) if_Shear_tissue
     read(112,*) if_Sudden_Shearing
     read(112,*) sudden_shearStrength
     read(112,*) sudden_shearWhen
     read(112,*) if_Oscillatory_Shearing
     read(112,*) Oscl_shearStrength
     read(112,*) Oscl_shearWhen
     read(112,*) Oscl_freq_wo
     read(112,*) if_Perturb_tissue
     read(112,*) if_sin_perturb
     read(112,*) sin_perturb_when
     read(112,*) sin_perturb_strength
     read(112,*) sin_perturb_waveNumber
     read(112,*) if_dirac_comb
     read(112,*) comb_onPeriod
     read(112,*) comb_offPeriod
     read(112,*) if_limb_force
     read(112,*) limb_force_strength
     read(112,*) if_cell_division
     read(112,*) area_0
     read(112,*) if_squeeze_tissue
     read(112,*) squeeze_when
     read(112,*) percent_squeeze
     read(112,*) if_active_contractility
     read(112,*) tau_contr
     read(112,*) active_contr_strength
     read(112,*) if_ABP
     read(112,*) vo
     read(112,*) Dr
     read(112,*) if_RhoROCK
     read(112,*) if_RK4
     read(112,*) nhill
     read(112,*) K_hill 
     read(112,*) A_Rho
     read(112,*) A_ROCK 
     read(112,*) A_Myosin
     read(112,*) D_rho
     read(112,*) D_ROCK
     read(112,*) D_Myosin
     read(112,*) Myosin_Coupling_Strength
     read(112,*) if_myosin_noise
     read(112,*) if_gaussian_noise
     read(112,*) Myosin_noise_strength
     read(112,*) if_coupling_noise
     read(112,*) coupling_noise_strength
     read(112,*) if_polar_motility
     read(112,*) polar_motility_strength
     read(112,*) if_free_edge_adhesion
     read(112,*) gamm_free
     read(112,*) if_shrink_edge
     read(112,*) shrink_v1
     read(112,*) shrink_v2
     read(112,*) shrink_nstages
     read(112,*) shrink_l_min
     read(112,*) shrink_stage_steps
     read(112,*) shrink_freeze_beta
     read(112,*) shrink_relax_steps
     ! Excluded volume (branch no_overlap): optional trailing lines -- a para
     ! file that ends at shrink_relax_steps (the comment block follows) leaves
     ! the feature off with its default range and stiffness.
     read(112,*,iostat=ios_ev) if_excluded_volume
     if (ios_ev /= 0) then
       if_excluded_volume = .false.
     else
       read(112,*) ev_d0
       read(112,*) ev_k
     end if
     ! Revision switches : optional, in this order; reading
     ! stops at the first missing line and the rest keep their defaults.
     if (ios_ev == 0) read(112,*,iostat=ios_ev) series_stride
     if (ios_ev == 0) read(112,*,iostat=ios_ev) rng_seed
     if (ios_ev == 0) read(112,*,iostat=ios_ev) mr_noise_mode
     if (ios_ev == 0) read(112,*,iostat=ios_ev) mr_noise_tau
     if (ios_ev == 0) read(112,*,iostat=ios_ev) mr_noise_sd
     if (ios_ev == 0) read(112,*,iostat=ios_ev) abp_force_mode
     ! Variable periodic box (merged from branch pbc_variable_box):
     ! optional, after abp_force_mode.
     if (ios_ev == 0) read(112,*,iostat=ios_ev) if_variable_box
     if (ios_ev == 0) read(112,*,iostat=ios_ev) box_mobility
     if (ios_ev == 0) read(112,*,iostat=ios_ev) box_pressure
     if (ios_ev == 0) read(112,*,iostat=ios_ev) box_aniso
     if (ios_ev == 0) read(112,*,iostat=ios_ev) mr_feedback_off
     if (ios_ev == 0) read(112,*,iostat=ios_ev) mr_beta_mean
     if (ios_ev == 0) read(112,*,iostat=ios_ev) fc_clip_negative
     ! if_variable_box only makes sense for a periodic box (an open tissue's
     ! area is already free).
     if (if_variable_box .and. .not. if_PBC) then
       write(*,*) 'read_input: if_variable_box needs if_PBC (an open tissue has no box).'
       stop 1
     end if
     if (series_stride < 1) series_stride = 1
     if (if_shrink_edge .and. series_stride /= 1) then
       write(*,*) 'read_input: if_shrink_edge needs series_stride = 1 (ShrinkEdge.f90 reads Energy every step).'
       stop 1
     end if
     if (mr_noise_sd < 0.0d0) mr_noise_sd = coupling_noise_strength * sqrt(1.0d-4 / (24.0d0 * mr_noise_tau))





     beta_0 = beta

     ! BUGFIX : beta/gamm must be nondimensionalized exactly ONCE here.
     ! Previously this rescale was done inside Force_Calculation (Force.f90) and
     ! again inside Calculate_StressTensor (Stress.f90), both of which run every
     ! timestep and mutate the global beta/gamm in place -- causing beta/gamm to
     ! be divided by (lambda*Ao) again on every call (exponential decay/blow-up
     ! over the run). Doing it once here, right after read, fixes that.
     beta = beta/(lambda*Ao)
     gamm = gamm/(lambda*(Ao)**1.5)
     ! gamm_free plays the exact same physical role as gamm (a line-tension
     ! coefficient), just for a different edge population -- must be
     ! nondimensionalized identically or the two would be in inconsistent
     ! units whenever lambda*Ao^1.5 /= 1 (if_free_edge_adhesion).
     gamm_free = gamm_free/(lambda*(Ao)**1.5)

     totT = int(totTr)
     nrun2_initialTime = int(nrun2_initialTime_r)

     ! CHANGE (MATLAB live-plotting request): Energy/ShearStress/
     ! T1_count/T2_count are now rewritten periodically during the run (not
     ! just once at the end) so they can be plotted while a simulation is
     ! still in progress -- see write_output. Naively rewriting the whole
     ! growing (1:it) slice every it_dump would cost O(totT^2/it_dump) total
     ! disk I/O over the run: for this session's actual production config
     ! (totT=1e7, it_dump=100) that works out to ~4TB written for Energy.dat
     ! alone (~16TB across all four files) -- large enough to dominate the
     ! run's wall-clock time. Instead, rewrite them only every
     ! summary_dump_interval, a multiple of it_dump chosen so there are
     ! about 200 rewrites total over the whole run regardless of totT/
     ! it_dump: this bounds the cumulative I/O to a small, fixed multiple of
     ! the final file size (~100x), while still refreshing often enough to
     ! watch a long run's progress. Never coarser than it_dump itself (a run
     ! with very few total dumps just gets a summary write every dump, same
     ! as before).
     summary_dump_interval = max(it_dump, it_dump * nint(dble(totT) / dble(it_dump * 200)))

     pi = acos(-1.0d0)

     ! Further if_PBC validation  -- these flags are read later in
     ! this same subroutine than if_Fixed_boundary/if_bottom_borders_fixed/
     ! if_top_borders_fixed above, so they're checked here instead:
     !  - if_squeeze_tissue/if_limb_force both assume a bounded domain with
     !    a real bottom edge / real corners to act on -- meaningless under
     !    full periodicity (no edge, no corners).
     !  - if_motility_gradient's Eulerian mode (if_motility_Eulerian)
     !    recomputes mot = etas_max*exp(-y/(mot_Lc*Ly)) from each vertex's
     !    CURRENT y every step -- not periodic in y, so a vertex wrapping
     !    y=Ly->0 would see a discontinuous motility jump. Disallowed for
     !    now per explicit user decision; a periodic-safe gradient formula
     !    is a separate future feature, not attempted here.
     if (if_PBC) then
       if (if_squeeze_tissue) then
         write(*,*) 'read_input: if_PBC is incompatible with if_squeeze_tissue', &
           ' (assumes a real bottom edge to compress against).'
         stop 1
       end if
       if (if_limb_force) then
         write(*,*) 'read_input: if_PBC is incompatible with if_limb_force', &
           ' (assumes real tissue corners to push outward from).'
         stop 1
       end if
       if (if_motility_gradient .and. if_motility_Eulerian) then
         write(*,*) 'read_input: if_PBC with if_motility_gradient+if_motility_Eulerian', &
           ' is not supported yet (the gradient formula is not periodic in y).'
         stop 1
       end if
       ! if_PBC + if_Shear_tissue is Lees-Edwards shear : ShearTissue
       ! applies the affine increment to every vertex AND advances
       ! LE_offset, and Unwrap_Relative/Wrap_Position (Geometry.f90) use
       ! the sheared image lattice -- so the affine move is consistent
       ! across the y-wrap instead of tearing the tissue there.
       ! nrun=2 restarts restore LE_strain/LE_offset (and the box length)
       ! from the source run's data/ShearStrain.dat and data/BoxLength.dat --
       ! see Restore_Box_And_Strain (read_data).
       if (if_excluded_volume) then
         write(*,*) 'read_input: if_PBC with if_excluded_volume is not supported yet -- ', &
           'the vertex-edge distances are not minimum-image.'
         stop 1
       end if
       if (if_free_edge_adhesion) then
         write(*,*) 'read_input: if_PBC with if_free_edge_adhesion makes no sense -- ', &
           'a periodic mesh has no free edges (every edge is shared by exactly two ', &
           'cells), so gamm_free would never be used. Set if_PBC=.false. to use this.'
         stop 1
       end if
     end if

     close(121)
     close(112)

   end subroutine read_input

   subroutine allocate_arrays
     implicit none
     integer :: iic
     call read_input

     allocate(num(num_dim))
     allocate(v(v_dim1,v_dim2))
     allocate(inn(inn_dim1,inn_dim2))
     allocate(fxx(v_dim2), fyy(v_dim2))
     allocate(fxx_ran(v_dim2), fyy_ran(v_dim2))
     allocate(fxx_Polar(v_dim2), fyy_Polar(v_dim2))
     allocate(fxx_temp(v_dim2), fyy_temp(v_dim2))
     allocate(fxx_ABP(v_dim2), fyy_ABP(v_dim2), theta_ABP(num_dim), rot_noise(num_dim))
     allocate(edgelength(Lx*Ly*inn_dim1))
     allocate(edgelengthIn(Lx*Ly*inn_dim1))
     allocate(workingzone(Lx*Ly+2*Lx+2*Ly+4))
     allocate(mainarea(Lx*Ly),leftP(Ly),rightP(Ly),topP(Lx),bottomP(Lx),corners(4))
     ! BUGFIX : Get_Boundary_info actually builds
     ! boundary = [leftP, rightP, topP, bottomP, corners], i.e. 2*Lx+2*Ly+4
     ! elements (corners was missing from this initial sizing) -- harmless
     ! today only because "boundary" is allocatable and Fortran's automatic
     ! reallocation-on-assignment (gfortran default, no -fno-realloc-lhs in
     ! compile.sh) silently resizes it on that whole-array assignment.
     ! Sized correctly here so it isn't fragile if that assignment pattern
     ! ever changes to a partial/indexed one (as the border arrays already do).
     allocate(boundary(2*Lx + 2*Ly + 4))
     allocate(GleftP(Ly),GrightP(Ly),GtopP(Lx),GbottomP(Lx),Gcorners(4))
     allocate(Gboundary(2*Lx + 2*Ly + 4))
     allocate(inside1((Lx - 2) * (Ly - 2)))
     allocate(inside2((Lx - 4) * (Ly - 4)))
     allocate(coefficients(6))
     ! Per-step series, one entry per block of series_stride steps:
     ! length (totT-1)/series_stride + 1, = totT for stride 1.
     allocate(Energy((totT - 1) / series_stride + 1))
     allocate(msdt((totT - 1) / series_stride + 1))
     allocate(Total_T1_count((totT - 1) / series_stride + 1), Total_T2_count((totT - 1) / series_stride + 1))
     allocate(cellCentInit(Lx*Ly-4*Lx - 4*Ly +16))
     allocate(borderver(v_dim2)) 
     allocate(mot(v_dim2),mot0(v_dim2))
     ! Uniform baseline -- every vertex starts at the single eta_0 read
     ! from para_Simulation.dat. A future flag/ODE would evolve eta(:)
     ! away from this per-vertex, same lifecycle as mot's gradient/hotspot
     ! setup (vertexmain.f90) relative to mot0.
     allocate(eta(v_dim2))
     eta(:) = eta_0
     allocate(coordNum(Lx*Ly))
     allocate(bound(2*Lx + 2*Ly))
     allocate(cellcen(num_dim, 2))

     ! BUGFIX : d_val/ver_no/ver_no_next/area_val must all share one
     ! consistent worst-case capacity. Previously d_val was sized Lx*Ly*6 while
     ! ver_no/ver_no_next were sized only inn_dim2 (much smaller) -- if find_T1
     ! (T1_swap.f90) ever found more than inn_dim2 short edges in one pass,
     ! ver_no/ver_no_next overflowed. Also, sizing off the *initial* Lx*Ly (not
     ! num_dim, the reserved max-cell capacity) meant these silently became too
     ! small again once Nc grows via cell division. Now all sized from num_dim.
     ! BUGFIX : the earlier num_dim*6 factor assumed
     ! an "average hexagonal cell" (~6 vertices/cell) but is not a proven
     ! bound against the actually configured per-cell vertex cap (inn_dim1);
     ! the true worst case (every live cell at inn_dim1 vertices, every edge
     ! short) is num_dim*inn_dim1. Sized off that instead.
     allocate(d_val(num_dim*inn_dim1))
     allocate(area_val(num_dim))
     allocate(ver_no(num_dim*inn_dim1), ver_no_next(num_dim*inn_dim1))
     ! BUGFIX : left_border/right_border collect
     ! vertices along the VERTICAL (Ly) edges of the tissue, not the
     ! horizontal (Lx) ones -- they were sized Lx*6 like bottom/top_border,
     ! which only happened to be safe because the shipped mesh is square
     ! (Lx==Ly). Sized off Ly instead, matching the dimension they actually
     ! scale with (c.f. Get_Boundary_info's GleftP/GrightP, already Ly-sized).
     allocate(bottom_border(Lx*6), top_border(Lx*6), left_border(Ly*6), right_border(Ly*6))
     allocate(all_borders(Lx*Ly*6))
     allocate(TotalSigma(2,2))
     allocate(ShearStress((totT - 1) / series_stride + 1))
     ! BUGFIX : ShearStress is unconditionally written to
     ! data/ShearStress.dat in write_output regardless of if_Shear_tissue.
     ! Without shearing it was only ever allocated, never assigned, so the
     ! output file contained uninitialized memory. Zero it here.
     ShearStress = 0.0d0
     allocate(ShearStrain((totT - 1) / series_stride + 1))
     allocate(BoxLength((totT - 1) / series_stride + 1))
     allocate(BoxLengthY((totT - 1) / series_stride + 1))
     ! Like ShearStrain: entry it is set after step it's box move, so the
     ! last entry of the file (written at the top of step totT) stays 0.
     BoxLength = 0.0d0
     BoxLengthY = 0.0d0
     allocate(beta_cell(num_dim))
     beta_cell = beta
     ShearStrain = 0.0d0

     allocate(chosen_cell(v_dim2))

     allocate(Rho(num_dim), ROCK(num_dim), Myosin(num_dim))

     allocate(cell_identity(num_dim))

     allocate(cell_centers(num_dim, 2), inside_cells(num_dim), outside_cells(num_dim))


     fxx = 0.0d0; fyy = 0.0d0
     fxx_ran = 0.0d0; fyy_ran = 0.0d0
     fxx_ABP = 0.0d0; fyy_ABP = 0.0d0; theta_ABP = 0.0d0
     fxx_Polar = 0.0d0; fyy_Polar = 0.0d0
     rot_noise = 0.0d0

     Rho = 0.0d0; ROCK = 0.0d0; Myosin = 0.0d0

      
     cell_identity = 'cell_0'
     do iic = 1, Lx*Ly
      write(cell_identity(iic), '(A,I0)') 'cell_', iic
     end do
     next_new_cell_id = Lx*Ly + 1

     inside_cells = 0
     outside_cells = 0

    end subroutine allocate_arrays

    subroutine read_data

      if(nrun.eq.1)then
       open(1, file='mesh/num_in.dat',status='old')
       open(2, file='mesh/v_in.dat', status='old')
       open(3, file='mesh/inn_in.dat', status='old')

       
       
    
       read(1,*)num
       read(2,*)v
       read(3,*)inn

    
    
       close(1)
       close(2)
       close(3)

     elseif(nrun.eq.2)then


       
        iunit_inn = 533
        iunit_num = 623
        iunit_v = 959

        write(fname_inn2, '("data/inn_", I8.8,".dat")')(nrun2_initialTime)
        write(fname_v2, '("data/v_", I8.8,".dat")')(nrun2_initialTime)
        write(fname_num2, '("data/num_", I8.8,".dat")')(nrun2_initialTime)
        write(*,*)fname_inn2

        open(unit = iunit_inn, file=fname_inn2, form = 'unformatted', status='old')
        open(unit = iunit_v, file=fname_v2, form = 'unformatted', status='old')
        open(unit = iunit_num, file=fname_num2, form = 'unformatted', status='old')


        read(iunit_inn)inn
        read(iunit_v)v
        read(iunit_num)num
       
        close(iunit_inn)
        close(iunit_num)
        close(iunit_v)

        ! RhoROCK / FC restart : carry the chemical state over from
        ! the start snapshot (write_output stores Rho, ROCK, Myosin
        ! interleaved per cell) instead of re-initialising it. An all-zero
        ! record means the source run had neither on -- nothing to restore,
        ! so vertexmain.f90 initialises as before.
        if (if_RhoROCK .or. if_active_contractility) then
          write(fname_Myosin, '("data/Myosin_", I8.8,".dat")')(nrun2_initialTime)
          inquire(file=fname_Myosin, exist=myosin_restored)
          if (myosin_restored) then
            open(unit = iunit_v, file=fname_Myosin, form = 'unformatted', status='old')
            read(iunit_v)(Rho(i), ROCK(i), Myosin(i), i = 1, num_dim)
            close(iunit_v)
            myosin_restored = any(Rho /= 0.0d0) .or. any(ROCK /= 0.0d0) .or. any(Myosin /= 0.0d0)
            if (.not. myosin_restored) then
              Rho = 0.0d0; ROCK = 0.0d0; Myosin = 0.0d0
            end if
          end if
          if (myosin_restored) then
            write(*,*) 'Myosin state restored from ', trim(fname_Myosin)
          else
            write(*,*) 'Myosin: no nonzero state in ', trim(fname_Myosin), ' -- initialising as a fresh run'
          end if
        end if

        if (if_PBC) call Restore_Box_And_Strain

        ! ABP restart : carry each cell's heading over, so the
        ! persistence survives the restart. Missing file (ABP was off in
        ! the source run) -> vertexmain.f90 draws random headings.
        if (if_ABP) then
          write(fname_Myosin, '("data/ABP_theta_", I8.8,".dat")')(nrun2_initialTime)
          inquire(file=fname_Myosin, exist=ABP_theta_restored)
          if (ABP_theta_restored) then
            open(unit = iunit_v, file=fname_Myosin, form = 'unformatted', status='old')
            read(iunit_v)(theta_ABP(i), i = 1, num_dim)
            close(iunit_v)
            write(*,*) 'ABP headings restored from ', trim(fname_Myosin)
          end if
        end if





     end if

     coefficients = [lambda, beta, gamm, Ao, Co, eta_0]

     ! BUGFIX : cell_no and vertex_occurance_count
     ! were sized off maxval(inn) -- the largest vertex index present in the
     ! *initial* mesh only, since read_data runs once before the main loop.
     ! vertex_occurance_count is indexed DIRECTLY by vertex id
     ! (Find_boundary_dynamic, Geometry.f90), and cell_division legitimately
     ! creates vertices with indices beyond that initial snapshot (up to
     ! v_dim2, the real reserved capacity) -- confirmed live: this produced
     ! "Index '4568' ... above upper bound of 4564" once cell division plus
     ! a T2 event grew the mesh past the initial maxval(inn). cell_no shares
     ! the same fragile sizing convention (it's filled in lockstep with
     ! ver_no in find_T1/find_T2's candidate scan) even though its true
     ! required capacity is candidate-count-based, not vertex-id-based.
     ! Sized off v_dim2 (the true reserved vertex capacity) for both.
     allocate(cell_no(v_dim2))
     allocate(vertex_occurance_count(v_dim2))

     ! BUGFIX : Nc (live cell count) used to be hardcoded to Lx*Ly
     ! in vertexmain.f90 regardless of nrun. That's correct for nrun=1 (a
     ! freshly generated mesh always has exactly Lx*Ly live cells), but
     ! wrong for nrun=2: cell division (Proliferation.f90: Nc=Nc+1) and T2
     ! extrusion (T2_swap.f90: Nc=Nc-1) both change the true live-cell
     ! count over the course of a run, and every physics loop (Force.f90,
     ! T1_swap.f90, T2_swap.f90) is bounded by `do ic=1,Nc`. Restarting
     ! (nrun=2) from a snapshot taken *after* any division/T2 event had a
     ! true Nc different from Lx*Ly, yet Nc was silently reset to Lx*Ly on
     ! restart -- confirmed live: restarting from a snapshot with 608 live
     ! cells (576 initial + 32 divisions) came back reporting Nc=576, so
     ! every already-divided daughter cell would have been silently
     ! ignored by physics from that point on despite loading correctly
     ! into v/inn/num. Derived here instead, from the invariant (already
     ! relied on by T2_swap.f90 and every MATLAB analysis script) that
     ! cells 1:Nc are live and contiguous and num(k)=0 for every k>Nc --
     ! this gives the same Lx*Ly for nrun=1 (all of num(1:Lx*Ly) are
     ! nonzero, the rest zero, by construction) and the correct restored
     ! value for nrun=2.
     Nc = count(num /= 0)

     ! PBC sanity check: a mesh built with if_periodic=.true.
     ! (Generate_Initial_Mesh.f90) has every vertex touched by exactly 3
     ! cells -- zero free-boundary (degree <3) vertices, by construction
     ! (the same invariant that generator itself verifies before ever
     ! writing the mesh out). Running such a mesh with if_PBC=.false. feeds
     ! cells that genuinely straddle the periodic wrap into non-wrap-aware
     ! geometry (CalculateArea/Force_Calculation/etc treat their raw,
     ! far-apart stored coordinates as one ordinary polygon) -- this
     ! reliably corrupts area/perimeter for those cells into huge, wrong
     ! values, blowing up forces into NaN within a handful of steps.
     ! Detect the mismatch here, before the run starts, instead of letting
     ! it silently NaN partway through.
     block
       integer, allocatable :: deg_check(:)
       integer :: kk, ll
       logical :: mesh_looks_periodic

       allocate(deg_check(v_dim2))
       deg_check = 0
       do kk = 1, Nc
         do ll = 1, num(kk)
           deg_check(inn(ll,kk)) = deg_check(inn(ll,kk)) + 1
         end do
       end do

       mesh_looks_periodic = .true.
       do kk = 1, v_dim2
         if (deg_check(kk) > 0 .and. deg_check(kk) < 3) then
           mesh_looks_periodic = .false.
           exit
         end if
       end do
       deallocate(deg_check)

       if ((.not. if_PBC) .and. mesh_looks_periodic) then
         write(*,*) 'read_data: STOPPING -- this mesh has ZERO free-boundary vertices', &
           ' (every vertex has degree exactly 3) -- it looks like a periodic mesh', &
           ' (if_periodic=.true. in Generate_Initial_Mesh.f90) but if_PBC is .false. in', &
           ' para_Simulation.dat. Running non-PBC-aware physics on a wrap-connected mesh', &
           ' corrupts area/perimeter for any cell straddling the wrap and reliably blows', &
           ' up into NaN within a few steps. Set if_PBC=.true. to match this mesh, or', &
           ' generate a non-periodic mesh (if_periodic=.false.) to match if_PBC=.false.'
         stop 1
       end if

       if (if_PBC .and. (.not. mesh_looks_periodic)) then
         write(*,*) 'read_data: WARNING -- if_PBC is .true. but this mesh has free-boundary', &
           ' (degree < 3) vertices -- it does not look like a mesh built with', &
           ' if_periodic=.true. in Generate_Initial_Mesh.f90. This is not immediately', &
           ' dangerous (if_Fixed_boundary is already forced off under if_PBC, and no real', &
           ' vertex pair is far enough apart to trigger a spurious wrap), but the tissue''s', &
           ' true free edge is now just an ordinary unpinned edge, not a periodic one --', &
           ' if_PBC is not doing what you likely intend with this mesh.'
       end if
     end block

    end subroutine read_data


    ! PBC restart (if_variable_box): the start snapshot
    ! (written at the top of step nrun2_initialTime) holds the positions
    ! after step nrun2_initialTime-1, so the box length and the Lees-Edwards
    ! strain of that same step are restored from the source run's
    ! data/BoxLength.dat (written only by if_variable_box runs -- missing =
    ! the box was fixed at Lx,Ly) and data/ShearStrain.dat (LE_strain after
    ! every step; all zeros if the source run did not shear). LE_offset =
    ! LE_strain*Ly_box mod Lx_box is exactly what Advance_LE_Strain and
    ! Variable_Box_Step would have left (isotropic scaling commutes with the
    ! modulo). Without this a restart silently reset the box and LE_offset.
    subroutine Restore_Box_And_Strain

      integer :: k, ios
      logical :: ex
      real*8, allocatable :: tmp(:)
      character(100) :: fname_box

      ! Preferred : the box state saved with the start snapshot.
      write(fname_box, '("data/Box_", I8.8,".dat")') nrun2_initialTime
      inquire(file=fname_box, exist=ex)
      if (ex) then
        open(unit=961, file=fname_box, form='unformatted', status='old')
        read(961) Lx_box, Ly_box, LE_strain, LE_offset
        close(961)
        write(*,*) 'Box restored from ', trim(fname_box), ': Lx_box =', Lx_box, ' Ly_box =', Ly_box, &
          ' LE_strain =', LE_strain
        return
      end if

      ! Fallback for source runs without Box_ files (branch pbc_variable_box):
      ! the per-step series, which is per step only if the source run had
      ! series_stride = 1.
      k = nrun2_initialTime - 1
      if (k .lt. 1) return
      allocate(tmp(k))

      inquire(file='data/BoxLength.dat', exist=ex)
      if (ex) then
        open(unit=961, file='data/BoxLength.dat', form='unformatted', status='old')
        read(961, iostat=ios) tmp
        close(961)
        if (ios .ne. 0) then
          write(*,*) 'Restore_Box_And_Strain: data/BoxLength.dat is shorter than step ', k
          stop 1
        end if
        Ly_box = tmp(k) * dble(Ly) / dble(Lx)
        Lx_box = tmp(k)
        write(*,*) 'Box restored from data/BoxLength.dat: Lx_box =', Lx_box, ' Ly_box =', Ly_box
      end if

      inquire(file='data/ShearStrain.dat', exist=ex)
      if (ex) then
        open(unit=961, file='data/ShearStrain.dat', form='unformatted', status='old')
        read(961, iostat=ios) tmp
        close(961)
        if (ios .ne. 0) then
          write(*,*) 'Restore_Box_And_Strain: data/ShearStrain.dat is shorter than step ', k
          stop 1
        end if
        LE_strain = tmp(k)
        LE_offset = modulo(LE_strain * Ly_box, Lx_box)
        if (LE_strain .ne. 0.0d0) write(*,*) 'Lees-Edwards strain restored: LE_strain =', LE_strain
      end if

      deallocate(tmp)

    end subroutine Restore_Box_And_Strain

    subroutine Open_Event_Logs
      ! T1/T2 event log : opened ONCE before the main it-loop
      ! (unlike write_output's snapshot files, which are reopened every
      ! dump) -- events are sparse and irregular, so they're appended to
      ! these two long-lived file handles as they happen, closed once at
      ! the very end (Close_Event_Logs). Same nrun2_ prefix convention as
      ! every other output file: a restart run's log never overwrites the
      ! original run's.
      implicit none
      character(100) :: fname_T1events, fname_T2events, fname_Divisionevents

      if (nrun.eq.1) then
        fname_T1events = 'data/T1_events.dat'
        fname_T2events = 'data/T2_events.dat'
        fname_Divisionevents = 'data/Division_events.dat'
      else
        fname_T1events = 'data/nrun2_T1_events.dat'
        fname_T2events = 'data/nrun2_T2_events.dat'
        fname_Divisionevents = 'data/nrun2_Division_events.dat'
      end if

      ! BINARY : one Fortran unformatted record per event (each
      ! WRITE statement on Do_T1/Do_T2's already-open unit produces its own
      ! fixed-size record, leading/trailing 4-byte length markers included
      ! -- the same convention every other output file in this codebase
      ! uses, just one small record per event instead of one big record for
      ! the whole run). All-real*8, matching this codebase's own habit of
      ! using real*8 even for integer-valued quantities (e.g.
      ! Total_T1_count) -- avoids ever mixing types within one unformatted
      ! record. Cell identities are stored as the numeric suffix of their
      ! 'cell_<N>' string (see CellIdNum below), 0.0d0 for 'cell_0'/no-cell
      ! padding -- matching this codebase's existing "0 means empty" idiom
      ! (e.g. num(i)==0). Roughly half the size of an equivalent formatted
      ! (ASCII) log, and keeps the same "readable while still running"
      ! property every other incrementally-written file here has.
      open(unit=iunit_T1events, file=fname_T1events, form='unformatted', status='replace')
      open(unit=iunit_T2events, file=fname_T2events, form='unformatted', status='replace')
      open(unit=iunit_Divisionevents, file=fname_Divisionevents, form='unformatted', status='replace')

    end subroutine Open_Event_Logs

    function CellIdNum(id_str) result(num_out)
      ! T1/T2 event log : cell_identity is always exactly
      ! 'cell_'//<integer> (allocation.f90's read_data, T2_core's 'cell_0'
      ! placeholder) -- extract that integer directly rather than storing
      ! the whole padded string, for the compact binary event log above.
      implicit none
      character(len=*), intent(in) :: id_str
      real*8 :: num_out
      integer :: n, ios
      read(id_str(6:), *, iostat=ios) n
      if (ios /= 0) then
        num_out = 0.0d0
      else
        num_out = dble(n)
      end if
    end function CellIdNum

    subroutine Seed_RNG
      ! rng_seed : > 0 seeds the generator, spreading the
      ! single integer over the seed array as Generate_Initial_Mesh does; 0
      ! leaves gfortran's own (random) seeding. Either way the generator state
      ! right after seeding (random_seed get; putting it back reproduces the
      ! stream) is written to data/[nrun2_]rng_seed.dat for provenance.
      integer :: n_seed, i_seed
      integer, allocatable :: seed_arr(:)
      character(100) :: fname_seed
      call random_seed(size=n_seed)
      allocate(seed_arr(n_seed))
      if (rng_seed > 0) then
        do i_seed = 1, n_seed
          seed_arr(i_seed) = rng_seed + 104729 * (i_seed - 1)
        end do
        call random_seed(put=seed_arr)
      end if
      call random_seed(get=seed_arr)
      if (nrun.eq.2) then
        fname_seed = 'data/nrun2_rng_seed.dat'
      else
        fname_seed = 'data/rng_seed.dat'
      end if
      open(unit=731, file=fname_seed, status='replace')
      write(731,*) rng_seed, '  ! rng_seed (0 = not seeded)'
      write(731,*) n_seed, '  ! length of the generator state below'
      write(731,*) seed_arr
      close(731)
      deallocate(seed_arr)
    end subroutine Seed_RNG

    subroutine Close_Event_Logs
      implicit none
      close(iunit_T1events)
      close(iunit_T2events)
      close(iunit_Divisionevents)
    end subroutine Close_Event_Logs

    subroutine write_output

      integer :: iunit_inn, iunit_num, iunit_v, iunit_force
      integer :: iunit_Myosin
      integer :: iunit_cell_identity
      integer :: iunit_motility
      integer :: iunit_eta
      character(100) :: fname_Energy, fname_ShearStress, fname_T1count, &
        fname_T2count, fname_motility, fname_eta, fname_ShearStrain, fname_BoxLength


       iunit_inn = 532
       iunit_num = 621
       iunit_v = 958
       iunit_force = 961
       iunit_Myosin = 966
       iunit_cell_identity = 967
       iunit_motility = 968
       ! eta : same per-frame-dump convention as motility just
       ! above (not a one-time "write once" file) -- eta(:) is already a
       ! per-vertex array that can evolve once a future flag/ODE drives
       ! it, so this needs to be visible per-frame from the start, same
       ! reasoning as motility's own bugfix comment right above.
       iunit_eta = 969

       if(nrun.eq.1)then
         write(fname_inn, '("data/inn_", I8.8,".dat")')(it)
         write(fname_num, '("data/num_", I8.8,".dat")')(it)
         write(fname_v, '("data/v_", I8.8,".dat")')(it)
         write(fname_force, '("data/force_", I8.8,".dat")')(it)
         write(fname_Myosin, '("data/Myosin_", I8.8,".dat")')(it)
         write(fname_cell_identity, '("data/cell_identity_", I8.8,".dat")')(it)
         ! BUGFIX : motility_store.dat used to be written ONCE, at
         ! it==1, and never again -- a frozen snapshot of the INITIAL
         ! per-vertex motility field. Any vertex created afterward (cell
         ! division; Proliferation.f90's mot inheritance fix) simply had no
         ! entry at all (its slot in the file was whatever the zero-filled
         ! array held at it=1), so Movie_Code.m's colorBy='Motility' showed
         ! every post-division cell as motility exactly 0 regardless of its
         ! true (correctly inherited) value -- confirmed directly: 598/598
         ! division-created vertices had frozen-file value 0.0 vs. nonzero
         ! live mot, in a real if_motility_hotspot+if_cell_division run.
         ! This also affected if_motility_decay/if_motility_Eulerian, which
         ! evolve mot over time for EVERY vertex, not just new ones. Now
         ! written every dump, same cadence/naming convention as every
         ! other per-vertex/per-cell array above -- no more special-cased
         ! "write once" file.
         write(fname_motility, '("data/motility_", I8.8,".dat")')(it)
         write(fname_eta, '("data/eta_", I8.8,".dat")')(it)
         fname_Energy = 'data/Energy.dat'
         fname_ShearStress = 'data/ShearStress.dat'
         fname_ShearStrain = 'data/ShearStrain.dat'
         fname_BoxLength = 'data/BoxLength.dat'
         fname_T1count = 'data/T1_count.dat'
         fname_T2count = 'data/T2_count.dat'
       elseif(nrun.eq.2)then
         write(fname_inn, '("data/nrun2_inn_", I8.8,".dat")')(it)
         write(fname_num, '("data/nrun2_num_", I8.8,".dat")')(it)
         write(fname_v, '("data/nrun2_v_", I8.8,".dat")')(it)
         write(fname_force, '("data/nrun2_force_", I8.8,".dat")')(it)
         write(fname_Myosin, '("data/nrun2_Myosin_", I8.8,".dat")')(it)
         write(fname_cell_identity, '("data/nrun2_cell_identity_", I8.8,".dat")')(it)
         write(fname_motility, '("data/nrun2_motility_", I8.8,".dat")')(it)
         write(fname_eta, '("data/nrun2_eta_", I8.8,".dat")')(it)
         ! BUGFIX : these 4 were still hardcoded to the same
         ! filenames as nrun=1 regardless of nrun -- unlike every other
         ! output above, an nrun=2 restart run silently overwrote the
         ! original nrun=1 run's Energy/ShearStress/T1_count/T2_count
         ! files with its own (shorter, restarted) data, destroying that
         ! part of the original run's record. Matched to the same nrun2_
         ! prefix convention as everything else here.
         fname_Energy = 'data/nrun2_Energy.dat'
         fname_ShearStress = 'data/nrun2_ShearStress.dat'
         fname_ShearStrain = 'data/nrun2_ShearStrain.dat'
         fname_BoxLength = 'data/nrun2_BoxLength.dat'
         fname_T1count = 'data/nrun2_T1_count.dat'
         fname_T2count = 'data/nrun2_T2_count.dat'
       end if



       open(unit = iunit_inn, file=fname_inn, form = 'unformatted', status='unknown')
       open(unit = iunit_num,file=fname_num, form ='unformatted', status='unknown')
       open(unit = iunit_v, file=fname_v, form='unformatted', status='unknown')
       open(unit = iunit_force, file=fname_force, form='unformatted', status='unknown')
       open(unit = iunit_Myosin,file=fname_Myosin, form ='unformatted', status='unknown')
       ! open(unit = iunit_force, file=fname_force, status='unknown')
       open(unit = iunit_cell_identity,file=fname_cell_identity, form ='unformatted', status='unknown')
       open(unit = iunit_motility,file=fname_motility, form ='unformatted', status='unknown')
       open(unit = iunit_eta,file=fname_eta, form ='unformatted', status='unknown')


       write(iunit_inn)((inn(i,j),i=1,inn_dim1),j=1,inn_dim2)
       write(iunit_num)(num(i), i=1,num_dim)
       write(iunit_v)((v(i,j), i=1,v_dim1),j=1,v_dim2)
       write(iunit_force)(fxx(i), fyy(i), fxx_ran(i), fyy_ran(i), &
        fxx_ABP(i), fyy_ABP(i), &
        fxx_Polar(i), fyy_Polar(i), i = 1, v_dim2)
       write(iunit_Myosin)(Rho(i), ROCK(i), Myosin(i), i = 1, num_dim)
       write(iunit_cell_identity)(cell_identity(i),  i=1,num_dim)
       write(iunit_motility)(mot(i), i = 1, v_dim2)
       write(iunit_eta)(eta(i), i = 1, v_dim2)
       ! ABP per-cell headings  -- only with if_ABP, so every other
       ! run's output is unchanged; nrun=2 restarts read them back.
       if (if_ABP) then
         if (nrun.eq.1) then
           write(fname_motility, '("data/ABP_theta_", I8.8,".dat")')(it)
         else
           write(fname_motility, '("data/nrun2_ABP_theta_", I8.8,".dat")')(it)
         end if
         open(unit = 970, file=fname_motility, form ='unformatted', status='unknown')
         write(970)(theta_ABP(i), i = 1, num_dim)
         close(970)
       end if
       ! Periodic box state of this snapshot : box lengths and the
       ! Lees-Edwards strain/offset that go with the dumped positions (both are
       ! values after step it-1, like v). nrun=2 restarts read this file
       ! (Restore_Box_And_Strain), so they do not depend on series_stride.
       if (if_PBC) then
         if (nrun.eq.1) then
           write(fname_motility, '("data/Box_", I8.8,".dat")')(it)
         else
           write(fname_motility, '("data/nrun2_Box_", I8.8,".dat")')(it)
         end if
         open(unit = 970, file=fname_motility, form ='unformatted', status='unknown')
         write(970) Lx_box, Ly_box, LE_strain, LE_offset
         close(970)
       end if
      
      
 
       ! CHANGE (MATLAB live-plotting request): these four used to
       ! only be written once, at it==totT -- meaning Energy/ShearStress/
       ! T1_count/T2_count did not exist at all until a run fully finished,
       ! which blocked plotting them while a run is still in progress. Now
       ! rewritten periodically (every summary_dump_interval, computed in
       ! read_input -- capped to ~200 rewrites over the whole run regardless
       ! of totT/it_dump, to avoid an O(totT^2/it_dump) I/O blowup; see the
       ! comment at summary_dump_interval's computation) and always on the
       ! final timestep, each time writing only the (1:it) slice actually
       ! computed so far (entries it+1:totT are not yet valid). The file
       ! simply grows longer at each such write; MATLAB's
       ! fread(fid,Inf,'float64') pattern (LoadData.m and friends) already
       ! reads "however much is there" with no hardcoded length assumption,
       ! so this needs no MATLAB-side change.
       if (modulo(it, summary_dump_interval).eq.0 .or. it.eq.totT) then

         open(unit = 711, file=fname_Energy, form='unformatted',  status='unknown')
         write(711)Energy(1:it_s)
         close(711)

         open(unit = 715, file=fname_ShearStress, form='unformatted',  status='unknown')
         write(715)ShearStress(1:it_s)
         close(715)

         open(unit = 716, file=fname_ShearStrain, form='unformatted',  status='unknown')
         write(716)ShearStrain(1:it_s)
         close(716)

         if (if_variable_box) then
           open(unit = 718, file=fname_BoxLength, form='unformatted',  status='unknown')
           write(718)BoxLength(1:it_s)
           close(718)
           if (box_aniso) then   ! Ly_box series 
             if (nrun.eq.1) then
               open(unit = 718, file='data/BoxLengthY.dat', form='unformatted',  status='unknown')
             else
               open(unit = 718, file='data/nrun2_BoxLengthY.dat', form='unformatted',  status='unknown')
             end if
             write(718)BoxLengthY(1:it_s)
             close(718)
           end if
         end if

         open(unit = 717, file=fname_T1count, form='unformatted',  status='unknown')
         write(717)Total_T1_count(1:it_s)
         close(717)

         open(unit = 719, file=fname_T2count, form='unformatted',  status='unknown')
         write(719)Total_T2_count(1:it_s)
         close(719)

       end if



 
 
       close(iunit_inn)
       close(iunit_num)
       close(iunit_v)
       close(iunit_force)
       ! BUGFIX : these two were never closed at all -- unlike
       ! inn/num/v/force just above, iunit_Myosin/iunit_cell_identity were
       ! left open after every single write_output call, relying on the
       ! NEXT call's open() (same fixed unit number, new filename) to
       ! implicitly close them. That leaves each snapshot's Myosin/
       ! cell_identity file in a not-yet-finalized state for up to a full
       ! it_dump interval (until the next dump event reopens the unit) --
       ! a real, independent bug from the read-side race this was found
       ! alongside (see that entry for the MATLAB-side half of the fix).
       close(iunit_Myosin)
       close(iunit_cell_identity)
       close(iunit_motility)
       close(iunit_eta)

    end subroutine write_output





end module allocation
