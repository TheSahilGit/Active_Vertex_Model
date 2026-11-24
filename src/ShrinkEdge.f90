module ShrinkEdge

  ! ---------------------------------------------------------------
  ! Edge shrinking -- T1 energy-landscape protocol after Bi, Lopez,
  ! Schwarz & Manning, Nat. Phys. 11, 1074 (2015): an edge is contracted
  ! quasi-statically to zero length, T1-swapped there, and the new edge
  ! expanded back to the original length; energy vs signed edge length
  ! (negative before the swap, positive after) gives the T1 path, and
  ! E(l->0) - E(l0) the barrier. (if_shrink_edge.)
  !
  ! Runs INSIDE the ordinary overdamped time loop -- every force, flag and
  ! topological rule behaves exactly as usual (so activity can be on).
  ! The only additions:
  !   - the edge's length is held fixed at a target: after every step the
  !     two vertices are put back exactly that far apart, symmetrically
  !     about their (freely moving) midpoint, along their (freely
  !     rotating) direction. For overdamped dynamics this is the same as
  !     removing the relative force along the edge, so with activity off
  !     the tissue gradient-descends to the constrained minimum E(l),
  !     which is Bi's protocol (they minimise by gradient descent in
  !     Surface Evolver at each length).
  !   - the target is stepped through stages of shrink_stage_steps steps:
  !       phase 0  hold at l0 (the length at it=1)            1 stage
  !       phase 1  l0 -> shrink_l_min in shrink_nstages        N stages
  !                at its end: the T1 of this edge, done by Do_T1 itself
  !                (find_T1's shrink_force_T1 branch), new edge set to
  !                shrink_l_min
  !       phase 2  shrink_l_min -> l0                          N+1 stages
  !       phase 3  released: no constraint, ordinary dynamics to totT
  !     Bi swap at exactly l=0; a zero-length edge has no direction, and
  !     its orientation relaxes at a rate ~ force/l (explicit Euler goes
  !     unstable below l ~ 1e-2 with dt=1e-3), hence
  !     shrink_l_min.
  !   - phases 0-2: find_T1 skips this edge (shrink_protect_edge), since
  !     it is driven below min_d_T1 on purpose. Every OTHER edge still
  !     T1s normally if if_Do_T1 is on.
  !   - shrink_freeze_beta (with if_RhoROCK or if_active_contractility):
  !     during phases 0-2 Myosin is not updated (no RhoROCK ODEs, no FC
  !     Ornstein-Uhlenbeck step) and the RhoROCK coupling noise is off, so
  !     every cell keeps its beta from the start snapshot -- Bi's protocol
  !     on the energy landscape of that instant (W = energy change again).
  !     Everything resumes at release.
  !   - shrink_relax_steps > 0: phase -1 before everything above -- the
  !     tissue relaxes freely for that many steps (edge unconstrained, beta
  !     frozen if shrink_freeze_beta), so the protocol starts from the
  !     local minimum of the (frozen) landscape; l0 = the edge length at
  !     the end of it. No rows are written during phase -1.
  !   - if the edge stops separating the expected two cells for any other
  !     reason (a neighbouring T1/T2 involving its vertices), the protocol
  !     aborts: constraint and protection off, noted in the output, and
  !     the run carries on as phase 4 (like phase 3: rows and snapshots
  !     every shrink_stage_steps, so the end-of-run summary files are
  !     still written when totT is a stage end).
  !
  ! Output: data/ShrinkEdge.dat (nrun2_ShrinkEdge.dat for nrun=2), one
  ! row at the end of every stage (and every shrink_stage_steps after
  ! release), plus an ordinary write_output snapshot at that same `it`
  ! (so main_MovieCode.m can render exactly the stage ends -- itList =
  ! column 1). Columns:
  !   it phase stage l_target l l_signed E E_mean dE_dl fmax dE_dl_mean
  ! E = Energy(it) at the stage end; E_mean = mean of Energy over the
  ! stage's second half (= E for a converged passive stage; the useful one
  ! with activity on); dE_dl = -(relative mechanical force along the
  ! edge)/2 per vertex, the constraint force, = dE/dl at a converged
  ! passive stage; fmax = largest mechanical vertex force with that
  ! constraint component removed (-> 0 at a passive constrained minimum);
  ! dE_dl_mean = dE_dl averaged over the stage's second half -- the
  ! integrand of the work W = int dE_dl dl done on the edge, which equals
  ! the energy change only when the forces derive from a fixed energy
  ! (not with RhoROCK, whose beta changes along the path).
  ! Header and events are '%' comment lines (MATLAB load() skips them).
  !
  ! data/(nrun2_)ShrinkEdge_cellE.dat: per-cell energies averaged over the
  ! same second-half window as E_mean, two rows per stage end:
  !   it 1 E_1 ... E_num_dim   (each cell's own beta, as Energy.dat)
  !   it 2 E_1 ... E_num_dim   (passive beta for every cell)
  ! column 2+ic is cell index ic (0 beyond Nc). Lets the energy of any
  ! patch around the edge be measured when single snapshots are too noisy
  ! (fluctuating contractility, RhoROCK coupling noise). The free-edge
  ! line tension (if_free_edge_adhesion) is not included.
  ! ---------------------------------------------------------------

  use allocation
  use Geometry
  use T1_swap

  character(500) :: shrink_cell_id(2)
  integer :: shrink_phase, shrink_stage, shrink_stage_start_it
  real*8 :: shrink_l0, shrink_l_target, shrink_dE_dl, shrink_fmax
  real*8 :: shrink_dE_dl_sum
  integer :: shrink_dE_dl_n
  integer, parameter :: iunit_shrink = 973, iunit_shrink_cellE = 974
  real*8, allocatable :: shrink_cellE_sum(:,:)
  integer :: shrink_cellE_n

  contains

  subroutine Init_Shrink_Edge
    implicit none
    integer :: ca, cb
    character(100) :: fname_shrink

    call Pick_Shrink_Edge(shrink_va, shrink_vb)
    call Find_Edge_Cells(shrink_va, shrink_vb, ca, cb)
    shrink_cell_id(1) = cell_identity(ca)
    shrink_cell_id(2) = cell_identity(cb)

    shrink_l0 = Edge_Length(shrink_va, shrink_vb)
    shrink_l_target = shrink_l0
    shrink_phase = 0
    shrink_stage = 0
    shrink_stage_start_it = 1
    shrink_protect_edge = .true.
    shrink_chem_frozen = shrink_freeze_beta .and. (if_RhoROCK .or. if_active_contractility)
    shrink_dE_dl = 0.0d0
    shrink_fmax = 0.0d0
    shrink_dE_dl_sum = 0.0d0
    shrink_dE_dl_n = 0
    allocate(shrink_cellE_sum(num_dim, 2))
    shrink_cellE_sum = 0.0d0
    shrink_cellE_n = 0

    if (nrun.eq.1) then
      fname_shrink = 'data/ShrinkEdge_cellE.dat'
    else
      fname_shrink = 'data/nrun2_ShrinkEdge_cellE.dat'
    end if
    open(unit=iunit_shrink_cellE, file=fname_shrink, status='replace')
    write(iunit_shrink_cellE,'(A)') '% per-cell energy, mean over each stage''s second half (ShrinkEdge.f90)'
    write(iunit_shrink_cellE,'(A)') '% columns: it kind E(cell 1) ... E(cell num_dim); kind 1 = own beta, 2 = passive beta'

    if (nrun.eq.1) then
      fname_shrink = 'data/ShrinkEdge.dat'
    else
      fname_shrink = 'data/nrun2_ShrinkEdge.dat'
    end if
    open(unit=iunit_shrink, file=fname_shrink, status='replace')

    write(iunit_shrink,'(A)') '% if_shrink_edge output (ShrinkEdge.f90) -- T1 path after Bi et al. 2015'
    write(iunit_shrink,'(A,I0,1X,I0)') '% edge_vertices ', shrink_va, shrink_vb
    write(iunit_shrink,'(A,A,1X,A)') '% edge_cells ', trim(shrink_cell_id(1)), trim(shrink_cell_id(2))
    write(iunit_shrink,'(A,ES24.16)') '% l0 ', shrink_l0
    write(iunit_shrink,'(A,I0,A,ES12.4,A,I0)') '% nstages ', shrink_nstages, &
      '  l_min ', shrink_l_min, '  stage_steps ', shrink_stage_steps
    if (shrink_chem_frozen) write(iunit_shrink,'(A)') '% beta frozen (Myosin not updated, no RhoROCK coupling noise) until release'
    ! FC's Ornstein-Uhlenbeck myosin is unbounded below. A cell frozen
    ! at beta < 0 has no energy minimum (-|beta|*P^2 lets it lower its
    ! energy without limit by growing its perimeter) and runs away over
    ! a long frozen protocol, so frozen FC myosin is clipped at 0
    ! . RhoROCK's Myosin is never negative.
    if (shrink_chem_frozen .and. if_active_contractility) then
      write(iunit_shrink,'(A,I0,A,ES12.4)') '% FC myosin clipped at 0 for the frozen protocol: ', &
        count(Myosin(1:Nc) < 0.0d0), ' cells, min was ', minval(Myosin(1:Nc))
      Myosin(1:Nc) = max(Myosin(1:Nc), 0.0d0)
    end if
    write(iunit_shrink,'(A)') '% phase: 0 hold at l0, 1 shrinking, 2 expanding after T1, 3 released, 4 aborted'
    write(iunit_shrink,'(A)') '% columns: it phase stage l_target l l_signed E E_mean dE_dl fmax dE_dl_mean'

    write(*,'(A,I0,A,I0,A,A,A,A,A,F8.5)') 'Shrink edge: vertices ', shrink_va, ' - ', &
      shrink_vb, '  (cells ', trim(shrink_cell_id(1)), ', ', trim(shrink_cell_id(2)), &
      ')  l0 = ', shrink_l0

    if (shrink_relax_steps > 0) then
      shrink_phase = -1
      shrink_protect_edge = .false.
      write(iunit_shrink,'(A,I0,A)') '% relaxing freely for ', shrink_relax_steps, ' steps before the protocol (phase -1)'
    end if

    if (totT < shrink_relax_steps + (2 * shrink_nstages + 2) * shrink_stage_steps + 1) then
      write(*,'(A,I0,A)') 'Shrink edge: WARNING -- totT too short to finish the protocol (needs ', &
        (2 * shrink_nstages + 2) * shrink_stage_steps + 1, ' steps)'
    end if
  end subroutine Init_Shrink_Edge


  ! Stage bookkeeping. Called right after the energy of step `it` is
  ! computed, so the configuration, Energy(it) and the snapshot all
  ! describe the end of the stage that just finished.
  subroutine Shrink_Edge_Stage_Check
    implicit none
    integer :: t1_before, ca, cb

    if (shrink_phase == -1) then
      if (it - 1 < shrink_relax_steps) return
      ! Free relaxation done: the protocol starts here, from this length.
      shrink_l0 = Edge_Length(shrink_va, shrink_vb)
      shrink_l_target = shrink_l0
      shrink_phase = 0
      shrink_stage = 0
      shrink_stage_start_it = it
      shrink_protect_edge = .true.
      shrink_dE_dl_sum = 0.0d0
      shrink_dE_dl_n = 0
      shrink_cellE_sum = 0.0d0
      shrink_cellE_n = 0
      write(iunit_shrink,'(A,I0,A,ES24.16)') '% relaxed at it=', it, ': l0 ', shrink_l0
      write(*,'(/,A,I0,A,F8.5)') 'Shrink edge: relaxed at it=', it, ', l0 = ', shrink_l0
      return
    end if

    if (2 * (it - shrink_stage_start_it) > shrink_stage_steps) call Accumulate_Cell_Energy
    if (it - shrink_stage_start_it < shrink_stage_steps) return

    call Write_Stage_Row
    call write_output

    select case (shrink_phase)
    case (0)
      shrink_phase = 1
      shrink_stage = 1
    case (1)
      if (shrink_stage < shrink_nstages) then
        shrink_stage = shrink_stage + 1
      else
        ! The T1 of this edge, through the ordinary Do_T1.
        t1_before = int(Total_T1_count(it))
        shrink_force_T1 = .true.
        call Do_T1
        shrink_force_T1 = .false.
        if (int(Total_T1_count(it)) == t1_before) then
          write(iunit_shrink,'(A,I0,A)') '% ABORTED at it=', it, ': Do_T1 rejected the T1 of the shrunk edge'
          write(*,'(/,A,I0)') 'Shrink edge: ABORTED -- Do_T1 rejected the T1 at it=', it
          shrink_stage_start_it = it
          call Abort_Shrink
          return
        end if
        call Find_Edge_Cells(shrink_va, shrink_vb, ca, cb)
        shrink_cell_id(1) = cell_identity(ca)
        shrink_cell_id(2) = cell_identity(cb)
        write(iunit_shrink,'(A,I0,A,A,1X,A)') '% T1 at it=', it, ': edge now between cells ', &
          trim(shrink_cell_id(1)), trim(shrink_cell_id(2))
        write(*,'(/,A,I0)') 'Shrink edge: T1 done at it=', it
        shrink_phase = 2
        shrink_stage = 0
      end if
    case (2)
      if (shrink_stage < shrink_nstages) then
        shrink_stage = shrink_stage + 1
      else
        shrink_phase = 3
        shrink_stage = 0
        shrink_protect_edge = .false.
        shrink_chem_frozen = .false.
        write(iunit_shrink,'(A,I0,A)') '% released at it=', it, ': constraint off'
        write(*,'(/,A,I0)') 'Shrink edge: protocol done, edge released at it=', it
      end if
    case (3, 4)
      shrink_stage = shrink_stage + 1
    end select

    shrink_stage_start_it = it
    if (shrink_phase <= 2) then
      shrink_l_target = Stage_Length(shrink_phase, shrink_stage)
      call Set_Edge_Length(shrink_va, shrink_vb, shrink_l_target)
    end if
  end subroutine Shrink_Edge_Stage_Check


  real*8 function Stage_Length(phase, stage)
    implicit none
    integer, intent(in) :: phase, stage
    select case (phase)
    case (1)
      Stage_Length = shrink_l0 + (shrink_l_min - shrink_l0) * dble(stage) / dble(shrink_nstages)
    case (2)
      Stage_Length = shrink_l_min + (shrink_l0 - shrink_l_min) * dble(stage) / dble(shrink_nstages)
    case default
      Stage_Length = shrink_l0
    end select
  end function Stage_Length


  subroutine Write_Stage_Row
    implicit none
    integer :: i0
    real*8 :: l, l_signed, E_mean, dE_dl_mean

    l = Edge_Length(shrink_va, shrink_vb)
    l_signed = merge(-l, l, shrink_phase <= 1)
    i0 = it - shrink_stage_steps / 2 + 1
    E_mean = sum(Energy(i0:it)) / dble(it - i0 + 1)
    dE_dl_mean = 0.0d0
    if (shrink_dE_dl_n > 0) dE_dl_mean = shrink_dE_dl_sum / dble(shrink_dE_dl_n)
    shrink_dE_dl_sum = 0.0d0
    shrink_dE_dl_n = 0
    write(iunit_shrink,'(I10,1X,I1,1X,I4,6(1X,ES24.16),1X,ES12.4,1X,ES24.16)') it, shrink_phase, &
      shrink_stage, merge(shrink_l_target, l, shrink_phase <= 2), l, l_signed, &
      Energy(it), E_mean, shrink_dE_dl, shrink_fmax, dE_dl_mean
    flush(iunit_shrink)

    if (shrink_cellE_n > 0) shrink_cellE_sum = shrink_cellE_sum / dble(shrink_cellE_n)
    write(iunit_shrink_cellE,'(I10,1X,I1,*(1X,ES16.8))') it, 1, shrink_cellE_sum(:,1)
    write(iunit_shrink_cellE,'(I10,1X,I1,*(1X,ES16.8))') it, 2, shrink_cellE_sum(:,2)
    flush(iunit_shrink_cellE)
    shrink_cellE_sum = 0.0d0
    shrink_cellE_n = 0
  end subroutine Write_Stage_Row


  ! Per-cell energies of step it, added to the stage's running sums. Same
  ! per-cell beta (and precedence) as Calculate_Energy_CellBeta.
  subroutine Accumulate_Cell_Energy
    implicit none
    integer :: jc, nv
    real*8 :: a, p, beta_cell, beta_pass, e0
    real*8, allocatable :: x(:), y(:)

    beta_pass = beta_0 / (lambda * Ao)
    do jc = 1, Nc
      nv = num(jc)
      call Gather_Cell_Vertices_PBC(inn(1:nv,jc), nv, x, y)
      call CalculateArea(x, y, nv, a)
      call CalculatePerimeter(x, y, nv, p)
      a = abs(a)
      if (if_active_contractility) then
        beta_cell = Myosin(jc) / (lambda * Ao)
      elseif (if_RhoROCK) then
        beta_cell = Myosin_Coupling_Strength * Myosin(jc) / (lambda * Ao)
      else
        beta_cell = beta_pass
      end if
      e0 = lambda * (a - Ao)**2 + gamm * p
      shrink_cellE_sum(jc,1) = shrink_cellE_sum(jc,1) + e0 + beta_cell * (p - Co)**2
      shrink_cellE_sum(jc,2) = shrink_cellE_sum(jc,2) + e0 + beta_pass * (p - Co)**2
    end do
    shrink_cellE_n = shrink_cellE_n + 1
  end subroutine Accumulate_Cell_Energy


  ! Abort if the edge no longer separates the expected two cells
  ! (tracked by cell_identity, so T2 renumbering of indices can't fool
  ! it). Called after the step's T1/T2 block.
  subroutine Check_Shrink_Edge
    implicit none
    integer :: ca, cb
    character(500) :: id_a, id_b

    if (shrink_phase > 2) return

    call Find_Edge_Cells(shrink_va, shrink_vb, ca, cb)
    id_a = ''
    id_b = ''
    if (ca > 0) id_a = cell_identity(ca)
    if (cb > 0) id_b = cell_identity(cb)
    if ((id_a == shrink_cell_id(1) .and. id_b == shrink_cell_id(2)) .or. &
        (id_a == shrink_cell_id(2) .and. id_b == shrink_cell_id(1))) return

    write(iunit_shrink,'(A,I0,A)') '% ABORTED at it=', it, &
      ': edge no longer between the expected cells (another topological change)'
    write(*,'(/,A,I0)') 'Shrink edge: ABORTED -- edge changed by another topological event at it=', it
    call Abort_Shrink
  end subroutine Check_Shrink_Edge


  subroutine Abort_Shrink
    implicit none
    ! shrink_stage_start_it is left alone -- the run stays on the same
    ! stage-end grid, so totT (a stage end) still gets its final snapshot.
    shrink_phase = 4
    shrink_stage = 0
    shrink_protect_edge = .false.
    shrink_chem_frozen = .false.
    flush(iunit_shrink)
  end subroutine Abort_Shrink


  ! Diagnostics only (forces are not modified -- the constraint acts on
  ! positions, Enforce_Shrink_Length). Must come after Force_Calculation.
  subroutine Shrink_Edge_Force_Diag
    implicit none
    integer :: nv, k
    real*8 :: ex, ey, s, fax, fay, fbx, fby

    if (shrink_phase > 2 .or. shrink_phase < 0) return

    call Edge_Unit_Vector(shrink_va, shrink_vb, ex, ey)
    s = 0.5d0 * ((fxx(shrink_vb) - fxx(shrink_va)) * ex + (fyy(shrink_vb) - fyy(shrink_va)) * ey)
    shrink_dE_dl = -s
    ! Second half of the stage (the same window as E_mean, one step
    ! earlier: these forces move the configuration into step it+1).
    if (2 * (it - shrink_stage_start_it) >= shrink_stage_steps) then
      shrink_dE_dl_sum = shrink_dE_dl_sum + shrink_dE_dl
      shrink_dE_dl_n = shrink_dE_dl_n + 1
    end if

    fax = fxx(shrink_va) + s * ex; fay = fyy(shrink_va) + s * ey
    fbx = fxx(shrink_vb) - s * ex; fby = fyy(shrink_vb) - s * ey
    nv = maxval(inn)
    shrink_fmax = max(fax**2 + fay**2, fbx**2 + fby**2)
    do k = 1, nv
      if (k == shrink_va .or. k == shrink_vb) cycle
      shrink_fmax = max(shrink_fmax, fxx(k)**2 + fyy(k)**2)
    end do
    shrink_fmax = sqrt(shrink_fmax)
  end subroutine Shrink_Edge_Force_Diag


  ! The constraint. Called at the very end of the step, after every
  ! position update (and the PBC wrap).
  subroutine Enforce_Shrink_Length
    implicit none
    if (shrink_phase > 2 .or. shrink_phase < 0) return
    call Set_Edge_Length(shrink_va, shrink_vb, shrink_l_target)
  end subroutine Enforce_Shrink_Length


  subroutine Close_Shrink_Edge
    implicit none
    close(iunit_shrink)
    close(iunit_shrink_cellE)
  end subroutine Close_Shrink_Edge


  ! Unit vector from va to vb (minimum image under PBC).
  subroutine Edge_Unit_Vector(va, vb, ex, ey)
    implicit none
    integer, intent(in) :: va, vb
    real*8, intent(out) :: ex, ey
    real*8 :: xb, yb, l
    xb = v(1,vb); yb = v(2,vb)
    call Unwrap_Relative(v(1,va), v(2,va), xb, yb)
    l = dsqrt((xb - v(1,va))**2 + (yb - v(2,va))**2)
    ex = (xb - v(1,va)) / l
    ey = (yb - v(2,va)) / l
  end subroutine Edge_Unit_Vector


  real*8 function Edge_Length(va, vb)
    implicit none
    integer, intent(in) :: va, vb
    real*8 :: xb, yb
    xb = v(1,vb); yb = v(2,vb)
    call Unwrap_Relative(v(1,va), v(2,va), xb, yb)
    Edge_Length = dsqrt((xb - v(1,va))**2 + (yb - v(2,va))**2)
  end function Edge_Length


  ! Place va/vb symmetrically about their current midpoint, along their
  ! current direction, exactly l apart.
  subroutine Set_Edge_Length(va, vb, l)
    implicit none
    integer, intent(in) :: va, vb
    real*8, intent(in) :: l
    real*8 :: xb, yb, xm, ym, ex, ey

    xb = v(1,vb); yb = v(2,vb)
    call Unwrap_Relative(v(1,va), v(2,va), xb, yb)
    xm = 0.5d0 * (v(1,va) + xb)
    ym = 0.5d0 * (v(2,va) + yb)
    call Edge_Unit_Vector(va, vb, ex, ey)

    v(1,va) = xm - 0.5d0 * l * ex
    v(2,va) = ym - 0.5d0 * l * ey
    v(1,vb) = xm + 0.5d0 * l * ex
    v(2,vb) = ym + 0.5d0 * l * ey

    call Wrap_Position(v(1,va), v(2,va))
    call Wrap_Position(v(1,vb), v(2,vb))
  end subroutine Set_Edge_Length


  ! The edge to shrink: shrink_v1/shrink_v2 if both given, otherwise the
  ! first interior edge of the cell whose centroid is nearest the tissue
  ! centre (box centre under PBC). Either way it must be an interior edge
  ! -- shared by exactly two cells, both endpoints 3-fold -- the only kind
  ! a T1 acts on.
  subroutine Pick_Shrink_Edge(va, vb)
    implicit none
    integer, intent(out) :: va, vb

    integer :: icell, kk, best_cell
    real*8 :: cx, cy, tx, ty, d2, best_d2
    real*8, allocatable :: ccx(:), ccy(:)

    if (shrink_v1 > 0 .and. shrink_v2 > 0) then
      va = shrink_v1
      vb = shrink_v2
      if (.not. Is_Interior_Edge(va, vb)) then
        write(*,'(A,I0,A,I0,A)') 'Pick_Shrink_Edge: vertices ', va, ' and ', vb, &
          ' are not an interior edge (consecutive in two cells, both endpoints 3-fold).'
        stop 1
      end if
      return
    end if

    ! Tissue centre = mean of cell centroids (each centroid PBC-unwrapped
    ! per cell). Under PBC the box centre is used instead -- a mean of
    ! wrapped centroids is not meaningful on a torus.
    allocate(ccx(Nc), ccy(Nc))
    do icell = 1, Nc
      call Gather_Cell_Vertices_PBC(inn(1:num(icell),icell), num(icell), vx, vy)
      ccx(icell) = sum(vx) / dble(num(icell))
      ccy(icell) = sum(vy) / dble(num(icell))
    end do
    if (if_PBC) then
      tx = 0.5d0 * Lx_box
      ty = 0.5d0 * Ly_box
    else
      tx = sum(ccx) / dble(Nc)
      ty = sum(ccy) / dble(Nc)
    end if

    best_cell = 0
    best_d2 = huge(1.0d0)
    do icell = 1, Nc
      cx = ccx(icell) - tx
      cy = ccy(icell) - ty
      if (if_PBC) then
        cx = cx - Lx_box * dnint(cx / Lx_box)
        cy = cy - Ly_box * dnint(cy / Ly_box)
      end if
      d2 = cx**2 + cy**2
      if (d2 < best_d2) then
        best_d2 = d2
        best_cell = icell
      end if
    end do
    deallocate(ccx, ccy)

    do kk = 1, num(best_cell)
      va = inn(kk, best_cell)
      vb = inn(merge(1, kk + 1, kk == num(best_cell)), best_cell)
      if (Is_Interior_Edge(va, vb)) return
    end do

    write(*,'(A,I0,A)') 'Pick_Shrink_Edge: centre cell ', best_cell, &
      ' has no interior edge -- tissue too small? Give shrink_v1/shrink_v2 explicitly.'
    stop 1

  end subroutine Pick_Shrink_Edge


  logical function Is_Interior_Edge(va, vb)
    implicit none
    integer, intent(in) :: va, vb
    integer :: ca, cb
    call Find_Edge_Cells(va, vb, ca, cb)
    Is_Interior_Edge = (ca > 0 .and. cb > 0 .and. &
      Vertex_Degree(va) == 3 .and. Vertex_Degree(vb) == 3)
  end function Is_Interior_Edge


  ! The (up to two) cells in which va and vb are consecutive vertices;
  ! 0 where there is no such cell.
  subroutine Find_Edge_Cells(va, vb, ca, cb)
    implicit none
    integer, intent(in) :: va, vb
    integer, intent(out) :: ca, cb
    integer :: icell, kk, knext

    ca = 0
    cb = 0
    do icell = 1, Nc
      do kk = 1, num(icell)
        knext = merge(1, kk + 1, kk == num(icell))
        if ((inn(kk,icell) == va .and. inn(knext,icell) == vb) .or. &
            (inn(kk,icell) == vb .and. inn(knext,icell) == va)) then
          if (ca == 0) then
            ca = icell
          else
            cb = icell
          end if
          exit
        end if
      end do
    end do
  end subroutine Find_Edge_Cells


  integer function Vertex_Degree(vid)
    implicit none
    integer, intent(in) :: vid
    integer :: icell
    Vertex_Degree = 0
    do icell = 1, Nc
      if (any(inn(1:num(icell),icell) == vid)) Vertex_Degree = Vertex_Degree + 1
    end do
  end function Vertex_Degree

end module ShrinkEdge
