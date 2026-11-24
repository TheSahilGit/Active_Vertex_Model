module Force
  use allocation
  use Geometry

  contains


  subroutine Force_Calculation

    implicit none

    integer :: prev_idx, next_idx
    real*8 :: grad_perimeter_X, grad_perimeter_Y, grad_area_X, grad_area_Y
    real*8 :: rann
    ! OPTIMIZATION : edge_dx/edge_dy/edge_len(k) cache the vector and
    ! length of the edge from vertex k to vertex next(k), computed once per
    ! cell instead of twice (each edge was previously recomputed once as
    ! vertex jc's "next" edge and again as vertex jc+1's "prev" edge, halving
    ! the sqrt calls in this hot loop -- run every cell, every timestep).
    real*8 :: edge_dx(inn_dim1), edge_dy(inn_dim1), edge_len(inn_dim1)

    ! Differential line tension (if_free_edge_adhesion): per-edge
    ! gamm, not a single scalar -- is_free_edge is computed ONCE per call
    ! (whole-mesh, not per-cell) since Get_Free_Edges needs every cell's
    ! edges to classify any one of them. gamm_prev_edge/gamm_jc_edge pick
    ! gamm_free or gamm for the two specific edges touching vertex jc.
    ! When the flag is off, both are just gamm every time -- byte-identical
    ! to the old single-gamm force below.
    logical :: is_free_edge(inn_dim1, num_dim)
    real*8 :: gamm_prev_edge, gamm_jc_edge

    fxx = 0.0d0
    fyy = 0.0d0

    ! BUGFIX : beta/gamm nondimensionalization moved to read_input
    ! (allocation.f90), done once instead of every call -- for why
    ! repeating it here compounded beta/gamm every timestep.

    if (if_free_edge_adhesion) call Get_Free_Edges(is_free_edge)

    do ic = 1, Nc !Lx*Ly
      
      nn = num(ic)

      call Gather_Cell_Vertices_PBC(inn(1:nn,ic), nn, vx, vy)
      

      if(if_RhoROCK)then
        ! mr_feedback_off (noise-only control): the
        ! deterministic part alpha*Myosin is replaced by the constant
        ! mr_beta_mean; the coupling noise is unchanged. A separate branch,
        ! so the standard expressions below (and their floating-point
        ! rounding) are untouched.
        if(mr_feedback_off)then
          if(if_coupling_noise .and. mr_noise_mode == 1)then
            beta = (mr_beta_mean + mr_xi(ic))/(lambda*Ao)
          elseif(if_coupling_noise .and. .not. shrink_chem_frozen)then
            call random_number(rann)
            rann = rann - 0.5d0
            beta = (mr_beta_mean + coupling_noise_strength * rann)/(lambda*Ao)
          else
            beta = mr_beta_mean/(lambda*Ao)
          end if
        elseif(if_coupling_noise .and. mr_noise_mode == 1)then
          ! Ornstein-Uhlenbeck coupling noise , advanced
          ! once per step in Update_MR_Noise.
          beta = (Myosin_Coupling_Strength * Myosin(ic) + mr_xi(ic))/(lambda*Ao)
        elseif(if_coupling_noise .and. .not. shrink_chem_frozen)then
          call random_number(rann)
          rann = rann - 0.5d0
          ! BUGFIX : this branch omitted the /(lambda*Ao) factor that
          ! the else-branch below applies, silently changing beta's scale by
          ! a factor of (lambda*Ao) whenever coupling noise is enabled.
          beta = (Myosin_Coupling_Strength * Myosin(ic) + coupling_noise_strength * rann)/(lambda*Ao)
        else
          beta = Myosin_Coupling_Strength * Myosin(ic)/(lambda*Ao)
        end if
      end if

      if(if_active_contractility)then
        beta = Myosin(ic) / (lambda * Ao)
      end if
      beta_cell(ic) = beta   ! the beta this cell's force used (Calculate_StressTensor)



      call CalculateArea(vx,vy,nn,area)
      call CalculatePerimeter(vx,vy,nn,perimeter)
      ! NOTE : CalculateArea returns the signed shoelace area; only
      ! the magnitude is meaningful here (2D area has no physical sign), so
      ! area is intentionally take as abs(area). Confirmed with user: mesh
      ! winding is a fixed convention in this codebase, no per-cell sign
      ! handling needed in the gradient below.
      area = abs(area)

      ! edge_dx/dy/len(k) = vector/length of the edge from vertex k to
      ! vertex k's next neighbour (wraparound at nn).
      do jc = 1, nn
        next_idx = jc + 1
        if(jc == nn) next_idx = 1
        edge_dx(jc) = vx(next_idx) - vx(jc)
        edge_dy(jc) = vy(next_idx) - vy(jc)
        edge_len(jc) = dsqrt(edge_dx(jc)**2 + edge_dy(jc)**2)
      end do

      do jc = 1,num(ic)

       next_idx = jc + 1
       prev_idx = jc - 1
       if(jc == num(ic))then
         next_idx = 1
       elseif(jc == 1)then
         prev_idx = num(ic)
       end if

       if(prev_idx.eq.0)then
         write(*,*)'OOPS!!!', ic, jc, num(ic), inn(1:num(ic), ic)
       end if

        grad_area_X = 0.5d0 * (vy(prev_idx) - vy(next_idx))
        grad_area_Y = 0.5d0 * (vx(next_idx) - vx(prev_idx))

        ! dx/dy/len_d for the (prev_idx -> jc) edge = edge_dx/dy/len(prev_idx);
        ! for the (jc -> next_idx) edge = -edge_dx/dy(jc), edge_len(jc).
        grad_perimeter_X = edge_dx(prev_idx)/edge_len(prev_idx) - edge_dx(jc)/edge_len(jc)
        grad_perimeter_Y = edge_dy(prev_idx)/edge_len(prev_idx) - edge_dy(jc)/edge_len(jc)

        ! Differential line tension : grad_perimeter_X/Y above is
        ! really the SUM of two separate edge-length gradients (the
        ! prev_idx->jc edge and the jc->next_idx edge) -- when
        ! if_free_edge_adhesion is on, decompose back into its two terms so
        ! each can take its OWN gamm (free vs bulk) instead of one shared
        ! coefficient for both. Kept as an explicit branch (not folded into
        ! one formula with gamm_prev_edge=gamm_jc_edge=gamm when off) so the
        ! off-case is the ORIGINAL "-gamm*grad_perimeter_X" expression,
        ! bit-for-bit -- gamm*(A-B) is not guaranteed identical to
        ! gamm*A-gamm*B in floating point, which this codebase's
        ! byte-identical-regression testing convention depends on.
        if (if_free_edge_adhesion) then
          gamm_prev_edge = merge(gamm_free, gamm, is_free_edge(prev_idx, ic))
          gamm_jc_edge   = merge(gamm_free, gamm, is_free_edge(jc, ic))

          fxx(inn(jc,ic)) = fxx(inn(jc,ic)) - 2.0d0 * lambda * (area - Ao)* grad_area_X &
            - 2.0d0 * beta * (perimeter - Co)* grad_perimeter_X &
            - ( gamm_prev_edge * edge_dx(prev_idx)/edge_len(prev_idx) &
              - gamm_jc_edge * edge_dx(jc)/edge_len(jc) )

          fyy(inn(jc,ic)) = fyy(inn(jc,ic)) - 2.0d0 * lambda * (area - Ao)* grad_area_Y  &
            - 2.0d0 * beta * (perimeter - Co)* grad_perimeter_Y &
            - ( gamm_prev_edge * edge_dy(prev_idx)/edge_len(prev_idx) &
              - gamm_jc_edge * edge_dy(jc)/edge_len(jc) )
        else
          fxx(inn(jc,ic)) = fxx(inn(jc,ic)) - 2.0d0 * lambda * (area - Ao)* grad_area_X &
            - 2.0d0 * beta * (perimeter - Co)* grad_perimeter_X &
            - gamm * grad_perimeter_X

          fyy(inn(jc,ic)) = fyy(inn(jc,ic)) - 2.0d0 * lambda * (area - Ao)* grad_area_Y  &
            - 2.0d0 * beta * (perimeter - Co)* grad_perimeter_Y &
            - gamm * grad_perimeter_Y
        end if

      end do

   end do

   fxx = fxx 
   fyy = fyy 

  end subroutine Force_Calculation


  subroutine Motile_Force_Calculation

    implicit none

    integer :: nv
    ! OPTIMIZATION : fixed-size automatic (stack) arrays sized to
    ! the reserved capacity v_dim2 -- only the first nv elements are ever
    ! touched. Declared once per call, no heap allocate/deallocate.
    real*8 :: rann_x(v_dim2), rann_y(v_dim2), sigma_mot(v_dim2)

    ! OPTIMIZATION : mot/fxx_ran/fyy_ran are sized v_dim2, the
    ! RESERVED vertex capacity for future cell-division growth -- not the
    ! number of vertices actually in use right now. Vertices are only ever
    ! appended (Proliferation_Core writes new ones at maxval(inn)+1/+2,
    ! never compacted or reused), so maxval(inn) is a correct, cheap
    ! (O(inn_dim1*inn_dim2) integer scan) upper bound on which vertex
    ! indices are actually referenced by any live cell. Slots beyond it are
    ! never read by anything -- CalculateArea/CalculatePerimeter/
    ! Force_Calculation only ever iterate via inn(1:nn,ic).
    ! NOTE: with if_cell_division off (or growth far short of the reserved
    ! capacity), this trims only a small, fixed number of always-unused
    ! slots -- the much bigger win below is eliminating ~2*nv individual
    ! random_number() calls (measured as the dominant cost in profiling of
    ! a 10^7-step run) in favor of 2 whole-array calls, which removes the
    ! per-call subroutine overhead of tens of thousands of scalar calls per
    ! timestep, plus lets the compiler vectorize the arithmetic below.
    ! NOTE ON REPRODUCIBILITY: unlike the sqrt-caching/early-exit
    ! optimizations elsewhere in this codebase, this one is NOT guaranteed
    ! to reproduce the exact same trajectory as the old scalar-loop version
    ! for a given (unseeded) run -- both trimming to nv<v_dim2 and batching
    ! the random_number calls change how many/which calls consume the
    ! global RNG stream, which shifts everything drawn afterward (this
    ! codebase has no random_seed anywhere, so exact run-to-run
    ! reproducibility was never guaranteed regardless --). The
    ! physics is unaffected: every live vertex still gets an independent,
    ! identically-distributed noise term.
    nv = maxval(inn)

    if(if_motility_decay)then
     ! mot = mot0*exp(-it*dt/motility_decay_timeScale)

     mot(1:nv) = mot(1:nv) - dt * mot(1:nv)/motility_decay_timeScale
    end if

    call random_number(rann_x(1:nv))
    call random_number(rann_y(1:nv))

    ! sqrt(2*mot*eta) is identical for both components -- was computed
    ! twice per vertex (once for fxx_ran, once for fyy_ran); now once.
    ! Uses the LOCAL per-vertex eta(1:nv), not eta_0  -- the
    ! fluctuation-dissipation relation this implements requires the noise
    ! amplitude to match whatever friction actually acts at that vertex,
    ! not a global baseline, once eta varies in space.
    sigma_mot(1:nv) = sqrt(2.0d0 * mot(1:nv) * eta(1:nv))

    fxx_ran(1:nv) = sigma_mot(1:nv) * (2.0d0 * rann_x(1:nv) - 1.0d0)
    fyy_ran(1:nv) = sigma_mot(1:nv) * (2.0d0 * rann_y(1:nv) - 1.0d0)

  end subroutine Motile_Force_Calculation

  subroutine ABP_Force_Calculation
    implicit none

    real*8 :: rann, px, py
    integer :: kv
    integer, allocatable :: ncell_at(:)

    if (if_ABP)then

      ! Each cell is one self-propelled particle :
      ! one heading theta_ABP(ic) with its own rotational noise. Every
      ! vertex moves with the AVERAGE of vo*(cos,sin)(theta) over the cells
      ! it belongs to (as in Barton et al., PLoS Comput. Biol. 13, e1005569
      ! (2017)), so a cell whose neighbours share its heading moves at vo,
      ! and free-boundary vertices (1-2 cells) are treated alike.
      allocate(ncell_at(size(fxx_ABP)))
      fxx_ABP = 0.0d0
      fyy_ABP = 0.0d0
      ncell_at = 0

      do ic = 1, Nc

        call gaussian_random(rann, 0.0d0, 1.0d0)  ! number, mu, sigma
        rot_noise(ic) = sqrt(2.0d0 * Dr) * rann

        px = vo*dcos(theta_ABP(ic))
        py = vo*dsin(theta_ABP(ic))
        do jc = 1, num(ic)
          kv = inn(jc,ic)
          fxx_ABP(kv) = fxx_ABP(kv) + px
          fyy_ABP(kv) = fyy_ABP(kv) + py
          ncell_at(kv) = ncell_at(kv) + 1
        end do

      end do

      ! abp_force_mode 1: keep the SUM
      ! over the vertex's cells, i.e. the same distribution as the PM force.
      if (abp_force_mode == 0) then
        where (ncell_at > 0)
          fxx_ABP = fxx_ABP / dble(ncell_at)
          fyy_ABP = fyy_ABP / dble(ncell_at)
        end where
      end if
      deallocate(ncell_at)

    else
      fxx_ABP = 0.0d0
      fyy_ABP = 0.0d0
      theta_ABP = 0.0d0
      rot_noise = 0.0d0
    end if


  end subroutine ABP_Force_Calculation

  subroutine Polar_Motile_Force_Calculation
     implicit none
     real*8 :: rann_x, rann_y

     if(if_polar_motility)then

       ! BUGFIX : fxx_Polar/fyy_Polar were assigned with '=' instead
       ! of accumulated with '+', unlike every other per-vertex force in this
       ! module -- a vertex shared between cells simply got overwritten by
       ! whichever cell was processed last, discarding the others' contribution.
       fxx_Polar = 0.0d0
       fyy_Polar = 0.0d0

       do ic = 1, Nc
         call random_number(rann_x)
         call random_number(rann_y)

         do jc = 1, num(ic)
           fxx_Polar(inn(jc,ic)) = fxx_Polar(inn(jc,ic)) + & 
             polar_motility_strength * (rann_x - 0.50d0)
           fyy_Polar(inn(jc,ic)) = fyy_Polar(inn(jc,ic)) + & 
             polar_motility_strength * (rann_y - 0.50d0)
         end do
       end do


     else
       fxx_Polar = 0.0d0
       fyy_Polar = 0.0d0
     end if


  end subroutine Polar_Motile_Force_Calculation

  subroutine Give_Motility_Gradient

    implicit none

    ! OPTIMIZATION : the previous implementation binned referenced
    ! vertices into unit-width y "patches" (via verryCount/verry/veryyN/
    ! valid_indices) purely to decide which vertices to visit on each of
    ! numver outer-loop passes -- but the assignment formula inside always
    ! used the vertex's own continuous v(2,...), never anything
    ! patch-dependent, so every vertex ends up assigned the exact same
    ! value regardless of which patch pass caught it. The patches were
    ! therefore just an expensive (O(numver * total vertex-visits), i.e.
    ! roughly O(v_dim2^2) for a well-populated mesh) way of partitioning
    ! the vertex set into numver non-overlapping groups and then doing one
    ! ordinary full pass per group -- a single direct pass over every
    ! cell's vertices gives an identical result in O(total vertex-visits).
    ! Verified byte-identical against the old implementation (scratch
    ! test_motgrad.f90 harness, 3200-vertex/1024-cell mesh: 0 differing
    ! entries, maxdiff=0.0). This matters now because if_motility_Eulerian
    ! calls this every step (vertexmain.f90's main loop) instead of once
    ! before it -- the old O(v_dim2^2) cost per call would have been
    ! prohibitive at that frequency.
    mot = 0.0d0
    do ic = 1, Nc
        do jc = 1, num(ic)
            mot(inn(jc, ic)) = etas_max * exp(-v(2, inn(jc,ic)) / (mot_Lc*Ly))
        end do
    end do

end subroutine Give_Motility_Gradient


subroutine Give_Motility_Hotspot

  implicit none
  integer :: ic, jc, ip, jp, iv

!  integer, parameter :: number_of_hotspot = 4
  integer, dimension(:), allocatable :: hotspot_location  ! Cell index
  real*8, dimension(number_of_hotspot) :: xCM, yCM

!  real*8 :: sigma_hotspot
  real*8 :: xij, yij, rij
  logical, dimension(v_dim2) :: vertex_done

  print*,'number of hotspot, sigma_hotspot', number_of_hotspot, sigma_hotspot

  if(number_of_hotspot.eq.4)then

    allocate(hotspot_location(1:number_of_hotspot))
    hotspot_location = (/ 12, Ly-12, Lx*Ly-Ly+1+12,  Lx*Ly-12 /)

  elseif(number_of_hotspot.eq.2)then

    allocate(hotspot_location(1:number_of_hotspot))
    hotspot_location = (/ 12, Lx*Ly-12/)

  elseif(number_of_hotspot.eq.1)then

    allocate(hotspot_location(1:number_of_hotspot))
    hotspot_location = (/ int(dble(Ly)/2.0d0)/)

  else
    ! BUGFIX : without this branch, an unsupported number_of_hotspot
    ! left hotspot_location unallocated, and it is dereferenced unconditionally
    ! below (ip = hotspot_location(jp)) -- a crash. Fail loudly instead.
    write(*,*) 'ERROR: number_of_hotspot must be 1, 2, or 4. Got:', number_of_hotspot
    stop 1

  end if

!  print*, hotspot_location
!  print*, xCM
!  print*, yCM

  do jp = 1, number_of_hotspot
    ip = hotspot_location(jp)
    ! PBC : a plain average of raw v(:,inn(...)) is only correct
    ! if this hotspot cell doesn't straddle a periodic wrap -- use the
    ! same unwrap-relative-to-vertex-1 gather everything else uses.
    call Gather_Cell_Vertices_PBC(inn(1:num(ip),ip), num(ip), vx, vy)
    xCM(jp) = sum(vx(1:num(ip)))/dble(num(ip))
    yCM(jp) = sum(vy(1:num(ip)))/dble(num(ip))
  end do

 ! sigma_hotspot = 5.0d0

  ! BUGFIX : two issues found rechecking this against a user
  ! report of "no color in the Motility plot" (the actual cause turned
  ! out to be a MATLAB-side stale-file bug, but this recheck found real
  ! Fortran-side issues too, worth fixing regardless):
  !
  ! 1. mot was never zeroed here before accumulating into it, unlike its
  !    sibling Give_Motility_Gradient (which does `mot = 0.0` first).
  !    allocate(mot(v_dim2)) does NOT guarantee zero-initialized memory --
  !    this was silently adding the hotspot contribution on top of
  !    whatever undefined bytes happened to be there, working "by luck"
  !    whenever that memory happened to already be zero.
  !
  ! 2. The old loop visited each vertex once per INCIDENT CELL (do
  !    ic=1,Nc; do jc=1,num(ic)) and divided each contribution by a fixed
  !    3.0d0 -- which only gives the right total for an INTERIOR vertex
  !    (shared by exactly 3 cells: 3 visits x contribution/3 =
  !    contribution). A boundary vertex is shared by fewer than 3 cells
  !    (see compute_Circularity.m's border test), so it only got 1 or 2
  !    visits -- i.e. 1/3 or 2/3 of the intended value -- an artificial
  !    dip in motility right at the tissue edge that was purely a
  !    counting artifact, not physical. Fixed by visiting each unique
  !    vertex exactly once (vertex_done) and computing its full
  !    contribution directly, with no fudge-factor needed regardless of
  !    how many cells share it.
  mot = 0.0d0
  vertex_done = .false.

  do ic  = 1, Nc !Lx*Ly
    do jc = 1, num(ic)
      iv = inn(jc,ic)
      if (vertex_done(iv)) cycle
      vertex_done(iv) = .true.

      do ip = 1, number_of_hotspot
        xij = v(1,iv) - xCM(ip)
        yij = v(2,iv) - yCM(ip)
        ! PBC : minimum-image distance -- a vertex on the far
        ! side of a periodic wrap from this hotspot's center can still be
        ! physically close to it via the wrap-around direction.
        if (if_PBC) then
          xij = xij - Lx_box * dnint(xij / Lx_box)
          yij = yij - Ly_box * dnint(yij / Ly_box)
        end if
        rij = xij**2 + yij**2

        mot(iv) = mot(iv) + etas_max * exp(-(rij)/sigma_hotspot**2)
      end do

    end do
  end do

end subroutine Give_Motility_Hotspot

subroutine Apply_Fixed_Boundary

  implicit none

  integer :: im

  if(sum(Total_T2_count(1:it_s)).gt.0)then
     call Find_boundary_dynamic
   else 
     call Get_Boundary_info
   end if

  do ii = 1,size(boundary)
     im = boundary(ii)
     nn = num(im)
     fxx(inn(1:nn, im)) = 0.0d0
     fxx_ran(inn(1:nn, im)) = 0.0d0
     fyy(inn(1:nn, im)) = 0.0d0
     fyy_ran(inn(1:nn, im)) = 0.0d0
     fxx_ABP(inn(1:nn, im)) = 0.0d0
     fyy_ABP(inn(1:nn, im)) = 0.0d0
     fxx_Polar(inn(1:nn, im)) = 0.0d0
     fyy_Polar(inn(1:nn, im)) = 0.0d0

   end do


end subroutine Apply_Fixed_Boundary

subroutine Apply_bottom_border_Fixed
  implicit none

!  call Find_boundary_dynamic

!  print*, 'Bottom border fixed'

  fxx(bottom_border(1:bottom_border_count)) = 0.0d0
  fxx_ran(bottom_border(1:bottom_border_count)) = 0.0d0
  fyy(bottom_border(1:bottom_border_count)) = 0.0d0
  fyy_ran(bottom_border(1:bottom_border_count)) = 0.0d0
  fxx_ABP(bottom_border(1:bottom_border_count)) = 0.0d0
  fyy_ABP(bottom_border(1:bottom_border_count)) = 0.0d0
  fxx_Polar(bottom_border(1:bottom_border_count)) = 0.0d0
  fyy_Polar(bottom_border(1:bottom_border_count)) = 0.0d0


end subroutine  Apply_bottom_border_Fixed

subroutine Apply_top_border_Fixed
  implicit none

!  print*, 'Top border fixed'

  fxx(top_border(1:top_border_count)) = 0.0d0
  fxx_ran(top_border(1:top_border_count)) = 0.0d0
  fyy(top_border(1:top_border_count)) = 0.0d0
  fyy_ran(top_border(1:top_border_count)) = 0.0d0
  fxx_ABP(top_border(1:top_border_count)) = 0.0d0
  fyy_ABP(top_border(1:top_border_count)) = 0.0d0
  fxx_Polar(top_border(1:top_border_count)) = 0.0d0
  fyy_Polar(top_border(1:top_border_count)) = 0.0d0


end subroutine Apply_top_border_Fixed

subroutine Apply_perturbation
  implicit none
  real*8 :: comb


  if(if_sin_perturb)then
    if(it.ge.sin_perturb_when)then
       if(if_dirac_comb)then

          if(mod((it-sin_perturb_when), & 
            comb_onPeriod+comb_offPeriod) &
            < comb_onPeriod)then

            comb = 1.0d0
          else 
            comb = 0.0d0
          end if

          v(1,:) = v(1,:) + dt * sin_perturb_strength * &
            sin(2*pi*sin_perturb_waveNumber* v(2,:)/Ly) * comb


       end if
    end if
  end if




end subroutine Apply_perturbation


subroutine Apply_Limb_Force
! Limb force at four spots near the boundary, towards outward.  


  implicit none 

  integer :: ic, ip

  integer, dimension(4) :: applyforce_at


    applyforce_at = (/ 4, Ly-4, Lx*Ly-Ly+1+4,  Lx*Ly-4 /)

    do ip =  1, size(applyforce_at)
      ic = applyforce_at(ip)
      nn = num(ic)

      ! BUGFIX : every other force in the code converts to a
      ! displacement via dt*force/eta (see vertexmain.f90); this was missing
      ! the /eta mobility factor, making its effect inconsistent by a factor
      ! of eta relative to every other named force. Indexed eta(inn(1:nn,ic))
      ! (not eta_0), same as every other per-vertex use -- each affected
      ! vertex's own local mobility, consistent once eta varies in space.
      if(ic .lt. int(Lx*Ly/2))then   ! Cheap way
        v(1, inn(1:nn,ic)) = v(1, inn(1:nn,ic)) + dt * limb_force_strength * (-1.0d0) / eta(inn(1:nn,ic))
      else
        v(1, inn(1:nn,ic)) = v(1, inn(1:nn,ic)) + dt * limb_force_strength * (1.0d0) / eta(inn(1:nn,ic))
      end if

    end do


end subroutine Apply_Limb_Force

subroutine Update_MR_Noise
  ! MR coupling noise, mr_noise_mode 1 : per-cell
  ! Ornstein-Uhlenbeck process, Euler-Maruyama with Gaussian increments,
  ! d xi = -xi/tau dt + sd sqrt(2/tau) dW, stationary s.d. sd = mr_noise_sd.
  implicit none
  integer :: ic
  real*8 :: g
  do ic = 1, Nc
    call gaussian_random(g, 0.0d0, 1.0d0)
    mr_xi(ic) = mr_xi(ic) - dt * mr_xi(ic) / mr_noise_tau + mr_noise_sd * sqrt(2.0d0 * dt / mr_noise_tau) * g
  end do
end subroutine Update_MR_Noise

subroutine gaussian_random(rnumber, mu, sigma)
  implicit none
  real(8), intent(out) :: rnumber
  real(8), intent(in), optional :: mu, sigma

  real(8) :: u1, u2, r, theta, z
  real(8), parameter :: pi = 4.0d0 * datan(1.0d0)
  real(8) :: mean, std

  ! Set defaults
  mean = 0.0d0
  std  = 1.0d0
  if (present(mu)) mean = mu
  if (present(sigma)) std = sigma

  ! Generate two uniform random numbers in (0,1)
  call random_number(u1)
  call random_number(u2)

  ! Ensure u1 is not zero (log(0) issue)
  if (u1 < 1.0d-12) u1 = 1.0d-12

  ! Box–Muller transform
  r = dsqrt(-2.0d0 * dlog(u1))
  theta = 2.0d0 * pi * u2
  z = r * dcos(theta)

  ! Scale and shift
  rnumber = mean + std * z

end subroutine gaussian_random


subroutine Initialize_RhoROCK
  implicit none

  real*8 :: rann

  do ic = 1, Nc

    call random_number(rann)
    !Rho(ic) = rann
    !Rho(ic) = rann*1.0d-1
    Rho(ic) = rann*beta

    call random_number(rann)
    !ROCK(ic) = rann
    !ROCK(ic) = rann*1.0d-1
    ROCK(ic) = rann*beta

    call random_number(rann)
    !Myosin(ic) = rann
    !Myosin(ic) = rann*1.0d-1
    Myosin(ic) = rann*beta

  end do

!  Rho = 0.0d0
!  ROCK = 0.0d0
!  Myosin = 0.0d0

end subroutine Initialize_RhoROCK


subroutine Solve_RhoROCK_Euler
  implicit none

  real*8 :: Rho_RHS(size(num)), ROCK_RHS(size(num)), Myosin_RHS(size(num))
  real*8 :: area_diff
  real*8 :: f_Rho, f_ROCK, f_Myosin
  real*8 :: Myosin_noise(Nc)


  Myosin_noise = 0.0d0
  if(if_myosin_noise)then
    ! zero-mean noise, as in Solve_RhoROCK_RK4: Gaussian N(0,1) or uniform on (-1/2, 1/2)
    if(if_gaussian_noise)then
      do ic = 1, Nc
        call gaussian_random(Myosin_noise(ic), 0.0d0, 1.0d0)
      end do
    else
      call random_number(Myosin_noise)
      Myosin_noise = Myosin_noise - 0.5d0
    end if
  end if

  
  do ic = 1, Nc

    nn = num(ic)
    call Gather_Cell_Vertices_PBC(inn(1:nn,ic), nn, vx, vy)
    call CalculateArea(vx,vy,nn,area)
    area = abs(area)
    area_diff = area - Ao

    if (area_diff .lt. 0)then
      f_Rho = 0.0d0
    else
      f_Rho = A_Rho * area_diff**nhill /  (K_hill**nhill + area_diff**nhill) 
    end if


    f_ROCK = A_ROCK * Rho(ic)
    f_Myosin = A_Myosin * ROCK(ic)

    
    Rho_RHS(ic) = f_Rho * (1 - Rho(ic)) - D_Rho * Rho(ic)
    ROCK_RHS(ic) = f_ROCK * (1 - ROCK(ic)) - D_ROCK * ROCK(ic)
    Myosin_RHS(ic) = f_Myosin * (1 - Myosin(ic)) - D_Myosin * Myosin(ic)

  end do


  Rho(1:Nc) = Rho(1:Nc) +  dt * Rho_RHS(1:Nc)
  ROCK(1:Nc) = ROCK(1:Nc) + dt * ROCK_RHS(1:Nc)
  Myosin(1:Nc) = Myosin(1:Nc) + dt * Myosin_RHS(1:Nc) + & 
    sqrt(dt) * Myosin_noise_strength * Myosin_noise(1:Nc)



end subroutine Solve_RhoROCK_Euler

subroutine Solve_RhoROCK_RK4
  implicit none

  integer :: ic
  real*8 :: area, area_diff
  real*8 :: f_Rho, f_ROCK, f_Myosin

  ! RK4 intermediate arrays
  real*8 :: k1_Rho(Nc), k2_Rho(Nc), k3_Rho(Nc), k4_Rho(Nc)
  real*8 :: k1_ROCK(Nc), k2_ROCK(Nc), k3_ROCK(Nc), k4_ROCK(Nc)
  real*8 :: k1_M(Nc), k2_M(Nc), k3_M(Nc), k4_M(Nc)

  real*8 :: Rho_tmp(Nc), ROCK_tmp(Nc), M_tmp(Nc)

  real*8 :: Myosin_noise1(Nc), Myosin_noise2(Nc), propen_Myosin1(Nc), propen_Myosin2(Nc)
  real*8, parameter :: eps = 1.0d-3


  Myosin_noise1 = 0.0d0
  propen_Myosin1 = 0.0d0
  propen_Myosin2 = 0.0d0

  if(if_myosin_noise)then
    
    if(if_gaussian_noise)then
      do ic = 1, Nc
        call gaussian_random(Myosin_noise1(ic), 0.0d0, 1.0d0)
        call gaussian_random(Myosin_noise2(ic), 0.0d0, 1.0d0)
      end do
    else
      call random_number(Myosin_noise1)
      call random_number(Myosin_noise2)
      Myosin_noise1 = Myosin_noise1 - 0.5d0
      Myosin_noise2 = Myosin_noise2 - 0.5d0
    end if

    propen_Myosin1(1:Nc) = A_Myosin * ROCK(1:Nc) * (1.0d0 - Myosin(1:Nc))
    propen_Myosin2(1:Nc) = D_Myosin * Myosin(1:Nc)
  end if


  ! ============
  ! RK4 steps
  ! ============

  ! k1
  call compute_RHS(Rho, ROCK, Myosin, k1_Rho, k1_ROCK, k1_M)

  ! k2
  Rho_tmp  = Rho(1:Nc)   + 0.5d0 * dt * k1_Rho
  ROCK_tmp = ROCK(1:Nc)  + 0.5d0 * dt * k1_ROCK
  M_tmp    = Myosin(1:Nc) + 0.5d0 * dt * k1_M
  call compute_RHS(Rho_tmp, ROCK_tmp, M_tmp, k2_Rho, k2_ROCK, k2_M)

  ! k3
  Rho_tmp  = Rho(1:Nc)   + 0.5d0 * dt * k2_Rho
  ROCK_tmp = ROCK(1:Nc)  + 0.5d0 * dt * k2_ROCK
  M_tmp    = Myosin(1:Nc) + 0.5d0 * dt * k2_M
  call compute_RHS(Rho_tmp, ROCK_tmp, M_tmp, k3_Rho, k3_ROCK, k3_M)

  ! k4
  Rho_tmp  = Rho(1:Nc)   + dt * k3_Rho
  ROCK_tmp = ROCK(1:Nc)  + dt * k3_ROCK
  M_tmp    = Myosin(1:Nc) + dt * k3_M
  call compute_RHS(Rho_tmp, ROCK_tmp, M_tmp, k4_Rho, k4_ROCK, k4_M)

  ! Final RK4 update
  Rho(1:Nc) = Rho(1:Nc)   + (dt/6.0d0) * (k1_Rho + 2*k2_Rho + 2*k3_Rho + k4_Rho) 
  ROCK(1:Nc) = ROCK(1:Nc)  + (dt/6.0d0) * (k1_ROCK + 2*k2_ROCK + 2*k3_ROCK + k4_ROCK)
  Myosin(1:Nc) = Myosin(1:Nc) + (dt/6.0d0) * (k1_M + 2*k2_M + 2*k3_M + k4_M) &
    + sqrt(dt) * Myosin(1:Nc) * (1.0d0 - Myosin(1:Nc)) * & 
    Myosin_noise_strength * Myosin_noise1(1:Nc)


!! -----------------------------------------------------------------------------------------
!  Myosin(1:Nc) = Myosin(1:Nc) + (dt/6.0d0) * (k1_M + 2*k2_M + 2*k3_M + k4_M) &
!    + sqrt(dt) * sqrt(Propen_Myosin1(1:Nc)) * Myosin_noise_strength * Myosin_noise1(1:Nc) &
!    - sqrt(dt) * sqrt(Propen_Myosin2(1:Nc)) * Myosin_noise_strength * Myosin_noise2(1:Nc)

!do ic = 1, Nc
!  if(Myosin(ic) .lt. eps) then
!    Myosin(ic) = Myosin(ic) + (dt/6.0d0) * (k1_M(ic) + 2*k2_M(ic) + 2*k3_M(ic) + k4_M(ic)) &
!     + sqrt(dt) * sqrt(Propen_Myosin1(ic)) * Myosin_noise_strength * Myosin_noise1(ic) 

! elseif(Myosin(ic) .gt. 1.0d0 - eps)then
!    Myosin(ic) = Myosin(ic) + (dt/6.0d0) * (k1_M(ic) + 2*k2_M(ic) + 2*k3_M(ic) + k4_M(ic)) &
!      - sqrt(dt) * sqrt(Propen_Myosin2(ic)) * Myosin_noise_strength * Myosin_noise2(ic)

!  else
!    Myosin(ic) = Myosin(ic) + (dt/6.0d0) * (k1_M(ic) + 2*k2_M(ic) + 2*k3_M(ic) + k4_M(ic)) &
!     + sqrt(dt) * sqrt(Propen_Myosin1(ic)) * Myosin_noise_strength * Myosin_noise1(ic) &
!      - sqrt(dt) * sqrt(Propen_Myosin2(ic)) * Myosin_noise_strength * Myosin_noise2(ic)
!  end if

!  if (Myosin(ic) .lt. 0.0d0 .or. Myosin(ic) .gt. 1.0d0) then
!    print*, 'out of bounds', ic, Myosin(ic)
!  end if

!end do

!!--- Applying reflecting BC. 

!  do ic = 1, Nc
!    if(Myosin(ic) .gt. 1.0d0) then
!      Myosin(ic) = 2.0d0 - Myosin(ic)
!    end if
!    if(Myosin(ic) .lt. 0.0d0) then
!      Myosin(ic) =  - Myosin(ic)
!    end if
!  end do

!! -----------------------------------------------------------------------------------------

end subroutine Solve_RhoROCK_RK4

subroutine compute_RHS(Rho_in, ROCK_in, M_in, Rho_out, ROCK_out, M_out)
  ! BUGFIX : missing implicit none let f_Rho/f_ROCK/f_Myosin fall
  ! back to Fortran's default implicit typing (single-precision real) instead
  ! of real*8, silently degrading precision in the RK4 solver path relative
  ! to the equivalent Euler path (Solve_RhoROCK_Euler), which declares them
  ! explicitly as real*8.
  implicit none
  real*8, intent(in)  :: Rho_in(Nc), ROCK_in(Nc), M_in(Nc)
  real*8, intent(out) :: Rho_out(Nc), ROCK_out(Nc), M_out(Nc)
  real*8 :: area_diff
  real*8 :: f_Rho, f_ROCK, f_Myosin

  do ic = 1, Nc
    nn = num(ic)
    call Gather_Cell_Vertices_PBC(inn(1:nn,ic), nn, vx, vy)

    call CalculateArea(vx,vy,nn,area)
    area = abs(area)
    area_diff = area - Ao

    ! Rho production
    if (area_diff .lt. 0.0d0) then
        f_Rho = 0.0d0
    else
        f_Rho = A_Rho * area_diff**nhill / (K_hill**nhill + area_diff**nhill)
    end if

    ! ROCK and Myosin production
    f_ROCK   = A_ROCK   * Rho_in(ic)
    f_Myosin = A_Myosin * ROCK_in(ic)

    ! RHS
    Rho_out(ic)  = f_Rho   * (1.0d0 - Rho_in(ic))  - D_Rho   * Rho_in(ic)
    ROCK_out(ic) = f_ROCK  * (1.0d0 - ROCK_in(ic)) - D_ROCK  * ROCK_in(ic)
    M_out(ic)    = f_Myosin * (1.0d0 - M_in(ic))    - D_Myosin * M_in(ic)
  end do

end subroutine compute_RHS

subroutine Solve_Myosin_Act_Contr

  implicit none
  integer :: ic
  real*8 :: rann

  do ic = 1,Nc 

    call random_number(rann)
    rann = rann - 0.5d0

    Myosin(ic) = Myosin(ic) - dt * ( Myosin(ic) - beta_0)/tau_contr &
      + sqrt(dt) * active_contr_strength * rann
    ! fc_clip_negative: contractility cannot become negative (default off)
    if (fc_clip_negative .and. Myosin(ic) < 0.0d0) Myosin(ic) = 0.0d0

  end do


end subroutine Solve_Myosin_Act_Contr
  




  subroutine Excluded_Volume_Force
    ! Excluded volume (branch no_overlap). The vertex model has no
    ! contact force, so with free boundaries the boundary can fold onto itself
    ! and edges then cross (cells overlap): in the T2-off rounding runs every
    ! first crossing was between two free edges, with no T1 nearby, and the
    ! folds later got buried in the tissue. Here each vertex p is pushed away
    ! from every edge (a,b) that does not touch p,
    ! U = (ev_k/2) (ev_d0 - d)^2 for d < ev_d0, d = distance from p to the
    ! closest point c = a + s (b - a) of the edge. The force on p is
    ! ev_k (ev_d0 - d) (p - c)/d, and -(1-s) and -s of it go to a and b
    ! (equal and opposite, the exact gradient of U). For an edge that starts
    ! at a neighbour a of p, the closest point is usually a itself (the
    ! corner p-a-b, about 120 deg); that distance is the edge p-a, which T1
    ! handles, so such a pair counts only when the closest point lies inside
    ! the edge (a corner folding shut, or a collapsing triangle). Pairs are found with a uniform grid of bins of
    ! side ev_bin: each edge is listed in every bin its bounding box (grown by
    ! ev_d0) touches, each vertex looks only in its own bin. Open boundaries
    ! only (read_input stops for if_PBC). It does not enter Energy.dat.
    implicit none
    integer, parameter :: max_nb = 16
    real*8, parameter :: ev_bin = 0.5d0
    integer :: nv, ne, k, j, a, b, p, q, ix, iy, jx, jy, nbx, nby, ntot, ie
    integer, allocatable :: ea(:), eb(:), nnb(:), nb(:,:), bin_start(:), bin_fill(:), bin_edge(:)
    logical, allocatable :: used(:)
    real*8 :: xmin, ymin, xmax, ymax, abx, aby, l2, s, cx, cy, dx, dy, d, f, fx, fy
    logical :: near_a, near_b

    nv = 0
    do k = 1, Nc
      if (num(k) > 0) nv = max(nv, maxval(inn(1:num(k), k)))
    end do
    if (nv == 0) return
    ne = sum(num(1:Nc), mask=num(1:Nc) > 0)
    allocate(nnb(nv), nb(max_nb, nv), used(nv), ea(ne), eb(ne))
    nnb = 0; used = .false.; ne = 0

    ! unique edges and vertex neighbours, from the cell polygons
    do k = 1, Nc
      if (num(k) <= 0) cycle
      do j = 1, num(k)
        a = inn(j, k); b = inn(merge(1, j + 1, j == num(k)), k)
        used(a) = .true.; used(b) = .true.
        if (any(nb(1:nnb(a), a) == b)) cycle
        if (nnb(a) < max_nb) then
          nnb(a) = nnb(a) + 1; nb(nnb(a), a) = b
        end if
        if (nnb(b) < max_nb) then
          nnb(b) = nnb(b) + 1; nb(nnb(b), b) = a
        end if
        ne = ne + 1; ea(ne) = a; eb(ne) = b
      end do
    end do

    ! grid
    xmin = huge(1.0d0); ymin = huge(1.0d0); xmax = -huge(1.0d0); ymax = -huge(1.0d0)
    do p = 1, nv
      if (.not. used(p)) cycle
      xmin = min(xmin, v(1, p)); xmax = max(xmax, v(1, p)); ymin = min(ymin, v(2, p)); ymax = max(ymax, v(2, p))
    end do
    nbx = int((xmax - xmin) / ev_bin) + 1; nby = int((ymax - ymin) / ev_bin) + 1
    allocate(bin_start(nbx * nby + 1), bin_fill(nbx * nby))
    bin_fill = 0
    do ie = 1, ne                                   ! pass 1: entries per bin
      do jy = bin_of(min(v(2, ea(ie)), v(2, eb(ie))) - ev_d0, ymin, nby), &
              bin_of(max(v(2, ea(ie)), v(2, eb(ie))) + ev_d0, ymin, nby)
        do jx = bin_of(min(v(1, ea(ie)), v(1, eb(ie))) - ev_d0, xmin, nbx), &
                bin_of(max(v(1, ea(ie)), v(1, eb(ie))) + ev_d0, xmin, nbx)
          bin_fill((jy - 1) * nbx + jx) = bin_fill((jy - 1) * nbx + jx) + 1
        end do
      end do
    end do
    bin_start(1) = 1
    do k = 1, nbx * nby
      bin_start(k + 1) = bin_start(k) + bin_fill(k)
    end do
    ntot = bin_start(nbx * nby + 1) - 1
    allocate(bin_edge(max(ntot, 1)))
    bin_fill = 0
    do ie = 1, ne                                   ! pass 2: fill
      do jy = bin_of(min(v(2, ea(ie)), v(2, eb(ie))) - ev_d0, ymin, nby), &
              bin_of(max(v(2, ea(ie)), v(2, eb(ie))) + ev_d0, ymin, nby)
        do jx = bin_of(min(v(1, ea(ie)), v(1, eb(ie))) - ev_d0, xmin, nbx), &
                bin_of(max(v(1, ea(ie)), v(1, eb(ie))) + ev_d0, xmin, nbx)
          k = (jy - 1) * nbx + jx
          bin_edge(bin_start(k) + bin_fill(k)) = ie; bin_fill(k) = bin_fill(k) + 1
        end do
      end do
    end do

    ! vertex-edge repulsion
    do p = 1, nv
      if (.not. used(p)) cycle
      ix = bin_of(v(1, p), xmin, nbx); iy = bin_of(v(2, p), ymin, nby); k = (iy - 1) * nbx + ix
      do q = bin_start(k), bin_start(k + 1) - 1
        ie = bin_edge(q); a = ea(ie); b = eb(ie)
        if (a == p .or. b == p) cycle
        near_a = .false.; near_b = .false.
        do j = 1, nnb(p)
          if (nb(j, p) == a) near_a = .true.
          if (nb(j, p) == b) near_b = .true.
        end do
        abx = v(1, b) - v(1, a); aby = v(2, b) - v(2, a); l2 = abx * abx + aby * aby
        s = 0.0d0
        if (l2 > 0.0d0) s = min(1.0d0, max(0.0d0, ((v(1, p) - v(1, a)) * abx + (v(2, p) - v(2, a)) * aby) / l2))
        ! an edge that starts at a neighbour of p: its end at that neighbour is
        ! the corner p-a-b, whose distance is just the edge p-a (T1 territory);
        ! push only when the closest point is inside the edge (angle < 90 deg,
        ! a folding corner)
        if (near_a .and. s <= 0.0d0) cycle
        if (near_b .and. s >= 1.0d0) cycle
        cx = v(1, a) + s * abx; cy = v(2, a) + s * aby
        dx = v(1, p) - cx; dy = v(2, p) - cy; d = sqrt(dx * dx + dy * dy)
        if (d >= ev_d0 .or. d <= 1.0d-12) cycle
        f = ev_k * (ev_d0 - d) / d; fx = f * dx; fy = f * dy
        fxx(p) = fxx(p) + fx
        fyy(p) = fyy(p) + fy
        fxx(a) = fxx(a) - (1.0d0 - s) * fx
        fyy(a) = fyy(a) - (1.0d0 - s) * fy
        fxx(b) = fxx(b) - s * fx
        fyy(b) = fyy(b) - s * fy
      end do
    end do

    deallocate(nnb, nb, used, ea, eb, bin_start, bin_fill, bin_edge)

  contains

    integer function bin_of(x, x0, n)
      real*8, intent(in) :: x, x0
      integer, intent(in) :: n
      bin_of = min(n, max(1, int((x - x0) / ev_bin) + 1))
    end function bin_of

  end subroutine Excluded_Volume_Force

end module Force
