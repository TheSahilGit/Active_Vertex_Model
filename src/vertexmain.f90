! Here what I am not using any fortran library or subroutine to create the initial voronoi structure. 
! As, I need voronoi only to create the initial sturcture, I do it using matlab ( It is easier there) and write those
! data in a '*.dat' file, and then load those data in this fortran code. 
! It would not slow down the speed as we are gonna do it only once at the start and use the data file for the rest.
! (The time taken to create voronoi diagram with 100x100 points is approx 30 sec.) 
! REMEMBER : If we want to make changes like perfect honeycomb to randomly structured points, those we'll have to 
! do in the matlab code 'Main.m' . 



program vertexmain
  use allocation
  use array_info
  use System_Info
  use T1_Swap
  use T2_Swap
  use Force
  use Stress
  use Proliferation
  use ShrinkEdge


  call read_input
  call allocate_arrays
  call read_data   ! Initialization

  call Seed_RNG   ! rng_seed : before any random number is drawn

  call Open_Event_Logs   ! T1/T2 spatial/cell-identity event log 



  Total_T1_count = 0
  Total_T2_count = 0

  ! BUGFIX : Nc is now derived from the loaded num array at the
  ! end of read_data (allocation.f90), correctly for both nrun=1 and
  ! nrun=2 -- see that comment for why hardcoding it to Lx*Ly here broke
  ! nrun=2 restarts taken after any division/T2 event.

  if(if_motility)then
    print*, "Apolar Cell Motility, ", ' etas_max', etas_max
  end if

  if(if_motility_gradient)then
    if(if_motility_Eulerian)then
      print*, 'Motility Gradient (Eulerian -- recomputed from current position every step)'
    else
      print*, 'Motility Gradient (Lagrangian -- assigned once, carried per vertex)'
    end if
    call Give_Motility_Gradient
  elseif(if_motility_hotspot)then
    print*, 'Motility Hotspot'
    call Give_Motility_Hotspot
  else
     mot = etas_max
  end if
  mot0 = mot


  if(if_Fixed_boundary)then
    print*, "Fixed boundary"
  end if

   if(if_bottom_borders_fixed)then
     print*, "Bottom border fixed"
   end if

   if(if_top_borders_fixed)then
     print*, "Top border fixed"
   end if

   if(if_cell_division)then
     print*, "Cell division true. A_cut", area_0
   end if

   if(if_active_contractility)then
     print*, 'Active Contractility ; ', 'tau', tau_contr, 'strength:', active_contr_strength
     if(.not. myosin_restored)then
       Myosin = beta_0
     end if
   end if

   if(if_ABP)then
     print*, 'Active Brownian ; ', 'vo', vo, 'Dr', Dr
     ! Random initial orientations  -- with every theta_ABP at
     ! 0 the whole tissue starts self-propelling along +x, a collective
     ! drift that only decorrelates over ~1/Dr.
     ! One heading per cell ; kept if an nrun=2 restart restored it.
     if(.not. ABP_theta_restored)then
       call random_number(theta_ABP)
       theta_ABP = 2.0d0 * pi * theta_ABP
     end if
   end if


   if(if_RhoROCK)then
     print*, 'RhoROCK'
     if(.not. myosin_restored)then
       call  Initialize_RhoROCK
     end if
     ! MR coupling noise as an Ornstein-Uhlenbeck process (mr_noise_mode 1,
     !): one xi per cell slot, started from its stationary
     ! distribution (not carried over a restart).
     if (if_coupling_noise .and. mr_noise_mode == 1) then
       print*, 'MR coupling noise: Ornstein-Uhlenbeck, tau', mr_noise_tau, 'sd', mr_noise_sd
       allocate(mr_xi(num_dim))
       do ic = 1, num_dim
         call gaussian_random(mr_xi(ic), 0.0d0, mr_noise_sd)
       end do
     end if
   end if

  if(if_polar_motility)then
    print*, 'Polar Motility ; ' , 'Strength', polar_motility_strength
  end if

  if(if_shrink_edge)then
    call Init_Shrink_Edge
  end if




  do it = 1,totT

    it_s = (it - 1) / series_stride + 1   ! entry of the per-step series (series_stride)

    ! Differential line tension : branch at the call site, not
    ! inside Calculate_Energy, so that subroutine's compiled code is
    ! completely unaffected when if_free_edge_adhesion is off (see
    ! System_Info.f90's comment on Calculate_Energy_FreeEdge for why).
    ! With series_stride > 1 only the value at the end of each block (and at
    ! snapshot steps, which write the series) is kept, so the energy is only
    ! computed there ; stride 1 computes it every step.
    if (series_stride == 1 .or. modulo(it, series_stride) == 0 .or. modulo(it, it_dump) == 0 &
        .or. it == totT) then
    if (if_free_edge_adhesion) then
      call Calculate_Energy_FreeEdge
    elseif (if_active_contractility .or. if_RhoROCK) then
      call Calculate_Energy_CellBeta
    else
      call Calculate_Energy
    end if
    end if

    if(if_shrink_edge)then
      call Shrink_Edge_Stage_Check
    end if

!*********************************************************************
    if(modulo(it,it_dump).eq.0.or.it.eq.1.or.it.eq.2)then
      call write_output
    end if

!******************************************************************
   

   if(modulo(it,T2_time_interval).eq.0.0d0)then
     do iki= 1, 1
!       call Do_T2
     end do
   end if


    if(modulo(it,T1_time_interval).eq.0.0d0.or.it.eq.1)then
!    if(it==1)then
      do iki = 1,500

        if(if_Do_T1)then
          call Do_T1
        end if

        if(if_Do_T2)then
          call Do_T2
        end if

        ! OPTIMIZATION : find_T1/find_T2 (called at the top of
        ! Do_T1/Do_T2) do a full O(Nc) mesh scan every one of these 500
        ! iterations, even after no short edges/small triangles remain --
        ! count_T1/count_T2 are freshly recomputed by that scan on every
        ! call. Once both are simultaneously zero, nothing in the mesh can
        ! change for the rest of this loop (no candidates -> no possible
        ! mutation -> every remaining full rescan would deterministically
        ! find the same empty result), so stopping here skips zero
        ! legitimate T1/T2 events -- only the wasted rescans that would
        ! have found nothing anyway.
        if((.not.if_Do_T1 .or. count_T1==0) .and. (.not.if_Do_T2 .or. count_T2==0)) then
          exit
        end if

      end do
    end if

    if(if_shrink_edge)then
      call Check_Shrink_Edge
    end if





    if(if_RhoROCK .and. .not. shrink_chem_frozen)then
      if(if_RK4)then
        call Solve_RhoROCK_RK4
      else
        call Solve_RhoROCK_Euler
      end if
    end if

    if(if_active_contractility .and. .not. shrink_chem_frozen)then
      call Solve_Myosin_Act_Contr
    end if

    if (if_RhoROCK .and. if_coupling_noise .and. mr_noise_mode == 1 .and. .not. shrink_chem_frozen) then
      call Update_MR_Noise
    end if



    call Force_Calculation

    if(if_motility)then
      ! Eulerian motility : if_motility_gradient's
      ! initial assignment (above, before this loop) is a one-time
      ! Lagrangian tag -- correct per vertex ID forever after, but blind to
      ! where that vertex has since moved. if_motility_Eulerian re-derives
      ! mot from every vertex's CURRENT position each step instead, so the
      ! gradient stays anchored to physical location rather than drifting
      ! with the tissue. Cheap enough to call unconditionally every step
      ! after the Give_Motility_Gradient optimization (Force.f90) -- see
      ! that subroutine's comment. NOTE: if_motility_decay's effect is
      ! largely moot combined with this -- decay modifies mot in place
      ! after this call, but the very next step's recompute discards it,
      ! so decay never accumulates under Eulerian mode.
      if(if_motility_gradient .and. if_motility_Eulerian)then
        call Give_Motility_Gradient
      end if
      call Motile_Force_Calculation
    end if

     if(if_ABP)then
       call ABP_Force_Calculation
     end if

    if(if_polar_motility)then
      call Polar_Motile_Force_Calculation
    end if

    if(if_shrink_edge)then
      call Shrink_Edge_Force_Diag
    end if

    ! Excluded volume (branch no_overlap): vertex-edge repulsion,
    ! added to the deterministic force after every other force term.
    if(if_excluded_volume)then
      call Excluded_Volume_Force
    end if

    !! -- Apply Boundary conditions --

    if(if_Fixed_boundary)then
      call Apply_Fixed_Boundary
    end if

     if(if_bottom_borders_fixed)then
       call Apply_bottom_border_Fixed
     end if

     if(if_top_borders_fixed)then
       call Apply_top_border_Fixed
     end if
    
   !! ---------------------------------------



    v(1,:) = v(1,:) + dt * fxx(:)/eta(:) + sqrt(dt) * fxx_ran(:)/eta(:) + &
      dt * fxx_ABP(:) / eta(:) + &
      sqrt(dt) * fxx_Polar(:)/eta(:)


    v(2,:) = v(2,:) + dt * fyy(:)/eta(:) + sqrt(dt) * fyy_ran(:)/eta(:) + &
      dt * fyy_ABP(:) / eta(:) + &
      sqrt(dt) * fyy_Polar(:)/eta(:)

    if(if_ABP)then
      theta_ABP(:)  = theta_ABP(:) + dsqrt(dt) * rot_noise(:)
    end if


    
    if(if_Perturb_tissue)then
      call Apply_perturbation
    end if

    if(if_Shear_tissue)then
      call ShearTissue
      call Calculate_StressTensor
    end if

    ! Variable periodic box : after every force/shear move of the
    ! step, before the end-of-step wrap (Stress.f90, Variable_Box_Step).
    if(if_variable_box)then
      call Variable_Box_Step
      BoxLength(it_s) = Lx_box   ! value at the end of the series block (series_stride)
      BoxLengthY(it_s) = Ly_box
    end if



    if(if_limb_force)then
      call Apply_Limb_Force
    end if

    if(if_cell_division)then
      call Do_Proliferation
    end if

    if(if_squeeze_tissue)then
      if(it.eq.squeeze_when)then
        call Squeeze_Tissue
      end if
    end if

    ! PBC : wrap every vertex back into the canonical
    ! [0,Lx_box) x [0,Ly_box) box once per step, after every
    ! position-modifying step above has run. T1/T2/division already wrap
    ! their own newly-written/moved vertex positions individually at the
    ! point of writing (Wrap_Position calls in T1_swap.f90/T2_swap.f90/
    ! Proliferation.f90) -- wrapping again here is idempotent/harmless for
    ! those, and catches the ordinary force-integration position update
    ! (and Apply_perturbation/ShearTissue, if ever combined with if_PBC)
    ! in one place. A no-op when if_PBC is false.
    ! Lees-Edwards : a vertex crossing the y-boundary also
    ! crosses into an image box displaced by LE_offset in x -- same rule
    ! as Wrap_Position (Geometry.f90); LE_offset=0 leaves this the plain
    ! modulo wrap.
    if(if_PBC)then
      v(1,:) = modulo(v(1,:) - floor(v(2,:)/Ly_box) * LE_offset, Lx_box)
      v(2,:) = modulo(v(2,:), Ly_box)
    end if

    ! Edge-shrinking length constraint (ShrinkEdge.f90) -- last, after
    ! every position update of the step.
    if(if_shrink_edge)then
      call Enforce_Shrink_Length
    end if

    write(*,'(A,I0,A,I0)', advance='no') achar(13)//'Step: ', it, '/', totT

!   write(*,*)'it= ', it
!   write(*,*)'tttttttttttttttttttttttttttttttt'


  end do ! ends it

  call Close_Event_Logs

  if(if_shrink_edge)then
    call Close_Shrink_Edge
  end if

!  write(*,*)"Total T1 count", Total_T1_count
!  write(*,*)"Total T2 count", Total_T2_count

  write(*,*)"**********End******************"

end program vertexmain
