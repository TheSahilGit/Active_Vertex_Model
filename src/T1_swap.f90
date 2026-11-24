module T1_swap
  use array_info
  use Geometry
  use allocation

  integer, allocatable, dimension(:) :: Affected, Occurrences

  logical :: T1_pass
  ! No-crossing T1 (branch no_overlap; only with if_excluded_volume):
  ! length factor of the new edge, 1 except while Do_T1 retries a flip whose
  ! new edges crossed another edge.
  real*8 :: t1_len_fac = 1.0d0
  ! Flips skipped this step (every try crossed): find_T1 leaves these edges
  ! out until the next step, so the 500-iteration T1 loop does not retry them.
  integer :: rej_it = -1, n_rej = 0, rej_a(64), rej_b(64)


  contains


    subroutine find_T1
      implicit none
      integer :: n_tri, it_tri, tri_cell(Nc)   ! T2-off speed fix (see below)
      logical :: on_tri

      ! Edge shrinking (ShrinkEdge.f90): when the protocol asks
      ! for the T1 of its own edge, that edge is the only candidate --
      ! everything after find_T1 in Do_T1 (rejection checks, event log,
      ! T1_core, counting) then runs exactly as for any other T1.
      if (shrink_force_T1) then
        call find_T1_Shrink_Edge
        return
      end if

      cell_no = 0
      ver_no = 0
      d_val = 0


      ! T2-off speed fix (branch no_overlap): with T2 off, triangles
      ! collapse and their edges stay below min_d_T1. Do_T1 always refuses
      ! those flips (a cell may not drop below 3 vertices), so they never
      ! leave the candidate list and the 500-iteration T1 loop rescans the
      ! mesh on every step (MR fell to ~1 time unit per minute). With T2 off,
      ! edges that belong to a triangle are left out: the accepted T1s are the
      ! same, only the random-number order changes. T2-on runs are unchanged.
      n_tri = 0
      if (.not. if_Do_T2) then
        do ic = 1, Nc
          if(num(ic).eq.3)then
            n_tri = n_tri + 1
            tri_cell(n_tri) = ic
          end if
        end do
      end if

      count_T1 = 0
      do ic = 1, Nc !Lx*Ly

        if(num(ic).le.3)cycle


        call Gather_Cell_Vertices_PBC(inn(1:num(ic), ic), num(ic), vx, vy)

        do jc = 1, num(ic)

          next_idx = jc + 1
          prev_idx = jc - 1

          if(jc == num(ic))then
            next_idx = 1
          elseif(jc == 1)then
            prev_idx = num(ic)
          end if

          dx = vx(jc) - vx(next_idx)
          dy = vy(jc) - vy(next_idx)
          len_d_sq = dx*dx + dy*dy

          if(len_d_sq < min_d_T1*min_d_T1)then
            ! Edge shrinking : the protocol's own edge is driven
            ! below min_d_T1 on purpose and swapped by the protocol itself
            ! (at shrink_l_min), so the automatic T1 must leave it alone.
            if (shrink_protect_edge) then
              if ((inn(jc,ic) == shrink_va .and. inn(next_idx,ic) == shrink_vb) .or. &
                  (inn(jc,ic) == shrink_vb .and. inn(next_idx,ic) == shrink_va)) cycle
            end if
            if (if_excluded_volume .and. rej_it == it) then
              on_tri = .false.
              do it_tri = 1, n_rej
                if ((inn(jc,ic) == rej_a(it_tri) .and. inn(next_idx,ic) == rej_b(it_tri)) .or. &
                    (inn(jc,ic) == rej_b(it_tri) .and. inn(next_idx,ic) == rej_a(it_tri))) on_tri = .true.
              end do
              if(on_tri)cycle
            end if
            on_tri = .false.
            do it_tri = 1, n_tri
              if(any(inn(1:3,tri_cell(it_tri)) == inn(jc,ic)) .and. &
                 any(inn(1:3,tri_cell(it_tri)) == inn(next_idx,ic))) on_tri = .true.
            end do
            if(on_tri)cycle
            len_d = sqrt(len_d_sq)
            count_T1 = count_T1 + 1
            cell_no(count_T1) = ic
            ver_no(count_T1) = jc
            ver_no_next(count_T1) = next_idx
            d_val(count_T1) = len_d
          end if

        end do

      end do

      if(count_T1>0)then
        call random_number(rr)
        chosen_index = int(rr*count_T1 + 1)

         cellNoT1 = cell_no(chosen_index)
         verNoT1 = ver_no(chosen_index)
         verNoNextT1 = ver_no_next(chosen_index)
         len_d_T1 = d_val(chosen_index)
        
      !  cellNoT1 = cell_no(1)
      !  verNoT1 = ver_no(1)
      !  verNoNextT1 = ver_no_next(1)
      !  len_d_T1 = d_val(1)

      else
        cellNoT1 = 0
        verNoT1 = 0
        verNoNextT1 = 0
        len_d_T1 = 0.0d0
        T1_pass = .false.

      end if

!      write(*,*)'ver&cell', cellNoT1, verNoT1


      
   end subroutine find_T1

   ! find_T1's forced branch: select the shrink protocol's edge
   ! (shrink_va, shrink_vb) as the one T1 candidate.
   subroutine find_T1_Shrink_Edge
     implicit none
     integer :: kk, knext
     real*8 :: xb, yb

     count_T1 = 0
     T1_pass = .false.
     do ic = 1, Nc
       do kk = 1, num(ic)
         knext = merge(1, kk + 1, kk == num(ic))
         if ((inn(kk,ic) == shrink_va .and. inn(knext,ic) == shrink_vb) .or. &
             (inn(kk,ic) == shrink_vb .and. inn(knext,ic) == shrink_va)) then
           cellNoT1 = ic
           verNoT1 = kk
           verNoNextT1 = knext
           xb = v(1,inn(knext,ic)); yb = v(2,inn(knext,ic))
           call Unwrap_Relative(v(1,inn(kk,ic)), v(2,inn(kk,ic)), xb, yb)
           len_d_T1 = sqrt((xb - v(1,inn(kk,ic)))**2 + (yb - v(2,inn(kk,ic)))**2)
           count_T1 = 1
           T1_pass = .true.
           return
         end if
       end do
     end do
   end subroutine find_T1_Shrink_Edge

   subroutine find_T1_Affected

     implicit none

     integer :: ik, jcNext
     integer :: temp_array(Nc) !(Lx*Ly)
     

    ! write(*,*)inn(verNoT1,cellNoT1), inn(verNoNextT1, cellNoT1)
     

    ik = 0
     do ic = 1, Nc !Lx*Ly
       do jc = 1, num(ic)
         if((inn(jc,ic)==inn(verNoT1,cellNoT1)).or. & 
           inn(jc,ic)==inn(verNoNextT1,cellNoT1))then
        
       ! jcNext = jc + 1
       ! if(jc == num(ic))then
       !   jcNext = 1
       ! end if

       ! if((inn(jc,ic)==inn(verNoT1,cellNoT1)).and. &
       !    (inn(jcNext, ic) == inn(verNoNextT1,cellNoT1)) .or. &
       !    inn(jc,ic)==inn(verNoNextT1,cellNoT1) .and. &
       !    (inn(jcNext,ic)==inn(verNoT1,cellNoT1))) then

           ik = ik + 1
            
           temp_array(ik) = ic

         end if

       end do

     end do


     call find_unique_sorted(temp_array, ik, Affected, Occurrences)

!     write(*,*)size(Affected)
!     T1_pass = .true.
!     if(size(Affected).ne.4)then
!       T1_pass = .false.
!     end if


   end subroutine find_T1_Affected

   subroutine T1_core
     implicit none

     real*8 :: scaling_factor
     real*8 :: OldX(2), OldY(2), NewX(2), NewY(2)
     integer :: ik,il, cellNo, verNo_in1, verNo_in2, verNo_in1_indx, verNo_in2_indx
     integer :: inn_affected(inn_dim1), inn_temp(inn_dim1)

     inn_affected = 0
     


      !write(*,*)Affected
      !write(*,*)Occurrences


     scaling_factor = 2.0d0*min_d_T1/len_d_T1*t1_len_fac   ! t1_len_fac: 1 except in the no-crossing retries (Do_T1)

!     write(*,*)cellNoT1, verNoT1, verNoNextT1, len_d_T1, scaling_factor


    OldX(1) = v(1, inn(verNoT1, cellNoT1))
    OldX(2) = v(1, inn(verNoNextT1, cellNoT1))
    OldY(1) = v(2, inn(verNoT1, cellNoT1))
    OldY(2) = v(2, inn(verNoNextT1, cellNoT1))

    ! PBC : OldX(2)/OldY(2) must be in the same local unwrapped
    ! frame as OldX(1)/OldY(1) before the midpoint/bisector below, or a
    ! genuinely short edge that happens to straddle a periodic wrap would
    ! get a nonsensical (near-box-spanning) midpoint instead.
    call Unwrap_Relative(OldX(1), OldY(1), OldX(2), OldY(2))

    call PerpendicularBisector2(OldX,OldY,scaling_factor,NewX,NewY)

    ! PBC: NewX/NewY were computed in OldX(1)/OldY(1)'s local frame, which
    ! may itself sit outside the canonical box (or the flip may push a
    ! vertex just past a periodic edge) -- re-wrap before storing.
    call Wrap_Position(NewX(1), NewY(1))
    call Wrap_Position(NewX(2), NewY(2))

    v(1, inn(verNoT1, cellNoT1)) = NewX(2)
    v(1, inn(verNoNextT1, cellNoT1)) = NewX(1)
    v(2, inn(verNoT1, cellNoT1)) = NewY(2)
    v(2, inn(verNoNextT1, cellNoT1)) = NewY(1)
        
    verNo_in1 = inn(verNoT1, cellNoT1)
    verNo_in2 = inn(verNoNextT1, cellNoT1)
    
!    write(*,*)'verNo', verNo_in1, verNo_in2

    do ik = 1,size(Affected)

    
      if(Occurrences(ik).eq.2)then   ! These will lose.
        inn_temp = 0
        inn_affected = 0

        cellNo = Affected(ik)
        inn_affected(1:num(cellNo)) = inn(1:num(cellNo), cellNo)

!        write(*,*)'inn Affected', inn_affected(1:num(cellNo))

        do il = 1, num(cellNo)
          if((inn_affected(il)).eq.verNo_in1)then
            verNo_in1_indx = il
          end if
          if((inn_affected(il)).eq.verNo_in2)then
            verNo_in2_indx = il
          end if
        end do

        if(verNo_in1_indx.eq.1.and.verNo_in2_indx.eq.num(cellNo))then
          verNo_in1_indx = num(cellNo) + 1
        elseif(verNo_in2_indx.eq.1.and.verNo_in1_indx.eq.num(cellNo))then
          verNo_in2_indx = num(cellNo) + 1
        end if

       ! write(*,*)verNo_in1_indx, verNo_in2_indx

       if(verNo_in2_indx.lt.verNo_in1_indx)then
          
         inn_temp = 0

         do il = 1, verNo_in2_indx-1
            inn_temp(il) = inn_affected(il)
         end do


         do il = verNo_in2_indx, inn_dim1-1
           inn_temp(il) = inn_affected(il+1)
         end do

         do il = 1, num(cellNo)
           inn(il, cellNo) = inn_temp(il)
         end do

         num(cellNo) = num(cellNo) - 1

       elseif(verNo_in1_indx.lt.verNo_in2_indx)then
        

         inn_temp = 0

         do il = 1, verNo_in1_indx-1
            inn_temp(il) = inn_affected(il)
         end do

        
         do il = verNo_in1_indx, inn_dim1-1
           inn_temp(il) = inn_affected(il+1)
         end do

         do il = 1, num(cellNo)
           inn(il, cellNo) = inn_temp(il)
         end do

         num(cellNo) = num(cellNo) - 1

       end if
        

!       write(*,*)'inn updated', inn(1:num(cellNo), cellNo) 






      else if(Occurrences(ik).eq.1)then  ! These will gain. 
        
       inn_temp = 0
       inn_affected = 0

       verNo_in1_indx = 0
       verNo_in2_indx = 0

        cellNo = Affected(ik)
        inn_affected(1:num(cellNo)) = inn(1:num(cellNo), cellNo)

!       write(*,*)'inn_affected1',inn_affected

        do il = 1, num(cellNo)
          if((inn_affected(il)).eq.verNo_in1)then
             verNo_in1_indx  = il
           end if
           if((inn_affected(il)).eq.verNo_in2)then
             verNo_in2_indx  = il
           end if
        end do
       ! write(*,*)'index',verNo_in1_indx, verNo_in2_indx
        
       
       ! BUGFIX : guard against overflowing inn's first dimension
       ! (inn_dim1) when a cell that is already at max vertex capacity would
       ! gain a vertex from this T1 flip. Previously there was no check here,
       ! so num(cellNo) could exceed inn_dim1 and the write-back loop below
       ! would write past inn's declared bound.
       if(verNo_in1_indx.ne.0.and.verNo_in2_indx.eq.0)then

         if(num(cellNo) .ge. inn_dim1)then
           write(*,*)'T1_core: cell', cellNo, 'already at max vertex capacity (inn_dim1) -- skipping vertex gain'
         else

         inn_temp = 0

         do il = 1, verNo_in1_indx
           inn_temp(il) = inn_affected(il)
         end do

         inn_temp(verNo_in1_indx+1) = verNo_in2

         do il = verNo_in1_indx + 2, inn_dim1
           inn_temp(il) = inn_affected(il-1)
         end do

         num(cellNo) = num(CellNo) + 1

         do il = 1, num(cellNo)
            inn(il, cellNo) = inn_temp(il)
         end do

         end if

       end if


       if(verNo_in2_indx.ne.0.and.verNo_in1_indx.eq.0)then

         if(num(cellNo) .ge. inn_dim1)then
           write(*,*)'T1_core: cell', cellNo, 'already at max vertex capacity (inn_dim1) -- skipping vertex gain'
         else

         inn_temp = 0
         do il = 1, verNo_in2_indx
           inn_temp(il) = inn_affected(il)
         end do

         inn_temp(verNo_in2_indx+1) = verNo_in1

         do il = verNo_in2_indx + 2, inn_dim1
           inn_temp(il) = inn_affected(il-1)
         end do

         num(cellNo) = num(CellNo) + 1

         do il = 1, num(cellNo)
            inn(il, cellNo) = inn_temp(il)
         end do

         end if

       end if
      
!       write(*,*)'inn_updated1',inn(:, cellNo)



     end if  ! Occurrence if
    
   end do   ! ik do
     



    end subroutine T1_core



    subroutine Do_T1


      implicit none
      integer :: il, im
      logical :: affected_computed

      ! Event-log locals (T1/T2 spatial/cell-identity tracking). Only used right before T1_core, once
      ! T1_pass has survived every rejection check above. ev_ids holds the
      ! numeric suffix of each affected cell's persistent cell_identity
      ! (CellIdNum, allocation.f90), not the string itself -- see
      ! Open_Event_Logs' comment for why.
      real*8 :: ev_x1, ev_y1, ev_x2, ev_y2, ev_mx, ev_my
      real*8 :: ev_ids(4)
      integer :: ev_ik
      ! no-crossing T1 (branch no_overlap)
      real*8, allocatable :: v_t1(:,:)
      integer, allocatable :: inn_t1(:,:), num_t1(:)
      integer :: va_t1, vb_t1, try_t1

      T1_pass = .true.
      affected_computed = .false.

      call find_T1

      if(T1_pass)then
        if(if_Fixed_boundary)then
          call find_T1_Affected
          affected_computed = .true.


       if(sum(Total_T2_count(1:it_s)).gt.0)then
         call Find_boundary_dynamic
       else
         call Get_Boundary_info
       end if



          do im = 1, size(boundary)
            do il = 1, size(Affected)
              if(boundary(im).eq.Affected(il))then
                  T1_pass = .false.
                  write(*,*)'Boundary Ignored T1'
              end if
            end do
          end do
       end if
      end if

      
      if(T1_pass .and. (if_bottom_borders_fixed .or. if_top_borders_fixed)) then
        call Find_boundary_dynamic
      end if

      ! BUGFIX : the second endpoint of the shrinking edge is
      ! verNoNextT1 (computed with wraparound in find_T1), not the raw index
      ! verNoT1+1. When the short edge is the polygon-closing edge
      ! (verNoT1 == num(cellNoT1)), inn(verNoT1+1, cellNoT1) read a stale/zero
      ! slot instead of the real second endpoint, letting a T1 flip through a
      ! fixed border undetected (and risking an out-of-bounds read on inn if
      ! the cell is at max vertex capacity).
      if(T1_pass.and.if_bottom_borders_fixed)then
      !!call Find_boundary_dynamic
        do im = 1, bottom_border_count
          if(bottom_border(im).eq.inn(verNoT1, cellNoT1) &
            .or.bottom_border(im).eq.inn(verNoNextT1, cellNoT1))then
            T1_pass = .false.
            write(*,*)'Bottom Border Ignored T1'
          end if
        end do
      end if

      if(T1_pass.and.if_top_borders_fixed)then
      !!call Find_boundary_dynamic
        do im = 1, top_border_count
          if(top_border(im).eq.inn(verNoT1, cellNoT1) &
            .or.top_border(im).eq.inn(verNoNextT1, cellNoT1))then
            T1_pass = .false.
            write(*,*)'Top Border Ignored T1'
          end if
        end do
      end if

     
!      write(*,*)T1_pass,count_T1
      
      if(T1_pass)then
          ! OPTIMIZATION : find_T1_Affected (an O(Nc*max_verts)
          ! full-mesh scan) was called unconditionally here even when the
          ! if_Fixed_boundary branch above already called it with the same
          ! inputs (cellNoT1/verNoT1/verNoNextT1 unchanged since), doubling
          ! the scan for every accepted T1. Skip it if already computed.
          if (.not. affected_computed) call find_T1_Affected

          ! BUGFIX : reject a T1 flip that would reduce any affected
          ! "losing" cell (Occurrences==2, i.e. it shares the collapsing
          ! edge) below the minimum viable polygon vertex count (3).
          ! T1_core's Occurrences==2 branch had no lower-bound guard -- a
          ! cell already at 3 vertices could be flipped down to 2, and later
          ! down to 1, at which point num(ic)==1 makes prev_idx wrap to 0 in
          ! Force_Calculation (an out-of-bounds array access). Confirmed live
          ! under an aggressive min_d_T1 stress test. A cell this small
          ! should be removed by T2, not further reduced by T1 -- reject the
          ! whole flip (no partial mutation) rather than trying to patch just
          ! this cell's update.
          do il = 1, size(Affected)
            if (Occurrences(il) == 2 .and. num(Affected(il)) <= 3) then
              T1_pass = .false.
              write(*,*)'Degenerate cell Ignored T1'
            end if
            ! a cell that would gain a vertex (Occurrences == 1) must have a free
            ! vertex slot; otherwise reject the whole flip before T1_core moves
            ! any vertex (T1_core would skip only that cell's update).
            if (Occurrences(il) == 1 .and. num(Affected(il)) >= inn_dim1) then
              T1_pass = .false.
              write(*,*)'Cell at vertex capacity Ignored T1'
            end if
          end do

          if (T1_pass) then
            ! Event log : capture the flipping edge's midpoint and
            ! the (up to 4) affected cells' PERSISTENT identities *before*
            ! T1_core overwrites the edge's vertex positions in place --
            ! afterward, v(:,inn(verNoT1/verNoNextT1,cellNoT1)) holds the
            ! POST-flip positions, not the event location. PBC-unwrapped
            ! (Unwrap_Relative/Wrap_Position), same convention as every
            ! other geometric routine in this codebase.
            ev_x1 = v(1, inn(verNoT1,     cellNoT1)); ev_y1 = v(2, inn(verNoT1,     cellNoT1))
            ev_x2 = v(1, inn(verNoNextT1, cellNoT1)); ev_y2 = v(2, inn(verNoNextT1, cellNoT1))
            call Unwrap_Relative(ev_x1, ev_y1, ev_x2, ev_y2)
            ev_mx = 0.5d0 * (ev_x1 + ev_x2)
            ev_my = 0.5d0 * (ev_y1 + ev_y2)
            call Wrap_Position(ev_mx, ev_my)

            ev_ids = 0.0d0
            do ev_ik = 1, min(4, size(Affected))
              ev_ids(ev_ik) = CellIdNum(cell_identity(Affected(ev_ik)))
            end do
            if (if_excluded_volume) then
              ! No-crossing T1 (branch no_overlap): the new edge is
              ! put on the old edge's perpendicular bisector, 2*min_d_T1 long.
              ! Next to small or distorted cells its ends can land beyond a
              ! neighbouring edge, and the flip itself then makes edges cross
              ! (most crossings left in PM with the excluded-volume force on).
              ! Keep the flip only if no edge at the two moved vertices
              ! crosses another edge; otherwise undo it and retry with a new
              ! edge half and a quarter as long, and skip it if all three cross.
              v_t1 = v; inn_t1 = inn; num_t1 = num
              va_t1 = inn(verNoT1, cellNoT1); vb_t1 = inn(verNoNextT1, cellNoT1)
              do try_t1 = 1, 3
                t1_len_fac = 0.5d0**(try_t1 - 1)
                call T1_core
                if (.not. T1_makes_crossing(va_t1, vb_t1)) exit
                v = v_t1; inn = inn_t1; num = num_t1
                if (try_t1 == 3) then
                  T1_pass = .false.
                  if (rej_it /= it) then
                    rej_it = it; n_rej = 0
                  end if
                  if (n_rej < size(rej_a)) then
                    n_rej = n_rej + 1; rej_a(n_rej) = va_t1; rej_b(n_rej) = vb_t1
                  end if
                end if
              end do
              t1_len_fac = 1.0d0
              if (T1_pass) then
                write(iunit_T1events) dble(it), ev_mx, ev_my, &
                  ev_ids(1), ev_ids(2), ev_ids(3), ev_ids(4)
                Total_T1_count(it_s) = Total_T1_count(it_s) + 1
              end if
            else
            write(iunit_T1events) dble(it), ev_mx, ev_my, &
              ev_ids(1), ev_ids(2), ev_ids(3), ev_ids(4)

            call T1_core
           ! print*, "T1 Happened at it = ", it
            Total_T1_count(it_s) = Total_T1_count(it_s) + 1
            end if
          end if
        end if



    end subroutine Do_T1







    logical function T1_makes_crossing(va, vb)
      ! No-crossing T1 (branch no_overlap): .true. if an edge at vertex
      ! va or vb crosses (properly, in its interior) an edge sharing no vertex
      ! with it. Brute force over all cell edges; T1s are rare enough.
      implicit none
      integer, intent(in) :: va, vb
      integer :: k, j, kk, jj, a, b, c, d
      T1_makes_crossing = .false.
      do k = 1, Nc
        if (num(k) <= 0) cycle
        do j = 1, num(k)
          a = inn(j, k); b = inn(merge(1, j + 1, j == num(k)), k)
          if (a /= va .and. a /= vb .and. b /= va .and. b /= vb) cycle
          do kk = 1, Nc
            if (num(kk) <= 0) cycle
            do jj = 1, num(kk)
              c = inn(jj, kk); d = inn(merge(1, jj + 1, jj == num(kk)), kk)
              if (c == a .or. c == b .or. d == a .or. d == b) cycle
              if (seg_cross(v(1:2, a), v(1:2, b), v(1:2, c), v(1:2, d))) then
                T1_makes_crossing = .true.; return
              end if
            end do
          end do
        end do
      end do
    end function T1_makes_crossing

    logical function seg_cross(p1, p2, p3, p4)
      implicit none
      real*8, intent(in) :: p1(2), p2(2), p3(2), p4(2)
      seg_cross = (orient(p1, p2, p3) * orient(p1, p2, p4) < 0.0d0) .and. &
                  (orient(p3, p4, p1) * orient(p3, p4, p2) < 0.0d0)
    end function seg_cross

    real*8 function orient(p, q, r)
      implicit none
      real*8, intent(in) :: p(2), q(2), r(2)
      real*8 :: cr
      cr = (q(1) - p(1)) * (r(2) - p(2)) - (q(2) - p(2)) * (r(1) - p(1))
      orient = 0.0d0
      if (cr > 0.0d0) orient = 1.0d0
      if (cr < 0.0d0) orient = -1.0d0
    end function orient

end module T1_swap
