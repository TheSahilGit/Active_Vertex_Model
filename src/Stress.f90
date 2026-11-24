module Stress

  use allocation 
  use System_Info


  contains

   subroutine Calculate_StressTensor
 
 
       implicit none
       integer :: cellNo
       real*8 :: term1, term2, dx, dy, len_d
       ! Differential line tension (if_free_edge_adhesion): see
       ! Force.f90 for the same pattern. term2 above used to be one constant
       ! per cell; term2_edge below re-derives it per edge, splitting out the
       ! gamm-dependent piece (beta_term2 stays the same for every edge of a
       ! given cell -- only the line-tension part varies by edge).
       logical :: is_free_edge(inn_dim1, num_dim)
       real*8 :: beta_term2, gamm_this_edge

     !  call Get_Boundary_info

       radius_from_core = sqrt(dble(Lx*Lx) + dble(Ly*Ly))/5.0d0
       !radius_from_core = 4.0d0 * sqrt(Ao)

       call Get_Cells_Within_Radius

       if (if_free_edge_adhesion) call Get_Free_Edges(is_free_edge)

!       open(1211, file='inside.dat', status='unknown')
!         write(1211,*)n_inside, inside_cells(1:n_inside)
!       close(1211)
 
       ! BUGFIX : beta/gamm nondimensionalization moved to read_input
       ! (allocation.f90), done once instead of every call to this subroutine
       ! -- for why repeating it here compounded beta/gamm.

       totalarea = 0.0d0
       TotalSigma = 0.0d0

 
       do ic = 1,n_inside
          cellNo = inside_cells(ic)
          nn = num(cellNo)

          call Gather_Cell_Vertices_PBC(inn(1:nn,cellNo), nn, vx, vy)

          call CalculateArea(vx,vy,nn,area)
          area = abs(area)
          totalarea = totalarea + area
       end do
 
 
       do ic = 1, n_inside
         cellNo = inside_cells(ic)
         nn = num(cellNo)

         call Gather_Cell_Vertices_PBC(inn(1:nn,cellNo), nn, vx, vy)

         call CalculateArea(vx,vy,nn,area)
         area = abs(area)
         call CalculatePerimeter(vx,vy,nn,perimeter)
 
         term1 =  2.0d0 * lambda * (area - Ao)
         ! each cell's own beta (beta_cell, set in Force_Calculation), not the
         ! module scalar beta, which holds the last cell's value under
         ! if_RhoROCK / if_active_contractility (stress_relaxation_LE copy)
         term2 = (2.0d0 * beta_cell(cellNo)* perimeter + gamm)/(2.0d0*area)
         beta_term2 = (2.0d0 * beta_cell(cellNo) * perimeter) / (2.0d0 * area)

         ! BUGFIX : sigma is a global (allocation.f90) that was never
         ! reset per cell, so each cell's edge-sum accumulated on top of the
         ! *previous cell's already fully-scaled* tensor -- corrupting
         ! TotalSigma/ShearStress for every cell after the first, every call.
         sigma = 0.0d0

         ! Differential line tension : kept as an explicit branch
         ! (not one formula with gamm_this_edge=gamm when off) so the
         ! off-case is the ORIGINAL "term2 * (summed unscaled edge tensor)"
         ! expression, bit-for-bit -- term2*(a+b+c) is not guaranteed
         ! identical to term2*a+term2*b+term2*c in floating point, which
         ! this codebase's byte-identical-regression testing convention
         ! depends on.
         if (if_free_edge_adhesion) then
           do jc = 1,nn
             jp = jc+1
             if (jc.eq.nn)then
               jp = 1
             end if
             dx = vx(jp) - vx(jc)
             dy = vy(jp) - vy(jc)
             len_d = sqrt(dx*dx + dy*dy)

             gamm_this_edge = merge(gamm_free, gamm, is_free_edge(jc, cellNo))

             sigma(1,1) = sigma(1,1) + (beta_term2 + gamm_this_edge/(2.0d0*area)) * dx*dx/len_d
             sigma(1,2) = sigma(1,2) + (beta_term2 + gamm_this_edge/(2.0d0*area)) * dx*dy/len_d
             sigma(2,1) = sigma(2,1) + (beta_term2 + gamm_this_edge/(2.0d0*area)) * dy*dx/len_d
             sigma(2,2) = sigma(2,2) + (beta_term2 + gamm_this_edge/(2.0d0*area)) * dy*dy/len_d
           end do

           sigma(1,1) = term1  + sigma(1,1)
           sigma(1,2) = sigma(1,2)
           sigma(2,1) = sigma(2,1)
           sigma(2,2) = term1 + sigma(2,2)
         else
           do jc = 1,nn
             jp = jc+1
             if (jc.eq.nn)then
               jp = 1
             end if
             dx = vx(jp) - vx(jc)
             dy = vy(jp) - vy(jc)
             len_d = sqrt(dx*dx + dy*dy)

             sigma(1,1) = sigma(1,1) + dx*dx/len_d
             sigma(1,2) = sigma(1,2) + dx*dy/len_d
             sigma(2,1) = sigma(2,1) + dy*dx/len_d
             sigma(2,2) = sigma(2,2) + dy*dy/len_d
           end do

           sigma(1,1) = term1  + term2 * sigma(1,1)
           sigma(1,2) = term2 * sigma(1,2)
           sigma(2,1) = term2 * sigma(2,1)
           sigma(2,2) = term1 + term2 * sigma(2,2)
         end if
 
 
 
 
        TotalSigma = TotalSigma + sigma * area/totalarea
 
 
        end do
 
 
        ShearStress(it_s) = TotalSigma(1,2)
 
    end subroutine Calculate_StressTensor



    subroutine ShearTissue

    implicit none

    real*8 :: comb

     ! write(*,*)if_Shearing, shearStrength

     if(if_Sudden_Shearing)then
       if(it.eq.sudden_shearWhen)then
           
          v(1,:) = v(1,:) + sudden_shearStrength * v(2,:)
          v(2,:) = v(2,:)
          call Advance_LE_Strain(sudden_shearStrength)

          if(if_bottom_borders_fixed)then
            call Find_boundary_dynamic
            v(1, bottom_border(1:bottom_border_count)) = &
              v(1,bottom_border(1:bottom_border_count)) - &
              sudden_shearStrength * v(2,bottom_border(1:bottom_border_count))
          end if

       end if
     end if 

     if(if_Oscillatory_Shearing)then
       if(it.gt.Oscl_shearWhen)then


         strainRate = Oscl_shearStrength * Oscl_freq_wo* cos(Oscl_freq_wo * it * dt)

         v(1,:) = v(1,:) + dt * strainRate * v(2,:)      
         call Advance_LE_Strain(dt * strainRate)
    
       end if

     end if

     ShearStrain(it_s) = LE_strain   ! value at the end of the series block (series_stride)




    end subroutine ShearTissue



    ! Lees-Edwards : every affine increment x += dgamma*y applied
    ! to the vertices must shift the image lattice by the same strain, or
    ! the tissue tears at the y-wrap (a cell straddling it would see its
    ! top vertices moved by ~dgamma*Ly relative to its bottom ones). Under
    ! PBC this keeps the affine move consistent in every image frame;
    ! without PBC LE_offset is never read, and LE_strain just records the
    ! applied strain for data/ShearStrain.dat.
    subroutine Advance_LE_Strain(dgamma)

    implicit none
    real*8, intent(in) :: dgamma

     LE_strain = LE_strain + dgamma
     LE_offset = modulo(LE_offset + dgamma * Ly_box, Lx_box)

    end subroutine Advance_LE_Strain


    ! Variable periodic box (if_variable_box): one overdamped step of
    ! the isotropic box scale s. With every vertex position and the box scaled
    ! by s (x -> s*x, Lx_box -> s*Lx_box, LE_offset -> s*LE_offset), the box
    ! moves under the generalized force of everything the vertices feel:
    !
    !   d(ln s) = [ dt*( -dE/dlns - 2*box_pressure*A_box + V(f_ABP) )
    !               + sqrt(dt)*V(f_ran + f_Polar) ] / Gam
    !
    !  - dE/dlns = sum_c [4*lambda*(A_c-Ao)*A_c + (2*beta_c*(P_c-Co)+gamm)*P_c]
    !    (area ~ s^2, perimeter ~ s; beta_c = beta_cell of this step, so FC and
    !    RhoROCK enter here), the conservative virial.
    !  - V(f) = sum_c sum_{v in c} w_v f_v . (r_v - R_c): the virial of the
    !    active vertex forces, each taken relative to the centroid R_c of the
    !    cells it acts through (periodic-safe; a force pushing a vertex away
    !    from its cells' centres expands the box). w_v = 1/(number of cells
    !    sharing v). The forces enter with the same dt / sqrt(dt) factors as
    !    in the vertex update (vertexmain.f90): AM noise and PM kicks are
    !    sqrt(dt)-scaled, ABP is dt-scaled.
    !  - Gam = sum_c sum_{v in c} w_v eta_v |r_v - R_c|^2: the vertex friction
    !    of dilating every cell about its centre -- the box responds with the
    !    friction the vertices themselves feel, no free parameter.
    !    box_mobility > 0 replaces 1/Gam by box_mobility/Nc (a fixed box
    !    mobility; used for testing).
    ! Isotropic scaling leaves the Lees-Edwards strain LE_offset/Ly_box
    ! unchanged and keeps wrapped positions in the canonical box. All vertex
    ! slots are scaled (unused pool slots too -- their positions are never
    ! read). Called after the step's vertex move, with that step's forces.
    ! box_aniso : the same with Lx_box and Ly_box relaxing
    ! separately, d ln s_a = [dt*(-dE/dln s_a - box_pressure*A_box + V_aa(f_ABP))
    ! + sqrt(dt)*V_aa(f_ran + f_Polar)] / Gam_a (a = x, y), with
    ! dE/dln s_x = sum_c [2*lambda*(A_c-Ao)*A_c + T_c * sum_e l_x^2/l],
    ! V_xx(f) = sum_c sum_v w_v f_x (x_v - X_c), Gam_x = sum_c sum_v w_v eta_v
    ! (x_v - X_c)^2 (the x and y parts sum to the isotropic expressions).
    ! The Lees-Edwards strain is held: the vertices move by
    ! x' = s_x x + LE_strain (s_y - s_x) y, y' = s_y y, which maps the image
    ! vector (LE_offset, Ly_box) onto (LE_strain Ly_box', Ly_box') mod Lx_box'.
    subroutine Variable_Box_Step

     implicit none
     integer :: cellNo, jv, kv
     real*8 :: dEdlns, s, a_c, p_c, cx, cy, rx, ry, w, Gam, Vdet, Vnoise, Q
     real*8 :: dEx, dEy, Gx, Gy, Vdx, Vdy, Vnx, Vny, sx, sy, t_c, lxx, lyy
     integer, allocatable :: ncell_at(:)

     allocate(ncell_at(size(v, 2)))
     ncell_at = 0
     do cellNo = 1, Nc
       do jv = 1, num(cellNo)
         ncell_at(inn(jv,cellNo)) = ncell_at(inn(jv,cellNo)) + 1
       end do
     end do

     dEdlns = 0.0d0
     Gam = 0.0d0
     Vdet = 0.0d0
     Vnoise = 0.0d0
     dEx = 0.0d0; dEy = 0.0d0; Gx = 0.0d0; Gy = 0.0d0
     Vdx = 0.0d0; Vdy = 0.0d0; Vnx = 0.0d0; Vny = 0.0d0
     do cellNo = 1, Nc
       nn = num(cellNo)
       call Gather_Cell_Vertices_PBC(inn(1:nn,cellNo), nn, vx, vy)
       call CalculateArea(vx,vy,nn,a_c)
       a_c = abs(a_c)
       call CalculatePerimeter(vx,vy,nn,p_c)
       dEdlns = dEdlns + 4.0d0 * lambda * (a_c - Ao) * a_c &
                       + (2.0d0 * beta_cell(cellNo) * (p_c - Co) + gamm) * p_c
       if (box_aniso) then
         t_c = 2.0d0 * beta_cell(cellNo) * (p_c - Co) + gamm
         lxx = 0.0d0; lyy = 0.0d0
         do jv = 1, nn
           rx = vx(mod(jv, nn) + 1) - vx(jv); ry = vy(mod(jv, nn) + 1) - vy(jv)
           lxx = lxx + rx*rx / sqrt(rx*rx + ry*ry); lyy = lyy + ry*ry / sqrt(rx*rx + ry*ry)
         end do
         dEx = dEx + 2.0d0 * lambda * (a_c - Ao) * a_c + t_c * lxx
         dEy = dEy + 2.0d0 * lambda * (a_c - Ao) * a_c + t_c * lyy
       end if
       cx = sum(vx) / dble(nn)
       cy = sum(vy) / dble(nn)
       do jv = 1, nn
         kv = inn(jv,cellNo)
         w = 1.0d0 / dble(ncell_at(kv))
         rx = vx(jv) - cx
         ry = vy(jv) - cy
         Gam = Gam + w * eta(kv) * (rx*rx + ry*ry)
         Vdet = Vdet + w * (fxx_ABP(kv)*rx + fyy_ABP(kv)*ry)
         Vnoise = Vnoise + w * ((fxx_ran(kv) + fxx_Polar(kv))*rx + (fyy_ran(kv) + fyy_Polar(kv))*ry)
         Gx = Gx + w * eta(kv) * rx*rx; Gy = Gy + w * eta(kv) * ry*ry
         Vdx = Vdx + w * fxx_ABP(kv)*rx; Vdy = Vdy + w * fyy_ABP(kv)*ry
         Vnx = Vnx + w * (fxx_ran(kv) + fxx_Polar(kv))*rx; Vny = Vny + w * (fyy_ran(kv) + fyy_Polar(kv))*ry
       end do
     end do
     deallocate(ncell_at)
     box_dEdlns = dEdlns

     if (box_aniso) then
       Q = dt * (-dEx - box_pressure * Lx_box * Ly_box + Vdx) + sqrt(dt) * Vnx
       if (box_mobility .gt. 0.0d0) then
         sx = exp(box_mobility * Q / dble(Nc))
       else
         sx = exp(Q / Gx)
       end if
       Q = dt * (-dEy - box_pressure * Lx_box * Ly_box + Vdy) + sqrt(dt) * Vny
       if (box_mobility .gt. 0.0d0) then
         sy = exp(box_mobility * Q / dble(Nc))
       else
         sy = exp(Q / Gy)
       end if
       v(1,:) = sx * v(1,:) + LE_strain * (sy - sx) * v(2,:)
       v(2,:) = sy * v(2,:)
       Lx_box = sx * Lx_box
       Ly_box = sy * Ly_box
       LE_offset = modulo(LE_strain * Ly_box, Lx_box)
       return
     end if

     Q = dt * (-dEdlns - 2.0d0 * box_pressure * Lx_box * Ly_box + Vdet) + sqrt(dt) * Vnoise
     if (box_mobility .gt. 0.0d0) then
       s = exp(box_mobility * Q / dble(Nc))
     else
       s = exp(Q / Gam)
     end if

     v(1,:) = s * v(1,:)
     v(2,:) = s * v(2,:)
     Lx_box = s * Lx_box
     Ly_box = s * Ly_box
     LE_offset = s * LE_offset

    end subroutine Variable_Box_Step





end module Stress
