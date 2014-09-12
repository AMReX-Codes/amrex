module nodal_stencil_apply_module

  use bl_constants_module
  use bl_types
  use multifab_module
  use bc_functions_module
  use impose_neumann_bcs_module
  use stencil_types_module

  implicit none

contains

  subroutine stencil_apply_1d_nodal(sg, dd, uu, mm, ng_u, ng_d, diagonalize)
    integer, intent(in) :: ng_u, ng_d
    real (kind = dp_t), intent(in   ) :: sg(   0:)
    real (kind = dp_t), intent(inout) :: dd(1-ng_d:)
    real (kind = dp_t), intent(inout) :: uu(1-ng_u:)
    integer           , intent(in   ) :: mm(    :)
    logical           , intent(in   ) :: diagonalize

    integer :: i,lo(1),nx

    nx = size(sg,dim=1) - 2
 
    dd = ZERO

    lo = 1
    call impose_neumann_bcs_1d(uu,mm,lo,ng_u)
   
    i = 1
    if (.not. bc_dirichlet(mm(i),1,0)) then
      dd(i) = sg(i) * (uu(i+1) - uu(i)) + sg(i-1) * (uu(i-1) - uu(i))
      if (diagonalize) dd(i) = -dd(i) / (sg(i)+sg(i-1))
    end if

    do i = 2,nx-1
      dd(i) = sg(i) * (uu(i+1) - uu(i)) + sg(i-1) * (uu(i-1) - uu(i))
      if (diagonalize) dd(i) = -dd(i) / (sg(i)+sg(i-1))
    end do

    i = nx
    if (.not. bc_dirichlet(mm(i),1,0)) then
      dd(i) = sg(i) * (uu(i+1) - uu(i)) + sg(i-1) * (uu(i-1) - uu(i))
      if (diagonalize) dd(i) = -dd(i) / (sg(i)+sg(i-1))
    end if

  end subroutine stencil_apply_1d_nodal

  subroutine stencil_apply_2d_nodal(sg, dd, uu, mm, ng_u, ng_d, stencil_type, diagonalize)
    integer, intent(in) :: ng_u, ng_d
    real (kind = dp_t), intent(inout) :: uu(1-ng_u:,1-ng_u:)
    real (kind = dp_t), intent(inout) :: dd(1-ng_d:,1-ng_d:)
    real (kind = dp_t), intent(in   ) :: sg(   0:,   0:)
    integer           , intent(in   ) :: mm(    :,    :)
    integer           , intent(in   ) :: stencil_type 
    logical           , intent(in   ) :: diagonalize

    integer :: i,j,lo(2),nx,ny
    logical :: zeroit,iface,jface
    real (kind = dp_t) :: ss0

    lo = 1
    call impose_neumann_bcs_2d(uu,mm,lo,ng_u)
 
    nx = size(sg,dim=1) - 1
    ny = size(sg,dim=2) - 1

    if (stencil_type .eq. ND_CROSS_STENCIL) then

       do j = 1,ny
          jface = .false. ; if ( (j.eq.1).or.(j.eq.ny) ) jface = .true.
          do i = 1,nx
             iface = .false. ; if ( (i.eq.1).or.(i.eq.nx) ) iface = .true.

             zeroit = .false.

             if ( iface .or. jface ) then
                if (bc_dirichlet(mm(i,j),1,0)) zeroit = .true.
             end if

             if (zeroit) then
                dd(i,j) = ZERO
             else
                ss0 = -(sg(i-1,j-1)+sg(i,j-1)+sg(i-1,j)+sg(i,j))
                dd(i,j) = ss0 * uu(i,j) + &
                   HALF * ( (sg(i  ,j-1)+sg(i  ,j  )) * uu(i+1,j) + &
                            (sg(i-1,j-1)+sg(i-1,j  )) * uu(i-1,j) + &
                            (sg(i-1,j  )+sg(i  ,j  )) * uu(i,j+1) + &
                            (sg(i-1,j-1)+sg(i  ,j-1)) * uu(i,j-1) ) 
                if (diagonalize) then
                    dd(i,j) = dd(i,j) / ss0
                end if
             end if
          end do
       end do

    else if (stencil_type .eq. ND_DENSE_STENCIL) then

       do j = 1,ny
          jface = .false. ; if ( (j.eq.1).or.(j.eq.ny) ) jface = .true.
          do i = 1,nx
             iface = .false. ; if ( (i.eq.1).or.(i.eq.nx) ) iface = .true.

             zeroit = .false.

             if ( iface .or. jface ) then
                if (bc_dirichlet(mm(i,j),1,0)) zeroit = .true.
             end if

             if (zeroit) then
                dd(i,j) = ZERO
             else
                ss0 = -TWO * THIRD * (sg(i-1,j-1) + sg(i,j-1) + sg(i-1,j) + sg(i,j))
                dd(i,j) = ss0 * uu(i,j) + THIRD * ( &
                            sg(i-1,j-1) * uu(i-1,j-1) + &
                            sg(i  ,j-1) * uu(i+1,j-1) + &
                            sg(i-1,j  ) * uu(i-1,j+1) + &
                            sg(i  ,j  ) * uu(i+1,j+1) + &
                  HALF * ( (sg(i-1,j-1) + sg(i  ,j-1)) * uu(i,j-1) + &
                           (sg(i-1,j-1) + sg(i-1,j  )) * uu(i-1,j) + &
                           (sg(i  ,j-1) + sg(i  ,j  )) * uu(i+1,j) + &
                           (sg(i-1,j  ) + sg(i  ,j  )) * uu(i,j+1) ) )
                if (diagonalize) then
                    dd(i,j) = dd(i,j) / ss0
                end if
             end if
          end do
       end do

    else if (stencil_type .eq. ND_VATER_STENCIL) then

       do j = 1,ny
          jface = .false. ; if ( (j.eq.1).or.(j.eq.ny) ) jface = .true.
          do i = 1,nx
             iface = .false. ; if ( (i.eq.1).or.(i.eq.nx) ) iface = .true.

             zeroit = .false.

             if ( iface .or. jface ) then
                if (bc_dirichlet(mm(i,j),1,0)) zeroit = .true.
             end if

             if (zeroit) then
                dd(i,j) = ZERO
             else
                ss0 = -THREE4TH * (sg(i-1,j-1) + sg(i,j-1) + sg(i-1,j) + sg(i,j))
                dd(i,j) = ss0 * uu(i,j) + &
                FOURTH * ( &
                            sg(i-1,j-1) * uu(i-1,j-1) + &
                            sg(i  ,j-1) * uu(i+1,j-1) + &
                            sg(i-1,j  ) * uu(i-1,j+1) + &
                            sg(i  ,j  ) * uu(i+1,j+1) + &
                           (sg(i-1,j-1) + sg(i  ,j-1)) * uu(i,j-1) + &
                           (sg(i-1,j-1) + sg(i-1,j  )) * uu(i-1,j) + &
                           (sg(i  ,j-1) + sg(i  ,j  )) * uu(i+1,j) + &
                           (sg(i-1,j  ) + sg(i  ,j  )) * uu(i,j+1) )
                if (diagonalize) then
                    dd(i,j) = dd(i,j) / ss0
                end if
             end if
          end do
       end do

    else
       call bl_error("stencil_apply_2d_nodal: dont know this stencil_type")
    end if

  end subroutine stencil_apply_2d_nodal

  subroutine stencil_apply_3d_nodal(sg, dd, uu, mm, ng_u, ng_d, stencil_type, &
                                    uniform_dh, bottom_solver, diagonalize)
    integer           , intent(in   ) :: ng_u, ng_d
    real (kind = dp_t), intent(inout) :: uu(1-ng_u:,1-ng_u:,1-ng_u:)
    real (kind = dp_t), intent(inout) :: dd(1-ng_d:,1-ng_d:,1-ng_d:)
    real (kind = dp_t), intent(in   ) :: sg(   0:,   0:,   0:)
    integer           , intent(in   ) :: mm(    :,    :,    :)
    integer           , intent(in   ) :: stencil_type
    logical           , intent(in   ) :: uniform_dh, bottom_solver, diagonalize

    integer :: i,j,k,lo(3),nx,ny,nz
    logical :: jface,kface
    real (kind = dp_t) :: f0, fx, fy, fz, fxyz, f2y2zx, f2x2zy, f2x2yz, ss0

    lo = 1
    call impose_neumann_bcs_3d(uu,mm,lo,ng_u)

    nz = size(sg,dim=3) - 1
    ny = size(sg,dim=2) - 1
    nx = size(sg,dim=1) - 1

    if (stencil_type .eq. ND_CROSS_STENCIL) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,jface,kface,ss0) IF(.not.bottom_solver)
       do k = 1,nz
          kface = .false. ; if ( (k.eq.1) .or. (k.eq.nz) ) kface = .true.

          do j = 1,ny
             jface = .false. ; if ( (j.eq.1) .or. (j.eq.ny) ) jface = .true.

             do i = 1,nx

                if ( (jface .or. kface .or. (i.eq.1) .or. (i.eq.nx)) .and. bc_dirichlet(mm(i,j,k),1,0) ) then
                   dd(i,j,k) = ZERO
                else
                   ss0 = -( sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1) &
                           +sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1) &
                           +sg(i-1,j-1,k  ) + sg(i-1,j  ,k  ) &
                           +sg(i  ,j-1,k  ) + sg(i  ,j  ,k  )) * THREE
                   dd(i,j,k) = ss0 * uu(i,j,k) + &
                     (sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                     +sg(i-1,j  ,k-1) + sg(i-1,j  ,k  )) * uu(i-1,j  ,k  ) + &
                     (sg(i  ,j-1,k-1) + sg(i  ,j-1,k  ) &
                     +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  )) * uu(i+1,j  ,k  ) + &
                     (sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                     +sg(i  ,j-1,k-1) + sg(i  ,j-1,k  )) * uu(i  ,j-1,k  ) + &
                     (sg(i-1,j  ,k-1) + sg(i-1,j  ,k  ) &
                     +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  )) * uu(i  ,j+1,k  ) + &
                     (sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1) &
                     +sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1)) * uu(i  ,j  ,k-1) + &
                     (sg(i-1,j-1,k  ) + sg(i-1,j  ,k  ) &
                     +sg(i  ,j-1,k  ) + sg(i  ,j  ,k  )) * uu(i  ,j  ,k+1) 
                   if (diagonalize) then
                       dd(i,j,k) = dd(i,j,k) / ss0
                   end if
                end if
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else if (stencil_type .eq. ND_DENSE_STENCIL) then

       fx     = ONE/36._dp_t
       fy     = fx
       fz     = fx
       f0     = FOUR * (fx + fy + fz)
       fxyz   = (fx+fy+fz)
       f2y2zx = (TWO*fy+TWO*fz-fx)
       f2x2zy = (TWO*fx+TWO*fz-fy)
       f2x2yz = (TWO*fx+TWO*fy-fz)

       !$OMP PARALLEL DO PRIVATE(i,j,k,jface,kface,ss0) IF(.not.bottom_solver)
       do k = 1,nz
          kface = .false. ; if ( (k.eq.1) .or. (k.eq.nz) ) kface = .true.

          do j = 1,ny
             jface = .false. ; if ( (j.eq.1) .or. (j.eq.ny) ) jface = .true.

             do i = 1,nx

                if ( (jface .or. kface .or. (i.eq.1) .or. (i.eq.nx)) .and. bc_dirichlet(mm(i,j,k),1,0) ) then
                   dd(i,j,k) = ZERO
                else

                ss0 =  -( sg(i-1,j-1,k-1) + sg(i,j-1,k-1) &
                         +sg(i-1,j  ,k-1) + sg(i,j  ,k-1) &
                         +sg(i-1,j-1,k  ) + sg(i,j-1,k  ) &
                         +sg(i-1,j  ,k  ) + sg(i,j  ,k  ) ) * f0 
                dd(i,j,k) = fxyz * ( &   ! Corners
                     sg(i-1,j-1,k-1) * uu(i-1,j-1,k-1) + sg(i  ,j-1,k-1) * uu(i+1,j-1,k-1) + &
                     sg(i-1,j  ,k-1) * uu(i-1,j+1,k-1) + sg(i  ,j  ,k-1) * uu(i+1,j+1,k-1) + &
                     sg(i-1,j-1,k  ) * uu(i-1,j-1,k+1) + sg(i  ,j-1,k  ) * uu(i+1,j-1,k+1) + &
                     sg(i-1,j  ,k  ) * uu(i-1,j+1,k+1) + sg(i  ,j  ,k  ) * uu(i+1,j+1,k+1)) &
                     + f2y2zx * ( & ! Edges in x-direction
                     (sg(i  ,j-1,k-1) + sg(i-1,j-1,k-1)) * uu(i  ,j-1,k-1) + &
                     (sg(i  ,j  ,k-1) + sg(i-1,j  ,k-1)) * uu(i  ,j+1,k-1) + &
                     (sg(i  ,j-1,k  ) + sg(i-1,j-1,k  )) * uu(i  ,j-1,k+1) + &
                     (sg(i  ,j  ,k  ) + sg(i-1,j  ,k  )) * uu(i  ,j+1,k+1)) &
                     + f2x2zy * ( & ! Edges in y-direction
                     (sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1)) * uu(i-1,j  ,k-1) + &
                     (sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1)) * uu(i+1,j  ,k-1) + &
                     (sg(i-1,j-1,k  ) + sg(i-1,j  ,k  )) * uu(i-1,j  ,k+1) + &
                     (sg(i  ,j-1,k  ) + sg(i  ,j  ,k  )) * uu(i+1,j  ,k+1)) &
                     + f2x2yz * ( & ! Edges in z-direction
                     (sg(i-1,j-1,k-1) + sg(i-1,j-1,k  )) * uu(i-1,j-1,k  ) + &
                     (sg(i  ,j-1,k-1) + sg(i  ,j-1,k  )) * uu(i+1,j-1,k  ) + &
                     (sg(i-1,j  ,k-1) + sg(i-1,j  ,k  )) * uu(i-1,j+1,k  ) + &
                     (sg(i  ,j  ,k-1) + sg(i  ,j  ,k  )) * uu(i+1,j+1,k  )) & 
                     + ss0 * uu(i,j,k)

                   if (.not. uniform_dh) then
                      !
                      ! Add faces (only non-zero for non-uniform dx)
                      !
                      dd(i,j,k) = dd(i,j,k) + &
                           ( (FOUR*fx-TWO*fy-TWO*fz)*(sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                            +sg(i-1,j  ,k-1) + sg(i-1,j  ,k  )) ) * uu(i-1,j  ,k  ) + &
                           ( (FOUR*fx-TWO*fy-TWO*fz)*(sg(i  ,j-1,k-1) + sg(i  ,j-1,k  ) &
                            +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  )) ) * uu(i+1,j  ,k  ) + &
                           ( (FOUR*fy-TWO*fx-TWO*fz)*(sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
                            +sg(i  ,j-1,k-1) + sg(i  ,j-1,k  )) ) * uu(i  ,j-1,k  ) + &
                           ( (FOUR*fy-TWO*fx-TWO*fz)*(sg(i-1,j  ,k-1) + sg(i-1,j  ,k  ) &
                            +sg(i  ,j  ,k-1) + sg(i  ,j  ,k  )) ) * uu(i  ,j+1,k  ) + &
                           ( (FOUR*fz-TWO*fx-TWO*fy)*(sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1) &
                            +sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1)) ) * uu(i  ,j  ,k-1) + &
                           ( (FOUR*fz-TWO*fx-TWO*fy)*(sg(i-1,j-1,k  ) + sg(i-1,j  ,k  ) &
                            +sg(i  ,j-1,k  ) + sg(i  ,j  ,k  )) ) * uu(i  ,j  ,k+1)
                   end if

                   ! This accounts for the fact that fac = 1/(4*dx*dx) to be compatible with 
                   !      the cross stencil
                   if (diagonalize) then
                       dd(i,j,k) = dd(i,j,k) / ss0
                   else
                       dd(i,j,k) = FOUR * dd(i,j,k)
                   end if

                end if

             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else if (stencil_type .eq. ND_VATER_STENCIL) then
       call bl_error("stencil_apply_3d_nodal: ND_VATER_STENCIL not implemented in 3-d")
    else
       call bl_error("stencil_apply_3d_nodal: dont know this stencil_type")
    end if

  end subroutine stencil_apply_3d_nodal

end module nodal_stencil_apply_module
