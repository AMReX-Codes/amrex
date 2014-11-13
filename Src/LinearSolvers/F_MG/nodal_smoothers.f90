module nodal_smoothers_module

  use bl_constants_module
  use bc_functions_module
  use nodal_stencil_module

  implicit none

contains

  subroutine nodal_line_solve_1d(sg, uu, ff, mm, lo, ng)

    use tridiag_module, only: tridiag

    integer        , intent(in   ) :: lo(:),ng
    real(kind=dp_t), intent(in   ) :: ff(lo(1)- 1:)
    real(kind=dp_t), intent(inout) :: uu(lo(1)-ng:)
    real(kind=dp_t), intent(in   ) :: sg(lo(1)- 1:)
    integer        , intent(in   ) :: mm(lo(1)   :)

!   real(kind=dp_t)              :: dd
    real(kind=dp_t), allocatable :: a_ls(:), b_ls(:), c_ls(:), r_ls(:), u_ls(:)
    integer                      :: is, ie, ilen, i, hi(size(lo))

    hi = ubound(uu)-ng

    if (.not. bc_dirichlet(mm(lo(1)),1,0)) then
       is = lo(1)
    else
       is = lo(1)+1
    end if
    if (.not. bc_dirichlet(mm(hi(1)),1,0)) then
       ie = hi(1)
    else
       ie = hi(1)-1
    end if

    ilen = ie-is+1

    allocate(a_ls(0:ilen-1))
    allocate(b_ls(0:ilen-1))
    allocate(c_ls(0:ilen-1))
    allocate(r_ls(0:ilen-1))
    allocate(u_ls(0:ilen-1))

    do i = is,ie
      a_ls(i-is) = sg(i-1)
      b_ls(i-is) = -(sg(i)+sg(i-1))
      c_ls(i-is) = sg(i  )
      r_ls(i-is) = ff(i)
    end do

    ! Adjust low end for Neumann boundary condition
    if (bc_neumann(mm(is),1,-1)) then
      c_ls(0) = TWO*c_ls(0)
    end if

    ! Adjust high end for Neumann boundary condition
    if (bc_neumann(mm(ie),1,1)) then
      a_ls(0) = TWO*a_ls(0)
    end if

    r_ls(0)      = r_ls(0)      - sg(is-1) * uu(is-1)
    r_ls(ilen-1) = r_ls(ilen-1) - sg(ie  ) * uu(ie+1)

    call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,ilen)
 
    do i = is, ie
       uu(i) = u_ls(i-is)
    end do

    deallocate(a_ls,b_ls,c_ls,r_ls,u_ls)

  end subroutine nodal_line_solve_1d

  subroutine nodal_smoother_1d(sg, uu, ff, mm, lo, ng, red_black)
    integer, intent(in) :: lo(:)
    integer, intent(in) :: ng, red_black

    real (kind = dp_t), intent(in)    :: ff(lo(1)- 1:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:)
    real (kind = dp_t), intent(in)    :: sg(lo(1)- 1:)
    integer            ,intent(in)    :: mm(lo(1)   :)
    real (kind = dp_t) :: dd,ss0
    integer :: i, hi(size(lo))

    real (kind = dp_t), parameter :: omega = 1.33_dp_t

    hi = ubound(uu)-ng

    ! Red/black
    do i = lo(1)+red_black,hi(1),2
       if ( .not. bc_dirichlet(mm(i),1,0) ) then
          ss0 = -(sg(i) + sg(i-1))
          dd = sg(i) * uu(i+1) + sg(i-1) * uu(i-1) + ss0 * uu(i)
          uu(i) = uu(i) + (omega/ss0) * (ff(i) - dd)
       end if
    end do

  end subroutine nodal_smoother_1d

  subroutine nodal_smoother_2d(sg, uu, ff, mm, lo, ng, stencil_type, red_black)

    use impose_neumann_bcs_module

    integer, intent(in) :: ng
    integer, intent(in) :: lo(:)
    real (kind = dp_t), intent(in   ) :: ff(lo(1)- 1:,lo(2)- 1:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:)
    real (kind = dp_t), intent(in   ) :: sg(lo(1)- 1:,lo(2)- 1:)
    integer            ,intent(in   ) :: mm(lo(1)   :,lo(2)   :)
    integer            ,intent(in   ) :: stencil_type
    integer            ,intent(in   ) :: red_black

    integer            :: j, i, ipar, hi(size(lo))
    real (kind = dp_t) :: dd, ss0

    hi = ubound(uu)-ng

    call impose_neumann_bcs_2d(uu,mm,lo,ng)

    if ( stencil_type .eq. ND_DENSE_STENCIL ) then

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if ( .not. bc_dirichlet(mm(i,j),1,0) ) then
                ss0 = -TWO3RD * (  sg(i-1,j-1) + sg(i,j-1) + sg(i-1,j) + sg(i,j) )
                dd =  THIRD * ( & 
                          sg(i-1,j-1) * uu(i-1,j-1) + &
                          sg(i  ,j-1) * uu(i+1,j-1) + &
                          sg(i-1,j  ) * uu(i-1,j+1) + &
                          sg(i  ,j  ) * uu(i+1,j+1) + &
                HALF * ( (sg(i-1,j-1) + sg(i  ,j-1)) * uu(i,j-1) + &
                         (sg(i-1,j-1) + sg(i-1,j  )) * uu(i-1,j) + &
                         (sg(i  ,j-1) + sg(i  ,j  )) * uu(i+1,j) + &
                         (sg(i-1,j  ) + sg(i  ,j  )) * uu(i,j+1) ) ) + &
                          ss0                        * uu(i,j  )
                uu(i,j) = uu(i,j) + (one/ss0) * (ff(i,j) - dd)
             end if
          end do
       end do

    else if ( stencil_type .eq. ND_VATER_STENCIL ) then

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if ( .not. bc_dirichlet(mm(i,j),1,0) ) then
                ss0 = -THREE4TH * (  sg(i-1,j-1) + sg(i,j-1) + sg(i-1,j) + sg(i,j) )
                dd =  ss0 * uu(i,j) + &
                      FOURTH * ( & 
                          sg(i-1,j-1) * uu(i-1,j-1) + &
                          sg(i  ,j-1) * uu(i+1,j-1) + &
                          sg(i-1,j  ) * uu(i-1,j+1) + &
                          sg(i  ,j  ) * uu(i+1,j+1) + &
                         (sg(i-1,j-1) + sg(i  ,j-1)) * uu(i,j-1) + &
                         (sg(i-1,j-1) + sg(i-1,j  )) * uu(i-1,j) + &
                         (sg(i  ,j-1) + sg(i  ,j  )) * uu(i+1,j) + &
                         (sg(i-1,j  ) + sg(i  ,j  )) * uu(i,j+1) )
                uu(i,j) = uu(i,j) + (one/ss0) * (ff(i,j) - dd)
             end if
          end do
       end do

    else if ( stencil_type .eq. ND_CROSS_STENCIL ) then

       ipar = 1-red_black
       do j = lo(2),hi(2)
          ipar = 1 - ipar
          do i = lo(1)+ipar,hi(1),2
             if ( .not. bc_dirichlet(mm(i,j),1,0) ) then
                  ss0 = -(sg(i-1,j-1)+sg(i,j-1)+sg(i-1,j)+sg(i,j))
                   dd = HALF * ( (sg(i  ,j-1)+sg(i  ,j  )) * uu(i+1,j) + &
                                 (sg(i-1,j-1)+sg(i-1,j  )) * uu(i-1,j) + &
                                 (sg(i-1,j  )+sg(i  ,j  )) * uu(i,j+1) + &
                                 (sg(i-1,j-1)+sg(i  ,j-1)) * uu(i,j-1) ) + &
                                  ss0                      * uu(i,j)
                uu(i,j) = uu(i,j) + (one/ss0) * (ff(i,j) - dd) 
             end if
          end do
       end do

    else
      call bl_error('BAD STENCIL_TYPE IN NODAL_SMOOTHER ',stencil_type)
    end if

  end subroutine nodal_smoother_2d

  subroutine nodal_smoother_3d(sg, uu, ff, mm, lo, ng, uniform_dh, pmask, stencil_type, red_black)

    use impose_neumann_bcs_module

    integer,            intent(in   ) :: ng
    integer,            intent(in   ) :: lo(:)
    logical,            intent(in   ) :: pmask(:)
    real (kind = dp_t), intent(in   ) :: ff(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real (kind = dp_t), intent(in   ) :: sg(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    integer,            intent(in   ) :: mm(lo(1)   :,lo(2)   :,lo(3)   :)
    logical,            intent(in   ) :: uniform_dh
    integer            ,intent(in   ) :: stencil_type
    integer,            intent(in   ) :: red_black

    integer            :: i, j, k, ipar, hi(size(lo)), half_x, half_y
    logical            :: x_is_odd, y_is_odd, jface, kface, doit
    real (kind = dp_t) :: dd, ss0
    real (kind = dp_t) :: f0, fx, fy, fz, fxyz, f2y2zx, f2x2zy, f2x2yz
    real (kind = dp_t), allocatable :: wrk(:,:,:)

    hi(1) = lo(1) + size(mm,dim=1)-1
    hi(2) = lo(2) + size(mm,dim=2)-1
    hi(3) = lo(3) + size(mm,dim=3)-1

    call impose_neumann_bcs_3d(uu,mm,lo,ng)

    if ( stencil_type .eq. ND_CROSS_STENCIL ) then

      half_x = (hi(1)-lo(1))/2
      if ( 2*half_x .eq. ( hi(1)-lo(1) ) ) then
         x_is_odd = .false.
      else
         x_is_odd = .true.
      end if

      half_y = (hi(2)-lo(2))/2
      if ( 2*half_y .eq. ( hi(2)-lo(2) ) ) then
         y_is_odd = .false.
      else
         y_is_odd = .true.
      end if

      if ( (x_is_odd .and. pmask(1)) .or. (y_is_odd .and. pmask(2)) ) then
         !
         ! Use this for Jacobi iteration.
         !
         allocate(wrk(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

         !$OMP PARALLEL DO PRIVATE(i,j,k,dd,jface,kface,doit,ss0)
         do k = lo(3),hi(3)
            kface = .false. ; if ( (k.eq.lo(3)) .or. (k.eq.hi(3)) ) kface = .true.

            do j = lo(2),hi(2)
               jface = .false. ; if ( (j.eq.lo(2)) .or. (j.eq.hi(2)) ) jface = .true.

               do i = lo(1),hi(1)

                  doit = .true.

                  if ( jface .or. kface .or. (i.eq.lo(1)) .or. (i.eq.hi(1)) ) then
                     if ( bc_dirichlet(mm(i,j,k),1,0) ) doit = .false.
                  end if

                  if ( doit ) then
                     ss0 = -(sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1) &
                            +sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1) &
                            +sg(i-1,j-1,k  ) + sg(i-1,j  ,k  ) &
                            +sg(i  ,j-1,k  ) + sg(i  ,j  ,k  )) * THREE 
                     dd =  (sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
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
                           +sg(i  ,j-1,k  ) + sg(i  ,j  ,k  )) * uu(i  ,j  ,k+1) + &
                            ss0                                * uu(i,j,k)
                     wrk(i,j,k) = uu(i,j,k) + (one/ss0) * (ff(i,j,k) - dd) 
                  else
                     wrk(i,j,k) = uu(i,j,k)
                  end if
               end do
            end do
         end do
         !$OMP END PARALLEL DO

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  uu(i,j,k) = wrk(i,j,k)
               end do
            end do
         end do

         deallocate(wrk)

      else
         !
         ! Use this for Gauss-Seidel iteration.
         !
         !$OMP PARALLEL DO PRIVATE(k,ipar,j,i,dd,jface,kface,doit,ss0)
         do k = lo(3),hi(3)
            kface = .false. ; if ( (k.eq.lo(3)) .or. (k.eq.hi(3)) ) kface = .true.

            do j = lo(2),hi(2)
               jface = .false. ; if ( (j.eq.lo(2)) .or. (j.eq.hi(2)) ) jface = .true.

               ipar = MOD(j + k + red_black,2)

               do i = lo(1)+ipar,hi(1),2

                  doit = .true.

                  if ( jface .or. kface .or. (i.eq.lo(1)) .or. (i.eq.hi(1)) ) then
                     if ( bc_dirichlet(mm(i,j,k),1,0) ) doit = .false.
                  end if

                  if (doit) then
                     ss0 = -(sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1) &
                            +sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1) &
                            +sg(i-1,j-1,k  ) + sg(i-1,j  ,k  ) &
                            +sg(i  ,j-1,k  ) + sg(i  ,j  ,k  )) * THREE 
                     dd =  (sg(i-1,j-1,k-1) + sg(i-1,j-1,k  ) &
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
                           +sg(i  ,j-1,k  ) + sg(i  ,j  ,k  )) * uu(i  ,j  ,k+1) + &
                            ss0                                * uu(i,j,k)

                     uu(i,j,k) = uu(i,j,k) + (one/ss0) * (ff(i,j,k) - dd) 
                  end if
               end do
            end do
         end do
         !$OMP END PARALLEL DO

      end if

    else if ( stencil_type .eq. ND_DENSE_STENCIL ) then
       !
       ! Gauss-Seidel.
       !
 
       fx     = ONE / 36._dp_t
       fy     = fx
       fz     = fx
       f0     = FOUR * (fx + fy + fz)
       fxyz   = (fx+fy+fz)
       f2y2zx = (TWO*fy+TWO*fz-fx)
       f2x2zy = (TWO*fx+TWO*fz-fy)
       f2x2yz = (TWO*fx+TWO*fy-fz)             

       do k = lo(3),hi(3)
          kface = .false. ; if ( (k.eq.lo(3)) .or. (k.eq.hi(3)) ) kface = .true.

          do j = lo(2),hi(2)
             jface = .false. ; if ( (j.eq.lo(2)) .or. (j.eq.hi(2)) ) jface = .true.

             do i = lo(1),hi(1)

                doit = .true.

                if ( jface .or. kface .or. (i.eq.lo(1)) .or. (i.eq.hi(1)) ) then
                   if ( bc_dirichlet(mm(i,j,k),1,0) ) doit = .false.
                end if

                if ( doit ) then

                ss0 =  -( sg(i-1,j-1,k-1) + sg(i,j-1,k-1) &
                         +sg(i-1,j  ,k-1) + sg(i,j  ,k-1) &
                         +sg(i-1,j-1,k  ) + sg(i,j-1,k  ) &
                         +sg(i-1,j  ,k  ) + sg(i,j  ,k  ) ) * f0

                dd = fxyz * ( &   ! Corners
                     sg(i-1,j-1,k-1) * uu(i-1,j-1,k-1) + sg(i  ,j-1,k-1) * uu(i+1,j-1,k-1) + &
                     sg(i-1,j  ,k-1) * uu(i-1,j+1,k-1) + sg(i  ,j  ,k-1) * uu(i+1,j+1,k-1) + &
                     sg(i-1,j-1,k  ) * uu(i-1,j-1,k+1) + sg(i  ,j-1,k  ) * uu(i+1,j-1,k+1) + &
                     sg(i-1,j  ,k  ) * uu(i-1,j+1,k+1) + sg(i  ,j  ,k  ) * uu(i+1,j+1,k+1))  &
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
                      dd = dd + &
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
                   dd  = FOUR * dd
                   ss0 = FOUR * ss0

                   uu(i,j,k) = uu(i,j,k) + (one/ss0) * (ff(i,j,k) - dd) 
                end if
             end do
          end do
       end do

    else
      call bl_error('BAD STENCIL_TYPE IN NODAL_SMOOTHER ',stencil_type)
    end if

  end subroutine nodal_smoother_3d

end module nodal_smoothers_module
