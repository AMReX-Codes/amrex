module nodal_smoothers_module

  use bl_constants_module
  use bc_functions_module
  use nodal_stencil_module

  implicit none

contains

  subroutine nodal_line_solve_1d(ss, uu, ff, mm, lo, ng)

    use tridiag_module, only: tridiag

    integer        , intent(in   ) :: lo(:),ng
    real(kind=dp_t), intent(in   ) :: ff(lo(1)-1:)
    real(kind=dp_t), intent(inout) :: uu(lo(1)-ng:)
    real(kind=dp_t), intent(in   ) :: ss(0:,lo(1):)
    integer        , intent(in   ) :: mm(lo(1):)

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
      a_ls(i-is) = ss(2,i)
      b_ls(i-is) = ss(0,i)
      c_ls(i-is) = ss(1,i)
      r_ls(i-is) = ff(i)
    end do

    ! Adjust low end for Neumann boundary condition
    if (bc_neumann(mm(is),1,-1)) then
      c_ls(0) = 2.d0*c_ls(0)
    end if

    ! Adjust high end for Neumann boundary condition
    if (bc_neumann(mm(ie),1,1)) then
      a_ls(0) = 2.d0*a_ls(0)
    end if

    r_ls(0)      = r_ls(0)      - ss(2,is) * uu(is-1)
    r_ls(ilen-1) = r_ls(ilen-1) - ss(1,ie) * uu(ie+1)

    call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,ilen)
 
    do i = is, ie
       uu(i) = u_ls(i-is)
    end do

!   TESTING ONLY 
!   i = lo(1)
!   if (.not. bc_dirichlet(mm(i),1,0)) then
!      dd = ss(0,i)*uu(i) + ss(1,i)*uu(i+1)
!      print *,'RES AT ',i, ff(i) - dd
!   end if

!   do i = lo(1)+1,hi(1)-1
!      dd = ss(0,i)*uu(i) + ss(1,i)*uu(i+1) + ss(2,i)*uu(i-1)
!      print *,'RES AT ',i, ff(i) - dd
!   end do

!   i = hi(1)
!   if (.not. bc_dirichlet(mm(i),1,0)) then
!      dd = ss(0,i)*uu(i) + ss(1,i)*uu(i+1) + ss(2,i)*uu(i-1)
!      print *,'RES AT ',i, ff(i) - dd
!   end if

    deallocate(a_ls,b_ls,c_ls,r_ls,u_ls)

  end subroutine nodal_line_solve_1d

  subroutine nodal_smoother_1d(ss, uu, ff, mm, lo, ng, red_black)
    integer, intent(in) :: lo(:)
    integer, intent(in) :: ng, red_black

    real (kind = dp_t), intent(in)    :: ff(lo(1)-1:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:)
    real (kind = dp_t), intent(in)    :: ss(0:,lo(1):)
    integer            ,intent(in)    :: mm(lo(1):)
    real (kind = dp_t) :: dd
    integer :: i, hi(size(lo))

    real (kind = dp_t), parameter :: omega = 1.33_dp_t

    hi = ubound(uu)-ng

    ! Red/black
    do i = lo(1)+red_black,hi(1),2
       if ( .not. bc_dirichlet(mm(i),1,0) ) then
          dd = ss(0,i)*uu(i) + ss(1,i)*uu(i+1) + ss(2,i)*uu(i-1)
          uu(i) = uu(i) + (omega/ss(0,i)) * (ff(i) - dd)
       end if
    end do

  end subroutine nodal_smoother_1d

  subroutine nodal_smoother_2d(ss, uu, ff, mm, lo, ng, pmask, stencil_type, red_black)

    use bl_prof_module
    use impose_neumann_bcs_module

    integer, intent(in) :: ng
    integer, intent(in) :: lo(:)
    logical, intent(in) :: pmask(:)
    real (kind = dp_t), intent(in) :: ff(lo(1)-1:, lo(2)-1:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:, lo(2)-ng:)
    real (kind = dp_t), intent(in) :: ss(0:,lo(1):,lo(2):)
    integer            ,intent(in) :: mm(lo(1):,lo(2):)
    integer            ,intent(in) :: stencil_type
    integer            ,intent(in) :: red_black

    integer            :: j, i, ipar, hi(size(lo))
    real (kind = dp_t) :: dd
    type(bl_prof_timer), save :: bpt

    call build(bpt, "nodal_smoother_2d")

    hi = ubound(uu)-ng

    call impose_neumann_bcs_2d(uu,mm,lo,ng)

    if ( stencil_type .eq. ND_DENSE_STENCIL ) then

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if ( .not. bc_dirichlet(mm(i,j),1,0) ) then
                dd =   ss(0,i,j) * uu(i  ,j  ) &
                     + ss(1,i,j) * uu(i-1,j-1) &
                     + ss(2,i,j) * uu(i  ,j-1) &
                     + ss(3,i,j) * uu(i+1,j-1) &
                     + ss(4,i,j) * uu(i-1,j  ) &
                     + ss(5,i,j) * uu(i+1,j  ) &
                     + ss(6,i,j) * uu(i-1,j+1) &
                     + ss(7,i,j) * uu(i  ,j+1) &
                     + ss(8,i,j) * uu(i+1,j+1)
                uu(i,j) = uu(i,j) + (one/ss(0,i,j)) * (ff(i,j) - dd)
             end if
          end do
       end do

    else if ( stencil_type .eq. ND_CROSS_STENCIL ) then

       ipar = 1-red_black
       do j = lo(2),hi(2)
          ipar = 1 - ipar
          do i = lo(1)+ipar,hi(1),2
             if ( .not. bc_dirichlet(mm(i,j),1,0) ) then
                dd =   ss(0,i,j) * uu(i  ,j ) &
                     + ss(2,i,j) * uu(i-1,j  ) + ss(1,i,j) * uu(i+1,j  ) &
                     + ss(4,i,j) * uu(i  ,j-1) + ss(3,i,j) * uu(i  ,j+1) 
                uu(i,j) = uu(i,j) + (one/ss(0,i,j)) * (ff(i,j) - dd) 
             end if
          end do
       end do

    else
      call bl_error('BAD STENCIL_TYPE IN NODAL_SMOOTHER ',stencil_type)
    end if

    call destroy(bpt)

  end subroutine nodal_smoother_2d

  subroutine nodal_smoother_3d(ss, uu, ff, mm, lo, ng, uniform_dh, pmask, stencil_type, red_black)

    use bl_prof_module
    use impose_neumann_bcs_module

    integer,            intent(in   ) :: ng
    integer,            intent(in   ) :: lo(:)
    logical,            intent(in   ) :: pmask(:)
    real (kind = dp_t), intent(in   ) :: ff(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real (kind = dp_t), intent(in   ) :: ss(0:,lo(1):,lo(2):,lo(3):)
    integer,            intent(in   ) :: mm(lo(1):,lo(2):,lo(3):)
    logical,            intent(in   ) :: uniform_dh
    integer            ,intent(in   ) :: stencil_type
    integer,            intent(in   ) :: red_black

    integer            :: i, j, k, ipar, hi(size(lo))
    logical            :: jface, kface, doit
    real (kind = dp_t) :: dd
    type(bl_prof_timer), save :: bpt

    call build(bpt, "nodal_smoother_3d")

    hi(1) = lo(1) + size(mm,dim=1)-1
    hi(2) = lo(2) + size(mm,dim=2)-1
    hi(3) = lo(3) + size(mm,dim=3)-1

    call impose_neumann_bcs_3d(uu,mm,lo,ng)

    if ( stencil_type .eq. ND_CROSS_STENCIL ) then

       !$OMP PARALLEL DO PRIVATE(k,ipar,j,i,dd,jface,kface,doit)
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
                   dd =   ss(0,i,j,k) * uu(i  ,j  ,k  ) &
                        + ss(2,i,j,k) * uu(i-1,j  ,k  ) + ss(1,i,j,k) * uu(i+1,j  ,k  ) &
                        + ss(4,i,j,k) * uu(i  ,j-1,k  ) + ss(3,i,j,k) * uu(i  ,j+1,k  ) &
                        + ss(6,i,j,k) * uu(i  ,j  ,k-1) + ss(5,i,j,k) * uu(i  ,j  ,k+1)

                   uu(i,j,k) = uu(i,j,k) + (one/ss(0,i,j,k)) * (ff(i,j,k) - dd) 
                end if
             end do
          end do
       end do
     !$OMP END PARALLEL DO

    else if ( stencil_type .eq. ND_DENSE_STENCIL ) then
       !
       ! Gauss-Seidel.
       !
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
                   dd = ss(0,i,j,k)*uu(i,j,k) &
                        + ss( 1,i,j,k) * uu(i-1,j-1,k-1) + ss( 2,i,j,k) * uu(i  ,j-1,k-1) &
                        + ss( 3,i,j,k) * uu(i+1,j-1,k-1) + ss( 4,i,j,k) * uu(i-1,j  ,k-1) &
                        + ss( 5,i,j,k) * uu(i+1,j  ,k-1) + ss( 6,i,j,k) * uu(i-1,j+1,k-1) &
                        + ss( 7,i,j,k) * uu(i  ,j+1,k-1) + ss( 8,i,j,k) * uu(i+1,j+1,k-1) &
                        + ss( 9,i,j,k) * uu(i-1,j-1,k  ) + ss(10,i,j,k) * uu(i+1,j-1,k  ) &
                        + ss(11,i,j,k) * uu(i-1,j+1,k  ) + ss(12,i,j,k) * uu(i+1,j+1,k  ) &
                        + ss(13,i,j,k) * uu(i-1,j-1,k+1) + ss(14,i,j,k) * uu(i  ,j-1,k+1) &
                        + ss(15,i,j,k) * uu(i+1,j-1,k+1) + ss(16,i,j,k) * uu(i-1,j  ,k+1) &
                        + ss(17,i,j,k) * uu(i+1,j  ,k+1) + ss(18,i,j,k) * uu(i-1,j+1,k+1) &
                        + ss(19,i,j,k) * uu(i  ,j+1,k+1) + ss(20,i,j,k) * uu(i+1,j+1,k+1) 

                   if ( (size(ss,dim=1) .eq. 27) .and. (.not. uniform_dh) ) then
                      !
                      ! Add faces (only non-zero for non-uniform dx)
                      !
                      dd = dd + &
                           ss(21,i,j,k) * uu(i-1,j  ,k  ) + ss(22,i,j,k) * uu(i+1,j  ,k  ) &
                           + ss(23,i,j,k) * uu(i  ,j-1,k  ) + ss(24,i,j,k) * uu(i  ,j+1,k  ) &
                           + ss(25,i,j,k) * uu(i  ,j  ,k-1) + ss(26,i,j,k) * uu(i  ,j  ,k+1)
                   end if

                   uu(i,j,k) = uu(i,j,k) + (one/ss(0,i,j,k)) * (ff(i,j,k) - dd) 
                end if
             end do
          end do
       end do

    else
      call bl_error('BAD STENCIL_TYPE IN NODAL_SMOOTHER ',stencil_type)
    end if

    call destroy(bpt)

  end subroutine nodal_smoother_3d

  subroutine nodal_smoother_3d_opt(ss, uu, ff, mm, lo, ng, uniform_dh, pmask, &
                                   stencil_type, dx)

    use bl_prof_module
    use impose_neumann_bcs_module

    integer,            intent(in   ) :: ng
    integer,            intent(in   ) :: lo(:)
    logical,            intent(in   ) :: pmask(:)
    real (kind = dp_t), intent(in   ) :: ff(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real (kind = dp_t), intent(in   ) :: ss(0:,lo(1):,lo(2):,lo(3):)
    integer,            intent(in   ) :: mm(lo(1):,lo(2):,lo(3):)
    logical,            intent(in   ) :: uniform_dh
    integer            ,intent(in   ) :: stencil_type
    real (kind = dp_t), intent(in   ) :: dx(:)

    integer            :: i, j, k, hi(size(lo))
    logical            :: jface, kface, doit
    real (kind = dp_t) :: dd
    type(bl_prof_timer), save :: bpt

    call build(bpt, "nodal_smoother_3d")

    hi(1) = lo(1) + size(mm,dim=1)-1
    hi(2) = lo(2) + size(mm,dim=2)-1
    hi(3) = lo(3) + size(mm,dim=3)-1

    print *,'In nodal_smoother_3d_opt with dx ', dx(1)

    call impose_neumann_bcs_3d(uu,mm,lo,ng)

    if ( stencil_type .ne. ND_DENSE_STENCIL ) then
      call bl_error('Bad stencil_type in nodal_smoother_3d_opt: ',stencil_type)
    end if

    if ( .not. uniform_dh) then
      call bl_error('Want uniform_dh = true in nodal_smoother_3d_opt: ')
    end if

    !
    ! Gauss-Seidel.
    !
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
                   dd = ss(0,i,j,k)*uu(i,j,k) &
                        + ss( 1,i,j,k) * uu(i-1,j-1,k-1) + ss( 2,i,j,k) * uu(i  ,j-1,k-1) &
                        + ss( 3,i,j,k) * uu(i+1,j-1,k-1) + ss( 4,i,j,k) * uu(i-1,j  ,k-1) &
                        + ss( 5,i,j,k) * uu(i+1,j  ,k-1) + ss( 6,i,j,k) * uu(i-1,j+1,k-1) &
                        + ss( 7,i,j,k) * uu(i  ,j+1,k-1) + ss( 8,i,j,k) * uu(i+1,j+1,k-1) &
                        + ss( 9,i,j,k) * uu(i-1,j-1,k  ) + ss(10,i,j,k) * uu(i+1,j-1,k  ) &
                        + ss(11,i,j,k) * uu(i-1,j+1,k  ) + ss(12,i,j,k) * uu(i+1,j+1,k  ) &
                        + ss(13,i,j,k) * uu(i-1,j-1,k+1) + ss(14,i,j,k) * uu(i  ,j-1,k+1) &
                        + ss(15,i,j,k) * uu(i+1,j-1,k+1) + ss(16,i,j,k) * uu(i-1,j  ,k+1) &
                        + ss(17,i,j,k) * uu(i+1,j  ,k+1) + ss(18,i,j,k) * uu(i-1,j+1,k+1) &
                        + ss(19,i,j,k) * uu(i  ,j+1,k+1) + ss(20,i,j,k) * uu(i+1,j+1,k+1) 

                   uu(i,j,k) = uu(i,j,k) + (one/ss(0,i,j,k)) * (ff(i,j,k) - dd) 
                end if
          end do
       end do
    end do

    call destroy(bpt)

  end subroutine nodal_smoother_3d_opt

end module nodal_smoothers_module
