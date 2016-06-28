module cc_smoothers_module

  use bl_constants_module
  use cc_stencil_module
  use stencil_types_module
  use bc_functions_module

  implicit none

  private

  public :: gs_line_solve_1d, &
       fourth_order_smoother_2d, fourth_order_smoother_3d, &
       gs_rb_smoother_1d, gs_rb_smoother_2d, gs_rb_smoother_3d, &
       jac_smoother_1d, jac_smoother_2d, jac_smoother_3d, &
       gs_rb_smoother_ibc_2d, gs_rb_smoother_ibc_3d, &
       jac_smoother_ibc_2d, jac_smoother_ibc_3d

contains

  subroutine gs_line_solve_1d(ss, uu, ff, mm, lo, ng)

    use tridiag_module, only : tridiag

    integer, intent(in) :: lo(:)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in)    :: ff(lo(1):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:)
    real (kind = dp_t), intent(in)    :: ss(0:,lo(1):)
    integer            ,intent(in)    :: mm(lo(1):)

    real (kind = dp_t), allocatable :: a_ls(:), b_ls(:), c_ls(:), r_ls(:), u_ls(:)
    integer :: ilen, i, hi(size(lo))
    integer, parameter ::  XBC = 3

    hi = ubound(ff)

    ilen = hi(1)-lo(1)+1
    allocate(a_ls(0:hi(1)-lo(1)))
    allocate(b_ls(0:hi(1)-lo(1)))
    allocate(c_ls(0:hi(1)-lo(1)))
    allocate(r_ls(0:hi(1)-lo(1)))
    allocate(u_ls(0:hi(1)-lo(1)))

    do i = lo(1), hi(1)
      a_ls(i-lo(1)) = ss(2,i)
      b_ls(i-lo(1)) = ss(0,i)
      c_ls(i-lo(1)) = ss(1,i)
      r_ls(i-lo(1)) = ff(i)

      if ( hi(1) > lo(1) ) then
         if (bc_skewed(mm(i),1,+1)) then
            r_ls(i-lo(1)) = r_ls(i-lo(1)) - ss(XBC,i)*uu(i+2)
         else if (bc_skewed(mm(i),1,-1)) then
            r_ls(i-lo(1)) = r_ls(i-lo(1)) - ss(XBC,i)*uu(i-2)
         end if
      end if

    end do

    r_ls(0)           = r_ls(0)           - ss(2,lo(1)) * uu(lo(1)-1)
    r_ls(hi(1)-lo(1)) = r_ls(hi(1)-lo(1)) - ss(1,hi(1)) * uu(hi(1)+1)

    call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,ilen)
 
    do i = lo(1), hi(1)
       uu(i) = u_ls(i-lo(1))
    end do

  end subroutine gs_line_solve_1d

  subroutine gs_rb_smoother_1d(ss, uu, ff, mm, lo, ng, n, skwd)
    use bl_prof_module
    integer, intent(in) :: ng
    integer, intent(in) :: lo(:)
    integer, intent(in) :: n
    real (kind = dp_t), intent(in) :: ff(lo(1):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:)
    real (kind = dp_t), intent(in) :: ss(0:,lo(1):)
    integer            ,intent(in) :: mm(lo(1):)
    logical, intent(in), optional :: skwd
    integer :: i, hi(size(lo)), ioff
    integer, parameter ::  XBC = 3
    real (kind = dp_t) :: dd
    logical :: lskwd
    real(dp_t) :: lr(2)

    real (kind = dp_t), parameter :: omega = 0.6_dp_t

    type(bl_prof_timer), save :: bpt

    call build(bpt, "gs_rb_smoother_1d")

    lr = 0

    hi = ubound(ff)

    if (present(skwd) ) then
       lskwd = skwd
    else 
       lskwd = .false.
       if (bc_skewed(mm(lo(1)),1,+1)) lskwd = .true.
       if (bc_skewed(mm(hi(1)),1,-1)) lskwd = .true.
    end if

    !! assumption: ss(0,i) vanishes only for 1x1 problems
    if ( all(lo == hi) ) then
       i = lo(1)
       if ( mod(i,2) == n ) then
          if ( abs(ss(0,i)) .gt. zero ) then
             dd = ss(0,i)*uu(i) &
                + ss(1,i)*uu(i+1) + ss(2,i)*uu(i-1) 
             uu(i) = uu(i) + (omega/ss(0,i)) * (ff(i) - dd)
          end if
       end if
       call destroy(bpt)
       return
    end if

    !! Assumption; bc_skewed is never true if hi==lo
    if ( lskwd ) then

       if (bc_skewed(mm(lo(1)),1,+1)) lr(1) = uu(lo(1)+2)
       if (bc_skewed(mm(hi(1)),1,-1)) lr(2) = uu(hi(1)-2)

       ioff = 0; if ( mod(lo(1), 2) /= n ) ioff = 1
       do i = lo(1) + ioff, hi(1), 2

             dd = ss(0,i)*uu(i) &
                  + ss(1,i) * uu(i+1) + ss(2,i) * uu(i-1)
             if ( i == lo(1) .or. i == hi(1)) then
                if (bc_skewed(mm(i),1,+1)) then
                   dd = dd + ss(XBC,i)*lr(1)
                end if
                if (bc_skewed(mm(i),1,-1)) then
                   dd = dd + ss(XBC,i)*lr(2)
                end if
             end if
             uu(i) = uu(i) + (omega/ss(0,i)) * (ff(i) - dd)
          end do

    else

       ioff = 0; if ( mod(lo(1), 2) /= n ) ioff = 1
       do i = lo(1) + ioff, hi(1), 2
          if (abs(ss(0,i)) .gt. zero) then
            dd = ss(0,i)*uu(i) &
               + ss(1,i)*uu(i+1) + ss(2,i)*uu(i-1) 
            uu(i) = uu(i) + (omega/ss(0,i)) * (ff(i) - dd)
          end if
       end do

    end if

    call destroy(bpt)

  end subroutine gs_rb_smoother_1d

  subroutine gs_rb_smoother_2d(ss, uu, ff, mm, lo, ng, n, skwd)
    use bl_prof_module
    integer, intent(in) :: ng
    integer, intent(in) :: lo(:)
    integer, intent(in) :: n
    real (kind = dp_t), intent(in) :: ff(lo(1):, lo(2):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:, lo(2)-ng:)
    real (kind = dp_t), intent(in) :: ss(0:,lo(1):,lo(2):)
    integer            ,intent(in) :: mm(lo(1):,lo(2):)
    logical, intent(in), optional :: skwd
    integer :: j, i, hi(size(lo)), ioff
    integer, parameter ::  XBC = 5, YBC = 6
    real (kind = dp_t) :: dd
    logical :: lskwd
    real(dp_t) :: lr(lbound(ff,2):ubound(ff,2), 2)
    real(dp_t) :: tb(lbound(ff,1):ubound(ff,1), 2)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "gs_rb_smoother_2d")

    hi = ubound(ff)

    if (present(skwd) ) then
       lskwd = skwd
    else 
       lskwd = .false.
       do i = lo(1), hi(1)
          if (bc_skewed(mm(i,lo(2)),2,+1)) lskwd = .true.
          if (bc_skewed(mm(i,hi(2)),2,-1)) lskwd = .true.
       end do
       do j = lo(2), hi(2)
          if (bc_skewed(mm(lo(1),j),1,+1)) lskwd = .true.
          if (bc_skewed(mm(hi(1),j),1,-1)) lskwd = .true.
       end do
    end if

    !! assumption: ss(0,i,j) vanishes only for 1x1 problems
    if ( all(lo == hi) ) then
       i = lo(1); j = lo(2)
       if ( mod(i + j,2) == n ) then
          if ( abs(ss(0,i,j)) .gt. zero ) then
             dd = ss(0,i,j)*uu(i,j) &
                  + ss(1,i,j)*uu(i+1,j) + ss(2,i,j)*uu(i-1,j) &
                  + ss(3,i,j)*uu(i,j+1) + ss(4,i,j)*uu(i,j-1)
             uu(i,j) = uu(i,j) + (one/ss(0,i,j))*(ff(i,j) - dd) 
          end if
       end if
       call destroy(bpt)
       return
    end if
    !! Assumption; bc_skewed is never true if hi==lo

    if ( lskwd ) then

       do i = lo(1), hi(1)
          if (bc_skewed(mm(i,lo(2)),2,+1)) tb(i,1) = uu(i,lo(2)+2)
          if (bc_skewed(mm(i,hi(2)),2,-1)) tb(i,2) = uu(i,hi(2)-2)
       end do

       do j = lo(2), hi(2)
          if (bc_skewed(mm(lo(1),j),1,+1)) lr(j,1) = uu(lo(1)+2,j)
          if (bc_skewed(mm(hi(1),j),1,-1)) lr(j,2) = uu(hi(1)-2,j)
       end do

       do j = lo(2),hi(2)
          ioff = 0; if ( mod(lo(1) + j, 2) /= n ) ioff = 1
          do i = lo(1) + ioff, hi(1), 2

             dd = ss(0,i,j)*uu(i,j) &
                  + ss(1,i,j) * uu(i+1,j) + ss(2,i,j) * uu(i-1,j) &
                  + ss(3,i,j) * uu(i,j+1) + ss(4,i,j) * uu(i,j-1)
             if ( i == lo(1) .or. i == hi(1) .or. j == lo(2) .or. j == hi(2) ) then
                if (bc_skewed(mm(i,j),1,+1)) then
                   dd = dd + ss(XBC,i,j)*lr(j,1)
                end if
                if (bc_skewed(mm(i,j),1,-1)) then
                   dd = dd + ss(XBC,i,j)*lr(j,2)
                end if
                if (bc_skewed(mm(i,j),2,+1)) then
                   dd = dd + ss(YBC,i,j)*tb(i,1)
                end if
                if (bc_skewed(mm(i,j),2,-1)) then
                   dd = dd + ss(YBC,i,j)*tb(i,2)
                end if
             end if
             uu(i,j) = uu(i,j) + (one/ss(0,i,j))*(ff(i,j) - dd)
          end do
       end do
    else
       !
       ! USE THIS FOR GAUSS-SEIDEL
       !
       do j = lo(2),hi(2)
          ioff = 0; if ( mod(lo(1) + j, 2) /= n ) ioff = 1
          do i = lo(1) + ioff, hi(1), 2
             if (abs(ss(0,i,j)) .gt. zero) then
               dd = ss(0,i,j)*uu(i,j) &
                    + ss(1,i,j) * uu(i+1,j) + ss(2,i,j) * uu(i-1,j) &
                    + ss(3,i,j) * uu(i,j+1) + ss(4,i,j) * uu(i,j-1)
               uu(i,j) = uu(i,j) + (one/ss(0,i,j))*(ff(i,j) - dd) 
             end if
          end do
       end do

    end if

    call destroy(bpt)

  end subroutine gs_rb_smoother_2d

  subroutine gs_rb_smoother_3d(ss, uu, ff, mm, lo, tlo, thi, ng, n, skwd)
    use bl_prof_module
    use bl_error_module
    integer,    intent(in   ) :: ng, lo(:), tlo(:), thi(:), n
    real(dp_t), intent(in   ) :: ff(lo(1):,lo(2):,lo(3):)
    real(dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(dp_t), intent(in   ) :: ss(0:,lo(1):, lo(2):, lo(3):)
    integer,    intent(in   ) :: mm(lo(1):,lo(2):,lo(3):)

    logical, intent(in), optional :: skwd

    integer    :: i, j, k, ioff, hi(size(lo))
    logical    :: lskwd
    real(dp_t) :: dd
    !
    ! These are small so we'll put'm on the stack instead of the heap.
    !
    real(dp_t) lr(tlo(2):thi(2), tlo(3):thi(3), 2)
    real(dp_t) tb(tlo(1):thi(1), tlo(3):thi(3), 2)
    real(dp_t) fb(tlo(1):thi(1), tlo(2):thi(2), 2)

    integer, parameter ::  XBC = 7, YBC = 8, ZBC = 9

    real(dp_t), parameter :: omega = 1.15_dp_t

    type(bl_prof_timer), save :: bpt

    call bl_proffortfuncstart("gs_rb_smoother_3d")
    call build(bpt, "gs_rb_smoother_3d")

    hi = ubound(ff)

    if ( all(lo == hi) ) then
       k = lo(3); j = lo(2); i = lo(1)
       if ( mod(i + j + k, 2) == n ) then
          if (abs(ss(0,i,j,k)) .gt. zero) then
             dd = ss(0,i,j,k)*uu(i,j,k)   + &
                  ss(1,i,j,k)*uu(i+1,j,k) + ss(2,i,j,k)*uu(i-1,j,k) + &
                  ss(3,i,j,k)*uu(i,j+1,k) + ss(4,i,j,k)*uu(i,j-1,k) + &
                  ss(5,i,j,k)*uu(i,j,k+1) + ss(6,i,j,k)*uu(i,j,k-1)
             uu(i,j,k) = uu(i,j,k) + (omega/ss(0,i,j,k))*(ff(i,j,k) - dd)
          end if
       end if
       call destroy(bpt)
       call bl_proffortfuncstop("gs_rb_smoother_3d")
       return
    end if

    if ( present(skwd) ) then
       lskwd = skwd
    else
       lskwd = .false.
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             if ( bc_skewed(mm(lo(1),j,k),1,+1) .or. bc_skewed(mm(hi(1),j,k),1,-1) ) then
                lskwd = .true.
                goto 1234
             end if
          end do
       end do
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             if ( bc_skewed(mm(i,lo(2),k),2,+1) .or. bc_skewed(mm(i,hi(2),k),2,-1) ) then
                lskwd = .true.
                goto 1234
             end if
          end do
       end do
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if ( bc_skewed(mm(i,j,lo(3)),3,+1) .or. bc_skewed(mm(i,j,hi(3)),3,-1) ) then
                lskwd = .true.
                goto 1234
             end if
          end do
       end do
    end if

1234 if ( lskwd ) then

       call bl_assert(lo.eq.tlo .and. hi.eq.thi, "Tiling must be turned off for skewed stencil.")

       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             if (bc_skewed(mm(i,lo(2),k),2,+1)) tb(i,k,1) = uu(i,lo(2)+2,k)
             if (bc_skewed(mm(i,hi(2),k),2,-1)) tb(i,k,2) = uu(i,hi(2)-2,k)
          end do
       end do

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             if (bc_skewed(mm(lo(1),j,k),1,+1)) lr(j,k,1) = uu(lo(1)+2,j,k)
             if (bc_skewed(mm(hi(1),j,k),1,-1)) lr(j,k,2) = uu(hi(1)-2,j,k)
          end do
       end do

       do j = lo(2), hi(2) 
          do i = lo(1), hi(1)
             if (bc_skewed(mm(i,j,lo(3)),3,+1)) fb(i,j,1) = uu(i,j,lo(3)+2)
             if (bc_skewed(mm(i,j,hi(3)),3,-1)) fb(i,j,2) = uu(i,j,hi(3)-2)
          end do
       end do

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             ioff = 0; if ( mod(lo(1) + j + k, 2) /= n ) ioff = 1
             do i = lo(1)+ioff, hi(1), 2

                dd = ss(0,i,j,k)*uu(i,j,k)   + &
                     ss(1,i,j,k)*uu(i+1,j,k) + ss(2,i,j,k)*uu(i-1,j,k) + &
                     ss(3,i,j,k)*uu(i,j+1,k) + ss(4,i,j,k)*uu(i,j-1,k) + &
                     ss(5,i,j,k)*uu(i,j,k+1) + ss(6,i,j,k)*uu(i,j,k-1)

                if ( i == lo(1) ) then
                   if ( bc_skewed(mm(i,j,k),1,+1) ) dd = dd + ss(XBC,i,j,k)*lr(j,k,1)
                else if ( i == hi(1) ) then
                   if ( bc_skewed(mm(i,j,k),1,-1) ) dd = dd + ss(XBC,i,j,k)*lr(j,k,2)
                end if

                if ( j == lo(2) ) then
                   if ( bc_skewed(mm(i,j,k),2,+1) ) dd = dd + ss(YBC,i,j,k)*tb(i,k,1)
                else if ( j == hi(2) ) then
                   if ( bc_skewed(mm(i,j,k),2,-1) ) dd = dd + ss(YBC,i,j,k)*tb(i,k,2)
                end if

                if ( k == lo(3) ) then
                   if ( bc_skewed(mm(i,j,k),3,+1) ) dd = dd + ss(ZBC,i,j,k)*fb(i,j,1)
                else if ( k == hi(3) ) then
                   if ( bc_skewed(mm(i,j,k),3,-1) ) dd = dd + ss(ZBC,i,j,k)*fb(i,j,2)
                end if

                uu(i,j,k) = uu(i,j,k) + (omega/ss(0,i,j,k))*(ff(i,j,k) - dd)
             end do
          end do
       end do

    else
       !
       ! USE THIS FOR GAUSS-SEIDEL
       !
       do k = tlo(3), thi(3)
          do j = tlo(2), thi(2)
             ioff = 0; if ( mod (tlo(1) + j + k, 2) /= n ) ioff = 1
             do i = tlo(1)+ioff, thi(1), 2

                  dd = ss(0,i,j,k)*uu(i,j,k)   + &
                       ss(1,i,j,k)*uu(i+1,j,k) + ss(2,i,j,k)*uu(i-1,j,k) + &
                       ss(3,i,j,k)*uu(i,j+1,k) + ss(4,i,j,k)*uu(i,j-1,k) + &
                       ss(5,i,j,k)*uu(i,j,k+1) + ss(6,i,j,k)*uu(i,j,k-1)

                  uu(i,j,k) = uu(i,j,k) + (omega/ss(0,i,j,k))*(ff(i,j,k) - dd)
             end do
          end do
       end do

    end if

    call destroy(bpt)
    call bl_proffortfuncstop("gs_rb_smoother_3d")

  end subroutine gs_rb_smoother_3d

  subroutine fourth_order_smoother_2d(ss, uu, ff, lo, ng, stencil_type, n)
    use bl_prof_module
    integer           , intent(in) :: ng, n
    integer           , intent(in) :: lo(:)
    real (kind = dp_t), intent(in) :: ff(lo(1):, lo(2):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:, lo(2)-ng:)
    real (kind = dp_t), intent(in   ) :: ss(0:,lo(1):, lo(2):)
    integer           , intent(in   ) :: stencil_type

    integer            :: i, j, hi(size(lo)), ioff
    real (kind = dp_t) :: dd

    type(bl_prof_timer), save :: bpt

    call bl_proffortfuncstart("fourth_order_smoother_2d")
    call build(bpt, "fourth_order_smoother_2d")

    hi = ubound(ff)

    if (stencil_type .eq. HO_CROSS_STENCIL) then

       do j = lo(2),hi(2)
          ioff = 0; if ( mod(lo(1) + j, 2) /= n ) ioff = 1
          do i = lo(1) + ioff, hi(1), 2
             if (abs(ss(0,i,j)) .gt. zero) then
               dd =   ss(0,i,j) * uu(i,j) &
                    + ss(1,i,j) * uu(i-2,j) + ss(2,i,j) * uu(i-1,j) &
                    + ss(3,i,j) * uu(i+1,j) + ss(4,i,j) * uu(i+2,j) &
                    + ss(5,i,j) * uu(i,j-2) + ss(6,i,j) * uu(i,j-1) &
                    + ss(7,i,j) * uu(i,j+1) + ss(8,i,j) * uu(i,j+2)
               uu(i,j) = uu(i,j) + (one/ss(0,i,j)) * (ff(i,j) - dd)
             end if
          end do
       end do

    else if (stencil_type .eq. HO_DENSE_STENCIL) then

       do j = lo(2),hi(2)
          do i = lo(1), hi(1)
             if ( abs(ss(0,i,j)) .gt. zero ) then
               dd =   ss( 0,i,j) * uu(i,j) &
                    + ss( 1,i,j) * uu(i-2,j-2) + ss( 2,i,j) * uu(i-1,j-2) & ! AT J-2
                    + ss( 3,i,j) * uu(i  ,j-2) + ss( 4,i,j) * uu(i+1,j-2) & ! AT J-2
                    + ss( 5,i,j) * uu(i+2,j-2)                            & ! AT J-2
                    + ss( 6,i,j) * uu(i-2,j-1) + ss( 7,i,j) * uu(i-1,j-1) & ! AT J-1
                    + ss( 8,i,j) * uu(i  ,j-1) + ss( 9,i,j) * uu(i+1,j-1) & ! AT J-1
                    + ss(10,i,j) * uu(i+2,j-1)                            & ! AT J-1
                    + ss(11,i,j) * uu(i-2,j  ) + ss(12,i,j) * uu(i-1,j  ) & ! AT J
                    + ss(13,i,j) * uu(i+1,j  ) + ss(14,i,j) * uu(i+2,j  ) & ! AT J
                    + ss(15,i,j) * uu(i-2,j+1) + ss(16,i,j) * uu(i-1,j+1) & ! AT J+1
                    + ss(17,i,j) * uu(i  ,j+1) + ss(18,i,j) * uu(i+1,j+1) & ! AT J+1
                    + ss(19,i,j) * uu(i+2,j+1)                            & ! AT J+1
                    + ss(20,i,j) * uu(i-2,j+2) + ss(21,i,j) * uu(i-1,j+2) & ! AT J+2
                    + ss(22,i,j) * uu(i  ,j+2) + ss(23,i,j) * uu(i+1,j+2) & ! AT J+2
                    + ss(24,i,j) * uu(i+2,j+2)                              ! AT J+2

               uu(i,j) = uu(i,j) + (one/ss(0,i,j)) * (ff(i,j) - dd)

             end if
          end do
       end do

    end if

    call destroy(bpt)
    call bl_proffortfuncstop("fourth_order_smoother_2d")

  end subroutine fourth_order_smoother_2d

  subroutine fourth_order_smoother_3d(ss, uu, ff, lo, ng, stencil_type)
    use bl_prof_module
    integer           , intent(in   ) :: ng
    integer           , intent(in   ) :: lo(:)
    real (kind = dp_t), intent(in   ) :: ff(lo(1):, lo(2):, lo(3):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:, lo(2)-ng:,lo(3)-ng:)
    real (kind = dp_t), intent(in   ) :: ss(0:, lo(1):, lo(2):, lo(3):)
    integer           , intent(in   ) :: stencil_type

    integer            :: i, j, k, hi(size(lo))
    real (kind = dp_t) :: dd

    type(bl_prof_timer), save :: bpt

    call bl_proffortfuncstart("fourth_order_smoother_3d")
    call build(bpt, "fourth_order_smoother_3d")

    hi = ubound(ff)

    ! This is the fourth order stencil for constant coefficients.
    if (stencil_type .eq. HO_CROSS_STENCIL) then

       do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1), hi(1)
             if ( abs(ss(0,i,j,k)) .gt. zero ) then
               dd =   ss( 0,i,j,k) * uu(i,j,k) &
                    + ss( 1,i,j,k) * uu(i-2,j,k) + ss( 2,i,j,k) * uu(i-1,j,k) &
                    + ss( 3,i,j,k) * uu(i+1,j,k) + ss( 4,i,j,k) * uu(i+2,j,k) &
                    + ss( 5,i,j,k) * uu(i,j-2,k) + ss( 6,i,j,k) * uu(i,j-1,k) &
                    + ss( 7,i,j,k) * uu(i,j+1,k) + ss( 8,i,j,k) * uu(i,j+2,k) &
                    + ss( 9,i,j,k) * uu(i,j,k-2) + ss(10,i,j,k) * uu(i,j,k-1) &
                    + ss(11,i,j,k) * uu(i,j,k+1) + ss(12,i,j,k) * uu(i,j,k+2)
               uu(i,j,k) = uu(i,j,k) + (one/ss(0,i,j,k)) * (ff(i,j,k) - dd)
             end if
          end do
       end do
       end do

    ! This is the fourth order stencil for variable coefficients.
    else if (stencil_type .eq. HO_DENSE_STENCIL) then

       do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1), hi(1)
             if (abs(ss(0,i,j,k)) .gt. zero) then

                dd = &
                       ss( 0,i,j,k) * uu(i  ,j  ,k  ) &
                       ! Contributions from k-2
                     + ss( 1,i,j,k) * uu(i  ,j-2,k-2) + ss( 2,i,j,k) * uu(i  ,j-1,k-2) &
                     + ss( 3,i,j,k) * uu(i-2,j  ,k-2) + ss( 4,i,j,k) * uu(i-1,j  ,k-2) &
                     + ss( 5,i,j,k) * uu(i  ,j  ,k-2) + ss( 6,i,j,k) * uu(i+1,j  ,k-2) &
                     + ss( 7,i,j,k) * uu(i+2,j  ,k-2) + ss( 8,i,j,k) * uu(i  ,j+1,k-2) &
                     + ss( 9,i,j,k) * uu(i  ,j+2,k-2)                                  &
                       ! Contributions from k-1
                     + ss(10,i,j,k) * uu(i  ,j-2,k-1) + ss(11,i,j,k) * uu(i  ,j-1,k-1) &
                     + ss(12,i,j,k) * uu(i-2,j  ,k-1) + ss(13,i,j,k) * uu(i-1,j  ,k-1) &
                     + ss(14,i,j,k) * uu(i  ,j  ,k-1) + ss(15,i,j,k) * uu(i+1,j  ,k-1) &
                     + ss(16,i,j,k) * uu(i+2,j  ,k-1) + ss(17,i,j,k) * uu(i  ,j+1,k-1) &
                     + ss(18,i,j,k) * uu(i  ,j+2,k-1)                                  &
                       ! Contributions from j-2,k
                     + ss(19,i,j,k) * uu(i-2,j-2,k  ) + ss(20,i,j,k) * uu(i-1,j-2,k  ) &
                     + ss(21,i,j,k) * uu(i  ,j-2,k  ) + ss(22,i,j,k) * uu(i+1,j-2,k  ) &
                     + ss(23,i,j,k) * uu(i+2,j-2,k  )

                dd = dd &
                       ! Contributions from j-1,k
                     + ss(24,i,j,k) * uu(i-2,j-1,k  ) + ss(25,i,j,k) * uu(i-1,j-1,k  ) &
                     + ss(26,i,j,k) * uu(i  ,j-1,k  ) + ss(27,i,j,k) * uu(i+1,j-1,k  ) &
                     + ss(28,i,j,k) * uu(i+2,j-1,k  )                                  &
                       ! Contributions from j  ,k
                     + ss(29,i,j,k) * uu(i-2,j  ,k  ) + ss(30,i,j,k) * uu(i-1,j  ,k  ) &
                                                      + ss(31,i,j,k) * uu(i+1,j  ,k  ) &
                     + ss(32,i,j,k) * uu(i+2,j  ,k  )                                  &
                       ! Contributions from j+1,k
                     + ss(33,i,j,k) * uu(i-2,j+1,k  ) + ss(34,i,j,k) * uu(i-1,j+1,k  ) &
                     + ss(35,i,j,k) * uu(i  ,j+1,k  ) + ss(36,i,j,k) * uu(i+1,j+1,k  ) &
                     + ss(37,i,j,k) * uu(i+2,j+1,k  )                                  &
                       ! Contributions from j+2,k
                     + ss(38,i,j,k) * uu(i-2,j+2,k  ) + ss(39,i,j,k) * uu(i-1,j+2,k  ) &
                     + ss(40,i,j,k) * uu(i  ,j+2,k  ) + ss(41,i,j,k) * uu(i+1,j+2,k  ) &
                     + ss(42,i,j,k) * uu(i+2,j+2,k  )

                dd = dd &
                       ! Contributions from k+1
                     + ss(43,i,j,k) * uu(i  ,j-2,k+1) + ss(44,i,j,k) * uu(i  ,j-1,k+1) &
                     + ss(45,i,j,k) * uu(i-2,j  ,k+1) + ss(46,i,j,k) * uu(i-1,j  ,k+1) &
                     + ss(47,i,j,k) * uu(i  ,j  ,k+1) + ss(48,i,j,k) * uu(i+1,j  ,k+1) &
                     + ss(49,i,j,k) * uu(i+2,j  ,k+1) + ss(50,i,j,k) * uu(i  ,j+1,k+1) &
                     + ss(51,i,j,k) * uu(i  ,j+2,k+1)                                  &
                       ! Contributions from k+2
                     + ss(52,i,j,k) * uu(i  ,j-2,k+2) + ss(53,i,j,k) * uu(i  ,j-1,k+2) &
                     + ss(54,i,j,k) * uu(i-2,j  ,k+2) + ss(55,i,j,k) * uu(i-1,j  ,k+2) &
                     + ss(56,i,j,k) * uu(i  ,j  ,k+2) + ss(57,i,j,k) * uu(i+1,j  ,k+2) &
                     + ss(58,i,j,k) * uu(i+2,j  ,k+2) + ss(59,i,j,k) * uu(i  ,j+1,k+2) &
                     + ss(60,i,j,k) * uu(i  ,j+2,k+2)

               uu(i,j,k) = uu(i,j,k) + (one/ss(0,i,j,k)) * (ff(i,j,k) - dd)
             end if
          end do
       end do
       end do

    end if

    call destroy(bpt)
    call bl_proffortfuncstop("fourth_order_smoother_3d")

  end subroutine fourth_order_smoother_3d

  subroutine jac_smoother_1d(ss, uu, ff, mm, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) ::  ss(0:,:)
    real (kind = dp_t), intent(in) ::  ff(:)
    integer           , intent(in) ::  mm(:)
    real (kind = dp_t), intent(inout) ::  uu(1-ng:)
    real (kind = dp_t):: wrk(size(ff,1))
    integer :: nx
    integer :: i
    real (kind = dp_t) :: dd
    integer, parameter ::  XBC = 3

    nx = size(ff,dim=1)
    do i = 1, nx
       dd = + ss(1,i)*uu(i+1) + ss(2,i)*uu(i-1) 
       if ( nx > 1 ) then
          if ( bc_skewed(mm(i),1,+1) ) then
             dd = dd + ss(XBC,i)*uu(i+2)
          else if ( bc_skewed(mm(i),1,-1) ) then
             dd = dd + ss(XBC,i)*uu(i-2)
          end if
       end if
       wrk(i) = ff(i) - dd
    end do
    do i = 1, nx
       uu(i) = uu(i) + (wrk(i)/ss(0,i)-uu(i))
    end do

  end subroutine jac_smoother_1d

  ! subroutine jac_dense_smoother_1d(ss, uu, ff, ng)
  !   integer, intent(in) :: ng
  !   real (kind = dp_t), intent(in) ::  ss(:,:)
  !   real (kind = dp_t), intent(in) ::  ff(:)
  !   real (kind = dp_t), intent(inout) ::  uu(1-ng:)
  !   real (kind = dp_t):: wrk(size(ff,1))
  !   integer :: nx, i

  !   nx = size(ff,dim=1)
  !   do i = 1, nx
  !      wrk(i) = ( ff(i) - ss(1,i)*uu(i+1) - ss(3,i)*uu(i-1) )
  !   end do
  !   do i = 1, nx
  !      uu(i) = uu(i) + (wrk(i)/ss(2,i)-uu(i))
  !   end do

  ! end subroutine jac_dense_smoother_1d

  subroutine jac_smoother_2d(ss, uu, ff, mm, ng)
    use bl_prof_module
    integer, intent(in) :: ng

    real (kind = dp_t), intent(in) ::  ss(0:,:,:)
    real (kind = dp_t), intent(in) ::  ff(:,:)
    integer           , intent(in) ::  mm(:,:)
    real (kind = dp_t), intent(inout) ::  uu(1-ng:,1-ng:)
    real (kind = dp_t) :: wrk(size(ff,1),size(ff,2))
    integer :: nx, ny
    integer :: i, j
    real (kind = dp_t) :: dd
    integer, parameter ::  XBC = 5, YBC = 6

    real (kind = dp_t), parameter :: omega = 4.0_dp_t/5.0_dp_t

    type(bl_prof_timer), save :: bpt

    call build(bpt, "jac_smoother_2d")

    nx = size(ff,dim=1)
    ny = size(ff,dim=2)

    do j = 1, ny
       do i = 1, nx
          dd =    + ss(1,i,j) * uu(i+1,j) + ss(2,i,j) * uu(i-1,j)
          dd = dd + ss(3,i,j) * uu(i,j+1) + ss(4,i,j) * uu(i,j-1)
          if ( nx > 1 ) then
             if (bc_skewed(mm(i,j),1,+1)) then
                dd = dd + ss(XBC,i,j)*uu(i+2,j)
             else if (bc_skewed(mm(i,j),1,-1)) then
                dd = dd + ss(XBC,i,j)*uu(i-2,j)
             end if
          end if
          if ( ny > 1 ) then
             if (bc_skewed(mm(i,j),2,+1)) then
                dd = dd + ss(YBC,i,j)*uu(i,j+2)
             else if (bc_skewed(mm(i,j),2,-1)) then
                dd = dd + ss(YBC,i,j)*uu(i,j-2)
             end if
          end if
          wrk(i,j) = ff(i,j) - dd
       end do
    end do

    do j = 1, ny
       do i = 1, nx
          uu(i,j) = uu(i,j) + omega*(wrk(i,j)/ss(0,i,j)-uu(i,j))
       end do
    end do

    call destroy(bpt)

  end subroutine jac_smoother_2d

  ! subroutine jac_dense_smoother_2d(ss, uu, ff, ng)
  !   use bl_prof_module
  !   integer, intent(in) :: ng
  !   real (kind = dp_t), intent(in) :: ss(0:,:,:)
  !   real (kind = dp_t), intent(in) :: ff(:,:)
  !   real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:)
  !   real (kind = dp_t) :: wrk(size(ff,1),size(ff,2))
  !   integer :: i, j
  !   integer :: nx, ny

  !   real (kind = dp_t), parameter :: omega = 4.0_dp_t/5.0_dp_t

  !   type(bl_prof_timer), save :: bpt

  !   call build(bpt, "jac_dense_smoother_2d")

  !   nx=size(ss,dim=2)
  !   ny=size(ss,dim=3)

  !   do j = 1, ny
  !      do i = 1, nx
  !         wrk(i,j) = ff(i,j) - (&
  !              + ss(1,i,j)*uu(i-1,j-1) &
  !              + ss(2,i,j)*uu(i  ,j-1) &
  !              + ss(3,i,j)*uu(i+1,j-1) &
               
  !              + ss(4,i,j)*uu(i-1,j  ) &
  !              + ss(5,i,j)*uu(i+1,j  ) &
               
  !              + ss(6,i,j)*uu(i-1,j+1) &
  !              + ss(7,i,j)*uu(i  ,j+1) &
  !              + ss(8,i,j)*uu(i+1,j+1) &
  !              )
  !      end do
  !   end do

  !   do j = 1, ny
  !      do i = 1, nx
  !         uu(i,j) = uu(i,j) + omega*(wrk(i,j)/ss(0,i,j)-uu(i,j))
  !      end do
  !   end do

  !   call destroy(bpt)

  ! end subroutine jac_dense_smoother_2d

  subroutine jac_smoother_3d(ss, uu, ff, mm, ng)
    use bl_prof_module
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) :: ss(0:,:,:,:)
    real (kind = dp_t), intent(in) :: ff(:,:,:)
    integer           , intent(in) :: mm(:,:,:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), allocatable :: wrk(:,:,:)
    real (kind = dp_t) :: dd
    integer :: nx, ny, nz
    integer :: i, j, k
    integer, parameter ::  XBC = 7, YBC = 8, ZBC = 9

    real (kind = dp_t), parameter :: omega = 6.0_dp_t/7.0_dp_t

    type(bl_prof_timer), save :: bpt

    call build(bpt, "jac_smoother_3d")

    nx = size(ff,dim=1)
    ny = size(ff,dim=2)
    nz = size(ff,dim=3)

    allocate(wrk(nx,ny,nz))

    !$OMP PARALLEL DO PRIVATE(j,i,k,dd) IF(nz.ge.4)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             dd = ss(1,i,j,k)*uu(i+1,j,k) + ss(2,i,j,k)*uu(i-1,j,k) + &
                  ss(3,i,j,k)*uu(i,j+1,k) + ss(4,i,j,k)*uu(i,j-1,k) + &
                  ss(5,i,j,k)*uu(i,j,k+1) + ss(6,i,j,k)*uu(i,j,k-1)

             if ( nx > 1 ) then
                if (bc_skewed(mm(i,j,k),1,+1)) then
                   dd = dd + ss(XBC,i,j,k)*uu(i+2,j,k)
                else if (bc_skewed(mm(i,j,k),1,-1)) then
                   dd = dd + ss(XBC,i,j,k)*uu(i-2,j,k)
                end if
             end if

             if ( ny > 1 ) then
                if (bc_skewed(mm(i,j,k),2,+1)) then
                   dd = dd + ss(YBC,i,j,k)*uu(i,j+2,k)
                else if (bc_skewed(mm(i,j,k),2,-1)) then
                   dd = dd + ss(YBC,i,j,k)*uu(i,j-2,k)
                end if
             end if

             if ( nz > 1 ) then
                if (bc_skewed(mm(i,j,k),3,+1)) then
                   dd = dd + ss(ZBC,i,j,k)*uu(i,j,k+2)
                else if (bc_skewed(mm(i,j,k),3,-1)) then
                   dd = dd + ss(ZBC,i,j,k)*uu(i,j,k-2)
                end if
             end if
             wrk(i,j,k) = ff(i,j,k) - dd
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(j,i,k) IF(nz.ge.4)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             uu(i,j,k) = uu(i,j,k) + omega*(wrk(i,j,k)/ss(0,i,j,k) - uu(i,j,k))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    call destroy(bpt)

  end subroutine jac_smoother_3d

  ! subroutine jac_dense_smoother_3d(ss, uu, ff, ng)
  !   use bl_prof_module
  !   integer, intent(in) :: ng
  !   real (kind = dp_t), intent(in) :: ss(0:,:,:,:)
  !   real (kind = dp_t), intent(in) :: ff(:,:,:)
  !   real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
  !   real (kind = dp_t), allocatable :: wrk(:,:,:)
  !   integer :: nx, ny, nz
  !   integer :: i, j, k

  !   real (kind = dp_t), parameter :: omega = 6.0_dp_t/7.0_dp_t

  !   type(bl_prof_timer), save :: bpt

  !   call build(bpt, "jac_dense_smoother_3d")

  !   nx = size(ff,dim=1)
  !   ny = size(ff,dim=2)
  !   nz = size(ff,dim=3)

  !   allocate(wrk(nx,ny,nz))

  !   !$OMP PARALLEL DO PRIVATE(j,i,k) IF(nz.ge.4)
  !   do k = 1, nz
  !      do j = 1, ny
  !         do i = 1, nx
  !            wrk(i,j,k) = ff(i,j,k) - ( &
  !                 + ss( 0,i,j,k)*uu(i  ,j  ,k)  &
  !                 + ss( 1,i,j,k)*uu(i-1,j-1,k-1) &
  !                 + ss( 2,i,j,k)*uu(i  ,j-1,k-1) &
  !                 + ss( 3,i,j,k)*uu(i+1,j-1,k-1) &
  !                 + ss( 4,i,j,k)*uu(i-1,j  ,k-1) &
  !                 + ss( 5,i,j,k)*uu(i  ,j  ,k-1) &
  !                 + ss( 6,i,j,k)*uu(i+1,j  ,k-1) &
  !                 + ss( 7,i,j,k)*uu(i-1,j+1,k-1) &
  !                 + ss( 8,i,j,k)*uu(i  ,j+1,k-1) &
  !                 + ss( 9,i,j,k)*uu(i+1,j+1,k-1) &
                  
  !                 + ss(10,i,j,k)*uu(i-1,j-1,k  ) &
  !                 + ss(11,i,j,k)*uu(i  ,j-1,k  ) &
  !                 + ss(12,i,j,k)*uu(i+1,j-1,k  ) &
  !                 + ss(13,i,j,k)*uu(i-1,j  ,k  ) &
  !                 + ss(14,i,j,k)*uu(i+1,j  ,k  ) &
  !                 + ss(15,i,j,k)*uu(i-1,j+1,k  ) &
  !                 + ss(16,i,j,k)*uu(i  ,j+1,k  ) &
  !                 + ss(17,i,j,k)*uu(i+1,j+1,k  ) &
                  
  !                 + ss(18,i,j,k)*uu(i-1,j-1,k+1) &
  !                 + ss(19,i,j,k)*uu(i  ,j-1,k+1) &
  !                 + ss(20,i,j,k)*uu(i+1,j-1,k+1) &
  !                 + ss(21,i,j,k)*uu(i-1,j  ,k+1) &
  !                 + ss(22,i,j,k)*uu(i  ,j  ,k+1) &
  !                 + ss(23,i,j,k)*uu(i+1,j  ,k+1) &
  !                 + ss(24,i,j,k)*uu(i-1,j+1,k+1) &
  !                 + ss(25,i,j,k)*uu(i  ,j+1,k+1) &
  !                 + ss(26,i,j,k)*uu(i+1,j+1,k+1) &
  !                 )
  !         end do
  !      end do
  !   end do
  !   !$OMP END PARALLEL DO

  !   !$OMP PARALLEL DO PRIVATE(j,i,k) IF(nz.ge.4)
  !   do k = 1, nz
  !      do j = 1, ny
  !         do i = 1, nx
  !            uu(i,j,k) = uu(i,j,k) + omega*(wrk(i,j,k)/ss(14,i,j,k)-uu(i,j,k))
  !         end do
  !      end do
  !   end do
  !   !$OMP END PARALLEL DO

  !   call destroy(bpt)

  ! end subroutine jac_dense_smoother_3d

  subroutine gs_rb_smoother_ibc_2d(ss, uu, ff, lo, ng, n)
    use bl_prof_module
    integer, intent(in) :: ng
    integer, intent(in) :: lo(:)
    integer, intent(in) :: n
    real (kind = dp_t), intent(in   ) :: ff(lo(1)   :, lo(2)   :)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:, lo(2)-ng:)
    real (kind = dp_t), intent(in   ) :: ss(0:)
    integer :: j, i, hi(size(lo)), ioff
    real (kind = dp_t) :: dd, ss0inv

    type(bl_prof_timer), save :: bpt

    call build(bpt, "gs_rb_smoother_ibc_2d")

    hi = ubound(ff)

    ss0inv = one/ss(0)

    do j = lo(2),hi(2)
       ioff = 0; if ( mod(lo(1) + j, 2) /= n ) ioff = 1
       do i = lo(1) + ioff, hi(1), 2
          dd =   ss(1) * (uu(i-1,j) + uu(i+1,j)) &
               + ss(2) * (uu(i,j-1) + uu(i,j+1))
          uu(i,j) = ss0inv*(ff(i,j) - dd)
       end do
    end do

    call destroy(bpt)

  end subroutine gs_rb_smoother_ibc_2d

  subroutine gs_rb_smoother_ibc_3d(ss, uu, ff, lo, tlo, thi, ng, n)
    use bl_prof_module
    integer, intent(in) :: ng
    integer, intent(in) :: lo(:), tlo(:), thi(:)
    integer, intent(in) :: n
    real (kind = dp_t), intent(in   ) :: ff(lo(1)   :, lo(2)   :, lo(3)   :)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:)
    real (kind = dp_t), intent(in   ) :: ss(0:)
    integer :: i, j, k, hi(size(lo)), ioff
    real (kind = dp_t) :: dd, ss0inv

    real(dp_t), parameter :: omega = 1.15_dp_t

    type(bl_prof_timer), save :: bpt

    call build(bpt, "gs_rb_smoother_ibc_3d")

    hi = ubound(ff)

    ss0inv = omega/ss(0)

    do k = tlo(3), thi(3)
       do j = tlo(2),thi(2)
          ioff = 0; if ( mod (tlo(1) + j + k, 2) /= n ) ioff = 1
          do i = tlo(1) + ioff, thi(1), 2
             dd =   ss(1) * (uu(i-1,j,k) + uu(i+1,j,k)) &
                  + ss(2) * (uu(i,j-1,k) + uu(i,j+1,k)) &
                  + ss(3) * (uu(i,j,k-1) + uu(i,j,k+1))
             uu(i,j,k) = uu(i,j,k)*(one-omega) + ss0inv*(ff(i,j,k) - dd)
          end do
       end do
    end do

    call destroy(bpt)

  end subroutine gs_rb_smoother_ibc_3d

  subroutine jac_smoother_ibc_2d(ss, uu, ff, ng)
    use bl_prof_module
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) ::  ss(0:)
    real (kind = dp_t), intent(in) ::  ff(:,:)
    real (kind = dp_t), intent(inout) ::  uu(1-ng:,1-ng:)
    integer :: nx, ny
    integer :: i, j
    real (kind = dp_t) :: dd, ss0inv

    real (kind = dp_t), parameter :: omega = 4.0_dp_t/5.0_dp_t

    type(bl_prof_timer), save :: bpt

    call build(bpt, "jac_smoother_ibc_2d")

    nx = size(ff,dim=1)
    ny = size(ff,dim=2)

    ss0inv = omega/ss(0)

    do j = 1, ny
       do i = 1, nx
          dd =   ss(1)*(uu(i-1,j) + uu(i+1,j)) &
               + ss(2)*(uu(i,j-1) + uu(i,j+1))
          uu(i,j) = uu(i,j)*(one-omega) + (ff(i,j)-dd)*ss0inv
       end do
    end do

    call destroy(bpt)

  end subroutine jac_smoother_ibc_2d

  subroutine jac_smoother_ibc_3d(ss, uu, ff, ng)
    use bl_prof_module
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) :: ss(0:)
    real (kind = dp_t), intent(in) :: ff(:,:,:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t) :: dd, ss0inv
    integer :: nx, ny, nz
    integer :: i, j, k

    real (kind = dp_t), parameter :: omega = 6.0_dp_t/7.0_dp_t

    type(bl_prof_timer), save :: bpt

    call build(bpt, "jac_smoother_ibc_3d")

    nx = size(ff,dim=1)
    ny = size(ff,dim=2)
    nz = size(ff,dim=3)

    ss0inv = omega/ss(0)

    !$OMP PARALLEL DO PRIVATE(i,j,k,dd) collapse(2)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             dd =   ss(1)*(uu(i-1,j,k)+uu(i+1,j,k))  &
                  + ss(2)*(uu(i,j-1,k)+uu(i,j+1,k))  &
                  + ss(3)*(uu(i,j,k-1)+uu(i,j,k+1))
             uu(i,j,k) = uu(i,j,k)*(one-omega) + (ff(i,j,k) - dd)*ss0inv
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    call destroy(bpt)

  end subroutine jac_smoother_ibc_3d

end module cc_smoothers_module
