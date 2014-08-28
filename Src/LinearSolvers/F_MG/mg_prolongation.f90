module mg_prolongation_module

  use bl_constants_module
  use bl_types
  use box_module

  implicit none

  interface cc_prolongation
     module procedure cc_prolongation_1d
     module procedure cc_prolongation_2d
     module procedure cc_prolongation_3d
  end interface

  interface pc_c_prolongation
     module procedure pc_c_prolongation_1d
     module procedure pc_c_prolongation_2d
     module procedure pc_c_prolongation_3d
  end interface

  interface lin_c_prolongation
     module procedure lin_c_prolongation_1d
     module procedure lin_c_prolongation_2d
     module procedure lin_c_prolongation_3d
  end interface

  interface nodal_prolongation
     module procedure nodal_prolongation_1d
     module procedure nodal_prolongation_2d
     module procedure nodal_prolongation_3d
  end interface

  interface nodal_linear_prolongation
     module procedure nodal_linear_prolongation_1d
     module procedure nodal_linear_prolongation_2d
     module procedure nodal_linear_prolongation_3d
  end interface

  interface nodal_cubic_prolongation
     module procedure nodal_cubic_prolongation_1d
     module procedure nodal_cubic_prolongation_2d
     module procedure nodal_cubic_prolongation_3d
  end interface

  private :: cubicInterpolate, bicubicInterpolate, tricubicInterpolate
  !
  ! Do you want to use [bi,tri]cubic nodal prolongation instead of linear?
  !
  logical, private, save :: Use_Nodal_Cubic = .false.

contains
  !
  ! Is nodal cubic enabled?
  !
  pure function using_nodal_cubic () result(r)
    logical r
    r = Use_Nodal_Cubic
  end function using_nodal_cubic
  !
  ! Enable/Disable cubic stuff returning previous setting.
  !
  function set_using_nodal_cubic (val) result(r)
    logical, intent(in) :: val
    logical r, old_val
    old_val = Use_Nodal_Cubic
    Use_Nodal_Cubic = val
    r = old_val
  end function set_using_nodal_cubic

  subroutine cc_prolongation_1d(ff, lof, cc, loc, lo, hi, ir, lininterp)
    integer,     intent(in   ) :: loc(:),lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):)
    real (dp_t), intent(in   ) :: cc(loc(1):)
    integer,     intent(in   ) :: ir(:)
    logical,     intent(in   ) :: lininterp

    if ( lininterp ) then
       call lin_c_prolongation_1d(ff, lof, cc, loc, lo, hi, ir)
    else
       call pc_c_prolongation_1d(ff, lof, cc, loc, lo, hi, ir)
    end if

  end subroutine cc_prolongation_1d

  subroutine cc_prolongation_2d(ff, lof, cc, loc, lo, hi, ir, lininterp, ptype)
    integer,     intent(in   ) :: loc(:), lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):)
    real (dp_t), intent(in   ) :: cc(loc(1):,loc(2):)
    integer,     intent(in   ) :: ir(:), ptype
    logical,     intent(in   ) :: lininterp

    if ( lininterp ) then
       call lin_c_prolongation_2d(ff, lof, cc, loc, lo, hi, ir, ptype)
    else
       call pc_c_prolongation_2d(ff, lof, cc, loc, lo, hi, ir)
    end if

  end subroutine cc_prolongation_2d

  subroutine cc_prolongation_3d(ff, lof, cc, loc, lo, hi, ir, lininterp, ptype)
    integer,     intent(in   ) :: loc(:),lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):,lof(3):)
    real (dp_t), intent(in   ) :: cc(loc(1):,loc(2):,loc(3):)
    integer,     intent(in   ) :: ir(:), ptype
    logical,     intent(in   ) :: lininterp

    if ( lininterp ) then
       call lin_c_prolongation_3d(ff, lof, cc, loc, lo, hi, ir, ptype)
    else
       call pc_c_prolongation_3d(ff, lof, cc, loc, lo, hi, ir)
    end if

  end subroutine cc_prolongation_3d

  subroutine pc_c_prolongation_1d(ff, lof, cc, loc, lo, hi, ir)
    integer,     intent(in   ) :: loc(:),lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):)
    real (dp_t), intent(in   ) :: cc(loc(1):)
    integer,     intent(in   ) :: ir(:)

    integer :: i, ic

    do i = lo(1),hi(1)
       ic = i / ir(1) 
       ff(i) = ff(i) + cc(ic)
    end do

  end subroutine pc_c_prolongation_1d

  subroutine pc_c_prolongation_2d(ff, lof, cc, loc, lo, hi, ir)
    integer,     intent(in   ) :: loc(:), lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):)
    real (dp_t), intent(in   ) :: cc(loc(1):,loc(2):)
    integer,     intent(in   ) :: ir(:)

    integer :: i, j, ic, jc, twoi, twoj, twoip1, twojp1, clo(2), chi(2)

    if ( ir(1) == 2 .and. ir(2) == 2 ) then
       !
       ! Specialized unrolled 2D version.
       !
       clo = int_coarsen(lo, 2)
       chi = int_coarsen(hi, 2)

       do j = clo(2),chi(2)
             twoj   = 2*j
             twojp1 = twoj+1
          do i = clo(1),chi(1)
             twoi   = 2*i
             twoip1 = twoi+1

             ff(twoi,   twoj  ) = ff(twoi,   twoj  ) + cc(i,j)
             ff(twoip1, twoj  ) = ff(twoip1, twoj  ) + cc(i,j)
             ff(twoi,   twojp1) = ff(twoi,   twojp1) + cc(i,j)
             ff(twoip1, twojp1) = ff(twoip1, twojp1) + cc(i,j)
          end do
       end do

    else
       !
       ! Generic ir/=2 case.
       !
       do j = lo(2),hi(2)
          jc = j / ir(2)
          do i = lo(1),hi(1)
             ic = i / ir(1) 
             ff(i,j) = ff(i,j) + cc(ic,jc)
          end do
       end do

    end if

  end subroutine pc_c_prolongation_2d

  subroutine pc_c_prolongation_3d(ff, lof, cc, loc, lo, hi, ir)
    integer,     intent(in   ) :: loc(:),lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):,lof(3):)
    real (dp_t), intent(in   ) :: cc(loc(1):,loc(2):,loc(3):)
    integer,     intent(in   ) :: ir(:)

    integer :: i, j, k, ic, jc, kc
    integer :: clo(3), chi(3), twoi, twoj, twok, twoip1, twojp1, twokp1

    if ( ir(1) == 2 .and. ir(2) == 2 .and. ir(3) == 2 ) then
       !
       ! Specialized unrolled 3D version.
       !
       clo = int_coarsen(lo, 2)
       chi = int_coarsen(hi, 2)

       do k = clo(3),chi(3)
          twok   = 2*k
          twokp1 = twok+1
          do j = clo(2),chi(2)
             twoj   = 2*j
             twojp1 = twoj+1
             do i = clo(1),chi(1)
                twoi   = 2*i
                twoip1 = twoi+1

                ff(twoip1, twojp1, twokp1) = ff(twoip1, twojp1, twokp1) + cc(i,j,k)
                ff(twoi,   twojp1, twokp1) = ff(twoi,   twojp1, twokp1) + cc(i,j,k)
                ff(twoip1, twoj,   twokp1) = ff(twoip1, twoj,   twokp1) + cc(i,j,k)
                ff(twoi,   twoj,   twokp1) = ff(twoi,   twoj,   twokp1) + cc(i,j,k)
                ff(twoip1, twojp1, twok  ) = ff(twoip1, twojp1, twok  ) + cc(i,j,k)
                ff(twoi,   twojp1, twok  ) = ff(twoi,   twojp1, twok  ) + cc(i,j,k)
                ff(twoip1, twoj,   twok  ) = ff(twoip1, twoj,   twok  ) + cc(i,j,k)
                ff(twoi,   twoj,   twok  ) = ff(twoi,   twoj,   twok  ) + cc(i,j,k)
             end do
          end do
       end do

    else
       !
       ! General ir/=2 case.
       !
       do k = lo(3),hi(3)
          kc = k / ir(3)
          do j = lo(2),hi(2)
             jc = j / ir(2)
             do i = lo(1),hi(1)
                ic = i / ir(1)
                ff(i,j,k) = ff(i,j,k) + cc(ic,jc,kc)
             end do
          end do
       end do

    end if

  end subroutine pc_c_prolongation_3d

  subroutine lin_c_prolongation_1d(ff, lof, cc, loc, lo, hi, ir)
    use bl_error_module
    integer,     intent(in   ) :: loc(:),lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):)
    real (dp_t), intent(in   ) :: cc(loc(1):)
    integer,     intent(in   ) :: ir(:)

    integer :: i, ng, clo(1), chi(1)

    if ( ir(1) == 2 ) then

       clo = int_coarsen(lo, 2)
       chi = int_coarsen(hi, 2)

       ng = (clo(1)-loc(1))

       if ( ng == 0 ) then
          !
          ! Don't have any grow cells.  Do piecewise constant at boundaries.
          !
          i = clo(1)
          ff(2*i  ) = ff(2*i  ) + cc(i)
          ff(2*i+1) = ff(2*i+1) + cc(i)

          i = chi(1)
          ff(2*i  ) = ff(2*i  ) + cc(i)
          ff(2*i+1) = ff(2*i+1) + cc(i)
          !
          ! Twiddle clo/chi so we linear interp only on the interior.
          !
          clo = clo + 1
          chi = chi - 1
       end if

       do i = clo(1),chi(1)
          ff(2*i  ) = ff(2*i  ) + FOURTH*( 3*cc(i) + cc(i-1) )
          ff(2*i+1) = ff(2*i+1) + FOURTH*( 3*cc(i) + cc(i+1) )
       end do

    else
       !
       ! For now don't bother with a generic lininterp.
       !

       call pc_c_prolongation(ff, lof, cc, loc, lo, hi, ir)
    end if

  end subroutine lin_c_prolongation_1d

  subroutine lin_c_prolongation_2d(ff, lof, cc, loc, lo, hi, ir, ptype)
    use bl_error_module
    integer,     intent(in   ) :: loc(:), lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):)
    real (dp_t), intent(in   ) :: cc(loc(1):,loc(2):)
    integer,     intent(in   ) :: ir(:), ptype

    integer :: i, j, ng, clo(2), chi(2), twoi, twoj, twoip1, twojp1
    logical :: interior_i, interior_j

    real (dp_t), parameter :: one16th = 1.0_dp_t /16.0_dp_t

    if ( ir(1) == 2 .and. ir(2) == 2 ) then

       clo = int_coarsen(lo, 2)
       chi = int_coarsen(hi, 2)

       ng = min(clo(1)-loc(1),clo(2)-loc(2))

       if ( ng == 0 ) then
          !
          ! Don't have any grow cells.  Do piecewise constant at boundaries.
          !
          do j = clo(2),chi(2)
             twoj   = 2*j
             twojp1 = twoj+1
             interior_j = ( j > clo(2) .and. j < chi(2) )

             do i = clo(1),chi(1)
                interior_i = ( i > clo(1) .and. i < chi(1) )

                if ( interior_i .and. interior_j ) cycle

                twoi   = 2*i
                twoip1 = twoi+1

                ff(twoi,   twoj  ) = ff(twoi,   twoj  ) + cc(i,j)
                ff(twoip1, twoj  ) = ff(twoip1, twoj  ) + cc(i,j)
                ff(twoi,   twojp1) = ff(twoi,   twojp1) + cc(i,j)
                ff(twoip1, twojp1) = ff(twoip1, twojp1) + cc(i,j)
             end do
          end do
          !
          ! Twiddle clo/chi so we linear interp only on the interior.
          !
          clo = clo + 1
          chi = chi - 1
       end if

       select case ( ptype )
       case ( 1 )
          !
          ! Type 1 - {2,1,1} weighted average of our neighbors.  Doesn't use corner cells.
          !
          do j = clo(2), chi(2)
             twoj   = 2*j
             twojp1 = twoj+1

             do i = clo(1), chi(1)
                twoi   = 2*i
                twoip1 = twoi+1

                ff(twoip1, twojp1) = ff(twoip1, twojp1) + FOURTH * ( 2*cc(i,j) + cc(i+1,j) + cc(i,j+1) )
                ff(twoi,   twojp1) = ff(twoi,   twojp1) + FOURTH * ( 2*cc(i,j) + cc(i-1,j) + cc(i,j+1) )
                ff(twoip1, twoj  ) = ff(twoip1, twoj  ) + FOURTH * ( 2*cc(i,j) + cc(i+1,j) + cc(i,j-1) )
                ff(twoi,   twoj  ) = ff(twoi,   twoj  ) + FOURTH * ( 2*cc(i,j) + cc(i-1,j) + cc(i,j-1) )
             end do
          end do
       case ( 2 )
          !
          ! Type 2 - {4,1,1} weighted average of our neighbors.  Doesn't use corner cells.
          !
          do j = clo(2), chi(2)
             twoj   = 2*j
             twojp1 = twoj+1

             do i = clo(1), chi(1)
                twoi   = 2*i
                twoip1 = twoi+1

                ff(twoip1, twojp1) = ff(twoip1, twojp1) + SIXTH * ( 4*cc(i,j) + cc(i+1,j) + cc(i,j+1) )
                ff(twoi,   twojp1) = ff(twoi,   twojp1) + SIXTH * ( 4*cc(i,j) + cc(i-1,j) + cc(i,j+1) )
                ff(twoip1, twoj  ) = ff(twoip1, twoj  ) + SIXTH * ( 4*cc(i,j) + cc(i+1,j) + cc(i,j-1) )
                ff(twoi,   twoj  ) = ff(twoi,   twoj  ) + SIXTH * ( 4*cc(i,j) + cc(i-1,j) + cc(i,j-1) )
             end do
          end do
       case ( 3 )
          !
          ! Type 3 - bilinear.  Uses corner cells.
          !
          do j = clo(2), chi(2)
             twoj   = 2*j
             twojp1 = twoj+1

             do i = clo(1), chi(1)
                twoi   = 2*i
                twoip1 = twoi+1

                ff(twoip1, twojp1) = ff(twoip1, twojp1) + one16th * &
                     ( 9*cc(i,j) + 3*cc(i+1,j) + 3*cc(i,j+1) + cc(i+1,j+1) )
                ff(twoi,   twojp1) = ff(twoi,   twojp1) + one16th * &
                     ( 9*cc(i,j) + 3*cc(i-1,j) + 3*cc(i,j+1) + cc(i-1,j+1) )
                ff(twoip1, twoj  ) = ff(twoip1, twoj  ) + one16th * &
                     ( 9*cc(i,j) + 3*cc(i,j-1) + 3*cc(i+1,j) + cc(i+1,j-1) )
                ff(twoi,   twoj  ) = ff(twoi,   twoj  ) + one16th * &
                     ( 9*cc(i,j) + 3*cc(i-1,j) + 3*cc(i,j-1) + cc(i-1,j-1) )
             end do
          end do
       case default
          call bl_error("lin_c_prolongation_2d: unknown ptype", ptype)
       end select

    else
       !
       ! For now don't bother with a generic lininterp.
       !

       call pc_c_prolongation(ff, lof, cc, loc, lo, hi, ir)
    end if

  end subroutine lin_c_prolongation_2d

  subroutine lin_c_prolongation_3d(ff, lof, cc, loc, lo, hi, ir, ptype)
    use bl_error_module
    integer,     intent(in   ) :: loc(:),lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):,lof(3):)
    real (dp_t), intent(in   ) :: cc(loc(1):,loc(2):,loc(3):)
    integer,     intent(in   ) :: ir(:), ptype

    integer :: i, j, k, ng, clo(3), chi(3), twoi, twoj, twoip1, twojp1, twok, twokp1
    logical :: interior_i, interior_j, interior_k

    real (dp_t), parameter ::   ONE64TH  = 1.0_dp_t / 64.0_dp_t
    real (dp_t), parameter :: THREE64THS = 3.0_dp_t / 64.0_dp_t

    if ( ir(1) == 2 .and. ir(2) == 2 .and. ir(3) == 2 ) then
       !
       ! First do all face points using piecewise-constant interpolation.
       !
       clo = int_coarsen(lo, 2)
       chi = int_coarsen(hi, 2)

       ng = min(clo(1)-loc(1),clo(2)-loc(2),clo(3)-loc(3))

       if ( ng == 0 ) then
          !
          ! Don't have any grow cells.  Do piecewise constant at boundaries.
          !
          do k = clo(3),chi(3)
             twok   = 2*k
             twokp1 = twok+1
             interior_k = ( k > clo(3) .and. k < chi(3) )

             do j = clo(2),chi(2)
                twoj   = 2*j
                twojp1 = twoj+1
                interior_j = ( j > clo(2) .and. j < chi(2) )

                do i = clo(1),chi(1)
                   interior_i = ( i > clo(1) .and. i < chi(1) )

                   if ( interior_i .and. interior_j .and. interior_k ) cycle

                   twoi   = 2*i
                   twoip1 = twoi+1

                   ff(twoip1, twojp1, twokp1) = ff(twoip1, twojp1, twokp1) + cc(i,j,k)
                   ff(twoi,   twojp1, twokp1) = ff(twoi,   twojp1, twokp1) + cc(i,j,k)
                   ff(twoip1, twoj,   twokp1) = ff(twoip1, twoj,   twokp1) + cc(i,j,k)
                   ff(twoi,   twoj,   twokp1) = ff(twoi,   twoj,   twokp1) + cc(i,j,k)
                   ff(twoip1, twojp1, twok  ) = ff(twoip1, twojp1, twok  ) + cc(i,j,k)
                   ff(twoi,   twojp1, twok  ) = ff(twoi,   twojp1, twok  ) + cc(i,j,k)
                   ff(twoip1, twoj,   twok  ) = ff(twoip1, twoj,   twok  ) + cc(i,j,k)
                   ff(twoi,   twoj,   twok  ) = ff(twoi,   twoj,   twok  ) + cc(i,j,k)
                end do
             end do
          end do
          !
          ! Twiddle clo/chi so we linear interp only on the interior.
          !
          clo = clo + 1
          chi = chi - 1
       end if

       select case ( ptype )
       case ( 1 )
          !
          ! Type 1 - {1,1,1,1} weighted average of our neighbors.  Doesn't use corners.
          !
          do k = clo(3),chi(3)
             twok   = 2*k
             twokp1 = twok+1

             do j = clo(2),chi(2)
                twoj   = 2*j
                twojp1 = twoj+1

                do i = clo(1),chi(1)
                   twoi   = 2*i
                   twoip1 = twoi+1

                   ff(twoip1, twojp1, twokp1) = ff(twoip1, twojp1, twokp1) + &
                        FOURTH * ( cc(i,j,k) + cc(i+1,j,k) + cc(i,j+1,k) + cc(i,j,k+1) )
                   ff(twoi,   twojp1, twokp1) = ff(twoi,   twojp1, twokp1) + &
                        FOURTH * ( cc(i,j,k) + cc(i-1,j,k) + cc(i,j+1,k) + cc(i,j,k+1) )
                   ff(twoip1, twoj,   twokp1) = ff(twoip1, twoj,   twokp1) + &
                        FOURTH * ( cc(i,j,k) + cc(i+1,j,k) + cc(i,j-1,k) + cc(i,j,k+1) )
                   ff(twoi,   twoj,   twokp1) = ff(twoi,   twoj,   twokp1) + &
                        FOURTH * ( cc(i,j,k) + cc(i-1,j,k) + cc(i,j-1,k) + cc(i,j,k+1) )
                   ff(twoip1, twojp1, twok  ) = ff(twoip1, twojp1, twok  ) + &
                        FOURTH * ( cc(i,j,k) + cc(i+1,j,k) + cc(i,j+1,k) + cc(i,j,k-1) )
                   ff(twoi,   twojp1, twok  ) = ff(twoi,   twojp1, twok  ) + &
                        FOURTH * ( cc(i,j,k) + cc(i-1,j,k) + cc(i,j+1,k) + cc(i,j,k-1) )
                   ff(twoip1, twoj,   twok  ) = ff(twoip1, twoj,   twok  ) + &
                        FOURTH * ( cc(i,j,k) + cc(i+1,j,k) + cc(i,j-1,k) + cc(i,j,k-1) )
                   ff(twoi,   twoj,   twok  ) = ff(twoi,   twoj,   twok  ) + &
                        FOURTH * ( cc(i,j,k) + cc(i-1,j,k) + cc(i,j-1,k) + cc(i,j,k-1) )
                end do
             end do
          end do
       case ( 2 )
          !
          ! Type 2 - {3,1,1,1} weighted average of our neighbors.  Doesn't use corners.
          !
          do k = clo(3),chi(3)
             twok   = 2*k
             twokp1 = twok+1

             do j = clo(2),chi(2)
                twoj   = 2*j
                twojp1 = twoj+1

                do i = clo(1),chi(1)
                   twoi   = 2*i
                   twoip1 = twoi+1

                   ff(twoip1, twojp1, twokp1) = ff(twoip1, twojp1, twokp1) + &
                        SIXTH * ( 3*cc(i,j,k) + cc(i+1,j,k) + cc(i,j+1,k) + cc(i,j,k+1) )
                   ff(twoi,   twojp1, twokp1) = ff(twoi,   twojp1, twokp1) + &
                        SIXTH * ( 3*cc(i,j,k) + cc(i-1,j,k) + cc(i,j+1,k) + cc(i,j,k+1) )
                   ff(twoip1, twoj,   twokp1) = ff(twoip1, twoj,   twokp1) + &
                        SIXTH * ( 3*cc(i,j,k) + cc(i+1,j,k) + cc(i,j-1,k) + cc(i,j,k+1) )
                   ff(twoi,   twoj,   twokp1) = ff(twoi,   twoj,   twokp1) + &
                        SIXTH * ( 3*cc(i,j,k) + cc(i-1,j,k) + cc(i,j-1,k) + cc(i,j,k+1) )
                   ff(twoip1, twojp1, twok  ) = ff(twoip1, twojp1, twok  ) + &
                        SIXTH * ( 3*cc(i,j,k) + cc(i+1,j,k) + cc(i,j+1,k) + cc(i,j,k-1) )
                   ff(twoi,   twojp1, twok  ) = ff(twoi,   twojp1, twok  ) + &
                        SIXTH * ( 3*cc(i,j,k) + cc(i-1,j,k) + cc(i,j+1,k) + cc(i,j,k-1) )
                   ff(twoip1, twoj,   twok  ) = ff(twoip1, twoj,   twok  ) + &
                        SIXTH * ( 3*cc(i,j,k) + cc(i+1,j,k) + cc(i,j-1,k) + cc(i,j,k-1) )
                   ff(twoi,   twoj,   twok  ) = ff(twoi,   twoj,   twok  ) + &
                        SIXTH * ( 3*cc(i,j,k) + cc(i-1,j,k) + cc(i,j-1,k) + cc(i,j,k-1) )
                end do
             end do
          end do
       case ( 3 )
          !
          ! Type 3 - trilinear.  Uses corners.
          !
          do k = clo(3),chi(3)
             twok   = 2*k
             twokp1 = twok+1

             do j = clo(2),chi(2)
                twoj   = 2*j
                twojp1 = twoj+1

                do i = clo(1),chi(1)
                   twoi   = 2*i
                   twoip1 = twoi+1

                   ff(twoip1, twojp1, twokp1) = ff(twoip1, twojp1, twokp1) + &
                        THREE64THS * ( 9*cc(i,j,k  ) + 3*cc(i+1,j,k  ) + 3*cc(i,j+1,k  ) + cc(i+1,j+1,k  ) ) + &
                        ONE64TH * ( 9*cc(i,j,k+1) + 3*cc(i+1,j,k+1) + 3*cc(i,j+1,k+1) + cc(i+1,j+1,k+1) )
                   ff(twoi,   twojp1, twokp1) = ff(twoi,   twojp1, twokp1) + &
                        THREE64THS * ( 9*cc(i,j,k  ) + 3*cc(i-1,j,k  ) + 3*cc(i,j+1,k  ) + cc(i-1,j+1,k  ) ) + &
                        ONE64TH * ( 9*cc(i,j,k+1) + 3*cc(i-1,j,k+1) + 3*cc(i,j+1,k+1) + cc(i-1,j+1,k+1) )
                   ff(twoip1, twoj,   twokp1) = ff(twoip1, twoj,   twokp1) + &
                        THREE64THS * ( 9*cc(i,j,k  ) + 3*cc(i+1,j,k  ) + 3*cc(i,j-1,k  ) + cc(i+1,j-1,k  ) ) + &
                        ONE64TH * ( 9*cc(i,j,k+1) + 3*cc(i+1,j,k+1) + 3*cc(i,j-1,k+1) + cc(i+1,j-1,k+1) )
                   ff(twoi,   twoj,   twokp1) = ff(twoi,   twoj,   twokp1) + &
                        THREE64THS * ( 9*cc(i,j,k  ) + 3*cc(i-1,j,k  ) + 3*cc(i,j-1,k  ) + cc(i-1,j-1,k  ) ) + &
                        ONE64TH * ( 9*cc(i,j,k+1) + 3*cc(i-1,j,k+1) + 3*cc(i,j-1,k+1) + cc(i-1,j-1,k+1) )
                   ff(twoip1, twojp1, twok) = ff(twoip1, twojp1, twok) + &
                        THREE64THS * ( 9*cc(i,j,k  ) + 3*cc(i+1,j,k  ) + 3*cc(i,j+1,k  ) + cc(i+1,j+1,k  ) ) + &
                        ONE64TH * ( 9*cc(i,j,k-1) + 3*cc(i+1,j,k-1) + 3*cc(i,j+1,k-1) + cc(i+1,j+1,k-1) )
                   ff(twoi,   twojp1, twok) = ff(twoi,   twojp1, twok) + &
                        THREE64THS * ( 9*cc(i,j,k  ) + 3*cc(i-1,j,k  ) + 3*cc(i,j+1,k  ) + cc(i-1,j+1,k  ) ) + &
                        ONE64TH * ( 9*cc(i,j,k-1) + 3*cc(i-1,j,k-1) + 3*cc(i,j+1,k-1) + cc(i-1,j+1,k-1) )
                   ff(twoip1, twoj,   twok) = ff(twoip1, twoj,   twok) + &
                        THREE64THS * ( 9*cc(i,j,k  ) + 3*cc(i+1,j,k  ) + 3*cc(i,j-1,k  ) + cc(i+1,j-1,k  ) ) + &
                        ONE64TH * ( 9*cc(i,j,k-1) + 3*cc(i+1,j,k-1) + 3*cc(i,j-1,k-1) + cc(i+1,j-1,k-1) )
                   ff(twoi,   twoj,   twok) = ff(twoi,   twoj,   twok) + &
                        THREE64THS * ( 9*cc(i,j,k  ) + 3*cc(i-1,j,k  ) + 3*cc(i,j-1,k  ) + cc(i-1,j-1,k  ) ) + &
                        ONE64TH * ( 9*cc(i,j,k-1) + 3*cc(i-1,j,k-1) + 3*cc(i,j-1,k-1) + cc(i-1,j-1,k-1) )
                end do
             end do
          end do
       case default
          call bl_error("lin_c_prolongation_3d: unknown ptype", ptype)
       end select

    else
       !
       ! For now don't bother with a generic lininterp.
       !

       call pc_c_prolongation(ff, lof, cc, loc, lo, hi, ir)
    end if

  end subroutine lin_c_prolongation_3d

  pure function cubicInterpolate (p, x) result (r)
    real (dp_t), intent(in) :: p(0:3), x
    real (dp_t) r
    r=p(1)+HALF*x*(p(2)-p(0)+x*(2*p(0)-5*p(1)+4*p(2)-p(3)+x*(3*(p(1)-p(2))+p(3)-p(0))))
  end function cubicInterpolate

  function bicubicInterpolate (p, x, y) result (r)
    real (dp_t), intent(in) :: p(0:15), x, y
    real (dp_t) r, arr(0:3)
    !
    ! First interpolate in x.
    !
    arr(0) = cubicInterpolate(p( 0:3 ), x) ! row 0
    arr(1) = cubicInterpolate(p( 4:7 ), x) ! row 1
    arr(2) = cubicInterpolate(p( 8:11), x) ! row 2
    arr(3) = cubicInterpolate(p(12:15), x) ! row 3
    !
    ! Then using those four points interpolate in y.
    !
    r = cubicInterpolate(arr, y)
  end function bicubicInterpolate

  function tricubicInterpolate (p, x, y, z) result (r)
    real (dp_t), intent(in) :: p(0:63), x, y, z
    real (dp_t) r, arr(0:3)
    arr(0) = bicubicInterpolate(p( 0:15), x, y)
    arr(1) = bicubicInterpolate(p(16:31), x, y)
    arr(2) = bicubicInterpolate(p(32:47), x, y)
    arr(3) = bicubicInterpolate(p(48:63), x, y)
    r = cubicInterpolate(arr, z)
  end function tricubicInterpolate

  subroutine nodal_prolongation_1d(ff, lof, cc, loc, lo, hi, ir)
    integer,     intent(in   ) :: loc(:), lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):)
    real (dp_t), intent(in   ) :: cc(loc(1):)
    integer,     intent(in   ) :: ir(:)

    if ( ir(1) == 2 .and. Use_Nodal_Cubic ) then
       call nodal_cubic_prolongation_1d(ff, lof, cc, loc, lo, hi, ir)
    else
       call nodal_linear_prolongation_1d(ff, lof, cc, loc, lo, hi, ir)
    end if

  end subroutine nodal_prolongation_1d

  subroutine nodal_linear_prolongation_1d(ff, lof, cc, loc, lo, hi, ir)
    integer,     intent(in   ) :: loc(:), lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):)
    real (dp_t), intent(in   ) :: cc(loc(1):)
    integer,     intent(in   ) :: ir(:)

    integer                :: i, ic, l
    real (dp_t)            :: fac_left, fac_rght
    !
    ! Direct injection for fine points overlaying coarse ones.
    !
    do i = lo(1),hi(1),ir(1)
       ic = i / ir(1) 
       ff(i) = ff(i) + cc(ic)
    end do

    do l = 1, ir(1)-1
       fac_rght = real(l,dp_t) / real(ir(1),dp_t)
       fac_left = ONE - fac_rght
       do i = lo(1), hi(1)-1, ir(1)
          ic = i / ir(1) 
          ff(i+l) = ff(i+1) + fac_left*cc(ic) + fac_rght*cc(ic+1)
       end do
    end do

  end subroutine nodal_linear_prolongation_1d

  subroutine nodal_cubic_prolongation_1d(ff, lof, cc, loc, lo, hi, ir)
    integer,     intent(in   ) :: loc(:), lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):)
    real (dp_t), intent(in   ) :: cc(loc(1):)
    integer,     intent(in   ) :: ir(:)

    integer               :: i, ic
    real (dp_t)           :: coeffs(0:3)
    !
    ! Direct injection for fine points overlaying coarse ones.
    !
    do i = lo(1),hi(1),ir(1)
       ic = i / ir(1) 
       ff(i) = ff(i) + cc(ic)
    end do

    if ( ir(1) == 2 ) then
       !
       ! Use linear interpolation at points one off from boundaries.
       !
       i  = lo(1)
       ic = i / 2
       ff(i+1) = HALF * (cc(ic) + cc(ic+1))

       i  = hi(1) -1
       ic = i / 2
       ff(i+1) = HALF * (cc(ic) + cc(ic+1))
       !
       ! Use cubic interpolation only on interior points away from boundary.
       ! This way we don't need to assume we have any grow cells.
       !
       do i = lo(1)+1, hi(1)-2, 2

          ic = i / 2

          coeffs(0) = cc(ic-1)
          coeffs(1) = cc(ic+0)
          coeffs(2) = cc(ic+1)
          coeffs(3) = cc(ic+2)

          ff(i+1) = ff(i+1) + cubicInterpolate(coeffs, HALF)
       end do
    else
       !
       ! Call lininterp in general ir != 2 case.
       !
       call nodal_linear_prolongation_1d(ff, lof, cc, loc, lo, hi, ir)
    end if

  end subroutine nodal_cubic_prolongation_1d

  subroutine nodal_prolongation_2d(ff, lof, cc, loc, lo, hi, ir)
    integer,     intent(in   ) :: loc(:),lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):)
    real (dp_t), intent(in   ) :: cc(loc(1):,loc(2):)
    integer,     intent(in   ) :: ir(:)

    if ( ir(1) == 2 .and. ir(2) == 2 .and. Use_Nodal_Cubic ) then
       call nodal_cubic_prolongation_2d(ff, lof, cc, loc, lo, hi, ir)
    else
       call nodal_linear_prolongation_2d(ff, lof, cc, loc, lo, hi, ir)
    end if

  end subroutine nodal_prolongation_2d

  subroutine nodal_linear_prolongation_2d(ff, lof, cc, loc, lo, hi, ir)
    integer,     intent(in   ) :: loc(:),lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):)
    real (dp_t), intent(in   ) :: cc(loc(1):,loc(2):)
    integer,     intent(in   ) :: ir(:)

    integer     :: i, j, ic, jc, l, m, clo(2), chi(2)
    real (dp_t) :: fac_left, fac_rght
    real (dp_t) :: temp(lo(1):hi(1),lo(2):hi(2))

    if ( ir(1) == 2 .and. ir(2) == 2 ) then
       !
       ! Do fast unrolled linear interp.
       !
       clo = int_coarsen(lo, 2)
       chi = int_coarsen(hi, 2)

       do jc = clo(2),chi(2)
          j = 2*jc
          do ic = clo(1),chi(1)
             i = 2*ic
             !
             ! Direct injection for fine points overlaying coarse ones.
             !
             ff(i,j) = ff(i,j) + cc(ic,jc)

             if ( i < hi(1) ) then
                ff(i+1,j) = ff(i+1,j) + HALF*( cc(ic,jc)+cc(ic+1,jc) )

                if ( j < hi(2) ) then
                   ff(i+1,j+1) = ff(i+1,j+1) + FOURTH*( cc(ic,jc)+cc(ic+1,jc)+cc(ic,jc+1)+cc(ic+1,jc+1) )
                end if
             end if

             if ( j < hi(2) ) then
                ff(i,j+1) = ff(i,j+1) + HALF*( cc(ic,jc)+cc(ic,jc+1) )
             end if
          end do
       end do

    else
       !
       ! General ir/=2 case.
       !
       ! Direct injection for fine points overlaying coarse ones.
       !
       do j = lo(2),hi(2),ir(2)
          jc = j / ir(2) 
          do i = lo(1),hi(1),ir(1)
             ic = i / ir(1) 
             temp(i,j) = cc(ic,jc)
          end do
       end do
       !
       ! Interpolate at fine nodes between coarse nodes in the i-direction only.
       !
       do j = lo(2),hi(2),ir(2)
          do l = 1, ir(1)-1
             fac_rght = real(l,dp_t) / real(ir(1),dp_t)
             fac_left = ONE - fac_rght
             do i = lo(1),hi(1)-1,ir(1)
                temp(i+l,j) = fac_left*temp(i,j) + fac_rght*temp(i+ir(1),j)
             end do
          end do
       end do
       !
       ! Interpolate in the j-direction using previously interpolated "temp".
       !
       do m = 1, ir(2)-1
          fac_rght = real(m,dp_t) / real(ir(2),dp_t)
          fac_left = ONE - fac_rght
          do j = lo(2),hi(2)-1,ir(2)
             do i = lo(1),hi(1)
                temp(i,j+m) = fac_left*temp(i,j)+fac_rght*temp(i,j+ir(2))
             end do
          end do
       end do

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             ff(i,j) = ff(i,j) + temp(i,j)
          end do
       end do

    end if

  end subroutine nodal_linear_prolongation_2d

  recursive subroutine nodal_cubic_prolongation_2d(ff, lof, cc, loc, lo, hi, ir)
    integer,     intent(in   ) :: loc(:),lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):)
    real (dp_t), intent(in   ) :: cc(loc(1):,loc(2):)
    integer,     intent(in   ) :: ir(:)

    integer     :: i, j, ic, jc, m, ng, clo(2), chi(2), istart, jstart, ilo(2), ihi(2)
    real (dp_t) :: coeffs(0:15), ipnt, jpnt

    if ( ir(1) == 2 .and. ir(2) == 2 ) then
       clo = int_coarsen(lo, 2)
       chi = int_coarsen(hi, 2)

       ng = min((clo(1)-loc(1)),(clo(2)-loc(2)))

       if ( ng > 0 ) then
          !
          ! Bicubic was requested and we have one ghost cell.
          !
          do jc = clo(2),chi(2)
             j      = 2*jc
             jpnt   = ZERO
             jstart = jc-1

             if ( j == hi(2) ) then
                jpnt   = ONE
                jstart = jc-2
             end if

             do ic = clo(1),chi(1)
                i      = 2*ic
                ipnt   = ZERO
                istart = ic-1

                if ( i == hi(1) ) then
                   ipnt   = ONE
                   istart = ic-2
                end if
                !
                ! Direct injection for fine points overlaying coarse ones.
                !
                ff(i,j) = ff(i,j) + cc(ic,jc)

                if ( i < hi(1) .or. j < hi(2) ) then
                   coeffs( 0 :  3) = (/ (cc(m,jstart+0), m=istart,istart+3) /)
                   coeffs( 4 :  7) = (/ (cc(m,jstart+1), m=istart,istart+3) /)
                   coeffs( 8 : 11) = (/ (cc(m,jstart+2), m=istart,istart+3) /)
                   coeffs(12 : 15) = (/ (cc(m,jstart+3), m=istart,istart+3) /)
                end if

                if ( i < hi(1) ) then
                   ff(i+1,j) = ff(i+1,j) + bicubicInterpolate(coeffs, HALF, jpnt)

                   if ( j < hi(2) ) then
                      ff(i+1,j+1) = ff(i+1,j+1) + bicubicInterpolate(coeffs, HALF, HALF)
                   end if

                end if

                if ( j < hi(2) ) then
                   ff(i,j+1) = ff(i,j+1) + bicubicInterpolate(coeffs, ipnt, HALF)
                end if

             end do
          end do

       else
          !
          ! First do linear interp
          !
          call nodal_linear_prolongation_2d(ff, lof, cc, loc, lo, hi, ir)
          !
          ! If the fine region is big enough do bicubic on interior.
          !
          ilo = lo+2
          ihi = hi-2

          if ( all( (ihi-ilo) > 0) ) then
             call nodal_cubic_prolongation_2d(ff, lof, cc, loc, ilo, ihi, ir)
          end if
       end if
    else
       !
       ! Call lininterp in general ir != 2 case.
       !
       call nodal_linear_prolongation_2d(ff, lof, cc, loc, lo, hi, ir)
    end if

  end subroutine nodal_cubic_prolongation_2d

  subroutine nodal_prolongation_3d(ff, lof, cc, loc, lo, hi, ir)
    integer,     intent(in   ) :: loc(:), lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):,lof(3):)
    real (dp_t), intent(in   ) :: cc(loc(1):,loc(2):,loc(3):)
    integer,     intent(in   ) :: ir(:)

    if ( ir(1) == 2 .and. ir(2) == 2 .and. ir(3) == 2 .and. Use_Nodal_Cubic ) then
       call nodal_cubic_prolongation_3d(ff, lof, cc, loc, lo, hi, ir)
    else
       call nodal_linear_prolongation_3d(ff, lof, cc, loc, lo, hi, ir)
    endif

  end subroutine nodal_prolongation_3d

  subroutine nodal_linear_prolongation_3d(ff, lof, cc, loc, lo, hi, ir)
    integer,     intent(in   ) :: loc(:), lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):,lof(3):)
    real (dp_t), intent(in   ) :: cc(loc(1):,loc(2):,loc(3):)
    integer,     intent(in   ) :: ir(:)

    integer     :: i, j, k, ic, jc, kc, l, m, n
    integer     :: clo(3), chi(3)
    real (dp_t) :: fac_left, fac_rght

    real (dp_t), allocatable :: temp(:,:,:)

    if ( ir(1) == 2 .and. ir(2) == 2 .and. ir(3) == 2 ) then
       !
       ! Do fast unrolled linear interp.
       !
       clo = int_coarsen(lo, 2)
       chi = int_coarsen(hi, 2)

       do kc = clo(3),chi(3)
          k = 2*kc
          do jc = clo(2),chi(2)
             j = 2*jc
             do ic = clo(1),chi(1)
                i = 2*ic
                !
                ! Direct injection for fine points overlaying coarse ones.
                !
                ff(i,j,k) = ff(i,j,k) + cc(ic,jc,kc)

                if ( i < hi(1) ) then
                   ff(i+1,j,k) = ff(i+1,j,k) + HALF*( cc(ic,jc,kc) + cc(ic+1,jc,kc) )

                   if ( j < hi(2) ) then
                      ff(i+1,j+1,k) = ff(i+1,j+1,k) + &
                           FOURTH*( cc(ic,jc,kc) + cc(ic+1,jc,kc) + cc(ic,jc+1,kc) + cc(ic+1,jc+1,kc) )
                   end if
                end if

                if ( j < hi(2) ) then
                   ff(i,j+1,k) = ff(i,j+1,k) + HALF*( cc(ic,jc,kc) + cc(ic,jc+1,kc) )
                end if

                if ( k < hi(3) ) then
                   ff(i,j,k+1) = ff(i,j,k+1) + HALF*( cc(ic,jc,kc)+cc(ic,jc,kc+1) )

                   if ( i < hi(1) ) then
                      ff(i+1,j,k+1) = ff(i+1,j,k+1) + &
                           FOURTH*( cc(ic,jc,kc) + cc(ic+1,jc,kc) + cc(ic,jc,kc+1) + cc(ic+1,jc,kc+1) )

                      if ( j < hi(2) ) then
                         ff(i+1,j+1,k+1) = ff(i+1,j+1,k+1) + &
                              EIGHTH*( ( cc(ic,jc,kc) + cc(ic+1,jc,kc) + cc(ic,jc+1,kc) + cc(ic+1,jc+1,kc) ) +  &
                              ( cc(ic,jc,kc+1) + cc(ic+1,jc,kc+1) + cc(ic,jc+1,kc+1) + cc(ic+1,jc+1,kc+1) ) )
                      end if
                   end if

                   if ( j < hi(2) ) then
                      ff(i,j+1,k+1) = ff(i,j+1,k+1) + &
                           FOURTH*( cc(ic,jc,kc) + cc(ic,jc+1,kc) + cc(ic,jc,kc+1) + cc(ic,jc+1,kc+1) )
                   end if
                end if

             end do
          end do
       end do

    else
       !
       ! General ir/=2 case.
       !
       allocate(temp(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
       !
       ! Direct injection for fine points overlaying coarse ones.
       !
       do k = lo(3),hi(3),ir(3)
          kc = k / ir(3) 
          do j = lo(2),hi(2),ir(2)
             jc = j / ir(2) 
             do i = lo(1),hi(1),ir(1)
                ic = i / ir(1) 
                temp(i,j,k) = cc(ic,jc,kc)
             end do
          end do
       end do
       !
       ! Interpolate at fine nodes between coarse nodes in the i-direction only.
       !
       do k = lo(3),hi(3),ir(3)
          do j = lo(2),hi(2),ir(2)
             do l = 1, ir(1)-1
                fac_rght = real(l,dp_t) / real(ir(1),dp_t)
                fac_left = ONE - fac_rght
                do i = lo(1),hi(1)-1,ir(1)
                   temp(i+l,j,k) = fac_left*temp(i,j,k) + fac_rght*temp(i+ir(1),j,k)
                end do
             end do
          end do
       end do
       !
       ! Interpolate in the j-direction.
       !
       do k = lo(3),hi(3),ir(3)
          do m = 1, ir(2)-1
             fac_rght = real(m,dp_t) / real(ir(2),dp_t)
             fac_left = ONE - fac_rght
             do j = lo(2),hi(2)-1,ir(2)
                do i = lo(1),hi(1)
                   temp(i,j+m,k) = fac_left*temp(i,j,k)+fac_rght*temp(i,j+ir(2),k)
                end do
             end do
          end do
       end do
       !
       ! Interpolate in the k-direction.
       !
       do n = 1, ir(3)-1
          fac_rght = real(n,dp_t) / real(ir(3),dp_t)
          fac_left = ONE - fac_rght
          do j = lo(2),hi(2)
             do k = lo(3),hi(3)-1,ir(3)
                do i = lo(1),hi(1)
                   temp(i,j,k+n) = fac_left*temp(i,j,k)+fac_rght*temp(i,j,k+ir(3))
                end do
             end do
          end do
       end do
       !
       ! Finally do the prolongation.
       !
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                ff(i,j,k) = ff(i,j,k) + temp(i,j,k)
             end do
          end do
       end do

    endif

  end subroutine nodal_linear_prolongation_3d

  recursive subroutine nodal_cubic_prolongation_3d(ff, lof, cc, loc, lo, hi, ir)
    integer,     intent(in   ) :: loc(:), lof(:)
    integer,     intent(in   ) :: lo(:), hi(:)
    real (dp_t), intent(inout) :: ff(lof(1):,lof(2):,lof(3):)
    real (dp_t), intent(in   ) :: cc(loc(1):,loc(2):,loc(3):)
    integer,     intent(in   ) :: ir(:)

    integer     :: i, j, k, ic, jc, kc, m, ng
    integer     :: clo(3), chi(3), istart, jstart, kstart , ilo(3), ihi(3)
    real (dp_t) :: coeffs(0:63), ipnt, jpnt, kpnt

    if ( ir(1) == 2 .and. ir(2) == 2 .and. ir(3) == 2 ) then
       clo = int_coarsen(lo, 2)
       chi = int_coarsen(hi, 2)

       ng = min((clo(1)-loc(1)), (clo(2)-loc(2)), (clo(3)-loc(3)))

       if ( ng > 0 ) then
          !
          ! Tricubic was requested and we have one ghost cell.
          !
          do kc = clo(3),chi(3)
             k      = 2*kc
             kpnt   = ZERO
             kstart = kc-1

             if ( k == hi(3) ) then
                kpnt   = ONE
                kstart = kc-2
             end if

             do jc = clo(2),chi(2)
                j      = 2*jc
                jpnt   = ZERO
                jstart = jc-1

                if ( j == hi(2) ) then
                   jpnt   = ONE
                   jstart = jc-2
                end if

                do ic = clo(1),chi(1)
                   i      = 2*ic
                   ipnt   = ZERO
                   istart = ic-1

                   if ( i == hi(1) ) then
                      ipnt   = ONE
                      istart = ic-2
                   end if
                   !
                   ! Direct injection for fine points overlaying coarse ones.
                   !
                   ff(i,j,k) = ff(i,j,k) + cc(ic,jc,kc)

                   if ( i < hi(1) .or. j < hi(2) .or. k < hi(3) ) then

                      coeffs( 0 :  3) = (/ (cc(m,jstart+0,kstart+0), m=istart,istart+3) /)
                      coeffs( 4 :  7) = (/ (cc(m,jstart+1,kstart+0), m=istart,istart+3) /)
                      coeffs( 8 : 11) = (/ (cc(m,jstart+2,kstart+0), m=istart,istart+3) /)
                      coeffs(12 : 15) = (/ (cc(m,jstart+3,kstart+0), m=istart,istart+3) /)

                      coeffs(16 : 19) = (/ (cc(m,jstart+0,kstart+1), m=istart,istart+3) /)
                      coeffs(20 : 23) = (/ (cc(m,jstart+1,kstart+1), m=istart,istart+3) /)
                      coeffs(24 : 27) = (/ (cc(m,jstart+2,kstart+1), m=istart,istart+3) /)
                      coeffs(28 : 31) = (/ (cc(m,jstart+3,kstart+1), m=istart,istart+3) /)

                      coeffs(32 : 35) = (/ (cc(m,jstart+0,kstart+2), m=istart,istart+3) /)
                      coeffs(36 : 39) = (/ (cc(m,jstart+1,kstart+2), m=istart,istart+3) /)
                      coeffs(40 : 43) = (/ (cc(m,jstart+2,kstart+2), m=istart,istart+3) /)
                      coeffs(44 : 47) = (/ (cc(m,jstart+3,kstart+2), m=istart,istart+3) /)

                      coeffs(48 : 51) = (/ (cc(m,jstart+0,kstart+3), m=istart,istart+3) /)
                      coeffs(52 : 55) = (/ (cc(m,jstart+1,kstart+3), m=istart,istart+3) /)
                      coeffs(56 : 59) = (/ (cc(m,jstart+2,kstart+3), m=istart,istart+3) /)
                      coeffs(60 : 63) = (/ (cc(m,jstart+3,kstart+3), m=istart,istart+3) /)

                   end if

                   if ( i < hi(1) ) then
                      ff(i+1,j,k) = ff(i+1,j,k) + tricubicInterpolate(coeffs,HALF,jpnt,kpnt)

                      if ( j < hi(2) ) then
                         ff(i+1,j+1,k) = ff(i+1,j+1,k) + tricubicInterpolate(coeffs,HALF,HALF,kpnt)
                      end if
                   end if

                   if ( j < hi(2) ) then
                      ff(i,j+1,k) = ff(i,j+1,k) + tricubicInterpolate(coeffs,ipnt,HALF,kpnt)
                   end if

                   if ( k < hi(3) ) then
                      ff(i,j,k+1) = ff(i,j,k+1) + tricubicInterpolate(coeffs,ipnt,jpnt,HALF)

                      if ( i < hi(1) ) then
                         ff(i+1,j,k+1) = ff(i+1,j,k+1) + tricubicInterpolate(coeffs,HALF,jpnt,HALF)

                         if ( j < hi(2) ) then
                            ff(i+1,j+1,k+1) = ff(i+1,j+1,k+1) + tricubicInterpolate(coeffs,HALF,HALF,HALF)
                         end if
                      end if

                      if ( j < hi(2) ) then
                         ff(i,j+1,k+1) = ff(i,j+1,k+1) + tricubicInterpolate(coeffs,ipnt,HALF,HALF)
                      end if
                   end if

                end do
             end do
          end do

       else
          !
          ! First do linear interp
          !
          call nodal_linear_prolongation_3d(ff, lof, cc, loc, lo, hi, ir)
          !
          ! If the fine region is big enough do tricubic on interior.
          !
          ilo = lo+2
          ihi = hi-2

          if ( all( (ihi-ilo) > 0) ) then
             call nodal_cubic_prolongation_3d(ff, lof, cc, loc, ilo, ihi, ir)
          end if
       end if
    else
       !
       ! Call lininterp in general ir != 2 case.
       !
       call nodal_linear_prolongation_3d(ff, lof, cc, loc, lo, hi, ir)
    endif

  end subroutine nodal_cubic_prolongation_3d

end module mg_prolongation_module
