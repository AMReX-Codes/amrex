module mg_prolongation_module

  use bl_types

  implicit none

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

  private :: cubicInterpolate, bicubicInterpolate, tricubicInterpolate

contains

  subroutine pc_c_prolongation_1d(ff, cc, ir)
    real (dp_t), intent(inout) :: ff(0:)
    real (dp_t), intent(in)    :: cc(0:)
    integer,     intent(in)    :: ir(:)
    integer                    :: nx, i, l

    nx = size(cc,dim=1)

    do l = 0, ir(1)-1
       do i = 0, nx - 1
          ff(ir(1)*i+l) = ff(ir(1)*i+l) + cc(i)
       end do
    end do

  end subroutine pc_c_prolongation_1d

  subroutine pc_c_prolongation_2d(ff, cc, ir)
    real (dp_t), intent(inout) :: ff(0:,0:)
    real (dp_t), intent(in)    :: cc(0:,0:)
    integer,     intent(in)    :: ir(:)

    integer :: nx, ny, i, j, l, m, twoi, twoj, twoip1, twojp1

    nx = size(cc,dim=1)
    ny = size(cc,dim=2)

    if ( ir(1) == 2 .and. ir(2) == 2 ) then

       do j = 0, ny - 1
             twoj   = 2*j
             twojp1 = 2*j+1

          do i = 0, nx-1
             twoi   = 2*i
             twoip1 = 2*i+1

             ff(twoi,   twoj  ) = ff(twoi,   twoj  ) + cc(i,j)
             ff(twoip1, twoj  ) = ff(twoip1, twoj  ) + cc(i,j)
             ff(twoi,   twojp1) = ff(twoi,   twojp1) + cc(i,j)
             ff(twoip1, twojp1) = ff(twoip1, twojp1) + cc(i,j)
          end do
       end do

    else

       do m = 0, ir(2)-1
          do l = 0, ir(1)-1
             do j = 0, ny - 1
                do i = 0, nx - 1
                   ff(ir(1)*i+l, ir(2)*j+m) = ff(ir(1)*i+l, ir(2)*j+m) + cc(i,j)
                end do
             end do
          end do
       end do

    end if

  end subroutine pc_c_prolongation_2d

  subroutine pc_c_prolongation_3d(ff, cc, ir)
    real (dp_t), intent(inout) :: ff(0:,0:,0:)
    real (dp_t), intent(in)    :: cc(0:,0:,0:)
    integer,     intent(in)    :: ir(:)

    integer :: nx, ny, nz, i, j, k, l, m, n, twoi, twoj, twoip1, twojp1, twok, twokp1

    nx = size(cc,dim=1)
    ny = size(cc,dim=2)
    nz = size(cc,dim=3)

    if ( ir(1) == 2 .and. ir(2) == 2 .and. ir(3) == 2 ) then

       do k = 0, nz-1
          twok   = 2*k
          twokp1 = 2*k+1

          do j = 0, ny-1
             twoj   = 2*j
             twojp1 = 2*j+1

             do i = 0, nx-1
                twoi   = 2*i
                twoip1 = 2*i+1

                ff(twoi,   twoj,   twok  ) = ff(twoi,   twoj,   twok  ) + cc(i,j,k)
                ff(twoip1, twoj,   twok  ) = ff(twoip1, twoj,   twok  ) + cc(i,j,k)
                ff(twoi,   twojp1, twok  ) = ff(twoi,   twojp1, twok  ) + cc(i,j,k)
                ff(twoip1, twojp1, twok  ) = ff(twoip1, twojp1, twok  ) + cc(i,j,k)
                ff(twoi,   twoj,   twokp1) = ff(twoi,   twoj,   twokp1) + cc(i,j,k)
                ff(twoip1, twoj,   twokp1) = ff(twoip1, twoj,   twokp1) + cc(i,j,k)
                ff(twoi,   twojp1, twokp1) = ff(twoi,   twojp1, twokp1) + cc(i,j,k)
                ff(twoip1, twojp1, twokp1) = ff(twoip1, twojp1, twokp1) + cc(i,j,k)
             end do
          end do
       end do

    else

       do k = 0, nz - 1
          do n = 0, ir(3)-1
             do j = 0, ny - 1
                do m = 0, ir(2)-1
                   do i = 0, nx - 1
                      do l = 0, ir(1)-1
                         ff(ir(1)*i+l, ir(2)*j+m, ir(3)*k+n) = ff(ir(1)*i+l, ir(2)*j+m, ir(3)*k+n) + cc(i,j,k)
                      end do
                   end do
                end do
             end do
          end do
       end do

    end if

  end subroutine pc_c_prolongation_3d

  subroutine lin_c_prolongation_1d(ff, cc, ir)
    use bl_error_module
    real (dp_t), intent(inout) :: ff(0:)
    real (dp_t), intent(in)    :: cc(0:)
    integer,     intent(in)    :: ir(:)
    integer                    :: nx, i

    call bl_assert(ir(1)==2, 'lin_c_prolongation_1d: ir==2')

    nx = size(cc,dim=1)

    i = 0
    ff(2*i  ) = ff(2*i  ) + cc(i)
    ff(2*i+1) = ff(2*i+1) + cc(i)

    i = nx-1
    ff(2*i  ) = ff(2*i  ) + cc(i)
    ff(2*i+1) = ff(2*i+1) + cc(i)

    do i = 1, nx-2
       ff(2*i  ) = ff(2*i  ) + 0.25d0*( 3*cc(i) + cc(i-1) )
       ff(2*i+1) = ff(2*i+1) + 0.25d0*( 3*cc(i) + cc(i+1) )
    end do

  end subroutine lin_c_prolongation_1d

  subroutine lin_c_prolongation_2d(ff, cc, ir, ptype)
    use bl_error_module
    real (dp_t), intent(inout) :: ff(0:,0:)
    real (dp_t), intent(in)    :: cc(0:,0:)
    integer,     intent(in)    :: ir(:), ptype
    integer                    :: nx, ny, i, j
    logical                    :: interior_i, interior_j
    real (dp_t), parameter     :: one6th  = 1.0d0/6.0d0
    real (dp_t), parameter     :: one16th = 1.0d0/16.0d0

    nx = size(cc,dim=1)
    ny = size(cc,dim=2)

    call bl_assert(ir(1)==2 .and. ir(2)==2, 'lin_c_prolongation_2d: ir==2')
    !
    ! First do all face points using piecewise-constant interpolation.
    !
    do j = 0, ny-1
       interior_j = ( (j > 0) .and. (j < (ny-1)) )

       do i = 0, nx-1
          interior_i = ( (i > 0) .and. (i < (nx-1)) )

          if ( interior_i .and. interior_j ) cycle

          ff(2*i,   2*j  ) = ff(2*i,   2*j  ) + cc(i,j)
          ff(2*i+1, 2*j  ) = ff(2*i+1, 2*j  ) + cc(i,j)
          ff(2*i,   2*j+1) = ff(2*i,   2*j+1) + cc(i,j)
          ff(2*i+1, 2*j+1) = ff(2*i+1, 2*j+1) + cc(i,j)
       end do
    end do

    select case ( ptype )
    case ( 1 )
       !
       ! Type 1 - {2,1,1} weighted average of our neighbors.
       !
       do j = 1, ny-2
          do i = 1, nx-2
             ff(2*i+1, 2*j+1) = ff(2*i+1, 2*j+1) + .25d0 * ( 2*cc(i,j) + cc(i+1,j) + cc(i,j+1) )
             ff(2*i,   2*j+1) = ff(2*i,   2*j+1) + .25d0 * ( 2*cc(i,j) + cc(i-1,j) + cc(i,j+1) )
             ff(2*i+1, 2*j  ) = ff(2*i+1, 2*j  ) + .25d0 * ( 2*cc(i,j) + cc(i+1,j) + cc(i,j-1) )
             ff(2*i,   2*j  ) = ff(2*i,   2*j  ) + .25d0 * ( 2*cc(i,j) + cc(i-1,j) + cc(i,j-1) )
          end do
       end do
    case ( 2 )
       !
       ! Type 2 - {4,1,1} weighted average of our neighbors.
       !
       do j = 1, ny-2
          do i = 1, nx-2
             ff(2*i+1, 2*j+1) = ff(2*i+1, 2*j+1) + one6th * ( 4*cc(i,j) + cc(i+1,j) + cc(i,j+1) )
             ff(2*i,   2*j+1) = ff(2*i,   2*j+1) + one6th * ( 4*cc(i,j) + cc(i-1,j) + cc(i,j+1) )
             ff(2*i+1, 2*j  ) = ff(2*i+1, 2*j  ) + one6th * ( 4*cc(i,j) + cc(i+1,j) + cc(i,j-1) )
             ff(2*i,   2*j  ) = ff(2*i,   2*j  ) + one6th * ( 4*cc(i,j) + cc(i-1,j) + cc(i,j-1) )
          end do
       end do
    case ( 3 )
       do j = 1, ny-2
          do i = 1, nx-2
             !
             ! Type 3 - bilinear.
             !
             ff(2*i+1, 2*j+1) = ff(2*i+1, 2*j+1) + one16th * &
                  ( 9*cc(i,j) + 3*cc(i+1,j) + 3*cc(i,j+1) + cc(i+1,j+1) )
             ff(2*i,   2*j+1) = ff(2*i,   2*j+1) + one16th * &
                  ( 9*cc(i,j) + 3*cc(i-1,j) + 3*cc(i,j+1) + cc(i-1,j+1) )
             ff(2*i+1, 2*j  ) = ff(2*i+1, 2*j  ) + one16th * &
                  ( 9*cc(i,j) + 3*cc(i,j-1) + 3*cc(i+1,j) + cc(i+1,j-1) )
             ff(2*i,   2*j  ) = ff(2*i,   2*j  ) + one16th * &
                  ( 9*cc(i,j) + 3*cc(i-1,j) + 3*cc(i,j-1) + cc(i-1,j-1) )
          end do
       end do
    case default
       call bl_error("lin_c_prolongation_2d: unknown ptype", ptype)
    end select

  end subroutine lin_c_prolongation_2d

  subroutine lin_c_prolongation_3d(ff, cc, ir, ptype)
    use bl_error_module
    real (dp_t), intent(inout) :: ff(0:,0:,0:)
    real (dp_t), intent(in)    :: cc(0:,0:,0:)
    integer,     intent(in)    :: ir(:), ptype
    integer                    :: nx, ny, nz, i, j, k
    logical                    :: interior_i, interior_j, interior_k
    real (dp_t), parameter     ::   one64ths = 1.0d0/64.0d0
    real (dp_t), parameter     :: three64ths = 3.0d0/64.0d0
    real (dp_t), parameter     ::     sixth  = 1.0d0/6.0d0

    call bl_assert(ir(1)==2 .and. ir(2)==2 .and. ir(3)==2, 'lin_c_prolongation_3d: ir==2')

    nx = size(cc,dim=1)
    ny = size(cc,dim=2)
    nz = size(cc,dim=3)
    !
    ! First do all face points using piecewise-constant interpolation.
    !
    do k = 0, nz-1
       interior_k = ( (k > 0) .and. (k < (nz-1)) )

       do j = 0, ny-1
          interior_j = ( (j > 0) .and. (j < (ny-1)) )

          do i = 0, nx-1
             interior_i = ( (i > 0) .and. (i < (nx-1)) )

             if ( interior_i .and. interior_j .and. interior_k ) cycle

             ff(2*i+1, 2*j+1, 2*k+1) = ff(2*i+1, 2*j+1, 2*k+1) + cc(i,j,k)
             ff(2*i,   2*j+1, 2*k+1) = ff(2*i,   2*j+1, 2*k+1) + cc(i,j,k)
             ff(2*i+1, 2*j,   2*k+1) = ff(2*i+1, 2*j,   2*k+1) + cc(i,j,k)
             ff(2*i,   2*j,   2*k+1) = ff(2*i,   2*j,   2*k+1) + cc(i,j,k)
             ff(2*i+1, 2*j+1, 2*k  ) = ff(2*i+1, 2*j+1, 2*k  ) + cc(i,j,k)
             ff(2*i,   2*j+1, 2*k  ) = ff(2*i,   2*j+1, 2*k  ) + cc(i,j,k)
             ff(2*i+1, 2*j,   2*k  ) = ff(2*i+1, 2*j,   2*k  ) + cc(i,j,k)
             ff(2*i,   2*j,   2*k  ) = ff(2*i,   2*j,   2*k  ) + cc(i,j,k)
          end do
       end do
    end do

    select case ( ptype )
    case ( 1 )
       !
       ! Type 1 - {1,1,1,1} weighted average of our neighbors.
       !
       do k = 1, nz-2
          do j = 1, ny-2
             do i = 1, nx-2
                ff(2*i+1, 2*j+1, 2*k+1) = ff(2*i+1, 2*j+1, 2*k+1) + &
                     .25d0 * ( cc(i,j,k) + cc(i+1,j,k) + cc(i,j+1,k) + cc(i,j,k+1) )
                ff(2*i,   2*j+1, 2*k+1) = ff(2*i,   2*j+1, 2*k+1) + &
                     .25d0 * ( cc(i,j,k) + cc(i-1,j,k) + cc(i,j+1,k) + cc(i,j,k+1) )
                ff(2*i+1, 2*j,   2*k+1) = ff(2*i+1, 2*j,   2*k+1) + &
                     .25d0 * ( cc(i,j,k) + cc(i+1,j,k) + cc(i,j-1,k) + cc(i,j,k+1) )
                ff(2*i,   2*j,   2*k+1) = ff(2*i,   2*j,   2*k+1) + &
                     .25d0 * ( cc(i,j,k) + cc(i-1,j,k) + cc(i,j-1,k) + cc(i,j,k+1) )
                ff(2*i+1, 2*j+1, 2*k  ) = ff(2*i+1, 2*j+1, 2*k  ) + &
                     .25d0 * ( cc(i,j,k) + cc(i+1,j,k) + cc(i,j+1,k) + cc(i,j,k-1) )
                ff(2*i,   2*j+1, 2*k  ) = ff(2*i,   2*j+1, 2*k  ) + &
                     .25d0 * ( cc(i,j,k) + cc(i-1,j,k) + cc(i,j+1,k) + cc(i,j,k-1) )
                ff(2*i+1, 2*j,   2*k  ) = ff(2*i+1, 2*j,   2*k  ) + &
                     .25d0 * ( cc(i,j,k) + cc(i+1,j,k) + cc(i,j-1,k) + cc(i,j,k-1) )
                ff(2*i,   2*j,   2*k  ) = ff(2*i,   2*j,   2*k  ) + &
                     .25d0 * ( cc(i,j,k) + cc(i-1,j,k) + cc(i,j-1,k) + cc(i,j,k-1) )
             end do
          end do
       end do
    case ( 2 )
       !
       ! Type 2 - {3,1,1,1} weighted average of our neighbors.
       !
       do k = 1, nz-2
          do j = 1, ny-2
             do i = 1, nx-2
                ff(2*i+1, 2*j+1, 2*k+1) = ff(2*i+1, 2*j+1, 2*k+1) + &
                     sixth * ( 3*cc(i,j,k) + cc(i+1,j,k) + cc(i,j+1,k) + cc(i,j,k+1) )
                ff(2*i,   2*j+1, 2*k+1) = ff(2*i,   2*j+1, 2*k+1) + &
                     sixth * ( 3*cc(i,j,k) + cc(i-1,j,k) + cc(i,j+1,k) + cc(i,j,k+1) )
                ff(2*i+1, 2*j,   2*k+1) = ff(2*i+1, 2*j,   2*k+1) + &
                     sixth * ( 3*cc(i,j,k) + cc(i+1,j,k) + cc(i,j-1,k) + cc(i,j,k+1) )
                ff(2*i,   2*j,   2*k+1) = ff(2*i,   2*j,   2*k+1) + &
                     sixth * ( 3*cc(i,j,k) + cc(i-1,j,k) + cc(i,j-1,k) + cc(i,j,k+1) )
                ff(2*i+1, 2*j+1, 2*k  ) = ff(2*i+1, 2*j+1, 2*k  ) + &
                     sixth * ( 3*cc(i,j,k) + cc(i+1,j,k) + cc(i,j+1,k) + cc(i,j,k-1) )
                ff(2*i,   2*j+1, 2*k  ) = ff(2*i,   2*j+1, 2*k  ) + &
                     sixth * ( 3*cc(i,j,k) + cc(i-1,j,k) + cc(i,j+1,k) + cc(i,j,k-1) )
                ff(2*i+1, 2*j,   2*k  ) = ff(2*i+1, 2*j,   2*k  ) + &
                     sixth * ( 3*cc(i,j,k) + cc(i+1,j,k) + cc(i,j-1,k) + cc(i,j,k-1) )
                ff(2*i,   2*j,   2*k  ) = ff(2*i,   2*j,   2*k  ) + &
                     sixth * ( 3*cc(i,j,k) + cc(i-1,j,k) + cc(i,j-1,k) + cc(i,j,k-1) )
             end do
          end do
       end do
    case ( 3 )
       !
       ! Type 3 - trilinear.
       !
       do k = 1, nz-2
          do j = 1, ny-2
             do i = 1, nx-2
                 ff(2*i+1, 2*j+1, 2*k+1) = ff(2*i+1, 2*j+1, 2*k+1) + &
                      three64ths * ( 9*cc(i,j,k  ) + 3*cc(i+1,j,k  ) + 3*cc(i,j+1,k  ) + cc(i+1,j+1,k  ) ) + &
                        one64ths * ( 9*cc(i,j,k+1) + 3*cc(i+1,j,k+1) + 3*cc(i,j+1,k+1) + cc(i+1,j+1,k+1) )
                 ff(2*i,   2*j+1, 2*k+1) = ff(2*i,   2*j+1, 2*k+1) + &
                       three64ths * ( 9*cc(i,j,k  ) + 3*cc(i-1,j,k  ) + 3*cc(i,j+1,k  ) + cc(i-1,j+1,k  ) ) + &
                         one64ths * ( 9*cc(i,j,k+1) + 3*cc(i-1,j,k+1) + 3*cc(i,j+1,k+1) + cc(i-1,j+1,k+1) )
                 ff(2*i+1, 2*j,   2*k+1) = ff(2*i+1, 2*j,   2*k+1) + &
                      three64ths * ( 9*cc(i,j,k  ) + 3*cc(i+1,j,k  ) + 3*cc(i,j-1,k  ) + cc(i+1,j-1,k  ) ) + &
                        one64ths * ( 9*cc(i,j,k+1) + 3*cc(i+1,j,k+1) + 3*cc(i,j-1,k+1) + cc(i+1,j-1,k+1) )
                 ff(2*i,   2*j,   2*k+1) = ff(2*i,   2*j,   2*k+1) + &
                      three64ths * ( 9*cc(i,j,k  ) + 3*cc(i-1,j,k  ) + 3*cc(i,j-1,k  ) + cc(i-1,j-1,k  ) ) + &
                        one64ths * ( 9*cc(i,j,k+1) + 3*cc(i-1,j,k+1) + 3*cc(i,j-1,k+1) + cc(i-1,j-1,k+1) )
                 ff(2*i+1, 2*j+1, 2*k) = ff(2*i+1, 2*j+1, 2*k) + &
                      three64ths * ( 9*cc(i,j,k  ) + 3*cc(i+1,j,k  ) + 3*cc(i,j+1,k  ) + cc(i+1,j+1,k  ) ) + &
                        one64ths * ( 9*cc(i,j,k-1) + 3*cc(i+1,j,k-1) + 3*cc(i,j+1,k-1) + cc(i+1,j+1,k-1) )
                 ff(2*i,   2*j+1, 2*k) = ff(2*i,   2*j+1, 2*k) + &
                       three64ths * ( 9*cc(i,j,k  ) + 3*cc(i-1,j,k  ) + 3*cc(i,j+1,k  ) + cc(i-1,j+1,k  ) ) + &
                         one64ths * ( 9*cc(i,j,k-1) + 3*cc(i-1,j,k-1) + 3*cc(i,j+1,k-1) + cc(i-1,j+1,k-1) )
                 ff(2*i+1, 2*j,   2*k) = ff(2*i+1, 2*j,   2*k) + &
                      three64ths * ( 9*cc(i,j,k  ) + 3*cc(i+1,j,k  ) + 3*cc(i,j-1,k  ) + cc(i+1,j-1,k  ) ) + &
                        one64ths * ( 9*cc(i,j,k-1) + 3*cc(i+1,j,k-1) + 3*cc(i,j-1,k-1) + cc(i+1,j-1,k-1) )
                 ff(2*i,   2*j,   2*k) = ff(2*i,   2*j,   2*k) + &
                      three64ths * ( 9*cc(i,j,k  ) + 3*cc(i-1,j,k  ) + 3*cc(i,j-1,k  ) + cc(i-1,j-1,k  ) ) + &
                        one64ths * ( 9*cc(i,j,k-1) + 3*cc(i-1,j,k-1) + 3*cc(i,j-1,k-1) + cc(i-1,j-1,k-1) )
             end do
          end do
       end do
    case default
       call bl_error("lin_c_prolongation_3d: unknown ptype", ptype)
    end select

  end subroutine lin_c_prolongation_3d

  pure function cubicInterpolate (p, x) result (r)
    real (dp_t), intent(in) :: p(0:3), x
    real (dp_t) r
    r=p(1)+0.5*x*(p(2)-p(0)+x*(2*p(0)-5*p(1)+4*p(2)-p(3)+x*(3*(p(1)-p(2))+p(3)-p(0))))
  end function cubicInterpolate

  function bicubicInterpolate (p, x, y) result (r)
    real (dp_t), intent(in) :: p(0:15), x, y
    real (dp_t) r, arr(0:3)
    !
    ! First interpolate in y.
    !
    arr(0) = cubicInterpolate(p( 0:3), x) ! row 0
    arr(1) = cubicInterpolate(p( 4:7), x) ! row 1
    arr(2) = cubicInterpolate(p( 8:11),x) ! row 2
    arr(3) = cubicInterpolate(p(12:15),x) ! row 3
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

  subroutine nodal_prolongation_1d(ff, cc, ir, ng)
    real (dp_t), intent(inout) :: ff(0:)
    real (dp_t), intent(in)    :: cc(-ng:)
    integer,     intent(in)    :: ir(:), ng
    integer                    :: nx, i, l
    real (dp_t)                :: fac_left, fac_rght, cc_avg, coeffs(0:15)

    nx = size(cc,dim=1)-2*ng-1

    do i = 0, nx
       ff(ir(1)*i) = ff(ir(1)*i) + cc(i)
    end do

    if ( ir(1) == 2 ) then
       !
       ! Use cubic interpolation on remaining points.
       !
       do i = 0, nx-1

          coeffs(0) = cc(i-1)
          coeffs(1) = cc(i+0)
          coeffs(2) = cc(i+1)
          coeffs(3) = cc(i+2)

          ff(ir(1)*i+1) = ff(ir(1)*i+1) + cubicInterpolate(coeffs, 0.5d0)
       end do
    else
       !
       ! Otherwise use linear.
       !
       do l = 1, ir(1)-1
          fac_left = real(l,kind=dp_t) / real(ir(1),kind=dp_t)
          fac_rght = 1.0_dp_t - fac_left
          do i = 0, nx - 1
             cc_avg = fac_left*cc(i) + fac_rght*cc(i+1)
             ff(ir(1)*i+l) = ff(ir(1)*i+l) + cc_avg
          end do
       end do
    end if

  end subroutine nodal_prolongation_1d

  recursive subroutine nodal_prolongation_2d(ff, cc, ir, ng, rtype)
    real (dp_t), intent(inout) :: ff(0:,0:)
    real (dp_t), intent(inout) :: cc(-ng:,-ng:)
    integer,     intent(in)    :: ir(:), ng, rtype
    integer                    :: nx, ny, i, j, l, m
    real (dp_t)                :: fac_left, fac_rght, coeffs(0:15)
    real (dp_t)                :: temp(0:size(ff,dim=1)-1,0:size(ff,dim=2)-1)

    nx = size(cc,dim=1)-2*ng-1
    ny = size(cc,dim=2)-2*ng-1

    if ( rtype == 1 .and. ir(1) == 2 .and. ir(2) == 2 ) then
       !
       ! Don't bother when the FABs are small.
       !
       if (nx <= 4 .or. ny <= 4) call nodal_prolongation_2d(ff,cc,ir,ng,rtype=0)
       !
       ! Use bicubic interpolation.
       !
       do j = 0, ny
          do i = 0, nx
             !
             ! Direct injection for fine points overlaying coarse ones.
             !
             ff(2*i, 2*j) = ff(2*i, 2*j) + cc(i,j)

             if ( i < nx ) then

                if ( j == ny ) then
                   coeffs( 0 :  3) = (/ (cc(m,j-2), m=i-1,i+2) /)
                   coeffs( 4 :  7) = (/ (cc(m,j-1), m=i-1,i+2) /)
                   coeffs( 8 : 11) = (/ (cc(m,j+0), m=i-1,i+2) /)
                   coeffs(12 : 15) = (/ (cc(m,j+1), m=i-1,i+2) /)
                else
                   coeffs( 0 :  3) = (/ (cc(m,j-1), m=i-1,i+2) /)
                   coeffs( 4 :  7) = (/ (cc(m,j+0), m=i-1,i+2) /)
                   coeffs( 8 : 11) = (/ (cc(m,j+1), m=i-1,i+2) /)
                   coeffs(12 : 15) = (/ (cc(m,j+2), m=i-1,i+2) /)
                end if

                ff(2*i+1, 2*j) = ff(2*i+1, 2*j) + bicubicInterpolate(coeffs, 0.5d0, 0.0d0) 

                if ( j < ny ) then

                   ff(2*i+1, 2*j+1) = ff(2*i+1, 2*j+1) + bicubicInterpolate(coeffs, 0.5d0, 0.5d0) 
                end if
             end if

             if ( j < ny ) then

                if ( i == nx ) then
                   coeffs( 0 :  3) = (/ (cc(m,j-1), m=i-2,i+1) /)
                   coeffs( 4 :  7) = (/ (cc(m,j+0), m=i-2,i+1) /)
                   coeffs( 8 : 11) = (/ (cc(m,j+1), m=i-2,i+1) /)
                   coeffs(12 : 15) = (/ (cc(m,j+2), m=i-2,i+1) /)
                end if

                ff(2*i, 2*j+1) = ff(2*i, 2*j+1) + bicubicInterpolate(coeffs, 0.0d0, 0.5d0) 
             end if

          end do
       end do

    else
       !
       ! Interpolate at fine nodes on top of coarse nodes.
       !
       do j = 0,ny
          do i = 0,nx
             temp(ir(1)*i, ir(2)*j) = cc(i,j)
          end do
       end do
       !
       ! Interpolate at fine nodes between coarse nodes in the i-direction only.
       !
       do j = 0,ny
          do l = 1, ir(1)-1
             fac_left = real(l,kind=dp_t) / real(ir(1),kind=dp_t)
             fac_rght = 1.0_dp_t - fac_left
             do i = 0,nx-1
                temp(ir(1)*i+l, ir(2)*j) = fac_left*cc(i,j) + fac_rght*cc(i+1,j)
             end do
          end do
       end do
       !
       ! Interpolate in the j-direction using previously interpolated "temp".
       !
       do m = 1, ir(2)-1
          fac_left = real(m,kind=dp_t) / real(ir(2),kind=dp_t)
          fac_rght = 1.0_dp_t - fac_left
          do j = 0,ny-1
             do i = 0,ir(1)*nx
                temp(i, ir(2)*j+m) = fac_left*temp(i,ir(2)*(j  )) + &
                     fac_rght*temp(i,ir(2)*(j+1))
             end do
          end do
       end do

       do j = 0,ir(2)*ny
          do i = 0,ir(1)*nx
             ff(i,j) = ff(i,j) + temp(i,j)
          end do
       end do

    end if

  end subroutine nodal_prolongation_2d

  subroutine nodal_prolongation_3d(ff, cc, ir)
    real (dp_t), intent(inout) :: ff(0:,0:,0:)
    real (dp_t), intent(in)    :: cc(0:,0:,0:)
    integer,     intent(in)    :: ir(:)
    integer                    :: nx, ny, nz, i, j, k, l, m, n
    real (dp_t)                :: fac_left, fac_rght
    real (dp_t), allocatable   :: temp(:,:,:)

    allocate(temp(0:size(ff,dim=1)-1,0:size(ff,dim=2)-1,0:size(ff,dim=3)-1))

    nx = size(cc,dim=1)-1
    ny = size(cc,dim=2)-1
    nz = size(cc,dim=3)-1
    !
    ! Interpolate at fine nodes on top of coarse nodes.
    !
    do k = 0,nz
       do j = 0,ny
          do i = 0,nx
             temp(ir(1)*i,ir(2)*j,ir(3)*k) = cc(i,j,k)
          end do
       end do
    end do
    !
    ! Interpolate at fine nodes between coarse nodes in the i-direction only.
    !
    do k = 0,nz
       do j = 0,ny
          do l = 1, ir(1)-1
             fac_left = real(l,kind=dp_t) / real(ir(1),kind=dp_t)
             fac_rght = 1.0_dp_t - fac_left
             do i = 0,nx-1
                temp(ir(1)*i+l,ir(2)*j,ir(3)*k) = fac_left*cc(i  ,j,k)+ &
                                                  fac_rght*cc(i+1,j,k)
             end do
          end do
       end do
    end do
    !
    ! Interpolate in the j-direction using previously interpolated "temp".
    !
    do m = 1, ir(2)-1
       fac_left = real(m,kind=dp_t) / real(ir(2),kind=dp_t)
       fac_rght = 1.0_dp_t - fac_left
       do k = 0,nz
          do j = 0,ny-1
             do i = 0,ir(1)*nx
                temp(i,ir(2)*j+m,ir(3)*k) = fac_left*temp(i,ir(2)*(j  ),ir(3)*k) + &
                                            fac_rght*temp(i,ir(2)*(j+1),ir(3)*k)
             end do
          end do
       end do
    end do
    !
    ! Interpolate in the k-direction using previously interpolated "temp".
    !
    do n = 1, ir(3)-1
       fac_left = real(n,kind=dp_t) / real(ir(3),kind=dp_t)
       fac_rght = 1.0_dp_t - fac_left
       do j = 0,ir(2)*ny
          do k = 0,nz-1
             do i = 0,ir(1)*nx
                temp(i, j, ir(3)*k+n) = fac_left*temp(i,j,ir(3)*(k  )) + &
                                        fac_rght*temp(i,j,ir(3)*(k+1))
             end do
          end do
       end do
    end do

    do k = 0,ir(3)*nz
       do j = 0,ir(2)*ny
          do i = 0,ir(1)*nx
             ff(i,j,k) = ff(i,j,k) + temp(i,j,k)
          end do
       end do
    end do

  end subroutine nodal_prolongation_3d

end module mg_prolongation_module
