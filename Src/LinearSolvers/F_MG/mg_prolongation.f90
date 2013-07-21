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

  subroutine nodal_prolongation_1d(ff, cc, ir)
    real (dp_t), intent(inout) :: ff(0:)
    real (dp_t), intent(in)    :: cc(0:)
    integer,     intent(in)    :: ir(:)
    integer                    :: nx, i, l
    real (dp_t)                :: fac_left, fac_rght, cc_avg

    nx = size(cc,dim=1)-1

    do i = 0, nx
       ff(ir(1)*i) = ff(ir(1)*i) + cc(i)
    end do

    do l = 1, ir(1)-1
       fac_left = real(l,kind=dp_t) / real(ir(1),kind=dp_t)
       fac_rght = 1.0_dp_t - fac_left
       do i = 0, nx - 1
          cc_avg = fac_left*cc(i) + fac_rght*cc(i+1)
          ff(ir(1)*i+l) = ff(ir(1)*i+l) + cc_avg
       end do
    end do

  end subroutine nodal_prolongation_1d

  subroutine nodal_prolongation_2d(ff, cc, ir)
    real (dp_t), intent(inout) :: ff(0:,0:)
    real (dp_t), intent(inout) :: cc(0:,0:)
    integer,     intent(in)    :: ir(:)
    integer                    :: nx, ny, i, j, l, m
    real (dp_t)                :: fac_left, fac_rght
    real (dp_t)                :: temp(0:size(ff,dim=1)-1,0:size(ff,dim=2)-1)

    nx = size(cc,dim=1)-1
    ny = size(cc,dim=2)-1
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

  end subroutine nodal_prolongation_2d

  subroutine nodal_prolongation_3d(ff, cc, ir)
    real (dp_t), intent(inout) :: ff(0:,0:,0:)
    real (dp_t), intent(in)    :: cc(0:,0:,0:)
    integer,     intent(in)    :: ir(:)
    integer                    ::  nx, ny, nz, i, j, k, l, m, n
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
