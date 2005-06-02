module mg_prolongation_module

  use bl_types

  implicit none

  interface pc_c_prolongation
     module procedure pc_c_prolongation_1d
     module procedure pc_c_prolongation_2d
     module procedure pc_c_prolongation_3d
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
    integer                    :: nx, ny, i, j, l, m

    nx = size(cc,dim=1)
    ny = size(cc,dim=2)

    do m = 0, ir(2)-1
       do l = 0, ir(1)-1
          !$OMP PARALLEL DO PRIVATE(j,i)
          do j = 0, ny - 1
             do i = 0, nx - 1
                ff(ir(1)*i+l, ir(2)*j+m) = ff(ir(1)*i+l, ir(2)*j+m) + cc(i,j)
             end do
          end do
          !$OMP END PARALLEL DO
       end do
    end do

  end subroutine pc_c_prolongation_2d

  subroutine pc_c_prolongation_3d(ff, cc, ir)
    real (dp_t), intent(inout) :: ff(0:,0:,0:)
    real (dp_t), intent(in)    :: cc(0:,0:,0:)
    integer,     intent(in)    :: ir(:)
    integer                    :: nx, ny, nz, i, j, k, l, m, n

    nx = size(cc,dim=1)
    ny = size(cc,dim=2)
    nz = size(cc,dim=3)

    do n = 0, ir(3)-1
       do m = 0, ir(2)-1
          do l = 0, ir(1)-1
             !$OMP PARALLEL DO PRIVATE(j,i,k)
             do k = 0, nz - 1
                do j = 0, ny - 1
                   do i = 0, nx - 1
                      ff(ir(1)*i+l, ir(2)*j+m, ir(3)*k+n) = ff(ir(1)*i+l, ir(2)*j+m, ir(3)*k+n) + cc(i,j,k)
                   end do
                end do
             end do
             !$OMP END PARALLEL DO
          end do
       end do
    end do

  end subroutine pc_c_prolongation_3d

  subroutine nodal_prolongation_1d(ff, cc, ir)
    real (dp_t), intent(inout) :: ff(0:)
    real (dp_t), intent(in)    :: cc(0:)
    integer,     intent(in)    :: ir(:)
    integer                    :: nx, i, l
    real (dp_t)                :: fac_left, fac_rght

    nx = size(cc,dim=1)-1

    do i = 0, nx
       ff(ir(1)*i) = ff(ir(1)*i) + cc(i)
    end do

    do l = 1, ir(1)-1
       fac_left = real(l,kind=dp_t) / real(ir(1),kind=dp_t)
       fac_rght = 1.0_dp_t - fac_left
       do i = 0, nx - 1
          ff(ir(1)*i+l) = ff(ir(1)*i+l) + fac_left*cc(i) + fac_rght*cc(i+1)
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
    real (dp_t)                :: temp(0:size(ff,dim=1)-1,0:size(ff,dim=2)-1,0:size(ff,dim=3)-1)

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
       do k = 0,nz-1
          do j = 0,ir(2)*ny
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
