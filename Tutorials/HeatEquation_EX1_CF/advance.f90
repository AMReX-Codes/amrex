
module advance_module

  use boxlib_module

  implicit none

  private

  public :: advance

contains

  subroutine advance (old_phi, new_phi, geom, dt)
    type(multifab) :: old_phi, new_phi
    type(geometry) :: geom
    double precision :: dt

    integer :: plo(4), phi(4)
    double precision :: dx
    type(box) :: bx
    type(mfiter) :: mfi
    double precision, dimension(:,:,:,:), pointer :: po, pn

    ! For simplicity, we assume dx(1) == dx(2) == dx(3)
    dx = geom%dx(1)

    ! This fills periodic ghost cells and ghost cells from neighboring grids
    call old_phi%fill_boundary(geom)

    !$omp parallel private(mfi,bx,po,pn,plo,phi)
    call mfiter_build(mfi, old_phi, tiling=.true.)

    do while (mfi%next())

       bx = mfi%tilebox()

       po => old_phi%dataptr(mfi)
       pn => new_phi%dataptr(mfi)

       plo = lbound(po)
       phi = ubound(pn)
       
       select case (bl_num_dims)
       case (2)
          call update_phi_2d(bx%lo, bx%hi, po, pn, plo, phi, dx, dt)
       case (3)
          call update_phi_3d(bx%lo, bx%hi, po, pn, plo, phi, dx, dt)
       end select
    end do
    !$omp end parallel

  end subroutine advance

  subroutine update_phi_2d (lo, hi, pold, pnew, plo, phi, dx, dt)
    integer, intent(in) :: lo(2), hi(2), plo(2), phi(2)
    double precision, intent(in   ) :: pold(plo(1):phi(1), plo(2):phi(2))
    double precision, intent(inout) :: pnew(plo(1):phi(1), plo(2):phi(2))
    double precision, intent(in) :: dx, dt

    integer :: i,j
    double precision :: dxinv, dtdx
    double precision :: fx(lo(1):hi(1)+1,lo(2):hi(2)  )
    double precision :: fy(lo(1):hi(1)  ,lo(2):hi(2)+1)
    
    dxinv = 1.d0/dx
    dtdx = dt*dxinv
    
    ! x-fluxes
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)+1
          fx(i,j) = ( pold(i-1,j) - pold(i,j) ) * dxinv
       end do
    end do
    
    ! y-fluxes
    do j=lo(2),hi(2)+1
       do i=lo(1),hi(1)
          fy(i,j) = ( pold(i,j-1) - pold(i,j) ) * dxinv
       end do
    end do

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          
          pnew(i,j) = pold(i,j) - dtdx * &
               ( fx(i+1,j)-fx(i,j) &
               + fy(i,j+1)-fy(i,j) )

       end do
    end do

  end subroutine update_phi_2d

  subroutine update_phi_3d (lo, hi, pold, pnew, plo, phi, dx, dt)
    integer, intent(in) :: lo(3), hi(3), plo(3), phi(3)
    double precision, intent(in   ) :: pold(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
    double precision, intent(inout) :: pnew(plo(1):phi(1), plo(2):phi(2), plo(3):phi(3))
    double precision, intent(in) :: dx, dt

    integer :: i,j,k
    double precision :: dxinv, dtdx
    double precision :: fx(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3))
    double precision :: fy(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3))
    double precision :: fz(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1)
    
    dxinv = 1.d0/dx
    dtdx = dt*dxinv
    
    ! x-fluxes
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             fx(i,j,k) = ( pold(i-1,j,k) - pold(i,j,k) ) * dxinv
          end do
       end do
    end do
    
    ! y-fluxes
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             fy(i,j,k) = ( pold(i,j-1,k) - pold(i,j,k) ) * dxinv
          end do
       end do
    end do
    
    ! z-fluxes
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             fz(i,j,k) = ( pold(i,j,k-1) - pold(i,j,k) ) * dxinv
          end do
       end do
    end do
    
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             
             pnew(i,j,k) = pold(i,j,k) - dtdx * &
                  ( fx(i+1,j,k)-fx(i,j,k) &
                  + fy(i,j+1,k)-fy(i,j,k) &
                  + fz(i,j,k+1)-fz(i,j,k) )
             
          end do
       end do
    end do

  end subroutine update_phi_3d

end module advance_module

