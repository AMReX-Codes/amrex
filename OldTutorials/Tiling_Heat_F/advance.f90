module advance_module

  use multifab_module

  implicit none

  private

  public :: advance

contains

  subroutine advance(phi_old,phi_new,dx,dt)

    type(multifab) , intent(inout) :: phi_old,phi_new
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(in   ) :: dt

    ! local variables
    integer i, ng_p

    real(kind=dp_t), pointer ::  pp_old(:,:,:,:), pp_new(:,:,:,:)

    type(multifab) :: temp
    
    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: lo(3), hi(3), tlo(3), thi(3), tsize(3)

    logical :: do_tiling

    ! timer for profiling - enable by building code with PROF=t
    type(bl_prof_timer), save :: bpt

    ! open the timer
    call build(bpt,"advance")

    ! Swap pointers to multifabs
    temp = phi_old
    phi_old = phi_new
    phi_new = temp

    ng_p = phi_old%ng

    tsize = (/ 128, 4, 4 /)

    do_tiling = .true.

    if(do_tiling) then

       !$omp parallel private(i,mfi,tilebox,tlo,thi,pp_old,pp_new,lo,hi)
    
       call mfiter_build(mfi, phi_old, tiling= .true., tilesize= tsize)
       
       do while(next_tile(mfi,i))
          
          tilebox = get_tilebox(mfi)
          tlo = lwb(tilebox)
          thi = upb(tilebox)
          
          pp_old => dataptr(phi_old,i)
          pp_new => dataptr(phi_new,i)
          lo = lwb(get_box(phi_old,i))
          hi = upb(get_box(phi_old,i))
          
          call advance_phi(pp_old(:,:,:,1), pp_new(:,:,:,1), &
               ng_p, lo, hi, dx, dt, tlo, thi)
          
       end do
       !$omp end parallel

    else 
       do i=1,nfabs(phi_old)
          pp_old => dataptr(phi_old,i)
          pp_new => dataptr(phi_new,i)
          lo = lwb(get_box(phi_old,i))
          hi = upb(get_box(phi_old,i))

          call advance_phi2(pp_old(:,:,:,1), pp_new(:,:,:,1), &
               ng_p, lo, hi, dx, dt)

       end do
    end if

    call destroy(bpt)

    call multifab_fill_boundary(phi_new)


  end subroutine advance


  subroutine advance_phi(phi_old, phi_new, ng_p, glo, ghi, dx, dt, tlo, thi)
 
    integer :: glo(3), ghi(3), ng_p, tlo(3), thi(3)
    double precision :: phi_old(glo(1)-ng_p:,glo(2)-ng_p:,glo(3)-ng_p:)
    double precision :: phi_new(glo(1)-ng_p:,glo(2)-ng_p:,glo(3)-ng_p:)
    double precision :: dx, dt

    ! local variables   
    integer :: i,j,k

    double precision :: dxinv, dtdx

    double precision :: fx(tlo(1):thi(1)+1,tlo(2):thi(2)  ,tlo(3):thi(3))
    double precision :: fy(tlo(1):thi(1)  ,tlo(2):thi(2)+1,tlo(3):thi(3))
    double precision :: fz(tlo(1):thi(1)  ,tlo(2):thi(2)  ,tlo(3):thi(3)+1)

    dxinv = 1.d0/dx
    dtdx = dt*dxinv

    ! x-fluxes
     do k=tlo(3),thi(3)
        do j=tlo(2),thi(2)
           do i=tlo(1),thi(1)+1
              fx(i,j,k) = ( phi_old(i,j,k) - phi_old(i-1,j,k) ) * dxinv
           end do
        end do
     end do

     ! y-fluxes
     do k=tlo(3),thi(3)
        do j=tlo(2),thi(2)+1
           do i=tlo(1),thi(1)
              fy(i,j,k) = ( phi_old(i,j,k) - phi_old(i,j-1,k) ) * dxinv
           end do
        end do
     end do
     
     ! z-fluxes
     do k=tlo(3),thi(3)+1
        do j=tlo(2),thi(2)
           do i=tlo(1),thi(1)
              fz(i,j,k) = ( phi_old(i,j,k) - phi_old(i,j,k-1) ) * dxinv
           end do
        end do
     end do
     
     do k=tlo(3),thi(3)
        do j=tlo(2),thi(2)
           do i=tlo(1),thi(1)
              phi_new(i,j,k) = phi_old(i,j,k) + dtdx * &
                   ( fx(i+1,j,k) - fx(i,j,k) &
                   + fy(i,j+1,k) - fy(i,j,k) &
                   + fz(i,j,k+1) - fz(i,j,k) )
           end do
        end do
     end do

   end subroutine advance_phi

   subroutine advance_phi2(phi_old, phi_new, ng_p, lo, hi, dx, dt)
     
     integer :: lo(3), hi(3), ng_p
     double precision :: phi_old(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
     double precision :: phi_new(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
     double precision :: dx, dt
   

     ! local variables
     integer :: i,j,k

     double precision :: dxinv, dtdx

     double precision :: fx(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3))
     double precision :: fy(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3))
     double precision :: fz(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1)

     dxinv = 1.d0/dx
     dtdx = dt*dxinv

    !$omp parallel private(i,j,k)

    ! x-fluxes
    !$omp do collapse(2)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             fx(i,j,k) = ( phi_old(i,j,k) - phi_old(i-1,j,k) ) * dxinv
          end do
       end do
    end do
    !$omp end do nowait

    ! y-fluxes
    !$omp do collapse(2)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             fy(i,j,k) = ( phi_old(i,j,k) - phi_old(i,j-1,k) ) * dxinv
          end do
       end do
    end do
    !$omp end do nowait

    ! z-fluxes
    !$omp do collapse(2)
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             fz(i,j,k) = ( phi_old(i,j,k) - phi_old(i,j,k-1) ) * dxinv
          end do
       end do
    end do
    !$omp end do

    !$omp do collapse(2)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             phi_new(i,j,k) = phi_old(i,j,k) + dtdx * &
                  ( fx(i+1,j,k) - fx(i,j,k) &
                   +fy(i,j+1,k) - fy(i,j,k) &
                   +fz(i,j,k+1) - fz(i,j,k) )
          end do
       end do
    end do
    !$omp end do

    !$omp end parallel

   end subroutine advance_phi2

end module advance_module
