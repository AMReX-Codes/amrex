
module cns_hydro_module

  use amrex_fort_module, only : rt=>amrex_real
  use cns_module, only : urho, umx, umy, umz, ueden, ueint, utemp, nvar, &
       qrho,qu,qv,qw,qp,qc,qeint,qtemp,qvar
  use mempool_module, only : amrex_allocate, amrex_deallocate
  implicit none
  private

  integer, parameter :: nghost_plm = 2  ! number of ghost cells needed for plm

  public:: cns_compute_hydro_flux, cns_eb_compute_hydro_flux

contains

  subroutine cns_compute_hydro_flux (lo,hi, dudt, utlo, uthi, &
       u,ulo,uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx) &
       bind(c,name='cns_compute_hydro_flux')
    use advection_module, only : hyp_mol_gam_3d
    use cns_eb_flux_module, only : compute_diffop
    integer, dimension(3), intent(in) :: lo,hi,utlo,uthi,ulo,uhi, &
         fxlo,fxhi,fylo,fyhi,fzlo,fzhi
    real(rt), intent(inout) :: dudt(utlo(1):uthi(1),utlo(2):uthi(2),utlo(3):uthi(3),nvar)
    real(rt), intent(in   ) :: u ( ulo(1): uhi(1), ulo(2): uhi(2), ulo(3): uhi(3),nvar)
    real(rt), intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),nvar)
    real(rt), intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),nvar)
    real(rt), intent(inout) :: fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),nvar)
    real(rt), intent(in) :: dx(3)

    integer :: qlo(3), qhi(3)
    real(rt), pointer, contiguous :: q(:,:,:,:)

    qlo = lo - nghost_plm
    qhi = hi + nghost_plm
    call amrex_allocate(q, qlo(1),qhi(1), qlo(2),qhi(2), qlo(3),qhi(3), 1,qvar)
    
    call ctoprim(qlo, qhi, u, ulo, uhi, q, qlo, qhi)
    
    call hyp_mol_gam_3d(q, qlo, qhi, lo, hi, dx, fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi)

    call compute_diffop (lo,hi,5,dx,dudt,utlo,uthi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi)

    dudt(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),6:nvar) = 0.d0

    call amrex_deallocate(q)
  end subroutine cns_compute_hydro_flux

  subroutine cns_eb_compute_hydro_flux (lo,hi, dudt, utlo, uthi, &
       u,ulo,uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,flag,fglo,fghi, &
       volfrac,vlo,vhi,apx,axlo,axhi,apy,aylo,ayhi,apz,azlo,azhi, &
       centx,cxlo,cxhi,centy,cylo,cyhi,centz,czlo,czhi, &
       dx) bind(c,name='cns_eb_compute_hydro_flux')
    use eb_advection_module, only : hyp_mol_gam_eb_3d, nextra_eb
    use cns_eb_flux_module, only : compute_eb_diffop
    integer, dimension(3), intent(in) :: lo,hi,utlo,uthi,ulo,uhi, &
         fxlo,fxhi,fylo,fyhi,fzlo,fzhi, &
         vlo,vhi,axlo,axhi,aylo,ayhi,azlo,azhi, &
         cxlo,cxhi,cylo,cyhi,czlo,czhi, &
         fglo,fghi
    real(rt), intent(inout) :: dudt(utlo(1):uthi(1),utlo(2):uthi(2),utlo(3):uthi(3),nvar)
    real(rt), intent(in   ) :: u ( ulo(1): uhi(1), ulo(2): uhi(2), ulo(3): uhi(3),nvar)
    real(rt), intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),nvar)
    real(rt), intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),nvar)
    real(rt), intent(inout) :: fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),nvar)
    integer , intent(in) ::  flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
    real(rt), intent(in) :: volfrac(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
    real(rt), intent(in) :: apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(rt), intent(in) :: apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(rt), intent(in) :: apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(rt), intent(in) :: centx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
    real(rt), intent(in) :: centy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2)
    real(rt), intent(in) :: centz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2)
    real(rt), intent(in) :: dx(3)

    integer :: qlo(3), qhi(3), dvlo(3), dvhi(3), dmlo(3), dmhi(3)
    real(rt), pointer, contiguous :: q(:,:,:,:), divc(:,:,:), dm(:,:,:,:)
    integer, parameter :: nghost = nextra_eb + nghost_plm ! 

    qlo = lo - nghost
    qhi = hi + nghost
    call amrex_allocate(q, qlo(1),qhi(1), qlo(2),qhi(2), qlo(3),qhi(3), 1,qvar)
    
    dvlo = lo-2
    dvhi = hi+2
    call amrex_allocate(divc, dvlo(1),dvhi(1), dvlo(2),dvhi(2), dvlo(3),dvhi(3))

    dmlo(1:3) = lo - 1
    dmhi(1:3) = hi + 1
    call amrex_allocate(dm, dmlo(1),dmhi(1), dmlo(2),dmhi(2), dmlo(3),dmhi(3), 1,5)

    call ctoprim(qlo, qhi, u, ulo, uhi, q, qlo, qhi)
    
    call hyp_mol_gam_eb_3d(q, qlo, qhi, lo, hi, dx, fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi,&
         flag, fglo, fghi)

    call compute_eb_diffop(lo,hi,5,dx,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,&
         dudt,utlo,uthi, q,qlo,qhi, &
         divc,dvlo,dvhi, dm,dmlo,dmhi, &
         volfrac,vlo,vhi,apx,axlo,axhi,apy,aylo,ayhi,apz,azlo,azhi, &
         centx(:,:,:,1),cxlo,cxhi, centx(:,:,:,2),cxlo,cxhi, &
         centy(:,:,:,1),cylo,cyhi, centy(:,:,:,2),cylo,cyhi, &
         centz(:,:,:,1),czlo,czhi, centz(:,:,:,2),czlo,czhi, &
         flag,fglo,fghi)

    dudt(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),6:nvar) = 0.d0

    call amrex_deallocate(q)
    call amrex_deallocate(divc)
    call amrex_deallocate(dm)
  end subroutine cns_eb_compute_hydro_flux

  subroutine ctoprim (lo,hi,u,ulo,uhi,q,qlo,qhi)
    use cns_physics_module, only : gamma, cv
    use cns_module, only : smallp, smallr
    integer, intent(in) :: lo(3),hi(3),ulo(3),uhi(3),qlo(3),qhi(3)
    real(rt), intent(in   ) :: u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),nvar)
    real(rt), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),qvar)
    
    integer :: i,j,k
    real(rt) :: rhoinv, kineng
    
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             q(i,j,k,qrho) = max(smallr,u(i,j,k,urho))
             rhoinv = 1.d0/q(i,j,k,qrho)
             q(i,j,k,qu) = u(i,j,k,umx)*rhoinv
             q(i,j,k,qv) = u(i,j,k,umy)*rhoinv
             q(i,j,k,qw) = u(i,j,k,umz)*rhoinv
             kineng = 0.5d0*q(i,j,k,qrho)*(q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2)
             q(i,j,k,qeint) = (u(i,j,k,ueden)-kineng) * rhoinv
             if (q(i,j,k,qeint) .le. 0.d0) then
                q(i,j,k,qeint) = u(i,j,k,ueint) * rhoinv
             end if
             q(i,j,k,qp) = max(smallp,(gamma-1.d0)*q(i,j,k,qeint)*q(i,j,k,qrho))
             q(i,j,k,qc) = sqrt(gamma*q(i,j,k,qp)*rhoinv)
             q(i,j,k,qtemp) = q(i,j,k,qeint) * (1.d0/cv)
          end do
       end do
    end do
  end subroutine ctoprim

end module cns_hydro_module
