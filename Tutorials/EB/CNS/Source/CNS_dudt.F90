
module cns_dudt_module

  use amrex_fort_module, only : rt=>amrex_real
  use cns_module, only : urho, umx, umy, umz, ueden, ueint, utemp, nvar, qvar
  use mempool_module, only : amrex_allocate, amrex_deallocate
  implicit none
  private

  integer, parameter :: nghost_plm = 2  ! number of ghost cells needed for plm

  public:: cns_compute_dudt, cns_eb_compute_dudt

contains

  subroutine cns_compute_dudt (lo,hi, dudt, utlo, uthi, &
       u,ulo,uhi,dx,dt) &
       bind(c,name='cns_compute_dudt')
    use cns_nd_module, only : ctoprim
    use advection_module, only : hyp_mol_gam_3d
    use diffusion_module, only : diff_mol_3d
    use cns_eb_flux_module, only : compute_diffop
    integer, dimension(3), intent(in) :: lo,hi,utlo,uthi,ulo,uhi
    real(rt), intent(inout) :: dudt(utlo(1):uthi(1),utlo(2):uthi(2),utlo(3):uthi(3),nvar)
    real(rt), intent(in   ) :: u ( ulo(1): uhi(1), ulo(2): uhi(2), ulo(3): uhi(3),nvar)
    real(rt), intent(in) :: dx(3), dt

    integer :: qlo(3), qhi(3)
    real(rt), dimension(:,:,:,:), pointer, contiguous :: q, fhx,fhy,fhz,fdx,fdy,fdz

    qlo = lo - nghost_plm
    qhi = hi + nghost_plm
    call amrex_allocate(q, qlo(1),qhi(1), qlo(2),qhi(2), qlo(3),qhi(3), 1,qvar)
    call amrex_allocate(fhx,lo(1),hi(1)+1,lo(2),hi(2)  ,lo(3),hi(3)  ,1,5)
    call amrex_allocate(fhy,lo(1),hi(1)  ,lo(2),hi(2)+1,lo(3),hi(3)  ,1,5)
    call amrex_allocate(fhz,lo(1),hi(1)  ,lo(2),hi(2)  ,lo(3),hi(3)+1,1,5)
    call amrex_allocate(fdx,lo(1),hi(1)+1,lo(2),hi(2)  ,lo(3),hi(3)  ,1,5)
    call amrex_allocate(fdy,lo(1),hi(1)  ,lo(2),hi(2)+1,lo(3),hi(3)  ,1,5)
    call amrex_allocate(fdz,lo(1),hi(1)  ,lo(2),hi(2)  ,lo(3),hi(3)+1,1,5)
    
    call ctoprim(qlo, qhi, u, ulo, uhi, q, qlo, qhi)
    
    call hyp_mol_gam_3d(q, qlo, qhi, lo, hi, dx, fhx, fhy, fhz)

    call diff_mol_3d(q, qlo, qhi, lo, hi, dx, fdx, fdy, fdz)

    fhx = fhx + fdx
    fhy = fhy + fdy
    fhz = fhz + fdz
    call compute_diffop (lo,hi,5,dx,dudt,utlo,uthi, &
         fhx, lo, [hi(1)+1,hi(2)  ,hi(3)  ], &
         fhy, lo, [hi(1)  ,hi(2)+1,hi(3)  ], &
         fhz, lo, [hi(1)  ,hi(2)  ,hi(3)+1])

    dudt(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),6:nvar) = 0.d0

    call amrex_deallocate(fdz)
    call amrex_deallocate(fdy)
    call amrex_deallocate(fdx)
    call amrex_deallocate(fhz)
    call amrex_deallocate(fhy)
    call amrex_deallocate(fhx)
    call amrex_deallocate(q)
  end subroutine cns_compute_dudt

  subroutine cns_eb_compute_dudt (lo,hi, dudt, utlo, uthi, &
       u,ulo,uhi,flag,fglo,fghi, &
       volfrac,vlo,vhi,apx,axlo,axhi,apy,aylo,ayhi,apz,azlo,azhi, &
       centx,cxlo,cxhi,centy,cylo,cyhi,centz,czlo,czhi, dx,dt) &
       bind(c,name='cns_eb_compute_dudt')
    use cns_nd_module, only : ctoprim
    use eb_advection_module, only : hyp_mol_gam_eb_3d, nextra_eb
    use eb_diffusion_module, only : eb_diff_mol_3d
    use cns_eb_flux_module, only : compute_eb_diffop
    integer, dimension(3), intent(in) :: lo,hi,utlo,uthi,ulo,uhi, &
         vlo,vhi,axlo,axhi,aylo,ayhi,azlo,azhi, &
         cxlo,cxhi,cylo,cyhi,czlo,czhi, &
         fglo,fghi
    real(rt), intent(inout) :: dudt(utlo(1):uthi(1),utlo(2):uthi(2),utlo(3):uthi(3),nvar)
    real(rt), intent(in   ) :: u ( ulo(1): uhi(1), ulo(2): uhi(2), ulo(3): uhi(3),nvar)
    integer , intent(in) ::  flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
    real(rt), intent(in) :: volfrac(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
    real(rt), intent(in) :: apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(rt), intent(in) :: apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(rt), intent(in) :: apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(rt), intent(in) :: centx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
    real(rt), intent(in) :: centy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2)
    real(rt), intent(in) :: centz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2)
    real(rt), intent(in) :: dx(3), dt

    integer :: qlo(3), qhi(3), dvlo(3), dvhi(3), dmlo(3), dmhi(3)
    integer :: fxlo(3), fylo(3), fzlo(3), fxhi(3), fyhi(3), fzhi(3)
    real(rt), pointer, contiguous :: q(:,:,:,:), divc(:,:,:), dm(:,:,:,:)
    real(rt), dimension(:,:,:,:), pointer :: fhx,fhy,fhz,fdx,fdy,fdz
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

    fxlo = lo - nextra_eb - 1;  fxlo(1) = lo(1)-nextra_eb
    fxhi = hi + nextra_eb + 1
    call amrex_allocate(fhx, fxlo(1),fxhi(1),fxlo(2),fxhi(2),fxlo(3),fxhi(3),1,5)
    call amrex_allocate(fdx, fxlo(1),fxhi(1),fxlo(2),fxhi(2),fxlo(3),fxhi(3),1,5)

    fylo = lo - nextra_eb - 1;  fylo(2) = lo(2)-nextra_eb
    fyhi = hi + nextra_eb + 1
    call amrex_allocate(fhy, fylo(1),fyhi(1),fylo(2),fyhi(2),fylo(3),fyhi(3),1,5)
    call amrex_allocate(fdy, fylo(1),fyhi(1),fylo(2),fyhi(2),fylo(3),fyhi(3),1,5)

    fzlo = lo - nextra_eb - 1;  fzlo(3) = lo(3)-nextra_eb
    fzhi = hi + nextra_eb + 1
    call amrex_allocate(fhz, fzlo(1),fzhi(1),fzlo(2),fzhi(2),fzlo(3),fzhi(3),1,5)
    call amrex_allocate(fdz, fzlo(1),fzhi(1),fzlo(2),fzhi(2),fzlo(3),fzhi(3),1,5)

    call ctoprim(qlo, qhi, u, ulo, uhi, q, qlo, qhi)
    
    call hyp_mol_gam_eb_3d(q, qlo, qhi, lo, hi, dx, fhx, fxlo, fxhi, fhy, fylo, fyhi, fhz, fzlo, fzhi,&
         flag, fglo, fghi)

    call eb_diff_mol_3d(q, qlo, qhi, lo, hi, dx, fdx, fxlo, fxhi, fdy, fylo, fyhi, fdz, fzlo, fzhi,&
         flag, fglo, fghi)

    fhx = fhx + fdx
    fhy = fhy + fdy
    fhz = fhz + fdz
    call compute_eb_diffop(lo,hi,5,dx,dt,fhx,fxlo,fxhi,fhy,fylo,fyhi,fhz,fzlo,fzhi,&
         dudt,utlo,uthi, q,qlo,qhi, &
         divc,dvlo,dvhi, dm,dmlo,dmhi, &
         volfrac,vlo,vhi,apx,axlo,axhi,apy,aylo,ayhi,apz,azlo,azhi, &
         centx(:,:,:,1),cxlo,cxhi, centx(:,:,:,2),cxlo,cxhi, &
         centy(:,:,:,1),cylo,cyhi, centy(:,:,:,2),cylo,cyhi, &
         centz(:,:,:,1),czlo,czhi, centz(:,:,:,2),czlo,czhi, &
         flag,fglo,fghi)

    dudt(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),6:nvar) = 0.d0

    call amrex_deallocate(fdz)
    call amrex_deallocate(fhz)
    call amrex_deallocate(fdy)
    call amrex_deallocate(fhy)
    call amrex_deallocate(fdx)
    call amrex_deallocate(fhx)
    call amrex_deallocate(dm)
    call amrex_deallocate(divc)
    call amrex_deallocate(q)
  end subroutine cns_eb_compute_dudt

end module cns_dudt_module
