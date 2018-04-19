
module cns_dudt_module

  use amrex_fort_module, only : rt=>amrex_real
  use cns_module, only : urho, umx, umy, umz, ueden, ueint, utemp, nvar, qvar, &
       actually_2D
  use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
  implicit none
  private

  integer, parameter :: nghost_plm = 2  ! number of ghost cells needed for plm

  public:: cns_compute_dudt, cns_eb_compute_dudt

contains

  subroutine cns_compute_dudt (lo,hi, dudt, utlo, uthi, &
       u,ulo,uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx,dt,level) &
       bind(c,name='cns_compute_dudt')
    use cns_nd_module, only : ctoprim
    use advection_module, only : hyp_mol_gam_3d
    use diff_coef_module, only : compute_diff_coef
    use diffusion_module, only : diff_mol_3d
    use cns_divop_module, only : compute_divop
    integer, dimension(3), intent(in) :: lo,hi,utlo,uthi,ulo,uhi,fxlo,fxhi,fylo,fyhi,fzlo,fzhi
    real(rt), intent(inout) :: dudt(utlo(1):uthi(1),utlo(2):uthi(2),utlo(3):uthi(3),nvar)
    real(rt), intent(in   ) :: u ( ulo(1): uhi(1), ulo(2): uhi(2), ulo(3): uhi(3),nvar)
    real(rt), intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),nvar)
    real(rt), intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),nvar)
    real(rt), intent(inout) :: fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),nvar)
    real(rt), intent(in) :: dx(3), dt
    integer, intent(in) :: level

    integer :: qlo(3), qhi(3)
    integer :: clo(3), chi(3)
    real(rt), dimension(:,:,:,:), pointer, contiguous :: q, fhx,fhy,fhz,fdx,fdy,fdz
    real(rt), dimension(:,:,:), pointer, contiguous :: lambda, mu, xi
    integer :: k,n

    qlo = lo - nghost_plm
    qhi = hi + nghost_plm
    call amrex_allocate(q, qlo(1),qhi(1), qlo(2),qhi(2), qlo(3),qhi(3), 1,qvar)
    call amrex_allocate(fhx,lo(1),hi(1)+1,lo(2),hi(2)  ,lo(3),hi(3)  ,1,5)
    call amrex_allocate(fhy,lo(1),hi(1)  ,lo(2),hi(2)+1,lo(3),hi(3)  ,1,5)
    call amrex_allocate(fhz,lo(1),hi(1)  ,lo(2),hi(2)  ,lo(3),hi(3)+1,1,5)
    call amrex_allocate(fdx,lo(1),hi(1)+1,lo(2),hi(2)  ,lo(3),hi(3)  ,1,5)
    call amrex_allocate(fdy,lo(1),hi(1)  ,lo(2),hi(2)+1,lo(3),hi(3)  ,1,5)
    call amrex_allocate(fdz,lo(1),hi(1)  ,lo(2),hi(2)  ,lo(3),hi(3)+1,1,5)

    clo = lo - 1
    chi = hi + 1
    call amrex_allocate(lambda, clo, chi)
    call amrex_allocate(mu, clo, chi)
    call amrex_allocate(xi, clo, chi)
    
    call ctoprim(qlo, qhi, u, ulo, uhi, q, qlo, qhi)
    
    call hyp_mol_gam_3d(q, qlo, qhi, lo, hi, dx, fhx, fhy, fhz)

    call compute_diff_coef(q, qlo, qhi, lambda, mu, xi, clo, chi)

    call diff_mol_3d(lo, hi, dx, q, qlo, qhi, lambda, mu, xi, clo, chi, fdx, fdy, fdz)

    fx(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3)  ,1:5) = fhx + fdx
    fy(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3)  ,1:5) = fhy + fdy
    fz(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1,1:5) = fhz + fdz

    fx(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3)  ,6:) = 0.d0
    fy(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3)  ,6:) = 0.d0
    fz(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1,6:) = 0.d0

    call compute_divop (lo,hi,5,dx,dudt,utlo,uthi, &
         fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi)

    dudt(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),6:nvar) = 0.d0

    if (actually_2D) then
       do n = 1, 5
          do k = lo(3), hi(3)
             dudt(lo(1):hi(1),lo(2):hi(2),k,n) = dudt(lo(1):hi(1),lo(2):hi(2),(lo(3)+hi(3))/2,n)
             fx(lo(1):hi(1)+1,lo(2):hi(2)  ,k,n) = fx(lo(1):hi(1)+1,lo(2):hi(2)  ,(lo(3)+hi(3))/2,n)
             fy(lo(1):hi(1)  ,lo(2):hi(2)+1,k,n) = fy(lo(1):hi(1)  ,lo(2):hi(2)+1,(lo(3)+hi(3))/2,n)
          end do
       end do
       fz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,:) = 0.d0
    end if

    call amrex_deallocate(xi)
    call amrex_deallocate(mu)
    call amrex_deallocate(lambda)
    call amrex_deallocate(fdz)
    call amrex_deallocate(fdy)
    call amrex_deallocate(fdx)
    call amrex_deallocate(fhz)
    call amrex_deallocate(fhy)
    call amrex_deallocate(fhx)
    call amrex_deallocate(q)
  end subroutine cns_compute_dudt

  subroutine cns_eb_compute_dudt (lo,hi, dudt, utlo, uthi, &
       u,ulo,uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,flag,fglo,fghi, &
       volfrac,vlo,vhi, bcent,blo,bhi, &
       apx,axlo,axhi,apy,aylo,ayhi,apz,azlo,azhi, &
       centx,cxlo,cxhi,centy,cylo,cyhi,centz,czlo,czhi, &
       as_crse_in, rr_drho_crse, rdclo, rdchi, rr_flag_crse, rfclo, rfchi, &
       as_fine_in, dm_as_fine, dflo, dfhi, &
       levmsk, lmlo, lmhi, &
       dx,dt,level) &
       bind(c,name='cns_eb_compute_dudt')
    use cns_nd_module, only : ctoprim
    use eb_advection_module, only : hyp_mol_gam_eb_3d, nextra_eb
    use diff_coef_module, only : compute_diff_coef
    use eb_diffusion_module, only : eb_diff_mol_3d
    use cns_divop_module, only : compute_eb_divop
    integer, dimension(3), intent(in) :: lo,hi,utlo,uthi,ulo,uhi, &
         vlo,vhi,axlo,axhi,aylo,ayhi,azlo,azhi, &
         cxlo,cxhi,cylo,cyhi,czlo,czhi, &
         fglo,fghi, blo, bhi, fxlo,fxhi, fylo,fyhi, fzlo,fzhi, &
         rdclo, rdchi, rfclo, rfchi, dflo, dfhi, lmlo, lmhi
    integer, intent(in) :: as_crse_in, as_fine_in
    real(rt), intent(inout) :: dudt(utlo(1):uthi(1),utlo(2):uthi(2),utlo(3):uthi(3),nvar)
    real(rt), intent(in   ) :: u ( ulo(1): uhi(1), ulo(2): uhi(2), ulo(3): uhi(3),nvar)
    real(rt), intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),nvar)
    real(rt), intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),nvar)
    real(rt), intent(inout) :: fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),nvar)
    integer , intent(in) ::  flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))
    real(rt), intent(in) :: volfrac(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
    real(rt), intent(in) :: bcent  (blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3)
    real(rt), intent(in) :: apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(rt), intent(in) :: apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(rt), intent(in) :: apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(rt), intent(in) :: centx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
    real(rt), intent(in) :: centy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2)
    real(rt), intent(in) :: centz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2)
    real(rt), intent(inout) :: rr_drho_crse(rdclo(1):rdchi(1),rdclo(2):rdchi(2),rdclo(3):rdchi(3),nvar)
    integer,  intent(in) ::  rr_flag_crse(rfclo(1):rfchi(1),rfclo(2):rfchi(2),rfclo(3):rfchi(3))
    real(rt), intent(out) :: dm_as_fine(dflo(1):dfhi(1),dflo(2):dfhi(2),dflo(3):dfhi(3),nvar)
    integer,  intent(in) ::  levmsk (lmlo(1):lmhi(1),lmlo(2):lmhi(2),lmlo(3):lmhi(3))
    real(rt), intent(in) :: dx(3), dt
    integer, intent(in) :: level

    integer :: qlo(3), qhi(3), dvlo(3), dvhi(3), dmlo(3), dmhi(3)
    integer :: lfxlo(3), lfylo(3), lfzlo(3), lfxhi(3), lfyhi(3), lfzhi(3)
    integer :: clo(3), chi(3)
    real(rt), pointer, contiguous :: q(:,:,:,:)
    real(rt), dimension(:,:,:), pointer, contiguous :: divc, dm, optmp, rediswgt
    real(rt), dimension(:,:,:,:), pointer, contiguous :: fhx,fhy,fhz,fdx,fdy,fdz
    real(rt), dimension(:,:,:), pointer, contiguous :: lambda, mu, xi
    integer, parameter :: nghost = nextra_eb + max(nghost_plm,3) ! 3 because of wall flux

    integer :: k,n
    logical :: as_crse, as_fine

    as_crse = as_crse_in .ne. 0
    as_fine = as_fine_in .ne. 0

    qlo = lo - nghost
    qhi = hi + nghost
    call amrex_allocate(q, qlo(1),qhi(1), qlo(2),qhi(2), qlo(3),qhi(3), 1,qvar)
    
    dvlo = lo-2
    dvhi = hi+2
    call amrex_allocate(divc, dvlo, dvhi)
    call amrex_allocate(optmp, dvlo, dvhi) 
    call amrex_allocate(rediswgt, dvlo, dvhi)

    dmlo(1:3) = lo - 1
    dmhi(1:3) = hi + 1
    call amrex_allocate(dm, dmlo, dmhi)

    lfxlo = lo - nextra_eb - 1;  lfxlo(1) = lo(1)-nextra_eb
    lfxhi = hi + nextra_eb + 1
    call amrex_allocate(fhx, lfxlo(1),lfxhi(1),lfxlo(2),lfxhi(2),lfxlo(3),lfxhi(3),1,5)
    call amrex_allocate(fdx, lfxlo(1),lfxhi(1),lfxlo(2),lfxhi(2),lfxlo(3),lfxhi(3),1,5)

    lfylo = lo - nextra_eb - 1;  lfylo(2) = lo(2)-nextra_eb
    lfyhi = hi + nextra_eb + 1
    call amrex_allocate(fhy, lfylo(1),lfyhi(1),lfylo(2),lfyhi(2),lfylo(3),lfyhi(3),1,5)
    call amrex_allocate(fdy, lfylo(1),lfyhi(1),lfylo(2),lfyhi(2),lfylo(3),lfyhi(3),1,5)

    lfzlo = lo - nextra_eb - 1;  lfzlo(3) = lo(3)-nextra_eb
    lfzhi = hi + nextra_eb + 1
    call amrex_allocate(fhz, lfzlo(1),lfzhi(1),lfzlo(2),lfzhi(2),lfzlo(3),lfzhi(3),1,5)
    call amrex_allocate(fdz, lfzlo(1),lfzhi(1),lfzlo(2),lfzhi(2),lfzlo(3),lfzhi(3),1,5)

    clo = lo - nextra_eb - 1
    chi = hi + nextra_eb + 1
    call amrex_allocate(lambda, clo, chi)
    call amrex_allocate(mu, clo, chi)
    call amrex_allocate(xi, clo, chi)

    call bl_proffortfuncstart_int(6)
    call ctoprim(qlo, qhi, u, ulo, uhi, q, qlo, qhi)
    call bl_proffortfuncstop_int(6)

    
    call bl_proffortfuncstart_int(7)
    call hyp_mol_gam_eb_3d(q, qlo, qhi, lo, hi, dx, &
         fhx, lfxlo, lfxhi, fhy, lfylo, lfyhi, fhz, lfzlo, lfzhi,&
         flag, fglo, fghi)
    call bl_proffortfuncstop_int(7)

    call bl_proffortfuncstart_int(8)
    call compute_diff_coef(q, qlo, qhi, lambda, mu, xi, clo, chi)
    call bl_proffortfuncstop_int(8)

    call bl_proffortfuncstart_int(9)
    call eb_diff_mol_3d(lo, hi, dx, q, qlo, qhi, lambda, mu, xi, clo, chi, &
         fdx, lfxlo, lfxhi, fdy, lfylo, lfyhi, fdz, lfzlo, lfzhi,&
         flag, fglo, fghi)
    call bl_proffortfuncstop_int(9)

    fhx = fhx + fdx
    fhy = fhy + fdy
    fhz = fhz + fdz

    fx(      fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),1:5) = &
         fhx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),1:5)
    fy(      fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),1:5) = &
         fhy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),1:5)
    fz(      fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),1:5) = &
         fhz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),1:5)

    fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),6:nvar) = 0.d0
    fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),6:nvar) = 0.d0
    fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),6:nvar) = 0.d0

    dm_as_fine = 0.d0

    call bl_proffortfuncstart_int(10)
    call compute_eb_divop(lo,hi,5,dx,dt,fhx,lfxlo,lfxhi,fhy,lfylo,lfyhi,fhz,lfzlo,lfzhi,&
         fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi, &
         dudt,utlo,uthi, q,qlo,qhi, lambda, mu, xi, clo, chi, &
         divc, optmp, rediswgt, dvlo,dvhi, &
         dm,dmlo,dmhi, &
         volfrac,vlo,vhi, bcent,blo,bhi, &
         apx,axlo,axhi,apy,aylo,ayhi,apz,azlo,azhi, &
         centx(:,:,:,1),cxlo,cxhi, centx(:,:,:,2),cxlo,cxhi, &
         centy(:,:,:,1),cylo,cyhi, centy(:,:,:,2),cylo,cyhi, &
         centz(:,:,:,1),czlo,czhi, centz(:,:,:,2),czlo,czhi, &
         flag,fglo,fghi, &
         as_crse, rr_drho_crse, rdclo, rdchi, rr_flag_crse, rfclo, rfchi, &
         as_fine, dm_as_fine, dflo, dfhi, &
         levmsk, lmlo, lmhi)
    
    dudt(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),6:nvar) = 0.d0
    call bl_proffortfuncstop_int(10)

    if (actually_2D) then
       do n = 1, 5
          do k = lo(3), hi(3)
             dudt(lo(1):hi(1),lo(2):hi(2),k,n) = dudt(lo(1):hi(1),lo(2):hi(2),(lo(3)+hi(3))/2,n)
             fx(lo(1):hi(1)+1,lo(2):hi(2)  ,k,n) = fx(lo(1):hi(1)+1,lo(2):hi(2)  ,(lo(3)+hi(3))/2,n)
             fy(lo(1):hi(1)  ,lo(2):hi(2)+1,k,n) = fy(lo(1):hi(1)  ,lo(2):hi(2)+1,(lo(3)+hi(3))/2,n)
          end do
       end do
       fz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1,:) = 0.d0
    end if

    call amrex_deallocate(xi)
    call amrex_deallocate(mu)
    call amrex_deallocate(lambda)
    call amrex_deallocate(fdz)
    call amrex_deallocate(fhz)
    call amrex_deallocate(fdy)
    call amrex_deallocate(fhy)
    call amrex_deallocate(fdx)
    call amrex_deallocate(fhx)
    call amrex_deallocate(dm)
    call amrex_deallocate(rediswgt)
    call amrex_deallocate(optmp)
    call amrex_deallocate(divc)
    call amrex_deallocate(q)
  end subroutine cns_eb_compute_dudt

end module cns_dudt_module
