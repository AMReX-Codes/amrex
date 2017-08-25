module cns_eb_flux_module
  use mempool_module, only : amrex_allocate, amrex_deallocate
  use amrex_fort_module, only : rt=>amrex_real
  use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell, is_single_valued_cell, &
       get_neighbor_cells
  implicit none
  private
  public :: compute_diffop, compute_eb_diffop
contains

  ! non-eb version
  subroutine compute_diffop (lo,hi,ncomp,dx,ut,ulo,uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi)
    integer, intent(in) :: lo(3),hi(3),ncomp,ulo(3),uhi(3),fxlo(3),fxhi(3), &
         fylo(3),fyhi(3),fzlo(3),fzhi(3)
    real(rt), intent(inout) :: ut( ulo(1): uhi(1), ulo(2): uhi(2), ulo(3): uhi(3),ncomp)
    real(rt), intent(in   ) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),ncomp)
    real(rt), intent(in   ) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),ncomp)
    real(rt), intent(in   ) :: fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),ncomp)
    real(rt), intent(in) :: dx(3)
    
    integer :: i,j,k,n
    real(rt) :: dxinv(3)
    
    dxinv = 1.d0/dx
    
    do n = 1, ncomp
       do       k = lo(3),hi(3)
          do    j = lo(2),hi(2)
             do i = lo(1),hi(1)
                ut(i,j,k,n) = (fx(i,j,k,n)-fx(i+1,j,k,n)) * dxinv(1) &
                     +        (fy(i,j,k,n)-fy(i,j+1,k,n)) * dxinv(2) &
                     +        (fz(i,j,k,n)-fz(i,j,k+1,n)) * dxinv(3)
             end do
          end do
       end do
    end do
  end subroutine compute_diffop

  subroutine compute_eb_diffop (lo,hi,ncomp, dx, &
       fluxx,fxlo,fxhi, &
       fluxy,fylo,fyhi, &
       fluxz,fzlo,fzhi, &
       ebdiffop, oplo, ophi, &
       q, qlo, qhi, &
       divc, dvlo, dvhi, &
       delm, dmlo, dmhi, &
       vfrac, vlo, vhi, &
       apx, axlo, axhi, &
       apy, aylo, ayhi, &
       apz, azlo, azhi, &
       centx_y, cxylo, cxyhi, &
       centx_z, cxzlo, cxzhi, &
       centy_x, cyxlo, cyxhi, &
       centy_z, cyzlo, cyzhi, &
       centz_x, czxlo, czxhi, &
       centz_y, czylo, czyhi, &
       cellflag, cflo, cfhi)
    use cns_module, only : qvar, qrho, qu, qv, qw, qp, umx, umy, umz
    use cns_eb_hyp_wall_module, only : compute_hyp_wallflux
    integer, intent(in), dimension(3) :: lo, hi, fxlo,fxhi,fylo,fyhi,fzlo,fzhi,oplo,ophi,&
         dvlo,dvhi,dmlo,dmhi,axlo,axhi,aylo,ayhi,azlo,azhi,cxylo,cxyhi,cxzlo,cxzhi,&
         cyxlo,cyxhi,cyzlo,cyzhi,czxlo,czxhi,czylo,czyhi,vlo,vhi,cflo,cfhi, qlo,qhi
    integer, intent(in) :: ncomp
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in) :: fluxx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),ncomp)
    real(rt), intent(in) :: fluxy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),ncomp)
    real(rt), intent(in) :: fluxz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),ncomp)
    real(rt), intent(inout) :: ebdiffop(oplo(1):ophi(1),oplo(2):ophi(2),oplo(3):ophi(3),ncomp)
    real(rt), intent(in) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),qvar)
    real(rt), intent(inout) :: divc(dvlo(1):dvhi(1),dvlo(2):dvhi(2),dvlo(3):dvhi(3))
    real(rt), intent(inout) :: delm(dmlo(1):dmhi(1),dmlo(2):dmhi(2),dmlo(3):dmhi(3),ncomp)
    real(rt), intent(in) :: vfrac(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
    real(rt), intent(in) :: apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(rt), intent(in) :: apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(rt), intent(in) :: apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(rt), intent(in) :: centx_y(cxylo(1):cxyhi(1),cxylo(2):cxyhi(2),cxylo(3):cxyhi(3))
    real(rt), intent(in) :: centx_z(cxzlo(1):cxzhi(1),cxzlo(2):cxzhi(2),cxzlo(3):cxzhi(3))
    real(rt), intent(in) :: centy_x(cyxlo(1):cyxhi(1),cyxlo(2):cyxhi(2),cyxlo(3):cyxhi(3))
    real(rt), intent(in) :: centy_z(cyzlo(1):cyzhi(1),cyzlo(2):cyzhi(2),cyzlo(3):cyzhi(3))
    real(rt), intent(in) :: centz_x(czxlo(1):czxhi(1),czxlo(2):czxhi(2),czxlo(3):czxhi(3))
    real(rt), intent(in) :: centz_y(czylo(1):czyhi(1),czylo(2):czyhi(2),czylo(3):czyhi(3))
    integer, intent(in) :: cellflag(cflo(1):cfhi(1),cflo(2):cfhi(2),cflo(3):cfhi(3))

    integer :: i,j,k,n,ii,jj,kk, nbr(-1:1,-1:1,-1:1)
    real(rt) :: fxp,fxm,fyp,fym,fzp,fzm,divnc, vtot,mtot, fracx,fracy,fracz,dxinv(3)
    real(rt) :: divw(5), divwn
    real(rt), pointer, contiguous :: optmp(:,:,:)

    !  centroid nondimensional  and zero at face center

    dxinv = 1.d0/dx

    call amrex_allocate(optmp, lo-2, hi+2)

    do n = 1, ncomp

       !
       ! First, we compute conservative divergence on (lo-2,hi+2)
       !
       do       k = lo(3)-2, hi(3)+2
          do    j = lo(2)-2, hi(2)+2
             do i = lo(1)-2, hi(1)+2
                divc(i,j,k) = (fluxx(i,j,k,n)-fluxx(i+1,j,k,n))*dxinv(1) &
                     +        (fluxy(i,j,k,n)-fluxy(i,j+1,k,n))*dxinv(2) &
                     +        (fluxz(i,j,k,n)-fluxz(i,j,k+1,n))*dxinv(3)
             end do

             do i = lo(1)-2, hi(1)+2
                if (is_covered_cell(cellflag(i,j,k))) then
                   divc(i,j,k) = 0.d0
                else if (is_single_valued_cell(cellflag(i,j,k))) then

                   call get_neighbor_cells(cellflag(i,j,k),nbr)

                   ! x-direction lo face
                   if (apx(i,j,k).lt.1.d0) then
                      if (centx_y(i,j,k).le.0.d0) then
                         fracy = -centx_y(i,j,k)*nbr(0,-1,0)
                         if(centx_z(i,j,k).le. 0.0d0)then
                            fracz = - centx_z(i,j,k)*nbr(0,0,-1)
                            fxm = (1.d0-fracz)*(     fracy *fluxx(i,j-1,k  ,n)  + &
                                 &             (1.d0-fracy)*fluxx(i,j  ,k  ,n)) + &
                                 &      fracz *(     fracy *fluxx(i,j-1,k-1,n)  + &
                                 &             (1.d0-fracy)*fluxx(i,j  ,k-1,n))
                         else
                            fracz =  centx_z(i,j,k)*nbr(0,0,1)
                            fxm = (1.d0-fracz)*(     fracy *fluxx(i,j-1,k  ,n)  + &
                                 &             (1.d0-fracy)*fluxx(i,j  ,k  ,n)) + &
                                 &      fracz *(     fracy *fluxx(i,j-1,k+1,n)  + &
                                 &             (1.d0-fracy)*fluxx(i,j  ,k+1,n))
                         endif
                      else
                         fracy = centx_y(i,j,k)*nbr(0,1,0)
                         if(centx_z(i,j,k).le. 0.0d0)then
                            fracz = -centx_z(i,j,k)*nbr(0,0,-1)
                            fxm = (1.d0-fracz)*(     fracy *fluxx(i,j+1,k  ,n)  + &
                                 &             (1.d0-fracy)*fluxx(i,j  ,k  ,n)) + &
                                 &      fracz *(     fracy *fluxx(i,j+1,k-1,n)  + &
                                 &             (1.d0-fracy)*fluxx(i,j  ,k-1,n))
                         else
                            fracz = centx_z(i,j,k)*nbr(0,0,1)
                            fxm = (1.d0-fracz)*(     fracy *fluxx(i,j+1,k  ,n)  + &
                                 &             (1.d0-fracy)*fluxx(i,j  ,k  ,n)) + &
                                 &      fracz *(     fracy *fluxx(i,j+1,k+1,n)  + &
                                 &             (1.d0-fracy)*fluxx(i,j  ,k+1,n))
                         endif
                      end if
                   else
                      fxm = fluxx(i,j,k,n)
                   end if

                   ! x-direction hi face
                   if (apx(i+1,j,k).lt.1.d0) then
                      if (centx_y(i+1,j,k).le.0.d0) then
                         fracy = -centx_y(i+1,j,k)*nbr(0,-1,0)
                         if(centx_z(i+1,j,k).le. 0.0d0)then
                            fracz = - centx_z(i+1,j,k)*nbr(0,0,-1)
                            fxp = (1.d0-fracz)*(     fracy *fluxx(i+1,j-1,k  ,n)  + &
                                 &             (1.d0-fracy)*fluxx(i+1,j  ,k  ,n)) + &
                                 &      fracz *(     fracy *fluxx(i+1,j-1,k-1,n)  + &
                                 &             (1.d0-fracy)*fluxx(i+1,j  ,k-1,n))
                         else
                            fracz =  centx_z(i+1,j,k)*nbr(0,0,1)
                            fxp = (1.d0-fracz)*(     fracy *fluxx(i+1,j-1,k  ,n)  + &
                                 &             (1.d0-fracy)*fluxx(i+1,j  ,k  ,n)) + &
                                 &      fracz *(     fracy *fluxx(i+1,j-1,k+1,n)  + &
                                 &             (1.d0-fracy)*fluxx(i+1,j  ,k+1,n))
                         endif
                      else
                         fracy = centx_y(i+1,j,k)*nbr(0,1,0)
                         if(centx_z(i+1,j,k).le. 0.0d0)then
                            fracz = -centx_z(i+1,j,k)*nbr(0,0,-1)
                            fxp = (1.d0-fracz)*(     fracy *fluxx(i+1,j+1,k  ,n)  + &
                                 &             (1.d0-fracy)*fluxx(i+1,j  ,k  ,n)) + &
                                 &      fracz *(     fracy *fluxx(i+1,j+1,k-1,n)  + &
                                 &             (1.d0-fracy)*fluxx(i+1,j  ,k-1,n))
                         else
                            fracz = centx_z(i+1,j,k)*nbr(0,0,1)
                            fxp = (1.d0-fracz)*(     fracy *fluxx(i+1,j+1,k  ,n)  + &
                                 &             (1.d0-fracy)*fluxx(i+1,j  ,k  ,n)) + &
                                 &      fracz *(     fracy *fluxx(i+1,j+1,k+1,n)  + &
                                 &             (1.d0-fracy)*fluxx(i+1,j  ,k+1,n))
                         endif
                      end if
                   else
                      fxp = fluxx(i+1,j,k,n)
                   end if

                   ! y-direction lo face
                   if (apy(i,j,k).lt.1.d0) then
                      if (centy_x(i,j,k).le.0.d0) then
                         fracx = - centy_x(i,j,k)*nbr(-1,0,0)
                         if(centy_z(i,j,k).le. 0.0d0)then
                            fracz = - centy_z(i,j,k)*nbr(0,0,-1)
                            fym = (1.d0-fracz)*(     fracx *fluxy(i-1,j,k  ,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j,k  ,n)) + &
                                 &      fracz *(     fracx *fluxy(i-1,j,k-1,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j,k-1,n))
                         else
                            fracz =  centy_z(i,j,k)*nbr(0,0,1)
                            fym = (1.d0-fracz)*(     fracx *fluxy(i-1,j,k  ,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j,k  ,n)) + &
                                 &      fracz *(     fracx *fluxy(i-1,j,k+1,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j,k+1,n))
                         endif
                      else
                         fracx =  centy_x(i,j,k)*nbr(1,0,0)
                         if(centy_z(i,j,k).le. 0.0d0)then
                            fracz = -centx_z(i,j,k)*nbr(0,0,-1)
                            fym = (1.d0-fracz)*(     fracx *fluxy(i+1,j,k  ,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j,k  ,n)) + &
                                 &      fracz *(     fracx *fluxy(i+1,j,k-1,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j,k-1,n))
                         else
                            fracz = centx_z(i,j,k)*nbr(0,0,1)
                            fym = (1.d0-fracz)*(     fracx *fluxy(i+1,j,k  ,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j,k  ,n)) + &
                                 &      fracz *(     fracx *fluxy(i+1,j,k+1,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j,k+1,n))
                         endif
                      endif
                   else
                      fym = fluxy(i,j,k,n)
                   end if

                   ! y-direction hi face
                   if (apy(i,j+1,k).lt.1d0) then
                      if (centy_x(i,j+1,k).le.0.d0) then
                         fracx = - centy_x(i,j+1,k)*nbr(-1,0,0)
                         if(centy_z(i,j+1,k).le. 0.0d0)then
                            fracz = - centy_z(i,j+1,k)*nbr(0,0,-1)
                            fyp = (1.d0-fracz)*(     fracx *fluxy(i-1,j+1,k  ,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j+1,k  ,n)) + &
                                 &      fracz *(     fracx *fluxy(i-1,j+1,k-1,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j+1,k-1,n))
                         else
                            fracz =  centy_z(i,j+1,k)*nbr(0,0,1)
                            fyp = (1.d0-fracz)*(     fracx *fluxy(i-1,j+1,k  ,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j+1,k  ,n)) + &
                                 &      fracz *(     fracx *fluxy(i-1,j+1,k+1,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j+1,k+1,n))
                         endif
                      else
                         fracx =  centy_x(i,j+1,k)*nbr(1,0,0)
                         if(centy_z(i,j+1,k).le. 0.0d0)then
                            fracz = -centx_z(i,j+1,k)*nbr(0,0,-1)
                            fyp = (1.d0-fracz)*(     fracx *fluxy(i+1,j+1,k  ,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j+1,k  ,n)) + &
                                 &      fracz *(     fracx *fluxy(i+1,j+1,k-1,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j+1,k-1,n))
                         else
                            fracz = centx_z(i,j+1,k)*nbr(0,0,1)
                            fyp = (1.d0-fracz)*(     fracx *fluxy(i+1,j+1,k  ,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j+1,k  ,n)) + &
                                 &      fracz *(     fracx *fluxy(i+1,j+1,k+1,n)  + &
                                 &             (1.d0-fracx)*fluxy(i  ,j+1,k+1,n))
                         endif
                      endif
                   else
                      fyp = fluxy(i,j+1,k,n)
                   end if

                   ! z-direction lo face
                   if(apz(i,j,k).lt.1.d0)then
                      if(centz_x(i,j,k).le. 0.0d0)then
                         fracx = - centz_x(i,j,k)*nbr(-1,0,0)
                         if(centz_y(i,j,k).le. 0.0d0)then
                            fracy = - centz_y(i,j,k)*nbr(0,-1,0)
                            fzm = (1.d0-fracy)*(     fracx *fluxz(i-1,j  ,k,n)  + &
                                 &             (1.d0-fracx)*fluxz(i  ,j  ,k,n)) + &
                                 &      fracy* (     fracx *fluxz(i-1,j-1,k,n)  + &
                                 &             (1.d0-fracx)*fluxz(i  ,j-1,k,n))
                         else
                            fracy =  centz_y(i,j,k)*nbr(0,1,0)
                            fzm = (1.d0-fracy)*(     fracx *fluxz(i-1,j  ,k,n)  + &
                                 &             (1.d0-fracx)*fluxz(i  ,j  ,k,n)) + &
                                 &      fracy *(     fracx *fluxz(i-1,j+1,k,n)  + &
                                 &             (1.d0-fracx)*fluxz(i  ,j+1,k,n))
                         endif
                      else
                         fracx =  centz_x(i,j,k)*nbr(1,0,0)
                         if(centz_y(i,j,k).le. 0.0d0)then
                            fracy = -centz_y(i,j,k)*nbr(0,-1,0)
                            fzm = (1.d0-fracy)*(     fracx *fluxz(i+1,j  ,k,n)  + &
                                 &             (1.d0-fracx)*fluxz(i  ,j  ,k,n)) + &
                                 &      fracy *(     fracx *fluxz(i+1,j-1,k,n)  + &
                                 &             (1.d0-fracx)*fluxz(i  ,j-1,k,n))
                         else
                            fracy = centz_y(i,j,k)*nbr(0,1,0)
                            fzm = (1.d0-fracy)*(     fracx *fluxz(i+1,j  ,k,n)+ &
                                 &             (1.d0-fracx)*fluxz(i  ,j  ,k,n)) + &
                                 &      fracy* (     fracx *fluxz(i+1,j+1,k,n)+ &
                                 &             (1.d0-fracx)*fluxz(i  ,j+1,k,n))
                         endif
                      endif
                   else
                      fzm = fluxz(i,j,k,n)
                   endif

                   ! z-direction hi face
                   if(apz(i,j,k+1).lt.1.d0)then
                      if(centz_x(i,j,k+1).le. 0.0d0)then
                         fracx = - centz_x(i,j,k+1)*nbr(-1,0,0)
                         if(centz_y(i,j,k+1).le. 0.0d0)then
                            fracy = - centz_y(i,j,k+1)*nbr(0,-1,0)
                            fzp = (1.d0-fracy)*(     fracx *fluxz(i-1,j  ,k+1,n)  + &
                                 &             (1.d0-fracx)*fluxz(i  ,j  ,k+1,n)) + &
                                 &      fracy* (     fracx *fluxz(i-1,j-1,k+1,n)  + &
                                 &             (1.d0-fracx)*fluxz(i  ,j-1,k+1,n))
                         else
                            fracy =  centz_y(i,j,k+1)*nbr(0,1,0)
                            fzp = (1.d0-fracy)*(     fracx *fluxz(i-1,j  ,k+1,n)  + &
                                 &             (1.d0-fracx)*fluxz(i  ,j  ,k+1,n)) + &
                                 &      fracy *(     fracx *fluxz(i-1,j+1,k+1,n)  + &
                                 &             (1.d0-fracx)*fluxz(i  ,j+1,k+1,n))
                         endif
                      else
                         fracx =  centz_x(i,j,k+1)*nbr(1,0,0)
                         if(centz_y(i,j,k+1).le. 0.0d0)then
                            fracy = -centz_y(i,j,k+1)*nbr(0,-1,0)
                            fzp = (1.d0-fracy)*(     fracx *fluxz(i+1,j  ,k+1,n)  + &
                                 &             (1.d0-fracx)*fluxz(i  ,j  ,k+1,n)) + &
                                 &      fracy *(     fracx *fluxz(i+1,j-1,k+1,n)  + &
                                 &             (1.d0-fracx)*fluxz(i  ,j-1,k+1,n))
                         else
                            fracy = centz_y(i,j,k+1)*nbr(0,1,0)
                            fzp = (1.d0-fracy)*(     fracx *fluxz(i+1,j  ,k+1,n)+ &
                                 &             (1.d0-fracx)*fluxz(i  ,j  ,k+1,n)) + &
                                 &      fracy* (     fracx *fluxz(i+1,j+1,k+1,n)+ &
                                 &             (1.d0-fracx)*fluxz(i  ,j+1,k+1,n))
                         endif
                      endif
                   else
                      fzp = fluxz(i,j,k+1,n)
                   endif

                   if (n .eq. umx .or. n .eq. umy .or. n .eq. umz) then
                      call compute_hyp_wallflux(divw, q(i,j,k,qrho), q(i,j,k,qu), &
                           q(i,j,k,qv), q(i,j,k,qw), q(i,j,k,qp), &
                           apx(i,j,k), apx(i+1,j,k), &
                           apy(i,j,k), apy(i,j+1,k), &
                           apz(i,j,k), apz(i,j,k+1))
                      divwn = divw(n)
                   else
                      divwn = 0.d0
                   end if

                   ! we assume dx == dy == dz
                   divc(i,j,k) = -((apx(i+1,j,k)*fxp - apx(i,j,k)*fxm) * dxinv(1) &
                        +          (apy(i,j+1,k)*fyp - apy(i,j,k)*fym) * dxinv(2) &
                        +          (apz(i,j,k+1)*fzp - apz(i,j,k)*fzm) * dxinv(3) &
                        +          divwn * dxinv(1)) / vfrac(i,j,k)
                end if
             end do
          end do
       end do

       optmp = 0.d0

       !
       ! Second, we compute delta M on (lo-1,hi+1)
       !
       do       k = lo(3)-1, hi(3)+1
          do    j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1
                if (is_single_valued_cell(cellflag(i,j,k))) then
                   vtot = 0.d0
                   divnc = 0.d0
                   call get_neighbor_cells(cellflag(i,j,k),nbr)
                   do kk = -1,1
                      do jj = -1,1
                         do ii = -1,1
                            if ((ii.ne. 0 .or. jj.ne.0 .or. kk.ne. 0) .and. nbr(ii,jj,kk).eq.1) then
                               vtot = vtot + vfrac(i+ii,j+jj,k+kk)
                               divnc = divnc + vfrac(i+ii,j+jj,k+kk)*divc(i+ii,j+jj,k+kk)
                            end if
                         end do
                      enddo
                   enddo
                   divnc = divnc / vtot
                   optmp(i,j,k) = (1.d0-vfrac(i,j,k))*(divnc-divc(i,j,k))
                   delm(i,j,k,n) = -vfrac(i,j,k)*optmp(i,j,k)
                else
                   delm(i,j,k,n) = 0.d0
                end if
             end do
          end do
       end do

       !
       ! Third, redistribution
       !
       do       k = lo(3)-1, hi(3)+1
          do    j = lo(2)-1, hi(2)+1
             do i = lo(1)-1, hi(1)+1
                if (is_single_valued_cell(cellflag(i,j,k))) then                   
                   mtot = 0.d0
                   call get_neighbor_cells(cellflag(i,j,k),nbr)
                   do kk = -1,1
                      do jj = -1,1
                         do ii = -1,1
                            if ((ii.ne. 0 .or. jj.ne.0 .or. kk.ne. 0) .and. nbr(ii,jj,kk).eq.1) then
                               mtot = mtot + vfrac(i+ii,j+jj,k+kk)*q(i+ii,j+jj,k+kk,qrho)
                            end if
                         end do
                      enddo
                   enddo
                   
                   mtot = 1.d0/mtot
                   do kk = -1,1
                      do jj = -1,1
                         do ii = -1,1
                            if((ii.ne. 0 .or. jj.ne.0 .or. kk.ne. 0) .and. nbr(ii,jj,kk).eq.1) then
                               optmp(i+ii,j+jj,k+kk) = optmp(i+ii,j+jj,k+kk) &
                                    + delm(i,j,k,n)*mtot*q(i+ii,j+jj,k+kk,qrho)
                            endif
                         enddo
                      enddo
                   end do
                end if
             end do
          end do
       end do

       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ebdiffop(i,j,k,n) = divc(i,j,k) + optmp(i,j,k)
             end do
          end do
       end do

    end do

    call amrex_deallocate(optmp)

  end subroutine compute_eb_diffop

end module cns_eb_flux_module
