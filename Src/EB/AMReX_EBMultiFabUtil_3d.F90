module amrex_eb_util_module
  use amrex_fort_module, only : amrex_real
  use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell, is_single_valued_cell
  use amrex_constants_module, only : half, zero, one
  implicit none
  
  private
  public :: amrex_eb_avgdown_sv, amrex_eb_avgdown, amrex_eb_avgdown_faces, &
       amrex_eb_avgdown_boundaries, amrex_compute_eb_divergence, &
       amrex_eb_avg_fc_to_cc, amrex_eb_set_covered_nodes, &
       amrex_eb_interpolate_to_face_centroid, &
       amrex_eb_interpolate_to_face_centroid_per_cell

contains

  subroutine amrex_eb_avgdown_sv (lo, hi, fine, flo, fhi, crse, clo, chi, &
       fv, fvlo, fvhi, vfrc, vflo, vfhi, lrat, ncomp) bind(c,name='amrex_eb_avgdown_sv')
    integer, intent(in) :: lo(3), hi(3), flo(3), fhi(3), clo(3), chi(3), &
         fvlo(3), fvhi(3), vflo(3), vfhi(3), lrat(3), ncomp
    real(amrex_real), intent(in   ) :: fine( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),ncomp)
    real(amrex_real), intent(inout) :: crse( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3),ncomp)
    real(amrex_real), intent(in   ) :: fv  (fvlo(1):fvhi(1),fvlo(2):fvhi(2),fvlo(3):fvhi(3))
    real(amrex_real), intent(in   ) :: vfrc(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3))
    
    integer :: i, j, k, ii, jj, kk, n, iref, jref, kref
    real(amrex_real) :: cv
    
    do n = 1, ncomp
       do k        = lo(3), hi(3)
          kk       = k * lrat(3)
          do j     = lo(2), hi(2)
             jj    = j * lrat(2)
             do i  = lo(1), hi(1)
                ii = i * lrat(1)
                crse(i,j,k,n) = 0.d0
                cv            = 0.d0
                do       kref = 0, lrat(3)-1
                   do    jref = 0, lrat(2)-1
                      do iref = 0, lrat(1)-1
                         cv = cv + (fv(ii+iref,jj+jref,kk+kref)*vfrc(ii+iref,jj+jref,kk+kref))
                         crse(i,j,k,n) = crse(i,j,k,n) + &
                              fine(ii+iref,jj+jref,kk+kref,n)*(fv(ii+iref,jj+jref,kk+kref)*vfrc(ii+iref,jj+jref,kk+kref))
                      end do
                   end do
                end do
                if (cv .gt. 1.d-30) then
                   crse(i,j,k,n) = crse(i,j,k,n) / cv
                else
                   crse(i,j,k,n) = fine(ii,jj,kk,n)  ! covered cell
                end if
             end do
          end do
       end do
    end do
  end subroutine amrex_eb_avgdown_sv


  subroutine amrex_eb_avgdown (lo, hi, fine, flo, fhi, crse, clo, chi, &
       vfrc, vflo, vfhi, lrat, ncomp) bind(c,name='amrex_eb_avgdown')
    integer, intent(in) :: lo(3), hi(3), flo(3), fhi(3), clo(3), chi(3), &
         vflo(3), vfhi(3), lrat(3), ncomp
    real(amrex_real), intent(in   ) :: fine( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),ncomp)
    real(amrex_real), intent(inout) :: crse( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3),ncomp)
    real(amrex_real), intent(in   ) :: vfrc(vflo(1):vfhi(1),vflo(2):vfhi(2),vflo(3):vfhi(3))
    
    integer :: i, j, k, ii, jj, kk, n, iref, jref, kref
    real(amrex_real) :: cv
    
    do n = 1, ncomp
       do k        = lo(3), hi(3)
          kk       = k * lrat(3)
          do j     = lo(2), hi(2)
             jj    = j * lrat(2)
             do i  = lo(1), hi(1)
                ii = i * lrat(1)
                crse(i,j,k,n) = 0.d0
                cv            = 0.d0
                do       kref = 0, lrat(3)-1
                   do    jref = 0, lrat(2)-1
                      do iref = 0, lrat(1)-1
                         cv = cv + vfrc(ii+iref,jj+jref,kk+kref)
                         crse(i,j,k,n) = crse(i,j,k,n) + &
                              fine(ii+iref,jj+jref,kk+kref,n)*vfrc(ii+iref,jj+jref,kk+kref)
                      end do
                   end do
                end do
                if (cv .gt. 1.d-30) then
                   crse(i,j,k,n) = crse(i,j,k,n) / cv
                else
                   crse(i,j,k,n) = fine(ii,jj,kk,n)  ! covered cell
                end if
             end do
          end do
       end do
    end do
  end subroutine amrex_eb_avgdown
  
  subroutine amrex_eb_avgdown_faces (lo, hi, fine, flo, fhi, crse, clo, chi, &
       ap, aplo, aphi, lrat, idir, ncomp) bind(c,name='amrex_eb_avgdown_faces')
    integer, dimension(3), intent(in) :: lo, hi, flo, fhi, clo,chi, aplo, aphi, lrat
    integer,               intent(in) :: idir, ncomp
    real(amrex_real),   intent(in   ) :: fine( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3),ncomp) 
    real(amrex_real),   intent(inout) :: crse( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3),ncomp)
    real(amrex_real),   intent(in   ) ::   ap(aplo(1):aphi(1),aplo(2):aphi(2),aplo(3):aphi(3))

    integer  :: i, j, k, ii, jj, kk, n, iref, jref, kref
    real(amrex_real) :: fa 
 
    if(idir.eq.0) then 
      do n              = 1, ncomp
         do k           = lo(3), hi(3)
            kk          = k*lrat(3)
            do j        = lo(2), hi(2)
               jj       = j*lrat(2)
               do i     = lo(1), hi(1)
                  ii    = i*lrat(1)
                  crse(i,j,k,n) = 0.d0
                  fa            = 0.d0
                  do    kref    = 0, lrat(3)-1
                    do  jref    = 0, lrat(2)-1
                        fa            = fa + ap(ii,jj+jref,kk+kref)
                        crse(i,j,k,n) = crse(i,j,k,n) + ap(ii,jj+jref,kk+kref)*fine(ii,jj+jref,kk+kref,n)
                    enddo
                  enddo
                  if(fa.gt.1.d-30) then 
                    crse(i,j,k,n) = crse(i,j,k,n)/fa
                  else
                    crse(i,j,k,n) = fine(ii,jj,kk,n) !covered face
                  endif
               enddo
            enddo
         enddo
      enddo 
    elseif(idir.eq.1) then 
      do n             = 1, ncomp   
         do k          = lo(3), hi(3)
            kk         = k*lrat(3)
            do j       = lo(2), hi(2)
               jj      = j*lrat(2)
               do i    = lo(1), hi(1)
                  ii   = i*lrat(1)
                  crse(i,j,k,n) = 0.d0
                  fa            = 0.d0
                  do    kref    = 0, lrat(3)-1
                    do  iref    = 0, lrat(1)-1
                        fa            = fa + ap(ii+iref, jj, kk+kref)
                        crse(i,j,k,n) = crse(i,j,k,n) + ap(ii+iref,jj,kk+kref)*fine(ii+iref,jj,kk+kref,n)
                    enddo
                  enddo
                  if(fa.gt.1.d-30) then
                    crse(i,j,k,n) = crse(i,j,k,n)/fa
                  else
                    crse(i,j,k,n) = fine(ii,jj,kk,n) !covered face
                  endif
               enddo
            enddo
         enddo
      enddo
    else
      do n            = 1, ncomp
         do k         = lo(3), hi(3)
            kk        = k*lrat(3)
            do j      = lo(2), hi(2)
               jj     = j*lrat(2)
               do i   = lo(1), hi(1)
                  ii  = i*lrat(1)
                  crse(i,j,k,n) = 0.d0
                  fa            = 0.d0
                  do    jref    = 0, lrat(2)-1
                    do  iref    = 0, lrat(1)-1
                        fa            = fa + ap(ii+iref,jj+jref,kk)
                        crse(i,j,k,n) = crse(i,j,k,n) + ap(ii+iref,jj+jref,kk)*fine(ii+iref,jj+jref,kk,n)
                    enddo
                  enddo
                  if(fa.gt.1.d-30) then
                    crse(i,j,k,n) = crse(i,j,k,n)/fa
                  else
                    crse(i,j,k,n) = fine(ii,jj,kk,n) !covered face
                  endif
               enddo
            enddo
         enddo
      enddo
    endif
  end subroutine amrex_eb_avgdown_faces

  
  subroutine amrex_eb_avgdown_boundaries (lo, hi, fine, flo, fhi, crse, clo, chi, &
       ba, blo, bhi, lrat, ncomp) bind(c,name='amrex_eb_avgdown_boundaries')
    integer, dimension(3), intent(in) :: lo, hi, flo, fhi, clo, chi, blo, bhi, lrat
    integer,               intent(in) :: ncomp
    real(amrex_real),   intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3),ncomp) 
    real(amrex_real),   intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),ncomp)
    real(amrex_real),   intent(in   ) ::  ba (blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))

    integer  :: i, j, k, ii, jj, kk, n, iref, jref, kref
    real(amrex_real) :: fa 
 
    do n              = 1, ncomp
       do k           = lo(3), hi(3)
          kk          = k*lrat(3)
          do j        = lo(2), hi(2)
             jj       = j*lrat(2)
             do i     = lo(1), hi(1)
                ii    = i*lrat(1)
                crse(i,j,k,n) = 0.d0
                fa            = 0.d0
                do    kref    = 0, lrat(3)-1
                   do  jref   = 0, lrat(2)-1
                      do iref = 0, lrat(1)-1
                         fa            = fa            + ba(ii+iref,jj+jref,kk+kref)
                         crse(i,j,k,n) = crse(i,j,k,n) + ba(ii+iref,jj+jref,kk+kref) &
                              &                       *fine(ii+iref,jj+jref,kk+kref,n)
                      enddo
                   enddo
                end do
                if(fa.gt.1.d-30) then 
                   crse(i,j,k,n) = crse(i,j,k,n)/fa
                else
                   crse(i,j,k,n) = zero
                endif
             enddo
          enddo
       enddo
    end do
  end subroutine amrex_eb_avgdown_boundaries

  subroutine amrex_compute_eb_divergence (lo, hi, divu, dlo, dhi, &
       u, ulo, uhi, v, vlo, vhi, w, wlo, whi, &
       ccm, cmlo, cmhi, flag, flo, fhi, vfrc, klo, khi, &
       apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi, &
       fcx, cxlo, cxhi, fcy, cylo, cyhi, fcz, czlo, czhi, &
       dxinv) bind(c,name='amrex_compute_eb_divergence')
    implicit none
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, ulo, uhi, vlo, vhi, wlo, whi, &
         cmlo, cmhi, flo, fhi, klo, khi, axlo, axhi, aylo, ayhi, azlo, azhi, &
         cxlo, cxhi, cylo, cyhi, czlo, czhi
    real(amrex_real), intent(inout) :: divu(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    real(amrex_real), intent(in   ) ::    u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3))
    real(amrex_real), intent(in   ) ::    v(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
    real(amrex_real), intent(in   ) ::    w(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3))
    integer         , intent(in   ) ::  ccm(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3)) 
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3)) 
    real(amrex_real), intent(in   ) :: vfrc( klo(1): khi(1), klo(2): khi(2), klo(3): khi(3)) 
    real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)) 
    real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(amrex_real), intent(in   ) ::  apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2),cxlo(3):cxhi(3),2)
    real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2),cylo(3):cyhi(3),2) 
    real(amrex_real), intent(in   ) ::  fcz(czlo(1):czhi(1),czlo(2):czhi(2),czlo(3):czhi(3),2) 
    real(amrex_real), intent(in) :: dxinv(3)

    integer :: i,j,k,ii,jj,kk
    real(amrex_real) :: fxm, fxp, fym, fyp, fzm, fzp, fracx, fracy, fracz

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (is_covered_cell(flag(i,j,k))) then
                divu(i,j,k) = zero
             else if (is_regular_cell(flag(i,j,k))) then
                divu(i,j,k) = dxinv(1) * (u(i+1,j,k)-u(i,j,k)) &
                     +        dxinv(2) * (v(i,j+1,k)-v(i,j,k)) &
                     +        dxinv(3) * (w(i,j,k+1)-w(i,j,k))
             else
                fxm = u(i,j,k)
                if (apx(i,j,k).ne.zero.and.apx(i,j,k).ne.one) then 
                    jj = j + int(sign(one, fcx(i,j,k,1)))
                    kk = k + int(sign(one, fcx(i,j,k,2)))
                    fracy = abs(fcx(i,j,k,1))*real(ior(ccm(i-1,jj,k),ccm(i,jj,k)),amrex_real)
                    fracz = abs(fcx(i,j,k,2))*real(ior(ccm(i-1,j,kk),ccm(i,j,kk)),amrex_real)
                    fxm = (one-fracy)*(one-fracz)*fxm + &
                         & fracy*(one-fracz)*u(i,jj,k ) + &
                         & fracz*(one-fracy)*u(i,j ,kk) + &
                         & fracy*     fracz *u(i,jj,kk)
                endif 

                fxp = u(i+1,j,k)
                if (apx(i+1,j,k).ne.zero.and.apx(i+1,j,k).ne.one) then 
                    jj = j + int(sign(one,fcx(i+1,j,k,1)))
                    kk = k + int(sign(one,fcx(i+1,j,k,2)))
                    fracy = abs(fcx(i+1,j,k,1))*real(ior(ccm(i,jj,k),ccm(i+1,jj,k)),amrex_real)
                    fracz = abs(fcx(i+1,j,k,2))*real(ior(ccm(i,j,kk),ccm(i+1,j,kk)),amrex_real)
                    fxp = (one-fracy)*(one-fracz)*fxp + &
                         & fracy*(one-fracz)*u(i+1,jj,k ) + & 
                         & fracz*(one-fracy)*u(i+1,j ,kk) + & 
                         & fracy*     fracz *u(i+1,jj,kk)
                endif 

                fym = v(i,j,k)
                if (apy(i,j,k).ne.zero.and.apy(i,j,k).ne.one) then 
                    ii = i + int(sign(one,fcy(i,j,k,1)))
                    kk = k + int(sign(one,fcy(i,j,k,2)))
                    fracx = abs(fcy(i,j,k,1))*real(ior(ccm(ii,j-1,k),ccm(ii,j,k)),amrex_real)
                    fracz = abs(fcy(i,j,k,2))*real(ior(ccm(i,j-1,kk),ccm(i,j,kk)),amrex_real)
                    fym = (one-fracx)*(one-fracz)*fym + &
                         & fracx*(one-fracz)*v(ii,j,k ) + & 
                         & fracz*(one-fracx)*v(i ,j,kk) + &
                         & fracx*     fracz *v(ii,j,kk)
                endif 

                fyp = v(i,j+1,k)
                if (apy(i,j+1,k).ne.zero.and.apy(i,j+1,k).ne.one) then 
                    ii = i + int(sign(one,fcy(i,j+1,k,1)))
                    kk = k + int(sign(one,fcy(i,j+1,k,2)))
                    fracx = abs(fcy(i,j+1,k,1))*real(ior(ccm(ii,j,k),ccm(ii,j+1,k)),amrex_real)
                    fracz = abs(fcy(i,j+1,k,2))*real(ior(ccm(i,j,kk),ccm(i,j+1,kk)),amrex_real)
                    fyp = (one-fracx)*(one-fracz)*fyp + &
                         & fracx*(one-fracz)*v(ii,j+1,k ) + &
                         & fracz*(one-fracx)*v(i ,j+1,kk) + & 
                         & fracx*     fracz *v(ii,j+1,kk)
                endif 

                fzm = w(i,j,k)
                if (apz(i,j,k).ne.zero.and.apz(i,j,k).ne.one) then 
                    ii = i + int(sign(one,fcz(i,j,k,1)))
                    jj = j + int(sign(one,fcz(i,j,k,2)))
                    fracx = abs(fcz(i,j,k,1))*real(ior(ccm(ii,j,k-1),ccm(ii,j,k)),amrex_real)
                    fracy = abs(fcz(i,j,k,2))*real(ior(ccm(i,jj,k-1),ccm(i,jj,k)),amrex_real)
                    fzm = (one-fracx)*(one-fracy)*fzm + &
                         & fracx*(one-fracy)*w(ii,j ,k) + & 
                         & fracy*(one-fracx)*w(i ,jj,k) + &
                         & fracx*     fracy *w(ii,jj,k)
                endif 

                fzp = w(i,j,k+1)
                if (apz(i,j,k+1).ne.zero.and.apz(i,j,k+1).ne.one) then 
                    ii = i + int(sign(one,fcz(i,j,k+1,1)))
                    jj = j + int(sign(one,fcz(i,j,k+1,2)))
                    fracx = abs(fcz(i,j,k+1,1))*real(ior(ccm(ii,j,k),ccm(ii,j,k+1)),amrex_real)
                    fracy = abs(fcz(i,j,k+1,2))*real(ior(ccm(i,jj,k),ccm(i,jj,k+1)),amrex_real)
                    fzp = (one-fracx)*(one-fracy)*fzp + & 
                         & fracx*(one-fracy)*w(ii,j ,k+1) + &
                         & fracy*(one-fracx)*w(i ,jj,k+1) + &
                         & fracx*     fracy *w(ii,jj,k+1)
                endif

                divu(i,j,k) = (one/vfrc(i,j,k)) * &
                     ( dxinv(1) * (apx(i+1,j,k)*fxp-apx(i,j,k)*fxm) &
                     + dxinv(2) * (apy(i,j+1,k)*fyp-apy(i,j,k)*fym) &
                     + dxinv(3) * (apz(i,j,k+1)*fzp-apz(i,j,k)*fzm) )
             end if
          end do
       end do
    end do
  end subroutine amrex_compute_eb_divergence


  subroutine amrex_eb_avg_fc_to_cc (lo, hi, cc, cclo, cchi, fx, fxlo, fxhi, fy, fylo, fyhi, &
       fz, fzlo, fzhi, ax, axlo, axhi, ay, aylo, ayhi, az, azlo, azhi, flag, flo, fhi) &
       bind(c,name='amrex_eb_avg_fc_to_cc')
    integer, dimension(3), intent(in) :: lo, hi, cclo, cchi, fxlo, fxhi, fylo, fyhi, fzlo, fzhi, &
         axlo, axhi, aylo, ayhi, azlo, azhi, flo, fhi
    real(amrex_real), intent(inout) :: cc  (cclo(1):cchi(1),cclo(2):cchi(2),cclo(3):cchi(3),3)
    real(amrex_real), intent(in   ) :: fx  (fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3))
    real(amrex_real), intent(in   ) :: fy  (fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3))
    real(amrex_real), intent(in   ) :: fz  (fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3))
    real(amrex_real), intent(in   ) :: ax  (axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(amrex_real), intent(in   ) :: ay  (aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(amrex_real), intent(in   ) :: az  (azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3))
    
    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (is_covered_cell(flag(i,j,k))) then
                cc(i,j,k,1) = zero
                cc(i,j,k,2) = zero
                cc(i,j,k,3) = zero
             else
                if (ax(i,j,k) .eq. zero) then
                   cc(i,j,k,1) = fx(i+1,j,k)
                else if (ax(i+1,j,k) .eq. zero) then
                   cc(i,j,k,1) = fx(i,j,k)
                else
                   cc(i,j,k,1) = half * (fx(i,j,k) + fx(i+1,j,k))
                end if
                
                if (ay(i,j,k) .eq. zero) then
                   cc(i,j,k,2) = fy(i,j+1,k)
                else if (ay(i,j+1,k) .eq. zero) then
                   cc(i,j,k,2) = fy(i,j,k)
                else
                   cc(i,j,k,2) = half * (fy(i,j,k) + fy(i,j+1,k))
                end if

                if (az(i,j,k) .eq. zero) then
                   cc(i,j,k,3) = fz(i,j,k+1)
                else if (az(i,j,k+1) .eq. zero) then
                   cc(i,j,k,3) = fz(i,j,k)
                else
                   cc(i,j,k,3) = half * (fz(i,j,k) + fz(i,j,k+1))
                end if
             end if
          end do
       end do
    end do
  end subroutine amrex_eb_avg_fc_to_cc

  subroutine amrex_eb_set_covered_nodes (lo, hi, d, dlo, dhi, f, flo, fhi, v, nc) &
       bind(c,name='amrex_eb_set_covered_nodes')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), flo(3), fhi(3), nc
    real(amrex_real), intent(inout) :: d(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nc)
    integer, intent(in) :: f(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(amrex_real), intent(in) :: v(nc)

    integer :: i, j, k, n

    do n = 1, nc
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (is_covered_cell(f(i-1,j-1,k-1)) .and. &
                    is_covered_cell(f(i  ,j-1,k-1)) .and. &
                    is_covered_cell(f(i-1,j  ,k-1)) .and. &
                    is_covered_cell(f(i  ,j  ,k-1)) .and. &
                    is_covered_cell(f(i-1,j-1,k  )) .and. &
                    is_covered_cell(f(i  ,j-1,k  )) .and. &
                    is_covered_cell(f(i-1,j  ,k  )) .and. &
                    is_covered_cell(f(i  ,j  ,k  ))) then
                   d(i,j,k,n) = v(n)
                end if
             end do
          end do
       end do
    end do

  end subroutine amrex_eb_set_covered_nodes

  ! Interpolate face-based variable from face center to face centroid -- this version 
  !   does one face on all grids
  subroutine amrex_eb_interpolate_to_face_centroid ( lo, hi, ivar, var, vlo, vhi, ncomp, &
        areafrac, alo, ahi, cent, clo, chi, flags, flo, fhi, face_type  ) &
       bind(c,name='amrex_eb_interpolate_to_face_centroid')

      use amrex_ebcellflag_module, only: is_covered_cell, get_neighbor_cells

      ! Tile bounds ( face centered )
      integer,  intent(in   ) :: lo(3),  hi(3)

      ! Array Bounds
      integer,  intent(in   ) :: vlo(3), vhi(3)
      integer,  intent(in   ) :: alo(3), ahi(3)
      integer,  intent(in   ) :: clo(3), chi(3)
      integer,  intent(in   ) :: flo(3), fhi(3)
      integer,  intent(in   ) :: ncomp

      ! Type of face (1=x, 2=y, 3=z)
      integer,  intent(in   ) :: face_type

      ! Arrays
      real(amrex_real), intent(inout) ::                            &
           & ivar(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),ncomp)  ! Interpolated Variable

      real(amrex_real), intent(inout) ::                            &
           & var(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),ncomp), &
           & areafrac(alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3)),  &
           & cent(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),2)

      integer, intent(in   ) :: &
           & flags(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))


      ! Local variables
      integer          :: i, j, k, n, nbr(-1:1,-1:1,-1:1)
      real(amrex_real) :: fracx, fracy, fracz

      select case ( face_type )
      case(1) ! >>>>>>>>>>>>>>>>>>>>>>  X-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !
         do n = 1, ncomp
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if ( ( areafrac(i,j,k) > zero ) .and. ( areafrac(i,j,k) < one ) ) then

                        call get_neighbor_cells( flags(i,j,k), nbr )

                        if ( cent(i,j,k,1) < zero ) then
                           fracy = - cent(i,j,k,1) * nbr(0,-1,0)
                           if ( cent(i,j,k,2) <= zero ) then
                              fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracy * var(i,j-1,k  ,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k  ,n)) + &
                                   &                fracz * (     fracy * var(i,j-1,k-1,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k-1,n))
                           else
                              fracz =  cent(i,j,k,2) * nbr(0,0,1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracy * var(i,j-1,k  ,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k  ,n)) + &
                                   &                fracz * (     fracy * var(i,j-1,k+1,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k+1,n))
                           endif
                        else
                           fracy = cent(i,j,k,1) * nbr(0,1,0)
                           if ( cent(i,j,k,2) <= zero ) then
                              fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracy * var(i,j+1,k  ,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k  ,n)) + &
                                   &                fracz * (     fracy * var(i,j+1,k-1,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k-1,n))
                           else
                              fracz =  cent(i,j,k,2) * nbr(0,0,1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracy * var(i,j+1,k  ,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k  ,n)) + &
                                   &                fracz * (     fracy * var(i,j+1,k+1,n)  + &
                                   &                        (one-fracy) * var(i,j  ,k+1,n))
                           endif
                        end if
                     else
                        ivar(i,j,k,n) = var(i,j,k,n)
                     end if
                  end do
               end do
            end do
         end do

      case(2)  ! >>>>>>>>>>>>>>>>>>>>>>  Y-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

         do n = 1, ncomp
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if ( ( areafrac(i,j,k) > zero ) .and. ( areafrac(i,j,k) < one ) ) then

                        call get_neighbor_cells( flags(i,j,k), nbr )

                        if ( cent(i,j,k,1) < zero ) then
                           fracx = - cent(i,j,k,1) * nbr(-1,0,0)
                           if ( cent(i,j,k,2) <= zero ) then
                              fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracx * var(i-1,j,k  ,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k  ,n)) + &
                                   &                fracz * (     fracx * var(i-1,j,k-1,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k-1,n))
                           else
                              fracz =  cent(i,j,k,2) * nbr(0,0,1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracx * var(i-1,j,k  ,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k  ,n)) + &
                                   &                fracz * (     fracx * var(i-1,j,k+1,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k+1,n))
                           endif
                        else
                           fracx = cent(i,j,k,1) * nbr(1,0,0)
                           if ( cent(i,j,k,2) <= zero ) then
                              fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracx * var(i+1,j,k  ,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k  ,n)) + &
                                   &                fracz * (     fracx * var(i+1,j,k-1,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k-1,n))
                           else
                              fracz =  cent(i,j,k,2) * nbr(0,0,1)
                              ivar(i,j,k,n) = (one-fracz) * (     fracx * var(i+1,j,k  ,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k  ,n)) + &
                                   &                fracz * (     fracx * var(i+1,j,k+1,n)  + &
                                   &                        (one-fracx) * var(i  ,j,k+1,n))
                           endif
                        end if
                     else
                        ivar(i,j,k,n) = var(i,j,k,n)
                     end if
                  end do
               end do
            end do
         end do

      case(3) ! >>>>>>>>>>>>>>>>>>>>>>  Z-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

         do n = 1, ncomp
            do k = lo(3), hi(3)
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     if ( ( areafrac(i,j,k) > zero ) .and. ( areafrac(i,j,k) < one ) ) then

                        call get_neighbor_cells( flags(i,j,k), nbr )

                        if ( cent(i,j,k,1) < zero ) then
                           fracx = - cent(i,j,k,1) * nbr(-1,0,0)
                           if ( cent(i,j,k,2) <= zero ) then
                              fracy = - cent(i,j,k,2) * nbr(0,-1,0)
                              ivar(i,j,k,n) = (one-fracy) * (     fracx * var(i-1,j  ,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j  ,k,n)) + &
                                   &                fracy * (     fracx * var(i-1,j-1,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j-1,k,n))
                           else
                              fracy =  cent(i,j,k,2) * nbr(0,1,0)
                              ivar(i,j,k,n) = (one-fracy) * (     fracx * var(i-1,j  ,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j  ,k,n)) + &
                                   &                fracy * (     fracx * var(i-1,j+1,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j+1,k,n))
                           endif
                        else
                           fracx = cent(i,j,k,1) * nbr(1,0,0)
                           if ( cent(i,j,k,2) <= zero ) then
                              fracy = - cent(i,j,k,2) * nbr(0,-1,0)
                              ivar(i,j,k,n) = (one-fracy) * (     fracx * var(i+1,j  ,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j  ,k,n)) + &
                                   &                fracy * (     fracx * var(i+1,j-1,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j-1,k,n))
                           else
                              fracy =  cent(i,j,k,2) * nbr(0,1,0)
                              ivar(i,j,k,n) = (one-fracy) * (     fracx * var(i+1,j  ,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j  ,k,n)) + &
                                   &                fracy * (     fracx * var(i+1,j+1,k,n)  + &
                                   &                        (one-fracx) * var(i  ,j+1,k,n))
                           endif
                        end if
                     else
                        ivar(i,j,k,n) = var(i,j,k,n)
                     end if
                  end do
               end do
            end do
         end do

      case default

         write(*,*) "amrex_eb_interpolate_to_face_centroid(): face_type = ", face_type, " but valid values are 1,2,3"
         stop

      end select

  end subroutine amrex_eb_interpolate_to_face_centroid

   !
   ! Returns flux at face centroid in direction dir for just cell (i,j,k) -- 
   !         note nbr is passed in 
   !
   function amrex_eb_interpolate_to_face_centroid_per_cell ( i, j, k, dir, var, vlo,  n,  &
        afrac, alo, cent, clo, nbr )  result(ivar)

      use amrex_ebcellflag_module, only: is_covered_cell
      use amrex_error_module,      only: amrex_abort

      ! Face indices: these must be consistent with a staggered indexing
      ! and therefore consistent with the value of dir
      integer,  intent(in   ) :: i, j, k

      ! Direction of staggering (1=x, 2=y, 3=z): this specify how (i,j,k) must
      ! be interpreted, i.e. which staggered numbering the indexing refer to
      integer,  intent(in   ) :: dir

      ! The component to interpolate
      integer,  intent(in   ) :: n

      ! Array Bounds ( only start index )
      integer,  intent(in   ) :: vlo(3), alo(3), clo(3)

      ! Arrays
      real(amrex_real), intent(in   ) ::           &
           &   var(vlo(1):, vlo(2):, vlo(3):,1:), &
           & afrac(alo(1):, alo(2):, alo(3):),    &
           &  cent(clo(1):, clo(2):, clo(3):,1:)

      ! Neighbors information
      integer,  intent(in   ) :: nbr(-1:1,-1:1,-1:1)

      ! Output: the interpolated value
      real(amrex_real)               :: ivar

      ! Local variables
      real(amrex_real)               :: fracx, fracy, fracz

      select case ( dir )
      case(1) ! >>>>>>>>>>>>>>>>>>>>>>  X-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

         if ( afrac(i,j,k) == zero ) then
            ivar = zero
         else if ( afrac(i,j,k) == one ) then
            ivar = var(i,j,k,n)
         else
            if ( cent(i,j,k,1) < zero ) then
               fracy = - cent(i,j,k,1) * nbr(0,-1,0)
               if ( cent(i,j,k,2) <= zero ) then
                  fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                  ivar = (one-fracz) * (     fracy * var(i,j-1,k  ,n)  + &
                       &               (one-fracy) * var(i,j  ,k  ,n)) + &
                       &       fracz * (     fracy * var(i,j-1,k-1,n)  + &
                       &               (one-fracy) * var(i,j  ,k-1,n))
               else
                  fracz =  cent(i,j,k,2) * nbr(0,0,1)
                  ivar = (one-fracz) * (     fracy * var(i,j-1,k  ,n)  + &
                       &               (one-fracy) * var(i,j  ,k  ,n)) + &
                       &       fracz * (     fracy * var(i,j-1,k+1,n)  + &
                       &               (one-fracy) * var(i,j  ,k+1,n))
               endif
            else
               fracy = cent(i,j,k,1) * nbr(0,1,0)
               if ( cent(i,j,k,2) <= zero ) then
                  fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                  ivar = (one-fracz) * (     fracy * var(i,j+1,k  ,n)  + &
                       &               (one-fracy) * var(i,j  ,k  ,n)) + &
                       &       fracz * (     fracy * var(i,j+1,k-1,n)  + &
                       &               (one-fracy) * var(i,j  ,k-1,n))
               else
                  fracz =  cent(i,j,k,2) * nbr(0,0,1)
                  ivar= (one-fracz) * (     fracy * var(i,j+1,k  ,n)  + &
                       &              (one-fracy) * var(i,j  ,k  ,n)) + &
                       &      fracz * (     fracy * var(i,j+1,k+1,n)  + &
                       &              (one-fracy) * var(i,j  ,k+1,n))
               endif
            end if
         end if


      case(2)  ! >>>>>>>>>>>>>>>>>>>>>>  Y-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

         if ( afrac(i,j,k) == zero ) then
            ivar = zero
         else if ( afrac(i,j,k) == one ) then
            ivar = var(i,j,k,n)
         else
            if ( cent(i,j,k,1) < zero ) then
               fracx = - cent(i,j,k,1) * nbr(-1,0,0)
               if ( cent(i,j,k,2) <= zero ) then
                  fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                  ivar = (one-fracz) * (     fracx * var(i-1,j,k  ,n)  + &
                       &               (one-fracx) * var(i  ,j,k  ,n)) + &
                       &       fracz * (     fracx * var(i-1,j,k-1,n)  + &
                       &               (one-fracx) * var(i  ,j,k-1,n))
               else
                  fracz =  cent(i,j,k,2) * nbr(0,0,1)
                  ivar = (one-fracz) * (     fracx * var(i-1,j,k  ,n)  + &
                       &               (one-fracx) * var(i  ,j,k  ,n)) + &
                       &       fracz * (     fracx * var(i-1,j,k+1,n)  + &
                       &               (one-fracx) * var(i  ,j,k+1,n))
               endif
            else
               fracx = cent(i,j,k,1) * nbr(1,0,0)
               if ( cent(i,j,k,2) <= zero ) then
                  fracz = - cent(i,j,k,2) * nbr(0,0,-1)
                  ivar = (one-fracz) * (     fracx * var(i+1,j,k  ,n)  + &
                       &               (one-fracx) * var(i  ,j,k  ,n)) + &
                       &       fracz * (     fracx * var(i+1,j,k-1,n)  + &
                       &               (one-fracx) * var(i  ,j,k-1,n))
               else
                  fracz =  cent(i,j,k,2) * nbr(0,0,1)
                  ivar = (one-fracz) * (     fracx * var(i+1,j,k  ,n)  + &
                       &               (one-fracx) * var(i  ,j,k  ,n)) + &
                       &       fracz * (     fracx * var(i+1,j,k+1,n)  + &
                       &               (one-fracx) * var(i  ,j,k+1,n))
               endif
            end if
         end if


      case(3) ! >>>>>>>>>>>>>>>>>>>>>>  Z-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

         if ( afrac(i,j,k) == zero ) then
            ivar = zero
         else if ( afrac(i,j,k) == one ) then
            ivar = var(i,j,k,n)
         else
            if ( cent(i,j,k,1) < zero ) then
               fracx = - cent(i,j,k,1) * nbr(-1,0,0)
               if ( cent(i,j,k,2) <= zero ) then
                  fracy = - cent(i,j,k,2) * nbr(0,-1,0)
                  ivar = (one-fracy) * (     fracx * var(i-1,j  ,k,n)  + &
                       &               (one-fracx) * var(i  ,j  ,k,n)) + &
                       &       fracy * (     fracx * var(i-1,j-1,k,n)  + &
                       &               (one-fracx) * var(i  ,j-1,k,n))
               else
                  fracy =  cent(i,j,k,2) * nbr(0,1,0)
                  ivar = (one-fracy) * (     fracx * var(i-1,j  ,k,n)  + &
                       &               (one-fracx) * var(i  ,j  ,k,n)) + &
                       &       fracy * (     fracx * var(i-1,j+1,k,n)  + &
                       &               (one-fracx) * var(i  ,j+1,k,n))
               endif
            else
               fracx = cent(i,j,k,1) * nbr(1,0,0)
               if ( cent(i,j,k,2) <= zero ) then
                  fracy = - cent(i,j,k,2) * nbr(0,-1,0)
                  ivar = (one-fracy) * (     fracx * var(i+1,j  ,k,n)  + &
                       &               (one-fracx) * var(i  ,j  ,k,n)) + &
                       &       fracy * (     fracx * var(i+1,j-1,k,n)  + &
                       &               (one-fracx) * var(i  ,j-1,k,n))
               else
                  fracy =  cent(i,j,k,2) * nbr(0,1,0)
                  ivar = (one-fracy) * (     fracx * var(i+1,j  ,k,n)  + &
                       &               (one-fracx) * var(i  ,j  ,k,n)) + &
                       &       fracy * (     fracx * var(i+1,j+1,k,n)  + &
                       &               (one-fracx) * var(i  ,j+1,k,n))
               endif
            end if
         end if

      case default

         call amrex_abort( "interpolate_to_face_centroid(): value of 'dir'"&
              //" is invalid. Must be either 1,2, or 3")

      end select

   end function amrex_eb_interpolate_to_face_centroid_per_cell

end module amrex_eb_util_module
