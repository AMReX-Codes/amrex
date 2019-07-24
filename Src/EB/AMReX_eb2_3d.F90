
module amrex_eb2_3d_module

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module, only : zero, one, two, three, four, five, six, seven, eight,&
       nine, ten, eleven, fifteen, sixteen, half, third, fourth, sixth, eighth, twelfth
  use amrex_ebcellflag_module, only : regular_cell => regular, covered_cell => covered, &
       is_regular_cell, is_single_valued_cell, is_covered_cell, get_cell_type, get_neighbor_cells, &
       set_regular_cell, set_single_valued_cell, set_covered_cell, &
       set_neighbor, clear_neighbor, clear_allneighbors, is_neighbor
  implicit none

  integer, parameter :: regular = 0
  integer, parameter :: covered = 1
  integer, parameter :: irregular = 2
  integer, parameter :: unknown = 3

  real(amrex_real), private, parameter :: small = 1.d-14

  private
  public :: amrex_eb2_coarsen_from_fine, amrex_eb2_build_cellflag_from_ap, amrex_eb2_check_mvmc

contains

  subroutine amrex_eb2_coarsen_from_fine (lo, hi, xlo, xhi, ylo, yhi, zlo, zhi, &
       cvol, cvlo, cvhi, fvol, fvlo, fvhi, ccent, cclo, cchi, fcent, fclo, fchi, &
       cba, cbalo, cbahi, fba, fbalo, fbahi, cbc, cbclo, cbchi, fbc, fbclo, fbchi, &
       cbn, cbnlo, cbnhi, fbn, fbnlo, fbnhi, capx, caxlo, caxhi, fapx, faxlo, faxhi, &
       capy, caylo, cayhi, fapy, faylo, fayhi, capz, cazlo, cazhi, fapz, fazlo, fazhi, &
       cfcx, cfxlo, cfxhi, ffcx, ffxlo, ffxhi, cfcy, cfylo, cfyhi, ffcy, ffylo, ffyhi, &
       cfcz, cfzlo, cfzhi, ffcz, ffzlo, ffzhi, cflag, cflo, cfhi, fflag, fflo, ffhi, ierr) &
       bind(c, name='amrex_eb2_coarsen_from_fine')
    integer, dimension(3), intent(in) :: lo, hi, xlo, xhi, ylo, yhi, zlo, zhi, &
         cvlo, cvhi,  fvlo, fvhi, cclo, cchi, fclo, fchi, &
         cbalo, cbahi, fbalo, fbahi, cbclo, cbchi, fbclo, fbchi, &
         cbnlo, cbnhi, fbnlo, fbnhi, caxlo, caxhi, faxlo, faxhi, &
         caylo, cayhi, faylo, fayhi, cazlo, cazhi, fazlo, fazhi, &
         cfxlo, cfxhi, ffxlo, ffxhi, cfylo, cfyhi, ffylo, ffyhi, &
         cfzlo, cfzhi, ffzlo, ffzhi, cflo, cfhi, fflo, ffhi
    integer, intent(inout) :: ierr
    real(amrex_real), intent(inout) :: cvol ( cvlo(1): cvhi(1), cvlo(2): cvhi(2), cvlo(3): cvhi(3))
    real(amrex_real), intent(in   ) :: fvol ( fvlo(1): fvhi(1), fvlo(2): fvhi(2), fvlo(3): fvhi(3))
    real(amrex_real), intent(inout) :: ccent( cclo(1): cchi(1), cclo(2): cchi(2), cclo(3): cchi(3),3)
    real(amrex_real), intent(in   ) :: fcent( fclo(1): fchi(1), fclo(2): fchi(2), fclo(3): fchi(3),3)
    real(amrex_real), intent(inout) :: cba  (cbalo(1):cbahi(1),cbalo(2):cbahi(2),cbalo(3):cbahi(3))
    real(amrex_real), intent(in   ) :: fba  (fbalo(1):fbahi(1),fbalo(2):fbahi(2),fbalo(3):fbahi(3))
    real(amrex_real), intent(inout) :: cbc  (cbclo(1):cbchi(1),cbclo(2):cbchi(2),cbclo(3):cbchi(3),3)
    real(amrex_real), intent(in   ) :: fbc  (fbclo(1):fbchi(1),fbclo(2):fbchi(2),fbclo(3):fbchi(3),3)
    real(amrex_real), intent(inout) :: cbn  (cbnlo(1):cbnhi(1),cbnlo(2):cbnhi(2),cbnlo(3):cbnhi(3),3)
    real(amrex_real), intent(in   ) :: fbn  (fbnlo(1):fbnhi(1),fbnlo(2):fbnhi(2),fbnlo(3):fbnhi(3),3)
    real(amrex_real), intent(inout) :: capx (caxlo(1):caxhi(1),caxlo(2):caxhi(2),caxlo(3):caxhi(3))
    real(amrex_real), intent(in   ) :: fapx (faxlo(1):faxhi(1),faxlo(2):faxhi(2),faxlo(3):faxhi(3))
    real(amrex_real), intent(inout) :: capy (caylo(1):cayhi(1),caylo(2):cayhi(2),caylo(3):cayhi(3))
    real(amrex_real), intent(in   ) :: fapy (faylo(1):fayhi(1),faylo(2):fayhi(2),faylo(3):fayhi(3))
    real(amrex_real), intent(inout) :: capz (cazlo(1):cazhi(1),cazlo(2):cazhi(2),cazlo(3):cazhi(3))
    real(amrex_real), intent(in   ) :: fapz (fazlo(1):fazhi(1),fazlo(2):fazhi(2),fazlo(3):fazhi(3))
    real(amrex_real), intent(inout) :: cfcx (cfxlo(1):cfxhi(1),cfxlo(2):cfxhi(2),cfxlo(3):cfxhi(3),2)
    real(amrex_real), intent(in   ) :: ffcx (ffxlo(1):ffxhi(1),ffxlo(2):ffxhi(2),ffxlo(3):ffxhi(3),2)
    real(amrex_real), intent(inout) :: cfcy (cfylo(1):cfyhi(1),cfylo(2):cfyhi(2),cfylo(3):cfyhi(3),2)
    real(amrex_real), intent(in   ) :: ffcy (ffylo(1):ffyhi(1),ffylo(2):ffyhi(2),ffylo(3):ffyhi(3),2)
    real(amrex_real), intent(inout) :: cfcz (cfzlo(1):cfzhi(1),cfzlo(2):cfzhi(2),cfzlo(3):cfzhi(3),2)
    real(amrex_real), intent(in   ) :: ffcz (ffzlo(1):ffzhi(1),ffzlo(2):ffzhi(2),ffzlo(3):ffzhi(3),2)
    integer         , intent(inout) :: cflag( cflo(1): cfhi(1), cflo(2): cfhi(2), cflo(3): cfhi(3))
    integer         , intent(in   ) :: fflag( fflo(1): ffhi(1), fflo(2): ffhi(2), fflo(3): ffhi(3))

    integer :: i,j,k, ii,jj,kk, ftype(2,2,2)
    real(amrex_real) :: cvolinv, cbainv, nx, ny, nz, nfac, apinv

    do       k = lo(3), hi(3)
       kk = k*2
       do    j = lo(2), hi(2)
          jj = j*2
          do i = lo(1), hi(1)
             ii = i*2

             ftype = get_cell_type(fflag(ii:ii+1,jj:jj+1,kk:kk+1))
             if (all(ftype.eq.regular_cell)) then
                ! nothing to do
             else if (all(ftype.eq.covered_cell)) then
                call set_covered_cell(cflag(i,j,k))
                cvol(i,j,k) = zero
             else

                call set_single_valued_cell(cflag(i,j,k))

                cvol(i,j,k) = eighth*sum(fvol(ii:ii+1,jj:jj+1,kk:kk+1))
                cvolinv = one/cvol(i,j,k)
                
                ccent(i,j,k,1) = eighth*cvolinv* &
                     ( fvol(ii  ,jj  ,kk  )*(half*fcent(ii  ,jj  ,kk  ,1)-fourth) &
                     + fvol(ii+1,jj  ,kk  )*(half*fcent(ii+1,jj  ,kk  ,1)+fourth) &
                     + fvol(ii  ,jj+1,kk  )*(half*fcent(ii  ,jj+1,kk  ,1)-fourth) &
                     + fvol(ii+1,jj+1,kk  )*(half*fcent(ii+1,jj+1,kk  ,1)+fourth) &
                     + fvol(ii  ,jj  ,kk+1)*(half*fcent(ii  ,jj  ,kk+1,1)-fourth) &
                     + fvol(ii+1,jj  ,kk+1)*(half*fcent(ii+1,jj  ,kk+1,1)+fourth) &
                     + fvol(ii  ,jj+1,kk+1)*(half*fcent(ii  ,jj+1,kk+1,1)-fourth) &
                     + fvol(ii+1,jj+1,kk+1)*(half*fcent(ii+1,jj+1,kk+1,1)+fourth) )
                ccent(i,j,k,2) = eighth*cvolinv* &
                     ( fvol(ii  ,jj  ,kk  )*(half*fcent(ii  ,jj  ,kk  ,2)-fourth) &
                     + fvol(ii+1,jj  ,kk  )*(half*fcent(ii+1,jj  ,kk  ,2)-fourth) &
                     + fvol(ii  ,jj+1,kk  )*(half*fcent(ii  ,jj+1,kk  ,2)+fourth) &
                     + fvol(ii+1,jj+1,kk  )*(half*fcent(ii+1,jj+1,kk  ,2)+fourth) &
                     + fvol(ii  ,jj  ,kk+1)*(half*fcent(ii  ,jj  ,kk+1,2)-fourth) &
                     + fvol(ii+1,jj  ,kk+1)*(half*fcent(ii+1,jj  ,kk+1,2)-fourth) &
                     + fvol(ii  ,jj+1,kk+1)*(half*fcent(ii  ,jj+1,kk+1,2)+fourth) &
                     + fvol(ii+1,jj+1,kk+1)*(half*fcent(ii+1,jj+1,kk+1,2)+fourth) )
                ccent(i,j,k,3) = eighth*cvolinv* &
                     ( fvol(ii  ,jj  ,kk  )*(half*fcent(ii  ,jj  ,kk  ,3)-fourth) &
                     + fvol(ii+1,jj  ,kk  )*(half*fcent(ii+1,jj  ,kk  ,3)-fourth) &
                     + fvol(ii  ,jj+1,kk  )*(half*fcent(ii  ,jj+1,kk  ,3)-fourth) &
                     + fvol(ii+1,jj+1,kk  )*(half*fcent(ii+1,jj+1,kk  ,3)-fourth) &
                     + fvol(ii  ,jj  ,kk+1)*(half*fcent(ii  ,jj  ,kk+1,3)+fourth) &
                     + fvol(ii+1,jj  ,kk+1)*(half*fcent(ii+1,jj  ,kk+1,3)+fourth) &
                     + fvol(ii  ,jj+1,kk+1)*(half*fcent(ii  ,jj+1,kk+1,3)+fourth) &
                     + fvol(ii+1,jj+1,kk+1)*(half*fcent(ii+1,jj+1,kk+1,3)+fourth) )
                
                cba(i,j,k) = fourth*sum(fba(ii:ii+1,jj:jj+1,kk:kk+1))
                cbainv = one/cba(i,j,k)
                
                cbc(i,j,k,1) = fourth*cbainv* &
                     ( fba(ii  ,jj  ,kk  )*(half*fbc(ii  ,jj  ,kk  ,1)-fourth) &
                     + fba(ii+1,jj  ,kk  )*(half*fbc(ii+1,jj  ,kk  ,1)+fourth) &
                     + fba(ii  ,jj+1,kk  )*(half*fbc(ii  ,jj+1,kk  ,1)-fourth) &
                     + fba(ii+1,jj+1,kk  )*(half*fbc(ii+1,jj+1,kk  ,1)+fourth) &
                     + fba(ii  ,jj  ,kk+1)*(half*fbc(ii  ,jj  ,kk+1,1)-fourth) &
                     + fba(ii+1,jj  ,kk+1)*(half*fbc(ii+1,jj  ,kk+1,1)+fourth) &
                     + fba(ii  ,jj+1,kk+1)*(half*fbc(ii  ,jj+1,kk+1,1)-fourth) &
                     + fba(ii+1,jj+1,kk+1)*(half*fbc(ii+1,jj+1,kk+1,1)+fourth) )
                cbc(i,j,k,2) = fourth*cbainv* &
                     ( fba(ii  ,jj  ,kk  )*(half*fbc(ii  ,jj  ,kk  ,2)-fourth) &
                     + fba(ii+1,jj  ,kk  )*(half*fbc(ii+1,jj  ,kk  ,2)-fourth) &
                     + fba(ii  ,jj+1,kk  )*(half*fbc(ii  ,jj+1,kk  ,2)+fourth) &
                     + fba(ii+1,jj+1,kk  )*(half*fbc(ii+1,jj+1,kk  ,2)+fourth) &
                     + fba(ii  ,jj  ,kk+1)*(half*fbc(ii  ,jj  ,kk+1,2)-fourth) &
                     + fba(ii+1,jj  ,kk+1)*(half*fbc(ii+1,jj  ,kk+1,2)-fourth) &
                     + fba(ii  ,jj+1,kk+1)*(half*fbc(ii  ,jj+1,kk+1,2)+fourth) &
                     + fba(ii+1,jj+1,kk+1)*(half*fbc(ii+1,jj+1,kk+1,2)+fourth) )
                cbc(i,j,k,3) = fourth*cbainv* &
                     ( fba(ii  ,jj  ,kk  )*(half*fbc(ii  ,jj  ,kk  ,3)-fourth) &
                     + fba(ii+1,jj  ,kk  )*(half*fbc(ii+1,jj  ,kk  ,3)-fourth) &
                     + fba(ii  ,jj+1,kk  )*(half*fbc(ii  ,jj+1,kk  ,3)-fourth) &
                     + fba(ii+1,jj+1,kk  )*(half*fbc(ii+1,jj+1,kk  ,3)-fourth) &
                     + fba(ii  ,jj  ,kk+1)*(half*fbc(ii  ,jj  ,kk+1,3)+fourth) &
                     + fba(ii+1,jj  ,kk+1)*(half*fbc(ii+1,jj  ,kk+1,3)+fourth) &
                     + fba(ii  ,jj+1,kk+1)*(half*fbc(ii  ,jj+1,kk+1,3)+fourth) &
                     + fba(ii+1,jj+1,kk+1)*(half*fbc(ii+1,jj+1,kk+1,3)+fourth) )
                
                nx =   fbn(ii  ,jj  ,kk  ,1)*fba(ii  ,jj  ,kk  ) &
                     + fbn(ii+1,jj  ,kk  ,1)*fba(ii+1,jj  ,kk  ) &
                     + fbn(ii  ,jj+1,kk  ,1)*fba(ii  ,jj+1,kk  ) &
                     + fbn(ii+1,jj+1,kk  ,1)*fba(ii+1,jj+1,kk  ) &
                     + fbn(ii  ,jj  ,kk+1,1)*fba(ii  ,jj  ,kk+1) &
                     + fbn(ii+1,jj  ,kk+1,1)*fba(ii+1,jj  ,kk+1) &
                     + fbn(ii  ,jj+1,kk+1,1)*fba(ii  ,jj+1,kk+1) &
                     + fbn(ii+1,jj+1,kk+1,1)*fba(ii+1,jj+1,kk+1)
                ny =   fbn(ii  ,jj  ,kk  ,2)*fba(ii  ,jj  ,kk  ) &
                     + fbn(ii+1,jj  ,kk  ,2)*fba(ii+1,jj  ,kk  ) &
                     + fbn(ii  ,jj+1,kk  ,2)*fba(ii  ,jj+1,kk  ) &
                     + fbn(ii+1,jj+1,kk  ,2)*fba(ii+1,jj+1,kk  ) &
                     + fbn(ii  ,jj  ,kk+1,2)*fba(ii  ,jj  ,kk+1) &
                     + fbn(ii+1,jj  ,kk+1,2)*fba(ii+1,jj  ,kk+1) &
                     + fbn(ii  ,jj+1,kk+1,2)*fba(ii  ,jj+1,kk+1) &
                     + fbn(ii+1,jj+1,kk+1,2)*fba(ii+1,jj+1,kk+1)
                nz =   fbn(ii  ,jj  ,kk  ,3)*fba(ii  ,jj  ,kk  ) &
                     + fbn(ii+1,jj  ,kk  ,3)*fba(ii+1,jj  ,kk  ) &
                     + fbn(ii  ,jj+1,kk  ,3)*fba(ii  ,jj+1,kk  ) &
                     + fbn(ii+1,jj+1,kk  ,3)*fba(ii+1,jj+1,kk  ) &
                     + fbn(ii  ,jj  ,kk+1,3)*fba(ii  ,jj  ,kk+1) &
                     + fbn(ii+1,jj  ,kk+1,3)*fba(ii+1,jj  ,kk+1) &
                     + fbn(ii  ,jj+1,kk+1,3)*fba(ii  ,jj+1,kk+1) &
                     + fbn(ii+1,jj+1,kk+1,3)*fba(ii+1,jj+1,kk+1)
                if (nx.eq.zero .and. ny.eq.zero .and. nz.eq.zero) then
                   ierr = 1
                   return
                end if
                nfac = one/sqrt(nx*nx+ny*ny+nz*nz)
                cbn(i,j,k,1) = nx*nfac
                cbn(i,j,k,2) = ny*nfac
                cbn(i,j,k,3) = nz*nfac
             end if
          end do
       end do
    end do

    do       k = xlo(3), xhi(3)
       kk = k*2
       do    j = xlo(2), xhi(2)
          jj = j*2
          do i = xlo(1), xhi(1)
             ii = i*2

             capx(i,j,k) = fourth*(fapx(ii,jj,kk)+fapx(ii,jj+1,kk)+fapx(ii,jj,kk+1)+fapx(ii,jj+1,kk+1))
             if (capx(i,j,k) .ne. zero) then
                apinv = one/capx(i,j,k)
                cfcx(i,j,k,1) = fourth*apinv* &
                     ( fapx(ii,jj  ,kk  )*(half*ffcx(ii,jj  ,kk  ,1)-fourth) &
                     + fapx(ii,jj+1,kk  )*(half*ffcx(ii,jj+1,kk  ,1)+fourth) &
                     + fapx(ii,jj  ,kk+1)*(half*ffcx(ii,jj  ,kk+1,1)-fourth) &
                     + fapx(ii,jj+1,kk+1)*(half*ffcx(ii,jj+1,kk+1,1)+fourth) )
                cfcx(i,j,k,2) = fourth*apinv* &
                     ( fapx(ii,jj  ,kk  )*(half*ffcx(ii,jj  ,kk  ,2)-fourth) &
                     + fapx(ii,jj+1,kk  )*(half*ffcx(ii,jj+1,kk  ,2)-fourth) &
                     + fapx(ii,jj  ,kk+1)*(half*ffcx(ii,jj  ,kk+1,2)+fourth) &
                     + fapx(ii,jj+1,kk+1)*(half*ffcx(ii,jj+1,kk+1,2)+fourth) )                
             end if
          end do
       end do
    end do

    do       k = ylo(3), yhi(3)
       kk = k*2
       do    j = ylo(2), yhi(2)
          jj = j*2
          do i = ylo(1), yhi(1)
             ii = i*2

             capy(i,j,k) = fourth*(fapy(ii,jj,kk)+fapy(ii+1,jj,kk)+fapy(ii,jj,kk+1)+fapy(ii+1,jj,kk+1))
             if (capy(i,j,k) .ne. zero) then
                apinv = one/capy(i,j,k)
                cfcy(i,j,k,1) = fourth*apinv* &
                     ( fapy(ii  ,jj,kk  )*(half*ffcy(ii  ,jj,kk  ,1)-fourth) &
                     + fapy(ii+1,jj,kk  )*(half*ffcy(ii+1,jj,kk  ,1)+fourth) &
                     + fapy(ii  ,jj,kk+1)*(half*ffcy(ii  ,jj,kk+1,1)-fourth) &
                     + fapy(ii+1,jj,kk+1)*(half*ffcy(ii+1,jj,kk+1,1)+fourth) )
                cfcy(i,j,k,2) = fourth*apinv* &
                     ( fapy(ii  ,jj,kk  )*(half*ffcy(ii  ,jj,kk  ,2)-fourth) &
                     + fapy(ii+1,jj,kk  )*(half*ffcy(ii+1,jj,kk  ,2)-fourth) &
                     + fapy(ii  ,jj,kk+1)*(half*ffcy(ii  ,jj,kk+1,2)+fourth) &
                     + fapy(ii+1,jj,kk+1)*(half*ffcy(ii+1,jj,kk+1,2)+fourth) )                
             end if
          end do
       end do
    end do

    do       k = zlo(3), zhi(3)
       kk = k*2
       do    j = zlo(2), zhi(2)
          jj = j*2
          do i = zlo(1), zhi(1)
             ii = i*2

             capz(i,j,k) = fourth*(fapz(ii,jj,kk)+fapz(ii+1,jj,kk)+fapz(ii,jj+1,kk)+fapz(ii+1,jj+1,kk))
             if (capz(i,j,k) .ne. zero) then
                apinv = one/capz(i,j,k)
                cfcz(i,j,k,1) = fourth*apinv* &
                     ( fapz(ii  ,jj  ,kk)*(half*ffcz(ii  ,jj  ,kk,1)-fourth) &
                     + fapz(ii+1,jj  ,kk)*(half*ffcz(ii+1,jj  ,kk,1)+fourth) &
                     + fapz(ii  ,jj+1,kk)*(half*ffcz(ii  ,jj+1,kk,1)-fourth) &
                     + fapz(ii+1,jj+1,kk)*(half*ffcz(ii+1,jj+1,kk,1)+fourth) )
                cfcz(i,j,k,2) = fourth*apinv* &
                     ( fapz(ii  ,jj  ,kk)*(half*ffcz(ii  ,jj  ,kk,2)-fourth) &
                     + fapz(ii+1,jj  ,kk)*(half*ffcz(ii+1,jj  ,kk,2)-fourth) &
                     + fapz(ii  ,jj+1,kk)*(half*ffcz(ii  ,jj+1,kk,2)+fourth) &
                     + fapz(ii+1,jj+1,kk)*(half*ffcz(ii+1,jj+1,kk,2)+fourth) )                
             end if
          end do
       end do
    end do

  end subroutine amrex_eb2_coarsen_from_fine


  subroutine amrex_eb2_build_cellflag_from_ap (lo, hi, &
       apx, xlo, xhi, apy, ylo, yhi, apz, zlo, zhi, cflag, flo, fhi) &
       bind(c,name='amrex_eb2_build_cellflag_from_ap')
    integer, dimension(3) :: lo, hi, xlo, xhi, ylo, yhi, zlo, zhi, flo, fhi
    real(amrex_real), intent(in   ) ::  apx (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(amrex_real), intent(in   ) ::  apy (ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    real(amrex_real), intent(in   ) ::  apz (zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3))
    integer         , intent(inout) :: cflag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

    integer :: i,j,k, flg

    ! By default, all neighbors are already set.
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             cflag(i,j,k) = clear_allneighbors(cflag(i,j,k))
                
             if (.not.is_covered_cell(cflag(i,j,k))) then
                flg = clear_allneighbors(cflag(i,j,k))
                
                if (apx(i  ,j,k).ne.zero) flg = set_neighbor(flg, -1,  0,  0)
                if (apx(i+1,j,k).ne.zero) flg = set_neighbor(flg,  1,  0,  0)
                if (apy(i,j  ,k).ne.zero) flg = set_neighbor(flg,  0, -1,  0)
                if (apy(i,j+1,k).ne.zero) flg = set_neighbor(flg,  0,  1,  0)
                if (apz(i,j,k  ).ne.zero) flg = set_neighbor(flg,  0,  0, -1)
                if (apz(i,j,k+1).ne.zero) flg = set_neighbor(flg,  0,  0,  1)
                
                if ( (apx(i,j,k).ne.zero .and. apy(i-1,j,k).ne.zero) .or. &
                     (apy(i,j,k).ne.zero .and. apx(i,j-1,k).ne.zero) ) then
                   flg = set_neighbor(flg, -1, -1, 0)
                   if (apz(i-1,j-1,k  ).ne.zero) flg = set_neighbor(flg,-1,-1,-1)
                   if (apz(i-1,j-1,k+1).ne.zero) flg = set_neighbor(flg,-1,-1, 1)
                end if
                
                if ( (apx(i+1,j,k).ne.zero .and. apy(i+1,j,k).ne.zero) .or. &
                     (apy(i,j,k).ne.zero .and. apx(i+1,j-1,k).ne.zero) ) then
                   flg = set_neighbor(flg, 1, -1, 0)
                   if (apz(i+1,j-1,k  ).ne.zero) flg = set_neighbor(flg, 1,-1,-1)
                   if (apz(i+1,j-1,k+1).ne.zero) flg = set_neighbor(flg, 1,-1, 1)
                end if
                
                if ( (apx(i,j,k).ne.zero .and. apy(i-1,j+1,k).ne.zero) .or. &
                     (apy(i,j+1,k).ne.zero .and. apx(i,j+1,k).ne.zero) ) then
                   flg = set_neighbor(flg, -1, 1, 0)
                   if (apz(i-1,j+1,k  ).ne.zero) flg = set_neighbor(flg,-1, 1,-1)
                   if (apz(i-1,j+1,k+1).ne.zero) flg = set_neighbor(flg,-1, 1, 1)
                end if
                
                if ( (apx(i+1,j,k).ne.zero .and. apy(i+1,j+1,k).ne.zero) .or. &
                     (apy(i,j+1,k).ne.zero .and. apx(i+1,j+1,k).ne.zero) ) then
                   flg = set_neighbor(flg, 1, 1, 0)
                   if (apz(i+1,j+1,k  ).ne.zero) flg = set_neighbor(flg, 1, 1,-1)
                   if (apz(i+1,j+1,k+1).ne.zero) flg = set_neighbor(flg, 1, 1, 1)
                end if
                
                if ( (apx(i,j,k).ne.zero .and. apz(i-1,j,k).ne.zero) .or. &
                     (apz(i,j,k).ne.zero .and. apx(i,j,k-1).ne.zero) ) then
                   flg = set_neighbor(flg, -1, 0, -1)
                   if (apy(i-1,j  ,k-1).ne.zero) flg = set_neighbor(flg,-1,-1,-1)
                   if (apy(i-1,j+1,k-1).ne.zero) flg = set_neighbor(flg,-1, 1,-1)
                end if
                
                if ( (apx(i+1,j,k).ne.zero .and. apz(i+1,j,k).ne.zero) .or. &
                     (apz(i,j,k).ne.zero .and. apx(i+1,j,k-1).ne.zero) ) then
                   flg = set_neighbor(flg, 1, 0, -1)
                   if (apy(i+1,j  ,k-1).ne.zero) flg = set_neighbor(flg, 1,-1,-1)
                   if (apy(i+1,j+1,k-1).ne.zero) flg = set_neighbor(flg, 1, 1,-1)
                end if
                
                if ( (apx(i,j,k).ne.zero .and. apz(i-1,j,k+1).ne.zero) .or. &
                     (apz(i,j,k+1).ne.zero .and. apx(i,j,k+1).ne.zero) ) then
                   flg = set_neighbor(flg, -1, 0, 1)
                   if (apy(i-1,j  ,k+1).ne.zero) flg = set_neighbor(flg,-1,-1, 1)
                   if (apy(i-1,j+1,k+1).ne.zero) flg = set_neighbor(flg,-1, 1, 1)
                end if
                
                if ( (apx(i+1,j,k).ne.zero .and. apz(i+1,j,k+1).ne.zero) .or. &
                     (apz(i,j,k+1).ne.zero .and. apx(i+1,j,k+1).ne.zero) ) then
                   flg = set_neighbor(flg, 1, 0, 1)
                   if (apy(i+1,j  ,k+1).ne.zero) flg = set_neighbor(flg, 1,-1, 1)
                   if (apy(i+1,j+1,k+1).ne.zero) flg = set_neighbor(flg, 1, 1, 1)
                end if
                
                if ( (apy(i,j,k).ne.zero .and. apz(i,j-1,k).ne.zero) .or. &
                     (apz(i,j,k).ne.zero .and. apy(i,j,k-1).ne.zero) ) then
                   flg = set_neighbor(flg, 0, -1, -1)
                   if (apx(i  ,j-1,k-1).ne.zero) flg = set_neighbor(flg,-1,-1,-1)
                   if (apx(i+1,j-1,k-1).ne.zero) flg = set_neighbor(flg, 1,-1,-1)
                end if
                
                if ( (apy(i,j+1,k).ne.zero .and. apz(i,j+1,k).ne.zero) .or. &
                     (apz(i,j,k).ne.zero .and. apy(i,j+1,k-1).ne.zero) ) then
                   flg = set_neighbor(flg, 0, 1, -1)
                   if (apx(i  ,j+1,k-1).ne.zero) flg = set_neighbor(flg,-1, 1,-1)
                   if (apx(i+1,j+1,k-1).ne.zero) flg = set_neighbor(flg, 1, 1,-1)
                end if
                
                if ( (apy(i,j,k).ne.zero .and. apz(i,j-1,k+1).ne.zero) .or. &
                     (apz(i,j,k+1).ne.zero .and. apy(i,j,k+1).ne.zero) ) then
                   flg = set_neighbor(flg, 0, -1, 1)
                   if (apx(i  ,j-1,k+1).ne.zero) flg = set_neighbor(flg,-1,-1, 1)
                   if (apx(i+1,j-1,k+1).ne.zero) flg = set_neighbor(flg, 1,-1, 1)
                end if
                
                if ( (apy(i,j+1,k).ne.zero .and. apz(i,j+1,k+1).ne.zero) .or. &
                     (apz(i,j,k+1).ne.zero .and. apy(i,j+1,k+1).ne.zero) ) then
                   flg = set_neighbor(flg, 0, 1, 1)
                   if (apx(i  ,j+1,k+1).ne.zero) flg = set_neighbor(flg,-1, 1, 1)
                   if (apx(i+1,j+1,k+1).ne.zero) flg = set_neighbor(flg, 1, 1, 1)
                end if
                   
                cflag(i,j,k) = flg
             end if
          end do
       end do
    end do
  end subroutine amrex_eb2_build_cellflag_from_ap


  subroutine amrex_eb2_check_mvmc (cclo, cchi, ndlo, ndhi, cls, clo, chi, fls, flo, fhi, &
       ncuts, tlo, thi, ierr) &
       bind(c,name='amrex_eb2_check_mvmc')
    integer, dimension(3), intent(in) :: cclo,cchi,ndlo,ndhi,clo,chi,flo,fhi,tlo,thi
    real(amrex_real), intent(inout) :: cls(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(amrex_real), intent(in   ) :: fls(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    integer, intent(inout) :: ncuts(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),6)
    integer, intent(inout) :: ierr

    integer :: i,j,k, ii,jj,kk, n, nopen

    do       k = ndlo(3), ndhi(3)
       kk = k*2
       do    j = ndlo(2), ndhi(2)
          jj = j*2
          do i = ndlo(1), ndhi(1)
             ii = i*2
             cls(i,j,k) = fls(ii,jj,kk)
          end do
       end do
    end do

    ! x-edges
    do       k = cclo(3), cchi(3)+1
       kk = k*2
       do    j = cclo(2), cchi(2)+1
          jj = j*2
          do i = cclo(1), cchi(1)
             ii = i*2
             ncuts(i,j,k,1) = 0
             if (has_cut(fls(ii,jj,kk),fls(ii+1,jj,kk))) then
                ncuts(i,j,k,1) = ncuts(i,j,k,1)+1
             end if
             if (has_cut(fls(ii+1,jj,kk),fls(ii+2,jj,kk))) then
                ncuts(i,j,k,1) = ncuts(i,j,k,1)+1
             end if
             if (ncuts(i,j,k,1) .eq. 2) then
                ierr = 1
                return
             end if
          end do
       end do
    end do

    ! y-edges
    do       k = cclo(3), cchi(3)+1
       kk = k*2
       do    j = cclo(2), cchi(2)
          jj = j*2
          do i = cclo(1), cchi(1)+1
             ii = i*2
             ncuts(i,j,k,2) = 0
             if (has_cut(fls(ii,jj,kk),fls(ii,jj+1,kk))) then
                ncuts(i,j,k,2) = ncuts(i,j,k,2)+1
             end if
             if (has_cut(fls(ii,jj+1,kk),fls(ii,jj+2,kk))) then
                ncuts(i,j,k,2) = ncuts(i,j,k,2)+1
             end if
             if (ncuts(i,j,k,2) .eq. 2) then
                ierr = 1
                return
             end if
          end do
       end do
    end do

    ! z-edges
    do       k = cclo(3), cchi(3)
       kk = k*2
       do    j = cclo(2), cchi(2)+1
          jj = j*2
          do i = cclo(1), cchi(1)+1
             ii = i*2
             ncuts(i,j,k,3) = 0
             if (has_cut(fls(ii,jj,kk),fls(ii,jj,kk+1))) then
                ncuts(i,j,k,3) = ncuts(i,j,k,3)+1
             end if
             if (has_cut(fls(ii,jj,kk+1),fls(ii,jj,kk+2))) then
                ncuts(i,j,k,3) = ncuts(i,j,k,3)+1
             end if
             if (ncuts(i,j,k,3) .eq. 2) then
                ierr = 1
                return
             end if
          end do
       end do
    end do

    ! x-faces
    do       k = cclo(3), cchi(3)
       do    j = cclo(2), cchi(2)
          do i = cclo(1), cchi(1)+1
             n =    ncuts(i,j,k,2) + ncuts(i,j,k+1,2) &
                  + ncuts(i,j,k,3) + ncuts(i,j+1,k,3)
             if (n .eq. 0) then
                ncuts(i,j,k,4) = 0
             else if (n .eq. 2) then
                ncuts(i,j,k,4) = 1
             else if (n .eq. 4) then
                ierr = 1
                return
             else
                call amrex_error("amrex_eb2_check_mvmc: how did this happen? wrong nubmer of cuts on x-face", n)
             end if
          end do
       end do
    end do
    
    ! y-faces
    do       k = cclo(3), cchi(3)
       do    j = cclo(2), cchi(2)+1
          do i = cclo(1), cchi(1)
             n =    ncuts(i,j,k,1) + ncuts(i,j,k+1,1) &
                  + ncuts(i,j,k,3) + ncuts(i+1,j,k,3)
             if (n .eq. 0) then
                ncuts(i,j,k,5) = 0
             else if (n .eq. 2) then
                ncuts(i,j,k,5) = 1
             else if (n .eq. 4) then
                ierr = 1
                return
             else
                call amrex_error("amrex_eb2_check_mvmc: how did this happen? wrong nubmer of cuts on y-face", n)
             end if
          end do
       end do
    end do

    ! z-faces
    do       k = cclo(3), cchi(3)+1
       do    j = cclo(2), cchi(2)
          do i = cclo(1), cchi(1)
             n =    ncuts(i,j,k,1) + ncuts(i,j+1,k,1) &
                  + ncuts(i,j,k,2) + ncuts(i+1,j,k,2)
             if (n .eq. 0) then
                ncuts(i,j,k,6) = 0
             else if (n .eq. 2) then
                ncuts(i,j,k,6) = 1
             else if (n .eq. 4) then
                ierr = 1
                return
             else
                call amrex_error("amrex_eb2_check_mvmc: how did this happen? wrong nubmer of cuts on x-face", n)
             end if
          end do
       end do
    end do

    do       k = cclo(3), cchi(3)
       do    j = cclo(2), cchi(2)
          do i = cclo(1), cchi(1)
             if ( ncuts(i  ,j  ,k  ,4).eq.1 .and. &
                  ncuts(i+1,j  ,k  ,4).eq.1 .and. &
                  ncuts(i  ,j  ,k  ,5).eq.1 .and. &
                  ncuts(i  ,j+1,k  ,5).eq.1 .and. &
                  ncuts(i  ,j  ,k  ,6).eq.1 .and. &
                  ncuts(i  ,j  ,k+1,6).eq.1 ) then
                ii = i*2
                jj = j*2
                kk = k*2
                nopen = 0
                if (fls(ii  ,jj  ,kk  ) .lt. zero) nopen = nopen + 1
                if (fls(ii+2,jj  ,kk  ) .lt. zero) nopen = nopen + 1
                if (fls(ii  ,jj+2,kk  ) .lt. zero) nopen = nopen + 1
                if (fls(ii+2,jj+2,kk  ) .lt. zero) nopen = nopen + 1
                if (fls(ii  ,jj  ,kk+2) .lt. zero) nopen = nopen + 1
                if (fls(ii+2,jj  ,kk+2) .lt. zero) nopen = nopen + 1
                if (fls(ii  ,jj+2,kk+2) .lt. zero) nopen = nopen + 1
                if (fls(ii+2,jj+2,kk+2) .lt. zero) nopen = nopen + 1
                if (nopen .eq. 2 .or. nopen .eq. 6) then
                   ierr = 1
                   return
                else if (nopen .ne. 4) then
                   call amrex_error("amrex_eb2_check_mvmc: how did this happen? nopen", nopen)
                end if
             end if
          end do
       end do
    end do
    
  contains
    pure function has_cut(a,b)
      real(amrex_real), intent(in) :: a,b
      logical has_cut
      has_cut = (a.ge.zero .and. b.lt.zero) .or. (b.ge.zero .and. a.lt.zero)
    end function has_cut

  end subroutine amrex_eb2_check_mvmc
  
end module amrex_eb2_3d_module
