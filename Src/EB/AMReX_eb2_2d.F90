
module amrex_eb2_2d_moudle

  use amrex_error_module
  use amrex_fort_module
  use amrex_constants_module, only : zero, one, half, fourth, sixth, eighth
  use amrex_ebcellflag_module, only : regular_cell => regular, covered_cell => covered, &
       is_regular_cell, is_single_valued_cell, is_covered_cell, get_cell_type, get_neighbor_cells, &
       set_regular_cell, set_single_valued_cell, set_covered_cell, &
       set_neighbor, clear_neighbor, is_neighbor
  implicit none

  integer, parameter :: regular = 0
  integer, parameter :: covered = 1
  integer, parameter :: irregular = 2
  integer, parameter :: unknown = 3

  real(amrex_real), private, parameter :: small = 1.d-14

  private
  public :: amrex_eb2_coarsen_from_fine, amrex_eb2_build_cellflag_from_ap, amrex_eb2_check_mvmc

contains

  subroutine amrex_eb2_coarsen_from_fine (lo, hi, xlo, xhi, ylo, yhi, &
       cvol, cvlo, cvhi, fvol, fvlo, fvhi, ccent, cclo, cchi, fcent, fclo, fchi, &
       cba, cbalo, cbahi, fba, fbalo, fbahi, cbc, cbclo, cbchi, fbc, fbclo, fbchi, &
       cbn, cbnlo, cbnhi, fbn, fbnlo, fbnhi, capx, caxlo, caxhi, fapx, faxlo, faxhi, &
       capy, caylo, cayhi, fapy, faylo, fayhi, &
       cfcx, cfxlo, cfxhi, ffcx, ffxlo, ffxhi, cfcy, cfylo, cfyhi, ffcy, ffylo, ffyhi, &
       cflag, cflo, cfhi, fflag, fflo, ffhi, ierr) &
       bind(c, name='amrex_eb2_coarsen_from_fine')
    integer, dimension(2), intent(in) :: lo, hi, xlo, xhi, ylo, yhi, &
         cvlo, cvhi,  fvlo, fvhi, cclo, cchi, fclo, fchi, &
         cbalo, cbahi, fbalo, fbahi, cbclo, cbchi, fbclo, fbchi, &
         cbnlo, cbnhi, fbnlo, fbnhi, caxlo, caxhi, faxlo, faxhi, &
         caylo, cayhi, faylo, fayhi, &
         cfxlo, cfxhi, ffxlo, ffxhi, cfylo, cfyhi, ffylo, ffyhi, &
         cflo, cfhi, fflo, ffhi
    integer, intent(inout) :: ierr
    real(amrex_real), intent(inout) :: cvol ( cvlo(1): cvhi(1), cvlo(2): cvhi(2))
    real(amrex_real), intent(in   ) :: fvol ( fvlo(1): fvhi(1), fvlo(2): fvhi(2))
    real(amrex_real), intent(inout) :: ccent( cclo(1): cchi(1), cclo(2): cchi(2),2)
    real(amrex_real), intent(in   ) :: fcent( fclo(1): fchi(1), fclo(2): fchi(2),2)
    real(amrex_real), intent(inout) :: cba  (cbalo(1):cbahi(1),cbalo(2):cbahi(2))
    real(amrex_real), intent(in   ) :: fba  (fbalo(1):fbahi(1),fbalo(2):fbahi(2))
    real(amrex_real), intent(inout) :: cbc  (cbclo(1):cbchi(1),cbclo(2):cbchi(2),2)
    real(amrex_real), intent(in   ) :: fbc  (fbclo(1):fbchi(1),fbclo(2):fbchi(2),2)
    real(amrex_real), intent(inout) :: cbn  (cbnlo(1):cbnhi(1),cbnlo(2):cbnhi(2),2)
    real(amrex_real), intent(in   ) :: fbn  (fbnlo(1):fbnhi(1),fbnlo(2):fbnhi(2),2)
    real(amrex_real), intent(inout) :: capx (caxlo(1):caxhi(1),caxlo(2):caxhi(2))
    real(amrex_real), intent(in   ) :: fapx (faxlo(1):faxhi(1),faxlo(2):faxhi(2))
    real(amrex_real), intent(inout) :: capy (caylo(1):cayhi(1),caylo(2):cayhi(2))
    real(amrex_real), intent(in   ) :: fapy (faylo(1):fayhi(1),faylo(2):fayhi(2))
    real(amrex_real), intent(inout) :: cfcx (cfxlo(1):cfxhi(1),cfxlo(2):cfxhi(2))
    real(amrex_real), intent(in   ) :: ffcx (ffxlo(1):ffxhi(1),ffxlo(2):ffxhi(2))
    real(amrex_real), intent(inout) :: cfcy (cfylo(1):cfyhi(1),cfylo(2):cfyhi(2))
    real(amrex_real), intent(in   ) :: ffcy (ffylo(1):ffyhi(1),ffylo(2):ffyhi(2))
    integer         , intent(inout) :: cflag( cflo(1): cfhi(1), cflo(2): cfhi(2))
    integer         , intent(in   ) :: fflag( fflo(1): ffhi(1), fflo(2): ffhi(2))

    integer :: i,j, ii,jj, ftype(2,2)
    real(amrex_real) :: cvolinv, cbainv, nx, ny, nfac, apinv

    do    j = lo(2), hi(2)
       jj = j*2
       do i = lo(1), hi(1)
          ii = i*2

          ftype = get_cell_type(fflag(ii:ii+1,jj:jj+1))
          if (all(ftype.eq.regular_cell)) then
             ! nothing to do
          else if (all(ftype.eq.covered_cell)) then
             call set_covered_cell(cflag(i,j))
             cvol(i,j) = zero
          else
             call set_single_valued_cell(cflag(i,j))

             cvol(i,j) = fourth*sum(fvol(ii:ii+1,jj:jj+1))
             cvolinv = one/cvol(i,j)
                
             ccent(i,j,1) = fourth*cvolinv* &
                  ( fvol(ii  ,jj  )*(half*fcent(ii  ,jj  ,1)-fourth) &
                  + fvol(ii+1,jj  )*(half*fcent(ii+1,jj  ,1)+fourth) &
                  + fvol(ii  ,jj+1)*(half*fcent(ii  ,jj+1,1)-fourth) &
                  + fvol(ii+1,jj+1)*(half*fcent(ii+1,jj+1,1)+fourth) )
             ccent(i,j,2) = fourth*cvolinv* &
                  ( fvol(ii  ,jj  )*(half*fcent(ii  ,jj  ,2)-fourth) &
                  + fvol(ii+1,jj  )*(half*fcent(ii+1,jj  ,2)-fourth) &
                  + fvol(ii  ,jj+1)*(half*fcent(ii  ,jj+1,2)+fourth) &
                  + fvol(ii+1,jj+1)*(half*fcent(ii+1,jj+1,2)+fourth) )
                
             cba(i,j) = half*sum(fba(ii:ii+1,jj:jj+1))
             cbainv = one/cba(i,j)
                
             cbc(i,j,1) = half*cbainv* &
                  ( fba(ii  ,jj  )*(half*fbc(ii  ,jj  ,1)-fourth) &
                  + fba(ii+1,jj  )*(half*fbc(ii+1,jj  ,1)+fourth) &
                  + fba(ii  ,jj+1)*(half*fbc(ii  ,jj+1,1)-fourth) &
                  + fba(ii+1,jj+1)*(half*fbc(ii+1,jj+1,1)+fourth) )
             cbc(i,j,2) = half*cbainv* &
                  ( fba(ii  ,jj  )*(half*fbc(ii  ,jj  ,2)-fourth) &
                  + fba(ii+1,jj  )*(half*fbc(ii+1,jj  ,2)-fourth) &
                  + fba(ii  ,jj+1)*(half*fbc(ii  ,jj+1,2)+fourth) &
                  + fba(ii+1,jj+1)*(half*fbc(ii+1,jj+1,2)+fourth) )
                
             nx =   fbn(ii  ,jj  ,1)*fba(ii  ,jj  ) &
                  + fbn(ii+1,jj  ,1)*fba(ii+1,jj  ) &
                  + fbn(ii  ,jj+1,1)*fba(ii  ,jj+1) &
                  + fbn(ii+1,jj+1,1)*fba(ii+1,jj+1)
             ny =   fbn(ii  ,jj  ,2)*fba(ii  ,jj  ) &
                  + fbn(ii+1,jj  ,2)*fba(ii+1,jj  ) &
                  + fbn(ii  ,jj+1,2)*fba(ii  ,jj+1) &
                  + fbn(ii+1,jj+1,2)*fba(ii+1,jj+1)
             if (nx.eq.zero .and. ny.eq.zero) then
                ierr = 1
                return
             end if
             nfac = one/sqrt(nx*nx+ny*ny)
             cbn(i,j,1) = nx*nfac
             cbn(i,j,2) = ny*nfac
          end if
       end do
    end do

    do    j = xlo(2), xhi(2)
       jj = j*2
       do i = xlo(1), xhi(1)
          ii = i*2
          
          capx(i,j) = half*(fapx(ii,jj)+fapx(ii,jj+1))
          if (capx(i,j) .ne. zero) then
             apinv = one/capx(i,j)
             cfcx(i,j) = half*apinv* &
                  ( fapx(ii,jj  )*(half*ffcx(ii,jj  )-fourth) &
                  + fapx(ii,jj+1)*(half*ffcx(ii,jj+1)+fourth) )
          end if
       end do
    end do

    do    j = ylo(2), yhi(2)
       jj = j*2
       do i = ylo(1), yhi(1)
          ii = i*2
          
          capy(i,j) = half*(fapy(ii,jj)+fapy(ii+1,jj))
          if (capy(i,j) .ne. zero) then
             apinv = one/capy(i,j)
             cfcy(i,j) = half*apinv* &
                  ( fapy(ii  ,jj)*(half*ffcy(ii  ,jj)-fourth) &
                  + fapy(ii+1,jj)*(half*ffcy(ii+1,jj)+fourth) )
          end if
       end do
    end do

  end subroutine amrex_eb2_coarsen_from_fine


  subroutine amrex_eb2_build_cellflag_from_ap (lo, hi, &
       apx, xlo, xhi, apy, ylo, yhi, cflag, flo, fhi) &
       bind(c,name='amrex_eb2_build_cellflag_from_ap')
    integer, dimension(2), intent(in) :: lo, hi, xlo, xhi, ylo, yhi, flo, fhi
    real(amrex_real), intent(in   ) ::  apx (xlo(1):xhi(1),xlo(2):xhi(2))
    real(amrex_real), intent(in   ) ::  apy (ylo(1):yhi(1),ylo(2):yhi(2))
    integer         , intent(inout) :: cflag(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i,j, flg

    ! By default, all neighbors are already set.
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          flg = cflag(i,j)

          if (apx(i  ,j  ).eq.zero) flg = clear_neighbor(flg, -1,  0)
          if (apx(i+1,j  ).eq.zero) flg = clear_neighbor(flg,  1,  0)
          if (apy(i  ,j  ).eq.zero) flg = clear_neighbor(flg,  0, -1)
          if (apy(i  ,j+1).eq.zero) flg = clear_neighbor(flg,  0,  1)

          if (apx(i,j).ne.zero .and. apy(i-1,j).ne.zero) then
          else if (apx(i,j-1).ne.zero .and. apy(i,j).ne.zero) then
          else
             flg = clear_neighbor(flg,-1,-1)
          end if

          if (apx(i+1,j).ne.zero .and. apy(i+1,j).ne.zero) then
          else if (apx(i+1,j-1).ne.zero .and. apy(i,j).ne.zero) then
          else
             flg = clear_neighbor(flg,1,-1)
          end if

          if (apx(i,j).ne.zero .and. apy(i-1,j+1).ne.zero) then
          else if (apx(i,j+1).ne.zero .and. apy(i,j+1).ne.zero) then
          else
             flg = clear_neighbor(flg,-1,1)
          end if

          if (apx(i+1,j).ne.zero .and. apy(i+1,j+1).ne.zero) then
          else if (apx(i+1,j+1).ne.zero .and. apy(i,j+1).ne.zero) then
          else
             flg = clear_neighbor(flg,1,1)
          end if

          cflag(i,j) = flg
       end do
    end do
  end subroutine amrex_eb2_build_cellflag_from_ap

  
  subroutine amrex_eb2_check_mvmc (cclo, cchi, ndlo, ndhi, cls, clo, chi, fls, flo, fhi, ierr) &
       bind(c,name='amrex_eb2_check_mvmc')
    integer, dimension(2), intent(in) :: cclo, cchi, ndlo, ndhi, clo, chi, flo, fhi
    real(amrex_real), intent(inout) :: cls(clo(1):chi(1),clo(2):chi(2))
    real(amrex_real), intent(in   ) :: fls(flo(1):fhi(1),flo(2):fhi(2))
    integer, intent(inout) :: ierr

    integer :: i,j, ii, jj, ncuts

    do    j = ndlo(2), ndhi(2)
       do i = ndlo(1), ndhi(1)
          cls(i,j) = fls(2*i,2*j)
       end do
    end do

    do    j = cclo(2), cchi(2)
       jj = 2*j
       do i = cclo(1), cchi(1)
          ii = 2*i
          ncuts = 0
          if (has_cut(fls(ii  ,jj  ),fls(ii+1,jj  ))) ncuts = ncuts+1
          if (has_cut(fls(ii+1,jj  ),fls(ii+2,jj  ))) ncuts = ncuts+1
          if (has_cut(fls(ii  ,jj+2),fls(ii+1,jj+2))) ncuts = ncuts+1
          if (has_cut(fls(ii+1,jj+2),fls(ii+2,jj+2))) ncuts = ncuts+1
          if (has_cut(fls(ii  ,jj  ),fls(ii  ,jj+1))) ncuts = ncuts+1
          if (has_cut(fls(ii  ,jj+1),fls(ii  ,jj+2))) ncuts = ncuts+1
          if (has_cut(fls(ii+2,jj  ),fls(ii+2,jj+1))) ncuts = ncuts+1
          if (has_cut(fls(ii+2,jj+1),fls(ii+2,jj+2))) ncuts = ncuts+1
          if (ncuts .ne. 0 .and. ncuts .ne. 2) then
             ierr = 2
             return
          end if
       end do
    end do
  contains
    pure function has_cut(a,b)
      real(amrex_real), intent(in) :: a,b
      logical has_cut
      has_cut = (a.ge.zero .and. b.lt.zero) .or. (b.ge.zero .and. a.lt.zero)
    end function has_cut
  end subroutine amrex_eb2_check_mvmc

end module amrex_eb2_2d_moudle
