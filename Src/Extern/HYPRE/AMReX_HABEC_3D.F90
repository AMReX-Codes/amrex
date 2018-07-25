module amrex_habec_module

  ! habec is Hypre abec, where abec is the form of the linear equation
  ! we are solving:
  ! 
  ! alpha*phi - div(beta*grad phi) + div(\vec{c}*phi) 

  use iso_c_binding
  use amrex_hypre_fort_module, only : hypre_int
  use amrex_fort_module, only : rt => amrex_real
  use amrex_lo_bctypes_module, only : amrex_lo_dirichlet, amrex_lo_neumann
  use amrex_error_module, only : amrex_error
  use amrex_constants_module, only : zero, one, half, three
  implicit none

contains

  subroutine amrex_hpacoef (lo, hi, mat, a, alo, ahi, sa) bind(c,name='amrex_hpacoef')
    integer, intent(in) :: lo(3), hi(3), alo(3), ahi(3)
    real(rt), intent(inout) :: mat(0:6, lo(1): hi(1), lo(2): hi(2), lo(3): hi(3))
    real(rt), intent(in   ) ::   a(    alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
    real(rt), intent(in) :: sa
    integer :: i, j, k
    if (sa .eq. zero) then
       !$omp parallel do private(i,j,k) collapse(2)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = zero
             enddo
          enddo
       enddo
       !$omp end parallel do
    else
       !$omp parallel do private(i,j,k) collapse(2)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = sa * a(i,j,k)
             enddo
          enddo
       enddo
       !$omp end parallel do
    end if
  end subroutine amrex_hpacoef

  subroutine amrex_hpbcoef (lo, hi, mat, b, blo, bhi, sb, dx, idim) &
       bind(c,name='amrex_hpbcoef')
    integer, intent(in) :: lo(3), hi(3), blo(3), bhi(3), idim
    real(rt), intent(inout) :: mat(0:6, lo(1): hi(1), lo(2): hi(2), lo(3): hi(3))
    real(rt), intent(in   ) ::   b(    blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))
    real(rt), intent(in) :: sb, dx(3)

    integer :: i, j, k
    real(rt) :: fac

    fac = sb / dx(idim+1)**2 

    if (idim .eq. 0) then
       !$omp parallel do private(i,j,k) collapse(2)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i+1,j,k))
                mat(1,i,j,k) = - fac * b(i,j,k)
                mat(2,i,j,k) = - fac * b(i+1,j,k)
             enddo
          enddo
       enddo
       !$omp end parallel do
    else if (idim .eq. 1) then
       !$omp parallel do private(i,j,k) collapse(2)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i,j+1,k))
                mat(3,i,j,k) = - fac * b(i,j,k)
                mat(4,i,j,k) = - fac * b(i,j+1,k)
             enddo
          enddo
       enddo
       !$omp end parallel do
    else
       !$omp parallel do private(i,j,k) collapse(2)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                mat(0,i,j,k) = mat(0,i,j,k) + fac * (b(i,j,k) + b(i,j,k+1))
                mat(5,i,j,k) = - fac * b(i,j,k)
                mat(6,i,j,k) = - fac * b(i,j,k+1)
             enddo
          enddo
       enddo
       !$omp end parallel do
    endif
    
  end subroutine amrex_hpbcoef

  
  subroutine amrex_hpmat (lo, hi, mat, b, blo, bhi, mask, mlo, mhi, &
       sb, dx, cdir, bct, bcl, bho) bind(c,name='amrex_hpmat')
    integer, intent(in) :: lo(3), hi(3), blo(3), bhi(3), mlo(3), mhi(3), cdir, bct, bho
    real(rt), intent(inout) ::  mat(0:6, lo(1): hi(1), lo(2): hi(2), lo(3): hi(3))
    real(rt), intent(in   ) ::    b(    blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3))
    integer , intent(in   ) :: mask(    mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    real(rt), intent(in) :: sb, dx(3), bcl

    integer :: i, j, k
    real(rt) :: fac, h, h2, h3, bf1, bf2

    if (cdir .eq. 0 .or. cdir .eq. 3) then
       h = dx(1)
    elseif (cdir .eq. 1 .or. cdir .eq. 4) then
       h = dx(2)
    else
       h = dx(3)
    endif
    fac = sb / (h**2)

    if (bct .eq. amrex_lo_dirichlet) then
       h2 = half * h
       if (bho.ge.1) then
          h3 = three * h2
          bf1 = fac * ((h3 - bcl) / (bcl + h2) - one)
          bf2 = fac * (bcl - h2) / (bcl + h3)
       else
          bf1 = fac * ( h / (bcl + h2) - one)          
          bf2 = zero
       end if
    else if (bct .eq. amrex_lo_neumann) then
       bf1 = -fac
       bf2 = zero
    else
       call amrex_error("hpmat: unsupported boundary type")       
    end if
    
    if (cdir .eq. 0) then
       i = lo(1)
       !$omp parallel do private(j,k)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             if (mask(i-1,j,k) .gt. 0) then
                mat(0,i,j,k) = mat(0,i,j,k) + bf1 * b(i,j,k)
                mat(1,i,j,k) = zero
                mat(2,i,j,k) = mat(2,i,j,k) + bf2 * b(i,j,k)
             endif
          enddo
       enddo
       !$omp end parallel do
    else if (cdir .eq. 3) then
       i = hi(1)
       !$omp parallel do private(j,k)
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             if (mask(i+1,j,k) .gt. 0) then
                mat(0,i,j,k) = mat(0,i,j,k) + bf1 * b(i+1,j,k)
                mat(2,i,j,k) = zero
                mat(1,i,j,k) = mat(1,i,j,k) + bf2 * b(i+1,j,k)
             endif
          enddo
       enddo
       !$omp end parallel do
    else if (cdir .eq. 1) then
       j = lo(2)
       !$omp parallel do private(i,k)
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             if (mask(i,j-1,k) .gt. 0) then
                mat(0,i,j,k) = mat(0,i,j,k) + bf1 * b(i,j,k)
                mat(3,i,j,k) = zero
                mat(4,i,j,k) = mat(4,i,j,k) + bf2 * b(i,j,k)
             endif
          enddo
       enddo
       !$omp end parallel do
    else if (cdir .eq. 4) then
       j = hi(2)
       !$omp parallel do private(i,k)
       do k = lo(3), hi(3)
          do i = lo(1), hi(1)
             if (mask(i,j+1,k) .gt. 0) then
                mat(0,i,j,k) = mat(0,i,j,k) + bf1 * b(i,j+1,k)
                mat(4,i,j,k) = zero
                mat(3,i,j,k) = mat(3,i,j,k) + bf2 * b(i,j+1,k)
             endif
          enddo
       enddo
       !$omp end parallel do
    else if (cdir .eq. 2) then
       k = lo(3)
       !$omp parallel do private(i,j)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (mask(i,j,k-1) .gt. 0) then
                mat(0,i,j,k) = mat(0,i,j,k) + bf1 * b(i,j,k)
                mat(5,i,j,k) = zero
                mat(6,i,j,k) = mat(6,i,j,k) + bf2 * b(i,j,k)
             endif
          enddo
       enddo
       !$omp end parallel do
    else if (cdir .eq. 5) then
       k = hi(3)
       !$omp parallel do private(i,j)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (mask(i,j,k+1) .gt. 0) then
                mat(0,i,j,k) = mat(0,i,j,k) + bf1 * b(i,j,k+1)
                mat(6,i,j,k) = zero
                mat(5,i,j,k) = mat(5,i,j,k) + bf2 * b(i,j,k+1)
             endif
          enddo
       enddo
       !$omp end parallel do
    else
       call amrex_error("hpmat: impossible face orientation")
    endif  
  end subroutine amrex_hpmat


  subroutine amrex_hpijmatrix (lo, hi, nrows, ncols, rows, cols, mat, &
       cell_id, clo, chi, cell_id_begin, a , alo, ahi, bx, xlo, xhi, &
       by, ylo, yhi, bz, zlo, zhi, sa, sb, dx, bct, bcl, bho) &
       bind(c,name='amrex_hpijmatrix')
    integer(hypre_int), intent(in) :: nrows, cell_id_begin;
    integer(hypre_int), dimension(0:nrows-1), intent(out) :: ncols, rows
    integer(hypre_int), dimension(0:nrows*7-1), intent(out) :: cols
    real(rt)          , dimension(0:nrows*7-1), intent(out) :: mat
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, alo, ahi, xlo, xhi, ylo, yhi, zlo, zhi
    integer(hypre_int), intent(in) :: cell_id(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt)          , intent(in) :: a      (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
    real(rt)          , intent(in) :: bx     (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(rt)          , intent(in) :: by     (ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    real(rt)          , intent(in) :: bz     (zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3))
    integer, intent(in) :: bct(0:5), bho
    real(rt), intent(in) :: sa, sb, dx(3), bcl(0:5)

    integer :: i,j,k, irow, imat, ic, cdir, idim
    integer(hypre_int) :: cols_tmp(0:6)
    real(rt) :: fac(3), mat_tmp(0:6)
    real(rt) :: bf1(0:5), bf2(0:5), h, h2, h3

    fac = sb/dx**2

    do cdir = 0, 5
       if (cdir .eq. 0 .or. cdir .eq. 3) then
          idim = 1
       else if (cdir .eq. 1 .or. cdir .eq. 4) then
          idim = 2
       else
          idim = 3
       end if
       h = dx(idim)
       if (bct(cdir) .eq. amrex_lo_dirichlet) then
          h2 = half * h
          if (bho.ge.1) then
             h3 = three * h2
             bf1(cdir) = fac(idim) * ((h3 - bcl(cdir)) / (bcl(cdir) + h2) - one)
             bf2(cdir) = fac(idim) * (bcl(cdir) - h2) / (bcl(cdir) + h3)
          else
             bf1(cdir) = fac(idim) * ( h / (bcl(cdir) + h2) - one)
             bf2(cdir) = zero
          end if
       else if (bct(cdir) .eq. amrex_lo_neumann) then
          bf1(cdir) = -fac(idim)
          bf2(cdir) = zero
       end if
    end do
    
    irow = 0
    imat = 0
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rows(irow) = cell_id(i,j,k)
             ncols(irow) = 0

             cols_tmp(0) = cell_id(i,j,k)
             mat_tmp(0) = sa * a(i,j,k) + fac(1)*(bx(i,j,k)+bx(i+1,j,k)) &
                  &                     + fac(2)*(by(i,j,k)+by(i,j+1,k)) &
                  &                     + fac(3)*(bz(i,j,k)+bz(i,j,k+1))

             cols_tmp(1) = cell_id(i-1,j,k)
             mat_tmp(1) = -fac(1)*bx(i,j,k)

             cols_tmp(2) = cell_id(i+1,j,k)
             mat_tmp(2) = -fac(1)*bx(i+1,j,k)

             cols_tmp(3) = cell_id(i,j-1,k)
             mat_tmp(3) = -fac(2)*by(i,j,k)

             cols_tmp(4) = cell_id(i,j+1,k)
             mat_tmp(4) = -fac(2)*by(i,j+1,k)

             cols_tmp(5) = cell_id(i,j,k-1)
             mat_tmp(5) = -fac(3)*bz(i,j,k)

             cols_tmp(6) = cell_id(i,j,k+1)
             mat_tmp(6) = -fac(3)*bz(i,j,k+1)

             if (i.eq.lo(1) .and. cell_id(i-1,j,k).lt.0) then
                cdir = 0
                mat_tmp(0) = mat_tmp(0) + bf1(cdir)*bx(i,j,k)
                mat_tmp(2) = mat_tmp(2) + bf2(cdir)*bx(i,j,k)
             end if

             if (i.eq.hi(1) .and. cell_id(i+1,j,k).lt.0) then
                cdir = 3
                mat_tmp(0) = mat_tmp(0) + bf1(cdir)*bx(i+1,j,k)
                mat_tmp(1) = mat_tmp(1) + bf2(cdir)*bx(i+1,j,k)
             end if

             if (j.eq.lo(2) .and. cell_id(i,j-1,k).lt.0) then
                cdir = 1
                mat_tmp(0) = mat_tmp(0) + bf1(cdir)*by(i,j,k)
                mat_tmp(4) = mat_tmp(4) + bf2(cdir)*by(i,j,k)
             end if

             if (j.eq.hi(2) .and. cell_id(i,j+1,k).lt.0) then
                cdir = 4
                mat_tmp(0) = mat_tmp(0) + bf1(cdir)*by(i,j+1,k)
                mat_tmp(3) = mat_tmp(3) + bf2(cdir)*by(i,j+1,k)
             end if

             if (k.eq.lo(3) .and. cell_id(i,j,k-1).lt.0) then
                cdir = 2
                mat_tmp(0) = mat_tmp(0) + bf1(cdir)*bz(i,j,k)
                mat_tmp(6) = mat_tmp(6) + bf2(cdir)*bz(i,j,k)
             end if

             if (k.eq.hi(3) .and. cell_id(i,j,k+1).lt.0) then
                cdir = 5
                mat_tmp(0) = mat_tmp(0) + bf1(cdir)*bz(i,j,k+1)
                mat_tmp(5) = mat_tmp(5) + bf2(cdir)*bz(i,j,k+1)
             end if
             
             do ic = 0, 6
                if (cols_tmp(ic) .ge. 0) then
                   ncols(irow) = ncols(irow) + 1
                   cols(imat) = cols_tmp(ic)
                   mat(imat) = mat_tmp(ic)
                   imat = imat + 1
                end if
             end do
             irow = irow+1

          end do
       end do
    end do
  end subroutine amrex_hpijmatrix


#ifdef AMREX_USE_EB

  subroutine amrex_hpeb_fill_cellid (lo, hi, nrows, cell_id, clo, chi, flag, flo, fhi) &
       bind(c,name='amrex_hpeb_fill_cellid')
    use amrex_ebcellflag_module, only : is_covered_cell
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, flo, fhi
    integer(hypre_int), intent(  out) :: nrows
    integer(hypre_int), intent(inout) :: cell_id(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    integer           , intent(in   ) :: flag   (flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

    integer :: i,j,k

    nrows = 0
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (.not.is_covered_cell(flag(i,j,k))) then
                cell_id(i,j,k) = nrows
                nrows = nrows+1
             end if
          end do
       end do
    end do
  end subroutine amrex_hpeb_fill_cellid

  subroutine amrex_hpeb_copy_from_vec (lo, hi, a, alo, ahi, v, nv, flag, flo, fhi) &
       bind(c,name='amrex_hpeb_copy_from_vec')
    use amrex_ebcellflag_module, only : is_covered_cell
    integer, dimension(3), intent(in) :: lo, hi, alo, ahi, flo, fhi
    integer(hypre_int), intent(in) :: nv
    real(rt), intent(inout) :: a   (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
    integer , intent(in   ) :: flag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(rt), intent(in) :: v(0:nv-1)
    
    integer :: i,j,k, nrows

    nrows = 0
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (.not.is_covered_cell(flag(i,j,k))) then
                a(i,j,k) = v(nrows)
                nrows = nrows+1
             end if
          end do
       end do
    end do
  end subroutine amrex_hpeb_copy_from_vec

  subroutine amrex_hpeb_copy_to_vec (lo, hi, a, alo, ahi, v, nv, flag, flo, fhi) &
       bind(c,name='amrex_hpeb_copy_to_vec')
    use amrex_ebcellflag_module, only : is_covered_cell
    integer, dimension(3), intent(in) :: lo, hi, alo, ahi, flo, fhi
    integer(hypre_int), intent(in) :: nv
    real(rt), intent(in) :: a   (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
    integer , intent(in) :: flag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(rt), intent(out) :: v(0:nv-1)
    
    integer :: i,j,k, nrows

    nrows = 0
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (.not.is_covered_cell(flag(i,j,k))) then
                v(nrows) = a(i,j,k)
                nrows = nrows+1
             end if
          end do
       end do
    end do
  end subroutine amrex_hpeb_copy_to_vec


  subroutine amrex_hpeb_ijmatrix (lo, hi, nrows, ncols, rows, cols, mat, &
       cell_id, clo, chi, cell_id_begin, a, alo, ahi, &
       bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, &
       flag, flo, fhi, vfrc, vlo, vhi, &
       apx, axlo, axhi, apy, aylo, ayhi, apz, azlo, azhi, &
       fcx, fxlo, fxhi, fcy, fylo, fyhi, fcz, fzlo, fzhi, &
       sa, sb, dx, bct, bcl, bho) &
       bind(c,name='amrex_hpeb_ijmatrix')
    use amrex_ebcellflag_module, only : is_covered_cell, is_regular_cell
    integer(hypre_int), intent(in) :: nrows, cell_id_begin;
    integer(hypre_int), dimension(0:nrows-1), intent(out) :: ncols, rows
    integer(hypre_int), dimension(0:nrows*27-1), intent(out) :: cols
    real(rt)          , dimension(0:nrows*27-1), intent(out) :: mat
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, alo, ahi, bxlo, bxhi, bylo, byhi, &
         bzlo, bzhi, flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi, azlo, azhi, &
         fxlo, fxhi, fylo, fyhi, fzlo, fzhi
    integer(hypre_int), intent(in) :: cell_id( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))
    real(rt)          , intent(in) :: a      ( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3))
    real(rt)          , intent(in) :: bx     (bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
    real(rt)          , intent(in) :: by     (bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
    real(rt)          , intent(in) :: bz     (bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
    integer           , intent(in) :: flag   ( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3))
    real(rt)          , intent(in) :: vfrc   ( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3))
    real(rt)          , intent(in) :: apx    (axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(rt)          , intent(in) :: apy    (aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(rt)          , intent(in) :: apz    (azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(rt)          , intent(in) :: fcx    (fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),2)
    real(rt)          , intent(in) :: fcy    (fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),2)
    real(rt)          , intent(in) :: fcz    (fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),2)
    integer, intent(in) :: bct(0:5), bho
    real(rt), intent(in) :: sa, sb, dx(3), bcl(0:5)

    integer :: i,j,k, irow, imat, cdir, idim, ii,jj,kk, ioff, joff, koff
    real(rt) :: fac(3), mat_tmp(-1:1,-1:1,-1:1)
    real(rt) :: bf1(0:5), bf2(0:5), h, h2, h3, bflo(0:5)
    real(rt) :: fracx, fracy, fracz, area, bc

    fac = sb/dx**2

    do cdir = 0, 5
       if (cdir .eq. 0 .or. cdir .eq. 3) then
          idim = 1
       else if (cdir .eq. 1 .or. cdir .eq. 4) then
          idim = 2
       else
          idim = 3
       end if
       h = dx(idim)
       if (bct(cdir) .eq. amrex_lo_dirichlet) then
          h2 = half * h
          bflo(cdir) = fac(idim) * ( h / (bcl(cdir) + h2) - one)
          if (bho.ge.1) then
             h3 = three * h2
             bf1(cdir) = fac(idim) * ((h3 - bcl(cdir)) / (bcl(cdir) + h2) - one)
             bf2(cdir) = fac(idim) * (bcl(cdir) - h2) / (bcl(cdir) + h3)
          else
             bf1(cdir) = bflo(cdir)
             bf2(cdir) = zero
          end if
       else if (bct(cdir) .eq. amrex_lo_neumann) then
          bflo(cdir) = -fac(idim)
          bf1(cdir) = -fac(idim)
          bf2(cdir) = zero
       end if
    end do

    irow = 0
    imat = 0
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (.not.is_covered_cell(flag(i,j,k))) then
                rows(irow) = cell_id(i,j,k)
                ncols(irow) = 0
                mat_tmp = zero

                if (is_regular_cell(flag(i,j,k))) then

                   mat_tmp(0,0,0) = sa*a(i,j,k) + fac(1)*(bx(i,j,k)+bx(i+1,j,k)) &
                        &                       + fac(2)*(by(i,j,k)+by(i,j+1,k)) &
                        &                       + fac(3)*(bz(i,j,k)+bz(i,j,k+1))
                   mat_tmp(-1, 0, 0) = -fac(1)*bx(i,j,k)
                   mat_tmp( 1, 0, 0) = -fac(1)*bx(i+1,j,k)
                   mat_tmp( 0,-1, 0) = -fac(2)*by(i,j,k)
                   mat_tmp( 0, 1, 0) = -fac(2)*by(i,j+1,k)
                   mat_tmp( 0, 0,-1) = -fac(3)*bz(i,j,k)
                   mat_tmp( 0, 0, 1) = -fac(3)*bz(i,j,k+1)

                   if (i.eq.lo(1) .and. cell_id(i-1,j,k).lt.0) then
                      cdir = 0
                      mat_tmp(0,0,0) = mat_tmp(0,0,0) + bf1(cdir)*bx(i,j,k)
                      mat_tmp(-1,0,0) = zero
                      mat_tmp(1,0,0) = mat_tmp(1,0,0) + bf2(cdir)*bx(i,j,k)
                   end if

                   if (i.eq.hi(1) .and. cell_id(i+1,j,k).lt.0) then
                      cdir = 3
                      mat_tmp(0,0,0) = mat_tmp(0,0,0) + bf1(cdir)*bx(i+1,j,k)
                      mat_tmp(1,0,0) = zero
                      mat_tmp(-1,0,0) = mat_tmp(-1,0,0) + bf2(cdir)*bx(i+1,j,k)
                   end if

                   if (j.eq.lo(2) .and. cell_id(i,j-1,k).lt.0) then
                      cdir = 1
                      mat_tmp(0,0,0) = mat_tmp(0,0,0) + bf1(cdir)*by(i,j,k)
                      mat_tmp(0,-1,0) = zero
                      mat_tmp(0,1,0) = mat_tmp(0,1,0) + bf2(cdir)*by(i,j,k)
                   end if

                   if (j.eq.hi(2) .and. cell_id(i,j+1,k).lt.0) then
                      cdir = 4
                      mat_tmp(0,0,0) = mat_tmp(0,0,0) + bf1(cdir)*by(i,j+1,k)
                      mat_tmp(0,1,0) = zero
                      mat_tmp(0,-1,0) = mat_tmp(0,-1,0) + bf2(cdir)*by(i,j+1,k)
                   end if

                   if (k.eq.lo(3) .and. cell_id(i,j,k-1).lt.0) then
                      cdir = 2
                      mat_tmp(0,0,0) = mat_tmp(0,0,0) + bf1(cdir)*bz(i,j,k)
                      mat_tmp(0,0,-1) = zero
                      mat_tmp(0,0,1) = mat_tmp(0,0,1) + bf2(cdir)*bz(i,j,k)
                   end if

                   if (k.eq.hi(3) .and. cell_id(i,j,k+1).lt.0) then
                      mat_tmp(0,0,0) = mat_tmp(0,0,0) + bf1(cdir)*bz(i,j,k+1)
                      mat_tmp(0,0,1) = zero
                      mat_tmp(0,0,-1) = mat_tmp(0,0,-1) + bf2(cdir)*bz(i,j,k+1)
                   end if
                   
                else
                end if
                
                do koff = -1, 1
                   do joff = -1, 1
                      do ioff = -1, 1
                         if (mat_tmp(ioff,joff,koff) .ne. zero) then
                            ncols(irow) = ncols(irow) + 1
                            cols(imat) = cell_id(i+ioff,j+joff,k+koff)
                            mat(imat) = mat_tmp(ioff,joff,koff)
                            imat = imat + 1
                         end if
                      end do
                end do
             end do
             irow = irow + 1
                
             end if
          end do
       end do
    end do
  end subroutine amrex_hpeb_ijmatrix
    
#endif

end module amrex_habec_module
