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
    integer, intent(in) :: lo(2), hi(2), alo(2), ahi(2)
    real(rt), intent(inout) :: mat(0:4, lo(1): hi(1), lo(2): hi(2))
    real(rt), intent(in   ) ::   a(    alo(1):ahi(1),alo(2):ahi(2))
    real(rt), intent(in) :: sa
    integer :: i, j
    if (sa .eq. zero) then
       !$omp parallel do private(i,j)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             mat(0,i,j) = zero
          enddo
       enddo
       !$omp end parallel do
    else
       !$omp parallel do private(i,j)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             mat(0,i,j) = sa * a(i,j)
          enddo
       enddo
       !$omp end parallel do
    end if
  end subroutine amrex_hpacoef


  subroutine amrex_hpbcoef (lo, hi, mat, b, blo, bhi, sb, dx, idim) &
       bind(c,name='amrex_hpbcoef')
    integer, intent(in) :: lo(2), hi(2), blo(2), bhi(2), idim
    real(rt), intent(inout) :: mat(0:4, lo(1): hi(1), lo(2): hi(2))
    real(rt), intent(in   ) ::   b(    blo(1):bhi(1),blo(2):bhi(2))
    real(rt), intent(in) :: sb, dx(2)

    integer :: i, j
    real(rt) :: fac

    fac = sb / dx(idim+1)**2 

    if (idim .eq. 0) then
       !$omp parallel do private(i,j)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             mat(0,i,j) = mat(0,i,j) + fac * (b(i,j) + b(i+1,j))
             mat(1,i,j) = - fac * b(i,j)
             mat(2,i,j) = - fac * b(i+1,j)
          enddo
       enddo
       !$omp end parallel do
    else
       !$omp parallel do private(i,j)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             mat(0,i,j) = mat(0,i,j) + fac * (b(i,j) + b(i,j+1))
             mat(3,i,j) = - fac * b(i,j)
             mat(4,i,j) = - fac * b(i,j+1)
          enddo
       enddo
       !$omp end parallel do
    endif
    
  end subroutine amrex_hpbcoef


  subroutine amrex_hpmat (lo, hi, mat, b, blo, bhi, mask, mlo, mhi, &
       sb, dx, cdir, bct, bcl, bho) bind(c,name='amrex_hpmat')
    integer, intent(in) :: lo(2), hi(2), blo(2), bhi(2), mlo(2), mhi(2), cdir, bct, bho
    real(rt), intent(inout) ::  mat(0:4, lo(1): hi(1), lo(2): hi(2))
    real(rt), intent(in   ) ::    b(    blo(1):bhi(1),blo(2):bhi(2))
    integer , intent(in   ) :: mask(    mlo(1):mhi(1),mlo(2):mhi(2))
    real(rt), intent(in) :: sb, dx(2), bcl

    integer :: i, j
    real(rt) :: fac, h, h2, h3, bf1, bf2

    if (cdir .eq. 0 .or. cdir .eq. 2) then
       h = dx(1)
    else
       h = dx(2)
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
       !$omp parallel do private(j)
       do j = lo(2), hi(2)
          if (mask(i-1,j) .gt. 0) then
             mat(0,i,j) = mat(0,i,j) + bf1 * b(i,j)
             mat(1,i,j) = zero
             mat(2,i,j) = mat(2,i,j) + bf2 * b(i,j)
          endif
       enddo
    else if (cdir .eq. 2) then
       i = hi(1)
       !$omp parallel do private(j)
       do j = lo(2), hi(2)
          if (mask(i+1,j) .gt. 0) then
             mat(0,i,j) = mat(0,i,j) + bf1 * b(i+1,j)
             mat(2,i,j) = zero
             mat(1,i,j) = mat(1,i,j) + bf2 * b(i+1,j)
          endif
       enddo
    else if (cdir .eq. 1) then
       j = lo(2)
       do i = lo(1), hi(1)
          if (mask(i,j-1) .gt. 0) then
             mat(0,i,j) = mat(0,i,j) + bf1 * b(i,j)
             mat(3,i,j) = zero
             mat(4,i,j) = mat(4,i,j) + bf2 * b(i,j)
          endif
       enddo
    else if (cdir .eq. 3) then
       j = hi(2)
       do i = lo(1), hi(1)
          if (mask(i,j+1) .gt. 0) then
             mat(0,i,j) = mat(0,i,j) + bf1 * b(i,j+1)
             mat(4,i,j) = zero
             mat(3,i,j) = mat(3,i,j) + bf2 * b(i,j+1)
          endif
       enddo
    else
       call amrex_error("hpmat: impossible face orientation")
    endif
    
  end subroutine amrex_hpmat


  subroutine amrex_hpijmatrix (lo, hi, nrows, ncols, rows, cols, mat, &
       cell_id, clo, chi, cell_id_begin, a , alo, ahi, bx, xlo, xhi, &
       by, ylo, yhi, sa, sb, dx, bct, bcl, bho) &
       bind(c,name='amrex_hpijmatrix')
    integer(hypre_int), intent(in) :: nrows, cell_id_begin;
    integer(hypre_int), dimension(0:nrows-1), intent(out) :: ncols, rows
    integer(hypre_int), dimension(0:nrows*5-1), intent(out) :: cols
    real(rt)          , dimension(0:nrows*5-1), intent(out) :: mat
    integer, dimension(2), intent(in) :: lo, hi, clo, chi, alo, ahi, xlo, xhi, ylo, yhi
    integer(hypre_int), intent(in) :: cell_id(clo(1):chi(1),clo(2):chi(2))
    real(rt)          , intent(in) :: a      (alo(1):ahi(1),alo(2):ahi(2))
    real(rt)          , intent(in) :: bx     (xlo(1):xhi(1),xlo(2):xhi(2))
    real(rt)          , intent(in) :: by     (ylo(1):yhi(1),ylo(2):yhi(2))
    integer, intent(in) :: bct(0:3), bho
    real(rt), intent(in) :: sa, sb, dx(2), bcl(0:3)

    integer :: i,j, irow, imat, ic, cdir, idim
    integer(hypre_int) :: cols_tmp(0:4)
    real(rt) :: fac(2), mat_tmp(0:4)
    real(rt) :: bf1(0:3), bf2(0:3), h, h2, h3

    fac = sb/dx**2

    do cdir = 0, 3
       if (cdir .eq. 0 .or. cdir .eq. 2) then
          idim = 1
       else
          idim = 2
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
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          rows(irow) = cell_id(i,j)
          ncols(irow) = 0
          
          cols_tmp(0) = cell_id(i,j)
          mat_tmp(0) = sa * a(i,j) + fac(1)*(bx(i,j)+bx(i+1,j)) &
               &                   + fac(2)*(by(i,j)+by(i,j+1))

          cols_tmp(1) = cell_id(i-1,j)
          mat_tmp(1) = -fac(1)*bx(i,j)

          cols_tmp(2) = cell_id(i+1,j)
          mat_tmp(2) = -fac(1)*bx(i+1,j)

          cols_tmp(3) = cell_id(i,j-1)
          mat_tmp(3) = -fac(2)*by(i,j)

          cols_tmp(4) = cell_id(i,j+1)
          mat_tmp(4) = -fac(2)*by(i,j+1)

          if (i.eq.lo(1) .and. cell_id(i-1,j).lt.0) then
             cdir = 0
             mat_tmp(0) = mat_tmp(0) + bf1(cdir)*bx(i,j)
             mat_tmp(2) = mat_tmp(2) + bf2(cdir)*bx(i,j)
          end if
          
          if (i.eq.hi(1) .and. cell_id(i+1,j).lt.0) then
             cdir = 2
             mat_tmp(0) = mat_tmp(0) + bf1(cdir)*bx(i+1,j)
             mat_tmp(1) = mat_tmp(1) + bf2(cdir)*bx(i+1,j)
          end if
          
          if (j.eq.lo(2) .and. cell_id(i,j-1).lt.0) then
             cdir = 1
             mat_tmp(0) = mat_tmp(0) + bf1(cdir)*by(i,j)
             mat_tmp(4) = mat_tmp(4) + bf2(cdir)*by(i,j)
          end if

          if (j.eq.hi(2) .and. cell_id(i,j+1).lt.0) then
             cdir = 3
             mat_tmp(0) = mat_tmp(0) + bf1(cdir)*by(i,j+1)
             mat_tmp(3) = mat_tmp(3) + bf2(cdir)*by(i,j+1)
          end if

          do ic = 0, 4
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
  end subroutine amrex_hpijmatrix


#ifdef AMREX_USE_EB

  subroutine amrex_hpeb_fill_cellid (lo, hi, nrows, cell_id, clo, chi, flag, flo, fhi) &
       bind(c,name='amrex_hpeb_fill_cellid')
    use amrex_ebcellflag_module, only : is_covered_cell
    integer, dimension(2), intent(in) :: lo, hi, clo, chi, flo, fhi
    integer(hypre_int), intent(  out) :: nrows
    integer(hypre_int), intent(inout) :: cell_id(clo(1):chi(1),clo(2):chi(2))
    integer           , intent(in   ) :: flag   (flo(1):fhi(1),flo(2):fhi(2))

    integer :: i,j

    nrows = 0
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (.not.is_covered_cell(flag(i,j))) then
             cell_id(i,j) = nrows
             nrows = nrows+1
          end if
       end do
    end do
  end subroutine amrex_hpeb_fill_cellid

  subroutine amrex_hpeb_copy_from_vec (lo, hi, a, alo, ahi, v, nv, flag, flo, fhi) &
       bind(c,name='amrex_hpeb_copy_from_vec')
    use amrex_ebcellflag_module, only : is_covered_cell
    integer, dimension(2), intent(in) :: lo, hi, alo, ahi, flo, fhi
    integer(hypre_int), intent(in) :: nv
    real(rt), intent(inout) :: a   (alo(1):ahi(1),alo(2):ahi(2))
    integer , intent(in   ) :: flag(flo(1):fhi(1),flo(2):fhi(2))
    real(rt), intent(in) :: v(0:nv-1)
    
    integer :: i,j, nrows

    nrows = 0
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (.not.is_covered_cell(flag(i,j))) then
             a(i,j) = v(nrows)
             nrows = nrows+1
          end if
       end do
    end do
  end subroutine amrex_hpeb_copy_from_vec

  subroutine amrex_hpeb_copy_to_vec (lo, hi, a, alo, ahi, v, nv, flag, flo, fhi) &
       bind(c,name='amrex_hpeb_copy_to_vec')
    use amrex_ebcellflag_module, only : is_covered_cell
    integer, dimension(2), intent(in) :: lo, hi, alo, ahi, flo, fhi
    integer(hypre_int), intent(in) :: nv
    real(rt), intent(in) :: a   (alo(1):ahi(1),alo(2):ahi(2))
    integer , intent(in) :: flag(flo(1):fhi(1),flo(2):fhi(2))
    real(rt), intent(out) :: v(0:nv-1)
    
    integer :: i,j, nrows

    nrows = 0
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (.not.is_covered_cell(flag(i,j))) then
             v(nrows) = a(i,j)
             nrows = nrows+1
          end if
       end do
    end do
  end subroutine amrex_hpeb_copy_to_vec

  subroutine amrex_hpeb_ijmatrix (lo, hi, nrows, ncols, rows, cols, mat, &
       cell_id, clo, chi, cell_id_begin, a, alo, ahi, bx, bxlo, bxhi, &
       by, bylo, byhi, flag, flo, fhi, vfrc, vlo, vhi, apx, axlo, axhi, apy, aylo, ayhi, &
       fcx, fxlo, fxhi, fcy, fylo, fyhi, &
       sa, sb, dx, bct, bcl, bho) &
       bind(c,name='amrex_hpeb_ijmatrix')
    use amrex_ebcellflag_module, only : is_covered_cell
    integer(hypre_int), intent(in) :: nrows, cell_id_begin;
    integer(hypre_int), dimension(0:nrows-1), intent(out) :: ncols, rows
    integer(hypre_int), dimension(0:nrows*9-1), intent(out) :: cols
    real(rt)          , dimension(0:nrows*9-1), intent(out) :: mat
    integer, dimension(2), intent(in) :: lo, hi, clo, chi, alo, ahi, bxlo, bxhi, bylo, byhi, &
         flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi, fxlo, fxhi, fylo, fyhi
    integer(hypre_int), intent(in) :: cell_id( clo(1): chi(1), clo(2): chi(2))
    real(rt)          , intent(in) :: a      ( alo(1): ahi(1), alo(2): ahi(2))
    real(rt)          , intent(in) :: bx     (bxlo(1):bxhi(1),bxlo(2):bxhi(2))
    real(rt)          , intent(in) :: by     (bylo(1):byhi(1),bylo(2):byhi(2))
    integer           , intent(in) :: flag   ( flo(1): fhi(1), flo(2): fhi(2))
    real(rt)          , intent(in) :: vfrc   ( vlo(1): vhi(1), vlo(2): vhi(2))
    real(rt)          , intent(in) :: apx    (axlo(1):axhi(1),axlo(2):axhi(2))
    real(rt)          , intent(in) :: apy    (aylo(1):ayhi(1),aylo(2):ayhi(2))
    real(rt)          , intent(in) :: fcx    (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    real(rt)          , intent(in) :: fcy    (fylo(1):fyhi(1),fylo(2):fyhi(2))
    integer, intent(in) :: bct(0:3), bho
    real(rt), intent(in) :: sa, sb, dx(2), bcl(0:3)

    integer :: i,j, irow, imat, ic, cdir, idim
    integer(hypre_int) :: cols_tmp(-1:1,-1:1)
    real(rt) :: fac(2), mat_tmp(-1:1,-1:1)
    real(rt) :: bf1(0:3), bf2(0:3), h, h2, h3

    fac = sb/dx**2

    do cdir = 0, 3
       if (cdir .eq. 0 .or. cdir .eq. 2) then
          idim = 1
       else
          idim = 2
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
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (.not.is_covered_cell(flag(i,j))) then
             rows(irow) = cell_id(i,j)
             ncols(irow) = 0
             cols_tmp = cell_id(i-1:i+1,j-1:j+1)
             mat_tmp = zero


          end if
       end do
    end do
  end subroutine amrex_hpeb_ijmatrix

#endif
  
end module amrex_habec_module
