module amrex_habec_module

  ! habec is Hypre abec, where abec is the form of the linear equation
  ! we are solving:
  ! 
  ! alpha*phi - div(beta*grad phi) + div(\vec{c}*phi) 

  use iso_c_binding
#ifdef AMREX_USE_PETSC
  use amrex_petsc_fort_module, only : it=>petsc_int
#else
  use amrex_hypre_fort_module, only : it=>hypre_int
#endif
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


  subroutine amrex_hpdiag (lo, hi, mat, diag, dlo, dhi) bind(c,name='amrex_hpdiag')
    integer, dimension(2), intent(in) :: lo, hi, dlo, dhi
    real(rt), intent(inout) :: mat (0:4, lo(1): hi(1), lo(2): hi(2))
    real(rt), intent(inout) :: diag(    dlo(1):dhi(1),dlo(2):dhi(2))
    integer :: i,j
    !$omp parallel do private(i,j)
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          diag(i,j) = one/mat(0,i,j)
          mat(0:4,i,j) = mat(0:4,i,j) * diag(i,j)
       enddo
    enddo
    !$omp end parallel do
  end subroutine amrex_hpdiag


  subroutine amrex_hpijmatrix (lo, hi, nrows, ncols, rows, cols, mat, &
       cell_id, clo, chi, cell_id_begin, diag, dlo, dhi, a, alo, ahi, &
       bx, xlo, xhi, by, ylo, yhi, sa, sb, dx, bct, bcl, bho) &
       bind(c,name='amrex_hpijmatrix')
    integer(it), intent(in) :: nrows, cell_id_begin;
    integer(it), dimension(0:nrows-1), intent(out) :: ncols, rows
    integer(it), dimension(0:nrows*5-1), intent(out) :: cols
    real(rt)   , dimension(0:nrows*5-1), intent(out) :: mat
    integer, dimension(2), intent(in) :: lo, hi, clo, chi, dlo, dhi, &
         alo, ahi, xlo, xhi, ylo, yhi
    integer(it), intent(in) :: cell_id(clo(1):chi(1),clo(2):chi(2))
    real(rt)   , intent(inout)::  diag(dlo(1):dhi(1),dlo(2):dhi(2))
    real(rt)   , intent(in) :: a      (alo(1):ahi(1),alo(2):ahi(2))
    real(rt)   , intent(in) :: bx     (xlo(1):xhi(1),xlo(2):xhi(2))
    real(rt)   , intent(in) :: by     (ylo(1):yhi(1),ylo(2):yhi(2))
    integer, intent(in) :: bct(0:3), bho
    real(rt), intent(in) :: sa, sb, dx(2), bcl(0:3)

    integer :: i,j, irow, imat, ic, cdir, idim
    integer(it) :: cols_tmp(0:4)
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

          diag(i,j) = one/mat_tmp(0)

          do ic = 0, 4
             if (cols_tmp(ic) .ge. 0) then
                ncols(irow) = ncols(irow) + 1
                cols(imat) = cols_tmp(ic)
                mat(imat) = mat_tmp(ic)*diag(i,j)
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
    integer(it), intent(  out) :: nrows
    integer(it), intent(inout) :: cell_id(clo(1):chi(1),clo(2):chi(2))
    integer    , intent(in   ) :: flag   (flo(1):fhi(1),flo(2):fhi(2))

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
    integer(it), intent(in) :: nv
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
    integer(it), intent(in) :: nv
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
       cell_id, clo, chi, cell_id_begin, diag, dlo, dhi, a, alo, ahi, &
       bx, bxlo, bxhi, by, bylo, byhi, flag, flo, fhi, vfrc, vlo, vhi, &
       apx, axlo, axhi, apy, aylo, ayhi, &
       fcx, fxlo, fxhi, fcy, fylo, fyhi, &
       ba, balo, bahi, bcen, bclo, bchi, beb, elo, ehi, is_eb_dirichlet, &
       sa, sb, dx, bct, bcl, bho) &
       bind(c,name='amrex_hpeb_ijmatrix')
    use amrex_ebcellflag_module, only : is_covered_cell, is_regular_cell
    use amrex_mlebabeclap_2d_module, only : amrex_get_dx_eb, blend_beta => amrex_blend_beta
    integer(it), intent(in) :: nrows, cell_id_begin
    integer(it), dimension(0:nrows-1), intent(out) :: ncols, rows
    integer(it), dimension(0:nrows*9-1), intent(out) :: cols
    real(rt)   , dimension(0:nrows*9-1), intent(out) :: mat
    integer, dimension(2), intent(in) :: lo, hi, clo, chi, dlo, dhi, &
         alo, ahi, bxlo, bxhi, bylo, byhi, &
         flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi, fxlo, fxhi, fylo, fyhi, &
         balo, bahi, bclo, bchi, elo, ehi
    integer(it), intent(in) :: cell_id( clo(1): chi(1), clo(2): chi(2))
    real(rt)   , intent(inout) :: diag( dlo(1): dhi(1), dlo(2): dhi(2))
    real(rt)   , intent(in) :: a      ( alo(1): ahi(1), alo(2): ahi(2))
    real(rt)   , intent(in) :: bx     (bxlo(1):bxhi(1),bxlo(2):bxhi(2))
    real(rt)   , intent(in) :: by     (bylo(1):byhi(1),bylo(2):byhi(2))
    integer    , intent(in) :: flag   ( flo(1): fhi(1), flo(2): fhi(2))
    real(rt)   , intent(in) :: vfrc   ( vlo(1): vhi(1), vlo(2): vhi(2))
    real(rt)   , intent(in) :: apx    (axlo(1):axhi(1),axlo(2):axhi(2))
    real(rt)   , intent(in) :: apy    (aylo(1):ayhi(1),aylo(2):ayhi(2))
    real(rt)   , intent(in) :: fcx    (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    real(rt)   , intent(in) :: fcy    (fylo(1):fyhi(1),fylo(2):fyhi(2))
    real(rt)   , intent(in) :: ba     (balo(1):bahi(1),balo(2):bahi(2))
    real(rt)   , intent(in) :: bcen   (bclo(1):bchi(1),bclo(2):bchi(2),2)
    real(rt)   , intent(in) :: beb    ( elo(1): ehi(1), elo(2): ehi(2))
    integer, intent(in) :: bct(0:3), bho
    real(rt), intent(in) :: sa, sb, dx(2), bcl(0:3)
    integer, intent(in) :: is_eb_dirichlet

    logical :: is_dirichlet
    integer :: i,j, irow, imat, cdir, idim, ii,jj, ioff, joff
    real(rt) :: fac(2), mat_tmp(-1:1,-1:1), phig1(4), phig2(4), feb(4)
    real(rt) :: bf1(0:3), bf2(0:3), h, h2, h3, bflo(0:3), c_0(-1:0,-1:0)
    real(rt) :: c_x(-1:0,-1:0), c_y(-1:0,-1:0), c_xy(-1:0,-1:0)
    real(rt) :: fracx, fracy, area, bc
    real(rt) :: gx, gy, anrmx, anrmy, anorm, anorminv, sx, sy
    real(rt) :: bctx, bcty, bsxinv, bsyinv
    real(rt) :: w1, w2, dg

    is_dirichlet = is_eb_dirichlet .ne. 0
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
          bflo(cdir) = fac(idim) * ( h / (bcl(cdir) + h2) - one) ! for low order
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
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (is_covered_cell(flag(i,j))) then
             diag(i,j) = zero
          else
             rows(irow) = cell_id(i,j)
             ncols(irow) = 0
             mat_tmp = zero

             if (is_regular_cell(flag(i,j))) then
                
                mat_tmp(0,0) = sa * a(i,j) + fac(1)*(bx(i,j)+bx(i+1,j)) &
                     &                     + fac(2)*(by(i,j)+by(i,j+1))
                mat_tmp(-1, 0) = -fac(1)*bx(i,j)
                mat_tmp( 1, 0) = -fac(1)*bx(i+1,j)
                mat_tmp( 0,-1) = -fac(2)*by(i,j)
                mat_tmp( 0, 1) = -fac(2)*by(i,j+1)
                
                if (i.eq.lo(1) .and. cell_id(i-1,j).lt.0) then
                   cdir = 0
                   mat_tmp(0,0) = mat_tmp(0,0) + bf1(cdir)*bx(i,j)
                   mat_tmp(-1,0) = zero
                   mat_tmp(1,0) = mat_tmp(1,0) + bf2(cdir)*bx(i,j)
                end if
                
                if (i.eq.hi(1) .and. cell_id(i+1,j).lt.0) then
                   cdir = 2
                   mat_tmp(0,0) = mat_tmp(0,0) + bf1(cdir)*bx(i+1,j)
                   mat_tmp(1,0) = zero
                   mat_tmp(-1,0) = mat_tmp(-1,0) + bf2(cdir)*bx(i+1,j)
                end if
                
                if (j.eq.lo(2) .and. cell_id(i,j-1).lt.0) then
                   cdir = 1
                   mat_tmp(0,0) = mat_tmp(0,0) + bf1(cdir)*by(i,j)
                   mat_tmp(0,-1) = zero
                   mat_tmp(0,1) = mat_tmp(0,1) + bf2(cdir)*by(i,j)
                end if
             
                if (j.eq.hi(2) .and. cell_id(i,j+1).lt.0) then
                   cdir = 3
                   mat_tmp(0,0) = mat_tmp(0,0) + bf1(cdir)*by(i,j+1)
                   mat_tmp(0,1) = zero
                   mat_tmp(0,-1) = mat_tmp(0,-1) + bf2(cdir)*by(i,j+1)
                end if

             else

                cdir = 0
                area = apx(i,j)
                bc = bx(i,j)
                if (area.ge.zero) then
                   if (area.ne.one) then
                      joff = int(sign(one,fcx(i,j)))
                      jj = j + joff
                      if (cell_id(i-1,jj).lt.0 .and. cell_id(i,jj).lt.0) then
                         fracy = zero
                      else
                         fracy = abs(fcx(i,j))
                      end if
                   else
                      joff = 0
                      jj = j
                      fracy = zero
                   end if

                   if (cell_id(i-1,j).ge.0) then
                      mat_tmp(0,0)  = mat_tmp(0,0)  + (one-fracy)*area*fac(1)*bc
                      mat_tmp(-1,0) = mat_tmp(-1,0) - (one-fracy)*area*fac(1)*bc
                   else if (cell_id(i+1,j).lt.0 .or. apx(i+1,j).eq.zero) then
                      mat_tmp(0,0)  = mat_tmp(0,0)  + (one-fracy)*area*(fac(1)+bflo(cdir))*bc
                   else
                      mat_tmp(0,0)  = mat_tmp(0,0)  + (one-fracy)*area*(fac(1)+bf1(cdir))*bc
                      mat_tmp(1,0)  = mat_tmp(1,0)  + (one-fracy)*area*        bf2(cdir) *bc
                   end if

                   if (fracy.gt.zero) then
                      if (cell_id(i-1,jj).ge.0 .and. cell_id(i,jj).ge.0) then
                         mat_tmp(-1,joff) = mat_tmp(-1,joff) - fracy*area*fac(1)*bx(i,jj)
                         mat_tmp( 0,joff) = mat_tmp( 0,joff) + fracy*area*fac(1)*bx(i,jj)
                      else if (cell_id(i+1,jj).lt.0 .or. apx(i+1,jj).eq.zero) then
                         mat_tmp(0,joff)  = mat_tmp(0,joff)  + (one-fracy)*area*(fac(1)+bflo(cdir))*bc
                      else
                         mat_tmp(0,joff)  = mat_tmp(0,joff)  + (one-fracy)*area*(fac(1)+bf1(cdir))*bc
                         mat_tmp(1,joff)  = mat_tmp(1,joff)  + (one-fracy)*area*        bf2(cdir) *bc
                      end if
                   end if
                end if

                cdir = 2
                area = apx(i+1,j)
                bc = bx(i+1,j)
                if (area.ge.zero) then
                   if (area.ne.one) then
                      joff = int(sign(one,fcx(i+1,j)))
                      jj = j + joff
                      if (cell_id(i,jj).lt.0 .and. cell_id(i+1,jj).lt.0) then
                         fracy = zero
                      else
                         fracy = abs(fcx(i+1,j))
                      end if
                   else
                      joff = 0
                      jj = j
                      fracy = zero
                   end if
                   
                   if (cell_id(i+1,j).ge.0) then
                      mat_tmp(0,0) = mat_tmp(0,0) + (one-fracy)*area*fac(1)*bc
                      mat_tmp(1,0) = mat_tmp(1,0) - (one-fracy)*area*fac(1)*bc
                   else if (cell_id(i-1,j).lt.0 .or. apx(i,j).eq.zero) then
                      mat_tmp(0,0) = mat_tmp(0,0) + (one-fracy)*area*(fac(1)+bflo(cdir))*bc
                   else
                      mat_tmp( 0,0) = mat_tmp( 0,0) + (one-fracy)*area*(fac(1)+bf1(cdir))*bc
                      mat_tmp(-1,0) = mat_tmp(-1,0) + (one-fracy)*area*        bf2(cdir) *bc
                   end if
                   
                   if (fracy.gt.zero) then
                      if (cell_id(i,jj).ge.0 .and. cell_id(i+1,jj).ge.0) then
                         mat_tmp(0,joff) = mat_tmp(0,joff) + fracy*area*fac(1)*bx(i+1,jj)
                         mat_tmp(1,joff) = mat_tmp(1,joff) - fracy*area*fac(1)*bx(i+1,jj)
                      else if (cell_id(i-1,jj).lt.0 .or. apx(i,jj).eq.zero) then
                         mat_tmp(0,joff) = mat_tmp(0,joff) + (one-fracy)*area*(fac(1)+bflo(cdir))*bc
                      else
                         mat_tmp( 0,joff) = mat_tmp( 0,joff) + (one-fracy)*area*(fac(1)+bf1(cdir))*bc
                         mat_tmp(-1,joff) = mat_tmp(-1,joff) + (one-fracy)*area*        bf2(cdir) *bc
                      end if
                   end if
                end if
                
                cdir = 1
                area = apy(i,j)
                bc = by(i,j)
                if (area.ge.zero) then
                   if (area.ne.one) then
                      ioff = int(sign(one,fcy(i,j)))
                      ii = i + ioff
                      if (cell_id(ii,j-1).lt.0 .and. cell_id(ii,j).lt.0) then
                         fracx = zero
                      else
                         fracx = abs(fcy(i,j))
                      end if
                   else
                      ioff = 0
                      ii = i
                      fracx = zero
                   end if

                   if (cell_id(i,j-1).ge.0) then
                      mat_tmp(0,0)  = mat_tmp(0,0)  + (one-fracx)*area*fac(2)*bc
                      mat_tmp(0,-1) = mat_tmp(0,-1) - (one-fracx)*area*fac(2)*bc
                   else if (cell_id(i,j+1).lt.0 .or. apy(i,j+1).eq.zero) then
                      mat_tmp(0,0)  = mat_tmp(0,0)  + (one-fracx)*area*(fac(2)+bflo(cdir))*bc
                   else
                      mat_tmp(0,0)  = mat_tmp(0,0)  + (one-fracx)*area*(fac(2)+bf1(cdir))*bc
                      mat_tmp(0,1)  = mat_tmp(0,1)  + (one-fracx)*area*        bf2(cdir) *bc
                   end if

                   if (fracx.gt.zero) then
                      if (cell_id(ii,j-1).ge.0 .and. cell_id(ii,j).ge.0) then
                         mat_tmp(ioff,-1) = mat_tmp(ioff,-1) - fracx*area*fac(2)*by(ii,j)
                         mat_tmp(ioff, 0) = mat_tmp(ioff, 0) + fracx*area*fac(2)*by(ii,j)
                      else if (cell_id(ii,j+1).lt.0 .or. apy(ii,j+1).eq.zero) then
                         mat_tmp(ioff,0)  = mat_tmp(ioff,0)  + (one-fracx)*area*(fac(2)+bflo(cdir))*bc
                      else
                         mat_tmp(ioff,0)  = mat_tmp(ioff,0)  + (one-fracx)*area*(fac(2)+bf1(cdir))*bc
                         mat_tmp(ioff,1)  = mat_tmp(ioff,1)  + (one-fracx)*area*        bf2(cdir) *bc
                      end if
                   end if
                end if

                cdir = 3
                area = apy(i,j+1)
                bc = by(i,j+1)
                if (area.ge.zero) then
                   if (area.ne.one) then
                      ioff = int(sign(one,fcy(i,j+1)))
                      ii = i + ioff
                      if (cell_id(ii,j).lt.0 .and. cell_id(ii,j+1).lt.0) then
                         fracx = zero
                      else
                         fracx = abs(fcy(i,j+1))
                      end if
                   else
                      ioff = 0
                      ii = i
                      fracx = zero
                   end if

                   if (cell_id(i,j+1).ge.0) then
                      mat_tmp(0,0) = mat_tmp(0,0) + (one-fracx)*area*fac(2)*bc
                      mat_tmp(0,1) = mat_tmp(0,1) - (one-fracx)*area*fac(2)*bc
                   else if (cell_id(i,j-1).lt.0 .or. apy(i,j).eq.zero) then
                      mat_tmp(0,0) = mat_tmp(0,0) + (one-fracx)*area*(fac(2)+bflo(cdir))*bc
                   else
                      mat_tmp(0, 0) = mat_tmp(0, 0) + (one-fracx)*area*(fac(2)+bf1(cdir))*bc
                      mat_tmp(0,-1) = mat_tmp(0,-1) + (one-fracx)*area*        bf2(cdir) *bc
                   end if

                   if (fracx.gt.zero) then
                      if (cell_id(ii,j).ge.0 .and. cell_id(ii,j+1).ge.0) then
                         mat_tmp(ioff,0) = mat_tmp(ioff,0) + fracx*area*fac(2)*bx(ii,j+1)
                         mat_tmp(ioff,1) = mat_tmp(ioff,1) - fracx*area*fac(2)*bx(ii,j+1)
                      else if (cell_id(ii,j-1).lt.0 .or. apy(ii,j).eq.zero) then
                         mat_tmp(ioff,0) = mat_tmp(ioff,0) + (one-fracx)*area*(fac(2)+bflo(cdir))*bc
                      else
                         mat_tmp(ioff, 0) = mat_tmp(ioff, 0) + (one-fracx)*area*(fac(2)+bf1(cdir))*bc
                         mat_tmp(ioff,-1) = mat_tmp(ioff,-1) + (one-fracx)*area*        bf2(cdir) *bc
                      end if
                   end if
                end if

                if (is_dirichlet) then
                   anorm = sqrt((apx(i,j) - apx(i+1,j))**2 + (apy(i,j) - apy(i,j+1))**2)
                   anorminv = one/anorm
                   anrmx = (apx(i,j) - apx(i+1,j))*anorminv
                   anrmy = (apy(i,j) - apy(i,j+1))*anorminv
                   bctx = bcen(i,j,1)
                   bcty = bcen(i,j,2)
                   if (abs(anrmx) .gt. abs(anrmy)) then
                      dg = amrex_get_dx_eb(vfrc(i,j)) / abs(anrmx)
                      gx = bctx - dg*anrmx
                      gy = bcty - dg*anrmy
                      sx = sign(one,anrmx)
                      sy = sign(one,anrmy)
                   else
                      dg = amrex_get_dx_eb(vfrc(i,j)) / abs(anrmy)
                      gx = bctx - dg*anrmx
                      gy = bcty - dg*anrmy
                      sx = sign(one,anrmx)
                      sy = sign(one,anrmy)
                   endif
                   ioff = -int(sx)
                   joff = -int(sy)

                   w1 = blend_beta(vfrc(i,j))
                   w2 = one - w1

                   if (w1.eq.zero) then 
                      phig1 = zero
                   else
                      phig1(1) = one + gx*sx + gy*sy + gx*gy*sx*sy
                      phig1(2) =     - gx*sx         - gx*gy*sx*sy
                      phig1(3) =             - gy*sy - gx*gy*sx*sy
                      phig1(4) =                     + gx*gy*sx*sy
                   endif

                   if (w2.eq.zero) then
                      phig2 = zero
                   else
                      bsxinv = one/(bctx+sx)
                      bsyinv = one/(bcty+sy)

                      ! c_0(0,0) = sx*sy*bsxinv*bsyinv
                      c_0(-1,0) = bctx*bsxinv
                      c_0(0,-1) = bcty*bsyinv
                      c_0(-1,-1) = -bctx*bcty*bsxinv*bsyinv

                      ! c_x(0,0) = sy*bsxinv*bsyinv
                      c_x(-1,0) = -bsxinv
                      c_x(0,-1) = sx*bcty*bsyinv
                      c_x(-1,-1) = -sx*bctx*bcty*bsxinv*bsyinv

                      ! c_y(0,0) = sx*bsxinv*bsyinv
                      c_y(-1,0) = sy*bctx*bsxinv
                      c_y(0,-1) = -bsyinv
                      c_y(-1,-1) = -sy*bctx*bcty*bsxinv*bsyinv

                      ! c_xy(0,0) = bsxinv*bsyinv
                      c_xy(-1,0) = -sy*bsxinv
                      c_xy(0,-1) = -sx*bsyinv
                      c_xy(-1,-1) = (one+sx*bctx+sy*bcty)*bsxinv*bsyinv

                      phig2(1) = zero
                      phig2(2) = (c_0(-1, 0) + gx*c_x(-1, 0) + gy*c_y(-1, 0) + gx*gy*c_xy(-1, 0))
                      phig2(3) = (c_0( 0,-1) + gx*c_x( 0,-1) + gy*c_y( 0,-1) + gx*gy*c_xy( 0,-1))
                      phig2(4) = (c_0(-1,-1) + gx*c_x(-1,-1) + gy*c_y(-1,-1) + gx*gy*c_xy(-1,-1))
                   endif

                   feb = -(w1*phig1 + w2*phig2) * (ba(i,j) * beb(i,j) / dg)
                   mat_tmp(0   , 0  ) = mat_tmp(0   , 0  ) - feb(1)*fac(1)
                   mat_tmp(ioff, 0  ) = mat_tmp(ioff, 0  ) - feb(2)*fac(1)
                   mat_tmp(0   ,joff) = mat_tmp(0   ,joff) - feb(3)*fac(1)
                   mat_tmp(ioff,joff) = mat_tmp(ioff,joff) - feb(4)*fac(1)
                endif

                mat_tmp = mat_tmp * (one/vfrc(i,j))
                mat_tmp(0,0) = mat_tmp(0,0) + sa*a(i,j)
             end if

             diag(i,j) = one/mat_tmp(0,0)

             do joff = -1, 1
                do ioff = -1, 1
                   if (mat_tmp(ioff,joff).ne.zero .and. cell_id(i+ioff,j+joff).ge.0) then
                      ncols(irow) = ncols(irow) + 1
                      cols(imat) = cell_id(i+ioff,j+joff)
                      mat(imat) = mat_tmp(ioff,joff)*diag(i,j)
                      imat = imat + 1
                   end if
                end do
             end do
             irow = irow + 1
             
          end if
       end do
    end do
  end subroutine amrex_hpeb_ijmatrix

#endif
  
end module amrex_habec_module
