module amrex_habec_module

  ! habec is Hypre abec, where abec is the form of the linear equation
  ! we are solving:
  ! 
  ! alpha*phi - div(beta*grad phi) + div(\vec{c}*phi) 

  use iso_c_binding
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
    
  

subroutine amrex_hmac(mat, a, &
                abox_l1,abox_l2,abox_l3,abox_h1,abox_h2,abox_h3, &
                reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
                alpha) bind(C, name="amrex_hmac")

  integer :: abox_l1,abox_l2,abox_l3,abox_h1,abox_h2,abox_h3
  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  real(rt)         :: a(abox_l1:abox_h1,abox_l2:abox_h2,abox_l3:abox_h3)
  real(rt)         :: mat(0:6, reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)
  real(rt)         :: alpha
  integer :: i, j, k

  if (alpha == 0.e0_rt) then
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = 0.e0_rt
           enddo
        enddo
     enddo
  else
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              mat(0,i,j,k) = alpha * a(i,j,k)
           enddo
        enddo
     enddo
  endif
  
end subroutine amrex_hmac


subroutine amrex_hmac_ij(a,abox_l1,abox_l2,abox_l3,abox_h1,abox_h2,abox_h3, &
     reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3,alpha, &
     Amat,Index, Gbox_l1,Gbox_l2,Gbox_l3,Gbox_h1,Gbox_h2,Gbox_h3) bind(C, name="amrex_hmac_ij")

  integer :: abox_l1,abox_l2,abox_l3,abox_h1,abox_h2,abox_h3
  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  real(rt)         :: a(abox_l1:abox_h1,abox_l2:abox_h2,abox_l3:abox_h3)
  real(rt)         :: alpha
  integer :: Gbox_l1,Gbox_l2,Gbox_l3,Gbox_h1,Gbox_h2,Gbox_h3
  integer :: Index(Gbox_l1:Gbox_h1,Gbox_l2:Gbox_h2,Gbox_l3:Gbox_h3)
  ! HYPRE objects
  integer(c_long) :: Amat
  integer :: i, j, k, ierr, count
  integer :: nrows
  integer,  dimension(:), pointer :: cols
  integer,  dimension(:), pointer :: rows
  integer,  dimension(:), pointer :: ncols
  real(rt), dimension(:), pointer :: values
  
  nrows = (reg_h3-reg_l3+1) * (reg_h2-reg_l2+1) * (reg_h1-reg_l1+1)
  
  allocate(rows   (nrows))
  allocate(cols   (nrows))
  allocate(values (nrows))
  allocate(ncols  (nrows))

  ncols(:) = 0
  count    = 0
  
  if (alpha == 0.e0_rt) then

     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              count = count + 1
              rows(count) = Index(i,j,k)
              
              ncols(count) = 1
              cols(count) = Index(i,j,k)
              values(count) = 0.e0_rt
           end do
        end do
     end do
     
  else

     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1
              count = count + 1
              rows(count) = Index(i,j,k)

              ncols(count) = 1
              cols(count) = Index(i,j,k)
              values(count) = alpha * a(i,j,k)
           end do
        end do
     end do
     
  endif

  call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)

  deallocate(rows)
  deallocate(cols)
  deallocate(values)
  deallocate(ncols)
 
end subroutine amrex_hmac_ij


subroutine amrexhmbc_ij(b, bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3, reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
                   beta, dx, n, &
                   Amat, Index, Gbox_l1,Gbox_l2,Gbox_l3,Gbox_h1,Gbox_h2,Gbox_h3) bind(C, name="amrex_hmbc_ij")

  integer :: bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3
  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  integer :: n
  real(rt)  :: b(bbox_l1:bbox_h1,bbox_l2:bbox_h2,bbox_l3:bbox_h3)
  real(rt)  :: beta, dx(3)
  integer :: Gbox_l1,Gbox_l2,Gbox_l3,Gbox_h1,Gbox_h2,Gbox_h3
  integer :: Index(Gbox_l1:Gbox_h1,Gbox_l2:Gbox_h2,Gbox_l3:Gbox_h3)
  ! HYPRE objects
  integer(c_long) :: Amat

  integer :: nrows, ncols, rows
  integer, dimension(3) :: cols
  real(rt), dimension(3) :: values
  real(rt) :: fac
  integer  :: i, j, k, ierr

  ncols = 3
  nrows = 1

  if (n == 0) then

     fac = beta / (dx(1)**2)
     
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           do i = reg_l1+1, reg_h1-1
              
              rows      = Index(i,j,k)
              
              cols(1)   = Index(i,j,k)
              values(1) = fac * (b(i,j,k) + b(i+1,j,k))
              
              cols(2)   = Index(i-1,j,k)
              values(2) = -fac * b(i,j,k)
              
              cols(3)   = Index(i+1,j,k)
              values(3) = -fac * b(i+1,j,k)
              
              call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)
              
           end do
        end do
     end do
     
  elseif (n == 1) then

     fac = beta / (dx(2)**2)
     do k = reg_l3, reg_h3
        do j = reg_l2+1, reg_h2-1
           do i = reg_l1, reg_h1
              
              rows = Index(i,j,k)

              cols(1)   = Index(i,j,k)
              values(1) = fac * (b(i,j,k) + b(i,j+1,k))

              cols(2)   = Index(i,j-1,k)
              values(2) = -fac * b(i,j,k)

              cols(3) = Index(i,j+1,k)
              values(3) = - fac * b(i,j+1,k)

              call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)
              
           end do
        end do
     end do
     
  else

     fac = beta / (dx(3)**2)

     do k = reg_l3+1, reg_h3-1
        do j = reg_l2, reg_h2
           do i = reg_l1, reg_h1

              rows = Index(i,j,k)

              cols(1) = Index(i,j,k)
              values(1) = fac * (b(i,j,k) + b(i,j,k+1))

              cols(2) = Index(i,j,k-1)
              values(2) = - fac * b(i,j,k)

              cols(3) = Index(i,j,k+1)
              values(3) = - fac * b(i,j,k+1)

              call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)
              
           end do
        end do
     end do
     
  endif
  
end subroutine amrexhmbc_ij

subroutine amrex_hmmat(mat, &
                 reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
                 cdir, bct, bho, bcl, &
                 mask, msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3, &
                 b, bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3, &
                 beta, dx) bind(C, name="amrex_hmmat")

  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  integer :: msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3
  integer :: bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3
  integer :: cdir, bct, bho
  real(rt)         :: bcl, beta, dx(3)
  real(rt)         :: mat(0:6, reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)
  integer :: mask(msk_l1:msk_h1,msk_l2:msk_h2,msk_l3:msk_h3)
  real(rt)         :: b(bbox_l1:bbox_h1,bbox_l2:bbox_h2,bbox_l3:bbox_h3)
  real(rt)         :: h, fac, bfm, bfv
  real(rt)         :: bfm2, h2, th2
  integer :: i, j, k

  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif

  fac = beta / (h**2)

  if (bct == amrex_lo_dirichlet) then
     if (bho >= 1) then
        h2 = 0.5e0_rt * h
        th2 = 3.e0_rt * h2
        bfm = fac * (th2 - bcl) / (bcl + h2) - fac
        bfm2 = fac * (bcl - h2) / (bcl + th2)
     else
        bfv = (beta / h) / (0.5e0_rt * h + bcl)
        bfm = bfv - fac
     endif
  else if (bct == amrex_lo_neumann) then
     bfm = -fac
     bfm2 = 0.e0_rt
  else
     print *, "hmmat: unsupported boundary type"
     stop
  endif
  
  if (cdir == 0) then
     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i-1,j,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
              mat(1,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(2,i,j,k) = mat(2,i,j,k) + bfm2 * b(i,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i+1,j,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i+1,j,k)
              mat(2,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(1,i,j,k) = mat(1,i,j,k) + bfm2 * b(i+1,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j-1,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
              mat(3,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(4,i,j,k) = mat(4,i,j,k) + bfm2 * b(i,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j+1,k) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j+1,k)
              mat(4,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(3,i,j,k) = mat(3,i,j,k) + bfm2 * b(i,j+1,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k-1) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k)
              mat(5,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(6,i,j,k) = mat(6,i,j,k) + bfm2 * b(i,j,k)
              endif
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k+1) > 0) then
              mat(0,i,j,k) = mat(0,i,j,k) + bfm * b(i,j,k+1)
              mat(6,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(5,i,j,k) = mat(5,i,j,k) + bfm2 * b(i,j,k+1)
              endif
           endif
        enddo
     enddo
  else
     print *, "hmmat: impossible face orientation"
  endif

end subroutine amrex_hmmat

subroutine amrex_hmmat_ij(reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
                 cdir, bct, bho, bcl, &
                 mask, msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3, &
                 b, bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3, &
                 beta, dx, Amat, Index, Gbox_l1,Gbox_l2,Gbox_l3,Gbox_h1,Gbox_h2,Gbox_h3) &
                 bind(C, name="amrex_hmmat_ij")

  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  integer :: msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3
  integer :: bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3
  integer :: cdir, bct, bho
  real(rt) :: bcl, beta, dx(3)
  integer :: mask(msk_l1:msk_h1,msk_l2:msk_h2,msk_l3:msk_h3)
  real(rt) :: b(bbox_l1:bbox_h1,bbox_l2:bbox_h2,bbox_l3:bbox_h3)
  integer :: Gbox_l1,Gbox_l2,Gbox_l3,Gbox_h1,Gbox_h2,Gbox_h3
  integer :: Index(Gbox_l1:Gbox_h1,Gbox_l2:Gbox_h2,Gbox_l3:Gbox_h3)
  ! HYPRE objects
  integer(c_long) :: Amat
  real(rt)         :: h, fac, bfm, bfv
  real(rt)         :: bfm2, h2, th2
  integer :: i, j, k, ierr

  integer, dimension(3)  :: cols
  real(rt), dimension(3) :: values
  integer :: ncols
  integer :: rows, nrows

  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  
  fac = beta / (h**2)
  if (bct == amrex_lo_dirichlet) then
     if (bho >= 1) then
        h2 = 0.5e0_rt * h
        th2 = 3.e0_rt * h2
        bfm = fac * (th2 - bcl) / (bcl + h2) - fac
        bfm2 = fac * (bcl - h2) / (bcl + th2)
     else
        bfv = (beta / h) / (0.5e0_rt * h + bcl)
        bfm = bfv - fac
     endif
  else if (bct == amrex_lo_neumann) then
     bfm = -fac
     bfm2 = 0.e0_rt
  else
     print *, "hmmat: unsupported boundary type"
     stop
  endif

  nrows = 1
  values(:) = 0.0
  
  if (cdir == 0) then

     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2

           rows = Index(i,j,k)
           ncols = 2
           
           cols(1) = Index(i,j,k)
           values(1) = fac*(b(i,j,k) + b(i+1,j,k)) 
           
           cols(2) = Index(i+1,j,k)
           values(2) = -fac*b(i+1,j,k)
           
           if (mask(i-1,j,k) > 0) then
              values(1) = values(1) + bfm * b(i,j,k)
              if (bho >= 1) then
                 values(2) = values(2) + bfm2 * b(i,j,k)
              endif
           else
              ncols = 3
              cols(3) = Index(i-1,j,k)
              values(3) = -fac*b(i,j,k)
           endif

           call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)
           
        enddo
     enddo

  else if (cdir == 3) then

     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2

           rows = Index(i,j,k)
           ncols = 2
           
           cols(1) = Index(i,j,k)
           values(1) = fac*(b(i,j,k) + b(i+1,j,k))
           
           cols(2) = Index(i-1,j,k)
           values(2) = -fac*b(i,j,k)

           if (mask(i+1,j,k) > 0) then
              values(1) = values(2) + bfm * b(i+1,j,k)
              if (bho >= 1) then
                 values(2) = values(2) + bfm2 * b(i+1,j,k)
              endif
           else
              ncols = 3
              cols(3) = Index(i+1,j,k)
              values(3) = -fac*b(i+1,j,k)
           end if

           call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)

        end do
     end do
   
  else if (cdir == 1) then

     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1

           rows = Index(i,j,k)
           ncols = 2
           
           cols(1) = Index(i,j,k)
           values(1) = fac * (b(i,j,k) + b(i,j+1,k))

           cols(2) = Index(i,j+1,k)
           values(2) = - fac * b(i,j+1,k)

           if (mask(i,j-1,k) > 0) then
              values(1) = values(1) + bfm * b(i,j,k)
              if (bho >= 1) then
                 values(2) = values(2) + bfm2 * b(i,j,k)
              endif
           else
              ncols = 3
              cols(3) = Index(i,j-1,k)
              values(3) = -fac*b(i,j,k)
           end if

           call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)
           
        end do
     end do
     
  else if (cdir == 4) then

     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1

           rows = Index(i,j,k)
           ncols = 2
           
           cols(1) = Index(i,j,k)
           values(1) = fac * (b(i,j,k) + b(i,j+1,k))

           cols(2) = Index(i,j-1,k)
           values(2) = - fac * b(i,j,k)

           if (mask(i,j+1,k) > 0) then
              values(1) = values(1) + bfm * b(i,j+1,k)
              if (bho >= 1) then
                 values(2) = values(2) + bfm2 * b(i,j+1,k)
              endif
           else
              ncols = 3
              cols(3) = Index(i,j+1,k)
              values(3) = -fac*b(i,j+1,k)
           endif
           
           call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)           

        end do
     end do
     
  else if (cdir == 2) then

     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1

           rows = Index(i,j,k)
           ncols = 2
           
           cols(1) = Index(i,j,k)
           values(1) = fac * (b(i,j,k) + b(i,j,k+1))
           
           cols(2) = Index(i,j,k+1)
           values(2) = - fac * b(i,j,k+1)
           
           if (mask(i,j,k-1) > 0) then
              values(1) = values(1) + bfm * b(i,j,k)
              if (bho >= 1) then
                 values(2) = values(2) + bfm2 * b(i,j,k)
              end if
           else
              ncols = 3
              cols(3) = Index(i,j,k-1)
              values(3) = -fac*b(i,j,k)
           end if

           call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)
           
        end do
     end do
     
  else if (cdir == 5) then

     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1

           rows = Index(i,j,k)
           ncols = 2
           
           cols(1) = Index(i,j,k)
           values(1) = fac * (b(i,j,k) + b(i,j,k+1))

           cols(2) = Index(i,j,k-1)
           values(2) = - fac * b(i,j,k)

           if (mask(i,j,k+1) > 0) then
              values(1) = values(1) + bfm * b(i,j,k+1)
              if (bho >= 1) then
                 values(2) = values(2) + bfm2 * b(i,j,k+1)
              end if
           else
              ncols = 3
              cols(3) = Index(i,j,k+1)
              values(3) = -fac*b(i,j,k+1)
           end if

           call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)
           
        end do
     end do

  else

     print *, "hmmat ij: impossible face orientation"

  endif
  
end subroutine amrex_hmmat_ij


subroutine amrex_hmmat3(mat, &
                  reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
                  cdir, bctype, bho, bcl, &
                  mask, msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3, &
                  b, bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3, &
                  beta, dx) bind(C, name="amrex_hmmat3")

  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  integer :: msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3
  integer :: bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3
  integer :: cdir, bctype, bho
  real(rt)         :: bcl, beta, dx(3)
  real(rt)         :: mat(0:6, reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)
  integer :: mask(msk_l1:msk_h1,msk_l2:msk_h2,msk_l3:msk_h3)
  real(rt)         :: b(bbox_l1:bbox_h1,bbox_l2:bbox_h2,bbox_l3:bbox_h3)
  real(rt)         :: h, fac, bfm, bfv
  real(rt)         :: bfm2, h2, th2
  integer :: i, j, k, bct
  ! The -fac * b(i,j,k) term applied to the matrix diagonal is the contribution
  ! from the interior stencil which must be removed at the boundary.
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  bct = bctype
  fac = beta / (h**2)
  if (cdir == 0) then
     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i-1,j,k) > 0) then
              if (bct == amrex_lo_dirichlet) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == amrex_lo_neumann) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
              mat(1,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(2,i,j,k) = mat(2,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 3) then
     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2
           if (mask(i+1,j,k) > 0) then
              if (bct == amrex_lo_dirichlet) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i+1,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i+1,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i+1,j,k)
                 endif
              else if (bct == amrex_lo_neumann) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i+1,j,k)
              mat(2,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(1,i,j,k) = mat(1,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 1) then
     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j-1,k) > 0) then
              if (bct == amrex_lo_dirichlet) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == amrex_lo_neumann) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
              mat(3,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(4,i,j,k) = mat(4,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 4) then
     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           if (mask(i,j+1,k) > 0) then
              if (bct == amrex_lo_dirichlet) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j+1,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j+1,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j+1,k)
                 endif
              else if (bct == amrex_lo_neumann) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j+1,k)
              mat(4,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(3,i,j,k) = mat(3,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 2) then
     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k-1) > 0) then
              if (bct == amrex_lo_dirichlet) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == amrex_lo_neumann) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k)
              mat(5,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(6,i,j,k) = mat(6,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else if (cdir == 5) then
     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           if (mask(i,j,k+1) > 0) then
              if (bct == amrex_lo_dirichlet) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k+1)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k+1)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k+1)
                 endif
              else if (bct == amrex_lo_neumann) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3: unsupported boundary type"
                 stop
              endif
              mat(0,i,j,k) = mat(0,i,j,k) + bfm - fac * b(i,j,k+1)
              mat(6,i,j,k) = 0.e0_rt
              if (bho >= 1) then
                 mat(5,i,j,k) = mat(5,i,j,k) + bfm2
              endif
           endif
        enddo
     enddo
  else
     print *, "hmmat3: impossible face orientation"

  endif
end subroutine amrex_hmmat3

subroutine amrex_hmmat3_ij(reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
                  cdir, bctype, bho, bcl, &
                  mask, msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3, &
                  b, bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3, &
                  beta, dx,      &
                  Amat, Index, Gbox_l1,Gbox_l2,Gbox_l3,Gbox_h1,Gbox_h2,Gbox_h3) &
                  bind(C, name="amrex_hmmat3_ij")
  
  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  integer :: msk_l1,msk_l2,msk_l3,msk_h1,msk_h2,msk_h3
  integer :: bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3
  integer  :: cdir, bctype, bho
  real(rt) :: bcl, beta, dx(3)
  integer  :: mask(msk_l1:msk_h1,msk_l2:msk_h2,msk_l3:msk_h3)
  real(rt) :: b(bbox_l1:bbox_h1,bbox_l2:bbox_h2,bbox_l3:bbox_h3)
  integer  :: Gbox_l1,Gbox_l2,Gbox_l3,Gbox_h1,Gbox_h2,Gbox_h3
  integer  :: Index(Gbox_l1:Gbox_h1,Gbox_l2:Gbox_h2,Gbox_l3:Gbox_h3)
  ! HYPRE objects
  integer(c_long) :: Amat
  real(rt)  :: h, fac, bfm, bfv
  real(rt)  :: bfm2, h2, th2
  integer :: i, j, k, bct, ierr
  integer,  dimension(3) :: cols
  real(rt), dimension(3) :: values
  integer :: ncols 
  integer :: rows, nrows

  nrows = 1
    
  ! The -fac * b(i,j,k) term applied to the matrix diagonal is the contribution
  ! from the interior stencil which must be removed at the boundary.
  if (cdir == 0 .OR. cdir == 3) then
     h = dx(1)
  elseif (cdir == 1 .OR. cdir == 4) then
     h = dx(2)
  else
     h = dx(3)
  endif
  
  bct = bctype
  fac = beta / (h**2)

  if (cdir == 0) then

     i = reg_l1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2

           rows = Index(i,j,k)
           ncols = 2

           cols(1) = Index(i,j,k)
           values(1) = fac*(b(i,j,k) + b(i+1,j,k)) 

           cols(2) = Index(i+1,j,k)
           values(2) = -fac*b(i+1,j,k)
           
           if (mask(i-1,j,k) > 0) then

              if (bct == amrex_lo_dirichlet) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == amrex_lo_neumann) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3 ij: unsupported boundary type"
                 stop
              endif

              values(1) = values(1) + bfm - fac * b(i,j,k)

              if (bho >= 1) then
                 values(2) = values(2) + bfm2
              endif
           else
              ncols = 3
              cols(3) = Index(i-1,j,k)
              values(3) = -fac*b(i,j,k)
           end if
           
           call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)
           
        end do
     end do
     
  else if (cdir == 3) then

     i = reg_h1
     do k = reg_l3, reg_h3
        do j = reg_l2, reg_h2

           rows = Index(i,j,k)
           ncols = 2
           
           cols(1) = Index(i,j,k)
           values(1) = fac*(b(i,j,k) + b(i+1,j,k))
           
           cols(2) = Index(i-1,j,k)
           values(2) = -fac*b(i,j,k)

           if (mask(i+1,j,k) > 0) then
              if (bct == amrex_lo_dirichlet) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i+1,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i+1,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i+1,j,k)
                 endif
              else if (bct == amrex_lo_neumann) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3 ij: unsupported boundary type"
                 stop
              endif
              values(1) = values(1) + bfm - fac * b(i+1,j,k)
              if (bho >= 1) then
                 values(2) = values(2) + bfm2
              endif
           else
              ncols = 3
              cols(3) = Index(i+1,j,k)
              values(3) = -fac*b(i+1,j,k)
           endif

           call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)
        end do
     end do
     
  else if (cdir == 1) then
     
     j = reg_l2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1
           
           rows = Index(i,j,k)
           ncols = 2
           
           cols(1) = Index(i,j,k)
           values(1) = fac * (b(i,j,k) + b(i,j+1,k))

           cols(2) = Index(i,j+1,k)
           values(2) = - fac * b(i,j+1,k)
           
           if (mask(i,j-1,k) > 0) then
              if (bct == amrex_lo_dirichlet) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == amrex_lo_neumann) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3 ij: unsupported boundary type"
                 stop
              endif

              values(1) = values(1) + bfm - fac * b(i,j,k)
              if (bho >= 1) then
                 values(2) = values(2) + bfm2
              end if
           else
              ncols = 3
              cols(3) = Index(i,j-1,k)
              values(3) = -fac*b(i,j,k)
           endif

           call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)

        end do
     end do
     
  else if (cdir == 4) then

     j = reg_h2
     do k = reg_l3, reg_h3
        do i = reg_l1, reg_h1

           rows = Index(i,j,k)
           ncols = 2
           
           cols(1) = Index(i,j,k)
           values(1) = fac * (b(i,j,k) + b(i,j+1,k))

           cols(2) = Index(i,j-1,k)
           values(2) = - fac * b(i,j,k)
           
           if (mask(i,j+1,k) > 0) then
              if (bct == amrex_lo_dirichlet) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j+1,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j+1,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j+1,k)
                 endif
              else if (bct == amrex_lo_neumann) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3 ij: unsupported boundary type"
                 stop
              endif
              values(1) = values(1) + bfm - fac * b(i,j+1,k)
              if (bho >= 1) then
                 values(2) = values(2) + bfm2
              endif
           else
              ncols = 3
              cols(3) = Index(i,j+1,k)
              values(3) = -fac*b(i,j+1,k)
           end if

           call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)
           
        end do
     end do
     
  else if (cdir == 2) then

     k = reg_l3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1

           rows = Index(i,j,k)
           ncols = 2
           
           cols(1) = Index(i,j,k)
           values(1) = fac * (b(i,j,k) + b(i,j,k+1))

           cols(2) = Index(i,j,k+1)
           values(2) = - fac * b(i,j,k+1)

           if (mask(i,j,k-1) > 0) then
              if (bct == amrex_lo_dirichlet) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k)
                 endif
              else if (bct == amrex_lo_neumann) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3 ij: unsupported boundary type"
                 stop
              endif
              values(1) = values(1) + bfm - fac * b(i,j,k)
              if (bho >= 1) then
                 values(2) = values(2) + bfm2
              endif
           else
              ncols = 3
              cols(3) = Index(i,j,k-1)
              values(3) = -fac*b(i,j,k)
           end if
           
           call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)
           
        end do
     end do
     
  else if (cdir == 5) then
     k = reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1

           rows = Index(i,j,k)
           ncols = 2
           
           cols(1) = Index(i,j,k)
           values(1) = fac * (b(i,j,k) + b(i,j,k+1))

           cols(2) = Index(i,j,k-1)
           values(2) = - fac * b(i,j,k)
           
           if (mask(i,j,k+1) > 0) then
              if (bct == amrex_lo_dirichlet) then
                 if (bho >= 1) then
                    h2 = 0.5e0_rt * h
                    th2 = 3.e0_rt * h2
                    bfm = fac * (th2 - bcl) / (bcl + h2)  * b(i,j,k+1)
                    bfm2 = fac * (bcl - h2) / (bcl + th2) * b(i,j,k+1)
                 else
                    bfv = (beta / h) / (0.5e0_rt * h + bcl)
                    bfm = bfv * b(i,j,k+1)
                 endif
              else if (bct == amrex_lo_neumann) then
                 bfm  = 0.e0_rt
                 bfm2 = 0.e0_rt
              else
                 print *, "hmmat3 ij: unsupported boundary type"
                 stop
              endif
              values(1) = values(1) + bfm - fac * b(i,j,k+1)
              if (bho >= 1) then
                 values(2) = values(2) + bfm2
              end if
              
           else
              ncols = 3
              cols(3) = Index(i,j,k+1)
              values(3) = -fac*b(i,j,k+1)
           end if
           
           call HYPRE_IJMatrixAddToValues(Amat, nrows, ncols, rows, cols, values, ierr)
           
        end do
     end do
  else
     print *, "hmmat3 ij: impossible face orientation"
  endif
  
end subroutine amrex_hmmat3_ij

subroutine amrex_buildglobalindex(Index, &
     bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3, &
     reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, &
     GIndex) bind(C, name="amrex_BuildGlobalIndex")

  integer :: bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3
  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  integer :: Index(bbox_l1:bbox_h1,bbox_l2:bbox_h2,bbox_l3:bbox_h3)
  integer :: GIndex
  integer :: i, j, k,count

  count = 0
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           Index(i,j,k) = Index(i,j,k) + GIndex + count
           count = count + 1
        end do
     end do
  end do

end subroutine amrex_buildglobalindex


subroutine amrex_conv_Vec_Local_Global(X, vec, nRows, &
     reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3, Index, bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3) &
     bind(C, name="amrex_conv_Vec_Local_Global")

  integer :: bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3
  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  integer :: Index(bbox_l1:bbox_h1,bbox_l2:bbox_h2,bbox_l3:bbox_h3)
  integer :: nRows
  real(rt) :: vec(reg_l1:reg_h1,reg_l2:reg_h2,reg_l3:reg_h3)
  integer(c_long) :: X
  real(rt) :: VecGB(0:(nRows-1))
  integer  :: Indices(0:(nRows-1))
  integer :: i, j, k, count, ierr

  count = 0

  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1

           VecGB(count) = vec(i,j,k)
           Indices(count) = Index(i,j,k)

           count = count + 1

        end do
     end do
  end do

  call HYPRE_IJVectorSetValues(X,nRows,Indices,VecGB, ierr)

end subroutine amrex_conv_Vec_Local_Global


subroutine amrex_conv_Vec_Global_Local(vec, bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3, VecGB, nRows, &
     reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3) bind(C, name="amrex_conv_Vec_Global_Local")

  integer :: reg_l1,reg_l2,reg_l3,reg_h1,reg_h2,reg_h3
  integer :: nRows
  real(rt) :: VecGB(0:(nRows-1))
  integer :: bbox_l1,bbox_l2,bbox_l3,bbox_h1,bbox_h2,bbox_h3
  real(rt) :: vec(bbox_l1:bbox_h1,bbox_l2:bbox_h2,bbox_l3:bbox_h3)
  integer :: i, j, k, count

  count = 0
  do k = reg_l3, reg_h3
     do j = reg_l2, reg_h2
        do i = reg_l1, reg_h1
           vec(i,j,k) = VecGB(count)
           count = count + 1
        end do
     end do
  end do

end subroutine amrex_conv_Vec_Global_Local


end module amrex_habec_module
