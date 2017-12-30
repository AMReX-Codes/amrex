module amrex_mlnodelap_3d_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use amrex_lo_bctypes_module, only : amrex_lo_dirichlet, amrex_lo_neumann, amrex_lo_inflow, amrex_lo_periodic
  implicit none

  ! external dirichlet at physical boundary or internal dirichlet at crse/fine boundary
  integer, parameter :: dirichlet = 1

  integer, parameter :: crse_cell = 0
  integer, parameter :: fine_cell = 1
  integer, parameter :: crse_node = 0
  integer, parameter :: crse_fine_node = 1
  integer, parameter :: fine_node = 2

  private
  public :: &
       ! masks
       amrex_mlndlap_set_nodal_mask, amrex_mlndlap_set_dirichlet_mask, &
       amrex_mlndlap_fixup_res_mask, amrex_mlndlap_any_fine_sync_cells, &
       ! coeffs
       amrex_mlndlap_avgdown_coeff, amrex_mlndlap_fillbc_cc, &
       ! bc
       amrex_mlndlap_applybc, &
       ! operator
       amrex_mlndlap_adotx_ha, amrex_mlndlap_adotx_aa, &
       amrex_mlndlap_jacobi_ha, amrex_mlndlap_jacobi_aa, &
       amrex_mlndlap_gauss_seidel_ha, amrex_mlndlap_gauss_seidel_aa, &
       ! restriction
       amrex_mlndlap_restriction, &
       ! interpolation
       amrex_mlndlap_interpolation_ha, amrex_mlndlap_interpolation_aa, &
       ! rhs & u
       amrex_mlndlap_divu, amrex_mlndlap_add_rhcc, amrex_mlndlap_mknewu, &
       amrex_mlndlap_divu_fine_contrib, amrex_mlndlap_divu_cf_contrib, &
       ! residual
       amrex_mlndlap_crse_resid, &
       amrex_mlndlap_res_fine_contrib, amrex_mlndlap_res_cf_contrib, &
       ! sync residual
       amrex_mlndlap_zero_fine

contains

  subroutine amrex_mlndlap_set_nodal_mask (lo, hi, nmsk, nlo, nhi, cmsk, clo, chi) &
       bind(c,name='amrex_mlndlap_set_nodal_mask')
    integer, dimension(3), intent(in) :: lo, hi, nlo, nhi, clo, chi
    integer, intent(inout) :: nmsk(nlo(1):nhi(1),nlo(2):nhi(2),nlo(3):nhi(3))
    integer, intent(in   ) :: cmsk(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    integer :: i,j,k
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (all(cmsk(i-1:i,j-1:j,k-1:k) .eq. crse_cell)) then
                nmsk(i,j,k) = crse_node
             else if (all(cmsk(i-1:i,j-1:j,k-1:k) .eq. fine_cell)) then
                nmsk(i,j,k) = fine_node
             else
                nmsk(i,j,k) = crse_fine_node
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_set_nodal_mask


  subroutine amrex_mlndlap_set_dirichlet_mask (dmsk, dlo, dhi, omsk, olo, ohi, &
       domlo, domhi, bclo, bchi) bind(c,name='amrex_mlndlap_set_dirichlet_mask')
    integer, dimension(3) :: dlo, dhi, olo, ohi, domlo, domhi, bclo, bchi
    integer, intent(inout) :: dmsk(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    integer, intent(in   ) :: omsk(olo(1):ohi(1),olo(2):ohi(2),olo(3):ohi(3))

    integer :: i,j,k
    
    do       k = dlo(3), dhi(3)
       do    j = dlo(2), dhi(2)
          do i = dlo(1), dhi(1)
             if (any(omsk(i-1:i,j-1:j,k-1:k).eq.1)) then
                dmsk(i,j,k) = dirichlet
             else
                dmsk(i,j,k) = 0
             end if
          end do
       end do
    end do

    if (dlo(1) .eq. domlo(1)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then
          dmsk(dlo(1),:,:) = 0
       end if
    end if

    if (dhi(1) .eq. domhi(1)) then
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          dmsk(dhi(1),:,:) = 0
       end if
    end if

    if (dlo(2) .eq. domlo(2)) then
       if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          dmsk(:,dlo(2),:) = 0
       end if
    end if

    if (dhi(2) .eq. domhi(2)) then
       if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          dmsk(:,dhi(2),:) = 0
       end if
    end if

    if (dlo(3) .eq. domlo(3)) then
       if (bclo(3) .eq. amrex_lo_neumann .or. bclo(3) .eq. amrex_lo_inflow) then
          dmsk(:,:,dlo(3)) = 0
       end if
    end if

    if (dhi(3) .eq. domhi(3)) then
       if (bchi(3) .eq. amrex_lo_neumann .or. bchi(3) .eq. amrex_lo_inflow) then
          dmsk(:,:,dhi(3)) = 0
       end if
    end if

    if (dlo(1) .eq. domlo(1)) then
       if (bclo(1) .eq. amrex_lo_dirichlet) then
          dmsk(dlo(1),:,:) = dirichlet
       end if
    end if

    if (dhi(1) .eq. domhi(1)) then
       if (bchi(1) .eq. amrex_lo_dirichlet) then
          dmsk(dhi(1),:,:) = dirichlet
       end if
    end if

    if (dlo(2) .eq. domlo(2)) then
       if (bclo(2) .eq. amrex_lo_dirichlet) then
          dmsk(:,dlo(2),:) = dirichlet
       end if
    end if

    if (dhi(2) .eq. domhi(2)) then
       if (bchi(2) .eq. amrex_lo_dirichlet) then
          dmsk(:,dhi(2),:) = dirichlet
       end if
    end if

    if (dlo(3) .eq. domlo(3)) then
       if (bclo(3) .eq. amrex_lo_dirichlet) then
          dmsk(:,:,dlo(3)) = dirichlet
       end if
    end if

    if (dhi(3) .eq. domhi(3)) then
       if (bchi(3) .eq. amrex_lo_dirichlet) then
          dmsk(:,:,dhi(3)) = dirichlet
       end if
    end if
    
  end subroutine amrex_mlndlap_set_dirichlet_mask


  subroutine amrex_mlndlap_fixup_res_mask (lo, hi, rmsk, rlo, rhi, fmsk, flo, fhi) &
       bind(c,name='amrex_mlndlap_fixup_res_mask')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, flo, fhi
    integer, intent(inout) :: rmsk(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    integer, intent(in   ) :: fmsk(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (fmsk(i,j,k) .eq. crse_fine_node) then
                rmsk(i,j,k) = 1
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_fixup_res_mask


  function amrex_mlndlap_any_fine_sync_cells (lo, hi, msk, mlo, mhi, fine_flag) result(r) &
       bind(c,name='amrex_mlndlap_any_fine_sync_cells')
    integer :: r
    integer, dimension(3), intent(in) :: lo, hi, mlo, mhi
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    integer, intent(in) :: fine_flag

    integer :: i,j,k

    r = 0

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .eq. fine_flag) then
                r = 1
                goto 100
             end if
          end do
       end do
    end do
100 continue
  end function amrex_mlndlap_any_fine_sync_cells


  subroutine amrex_mlndlap_avgdown_coeff (lo, hi, crse, clo, chi, fine, flo, fhi, idim) &
       bind(c,name='amrex_mlndlap_avgdown_coeff')
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, flo, fhi
    integer, intent(in) :: idim
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))

    integer :: i,j,k
    real(amrex_real) :: cl,cr

    if (idim .eq. 0) then
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                cl = 0.25d0*(fine(2*i  ,2*j,2*k  )+fine(2*i  ,2*j+1,2*k  ) &
                     &      +fine(2*i  ,2*j,2*k+1)+fine(2*i  ,2*j+1,2*k+1))
                cr = 0.25d0*(fine(2*i+1,2*j,2*k  )+fine(2*i+1,2*j+1,2*k  ) &
                     &      +fine(2*i+1,2*j,2*k+1)+fine(2*i+1,2*j+1,2*k+1))
                crse(i,j,k) = 2.d0*cl*cr/(cl+cr)
             end do
          end do
       end do
    else if (idim .eq. 1) then
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                cl = 0.25d0*(fine(2*i,2*j  ,2*k  )+fine(2*i+1,2*j  ,2*k  ) &
                     &      +fine(2*i,2*j  ,2*k+1)+fine(2*i+1,2*j  ,2*k+1))
                cl = 0.25d0*(fine(2*i,2*j+1,2*k  )+fine(2*i+1,2*j+1,2*k  ) &
                     &      +fine(2*i,2*j+1,2*k+1)+fine(2*i+1,2*j+1,2*k+1))
                crse(i,j,k) = 2.d0*cl*cr/(cl+cr)
             end do
          end do
       end do
    else
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                cl = 0.25d0*(fine(2*i,2*j  ,2*k  )+fine(2*i+1,2*j  ,2*k  ) &
                     &      +fine(2*i,2*j+1,2*k  )+fine(2*i+1,2*j+1,2*k  ))
                cr = 0.25d0*(fine(2*i,2*j  ,2*k+1)+fine(2*i+1,2*j  ,2*k+1) &
                     &      +fine(2*i,2*j+1,2*k+1)+fine(2*i+1,2*j+1,2*k+1))
                crse(i,j,k) = 2.d0*cl*cr/(cl+cr)
             end do
          end do
       end do
    end if
  end subroutine amrex_mlndlap_avgdown_coeff


  subroutine amrex_mlndlap_fillbc_cc (sigma, slo, shi, dlo, dhi, bclo, bchi) &
       bind(c, name='amrex_mlndlap_fillbc_cc')
    integer, dimension(3), intent(in) :: slo, shi, dlo, dhi, bclo, bchi
    real(amrex_real), intent(inout) :: sigma(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    
    integer :: ilo, ihi, jlo, jhi, klo, khi

    ilo = max(dlo(1), slo(1))
    ihi = min(dhi(1), shi(1))
    jlo = max(dlo(2), slo(2))
    jhi = min(dhi(2), shi(2))
    klo = max(dlo(3), slo(3))
    khi = min(dhi(3), shi(3))

    ! faces

    if (bclo(1) .ne. amrex_lo_periodic .and. slo(1) .lt. dlo(1)) then
       sigma(dlo(1)-1,jlo:jhi,klo:khi) = sigma(dlo(1),jlo:jhi,klo:khi)
    end if
    
    if (bchi(1) .ne. amrex_lo_periodic .and. shi(1) .gt. dhi(1)) then
       sigma(dhi(1)+1,jlo:jhi,klo:khi) = sigma(dhi(1),jlo:jhi,klo:khi)
    end if

    if (bclo(2) .ne. amrex_lo_periodic .and. slo(2) .lt. dlo(2)) then
       sigma(ilo:ihi,dlo(2)-1,klo:khi) = sigma(ilo:ihi,dlo(2),klo:khi)
    end if

    if (bchi(2) .ne. amrex_lo_periodic .and. shi(2) .gt. dhi(2)) then
       sigma(ilo:ihi,dhi(2)+1,klo:khi) = sigma(ilo:ihi,dhi(2),klo:khi)
    end if

    if (bclo(3) .ne. amrex_lo_periodic .and. slo(3) .lt. dlo(3)) then
       sigma(ilo:ihi,jlo:jhi,dlo(3)-1) = sigma(ilo:ihi,jlo:jhi,dlo(3))
    end if

    if (bchi(3) .ne. amrex_lo_periodic .and. shi(3) .gt. dhi(3)) then
       sigma(ilo:ihi,jlo:jhi,dhi(3)+1) = sigma(ilo:ihi,jlo:jhi,dhi(3))
    end if

    ! edges

    if (slo(1) .lt. dlo(1) .and. slo(2) .lt. dlo(2)) then
       if (bclo(1) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dlo(2)-1,klo:khi) = sigma(dlo(1),dlo(2)-1,klo:khi)
       else if (bclo(2) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dlo(2)-1,klo:khi) = sigma(dlo(1)-1,dlo(2),klo:khi)
       end if
    end if

    if (shi(1) .gt. dhi(1) .and. slo(2) .lt. dlo(2)) then
       if (bchi(1) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dlo(2)-1,klo:khi) = sigma(dhi(1),dlo(2)-1,klo:khi)
       else if (bclo(2) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dlo(2)-1,klo:khi) = sigma(dhi(1)+1,dlo(2),klo:khi)
       end if
    end if

    if (slo(1) .lt. dlo(1) .and. shi(2) .gt. dhi(2)) then
       if (bclo(1) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dhi(2)+1,klo:khi) = sigma(dlo(1),dhi(2)+1,klo:khi)
       else if (bchi(2) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dhi(2)+1,klo:khi) = sigma(dlo(1)-1,dhi(2),klo:khi)
       end if
    end if

    if (shi(1) .gt. dhi(1) .and. shi(2) .gt. dhi(2)) then
       if (bchi(1) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dhi(2)+1,klo:khi) = sigma(dhi(1),dhi(2)+1,klo:khi)
       else if (bchi(2) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dhi(2)+1,klo:khi) = sigma(dhi(1)+1,dhi(2),klo:khi)
       end if
    end if

    if (slo(1) .lt. dlo(1) .and. slo(3) .lt. dlo(3)) then
       if (bclo(1) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,jlo:jhi,dlo(3)-1) = sigma(dlo(1),jlo:jhi,dlo(3)-1)
       else if (bclo(3) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,jlo:jhi,dlo(3)-1) = sigma(dlo(1)-1,jlo:jhi,dlo(3))
       end if
    end if
    
    if (shi(1) .gt. dhi(1) .and. slo(3) .lt. dlo(3)) then
       if (bchi(1) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,jlo:jhi,dlo(3)-1) = sigma(dhi(1),jlo:jhi,dlo(3)-1)
       else if (bclo(3) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,jlo:jhi,dlo(3)-1) = sigma(dhi(1)+1,jlo:jhi,dlo(3))
       end if
    end if
    
    if (slo(1) .lt. dlo(1) .and. shi(3) .gt. dhi(3)) then
       if (bclo(1) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,jlo:jhi,dhi(3)+1) = sigma(dlo(1),jlo:jhi,dhi(3)+1)
       else if (bchi(3) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,jlo:jhi,dhi(3)+1) = sigma(dlo(1)-1,jlo:jhi,dhi(3))
       end if
    end if

    if (shi(1) .gt. dhi(1) .and. shi(3) .gt. dhi(3)) then
       if (bchi(1) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,jlo:jhi,dhi(3)+1) = sigma(dhi(1),jlo:jhi,dhi(3)+1)
       else if (bchi(3) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,jlo:jhi,dhi(3)+1) = sigma(dhi(1)+1,jlo:jhi,dhi(3))
       end if
    end if

    if (slo(2) .lt. dlo(2) .and. slo(3) .lt. dlo(3)) then
       if (bclo(2) .ne. amrex_lo_periodic) then
          sigma(ilo:ihi,dlo(2)-1,dlo(3)-1) = sigma(ilo:ihi,dlo(2),dlo(3)-1)
       else if (bclo(3) .ne. amrex_lo_periodic) then
          sigma(ilo:ihi,dlo(2)-1,dlo(3)-1) = sigma(ilo:ihi,dlo(2)-1,dlo(3))
       end if
    end if

    if (shi(2) .gt. dhi(2) .and. slo(3) .lt. dlo(3)) then
       if (bchi(2) .ne. amrex_lo_periodic) then
          sigma(ilo:ihi,dhi(2)+1,dlo(3)-1) = sigma(ilo:ihi,dhi(2),dlo(3)-1)
       else if (bclo(3) .ne. amrex_lo_periodic) then
          sigma(ilo:ihi,dhi(2)+1,dlo(3)-1) = sigma(ilo:ihi,dhi(2)+1,dlo(3))
       end if
    end if

    if (slo(2) .lt. dlo(2) .and. shi(3) .gt. dhi(3)) then
       if (bclo(2) .ne. amrex_lo_periodic) then
          sigma(ilo:ihi,dlo(2)-1,dhi(3)+1) = sigma(ilo:ihi,dlo(2),dhi(3)+1)
       else if (bchi(3) .ne. amrex_lo_periodic) then
          sigma(ilo:ihi,dlo(2)-1,dhi(3)+1) = sigma(ilo:ihi,dlo(2)-1,dhi(3))
       end if
    end if

    if (shi(2) .gt. dhi(2) .and. shi(3) .gt. dhi(3)) then
       if (bchi(2) .ne. amrex_lo_periodic) then
          sigma(ilo:ihi,dhi(2)+1,dhi(3)+1) = sigma(ilo:ihi,dhi(2),dhi(3)+1)
       else if (bchi(3) .ne. amrex_lo_periodic) then
          sigma(ilo:ihi,dhi(2)+1,dhi(3)+1) = sigma(ilo:ihi,dhi(2)+1,dhi(3))
       end if
    end if

    ! corners

    if (slo(1) .lt. dlo(1) .and. slo(2) .lt. dlo(2) .and. slo(3) .lt. dlo(3)) then
       if (bclo(1) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dlo(2)-1,dlo(3)-1) = sigma(dlo(1),dlo(2)-1,dlo(3)-1)
       else if (bclo(2) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dlo(2)-1,dlo(3)-1) = sigma(dlo(1)-1,dlo(2),dlo(3)-1)
       else if (bclo(3) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dlo(2)-1,dlo(3)-1) = sigma(dlo(1)-1,dlo(2)-1,dlo(3))
       end if
    end if

    if (shi(1) .gt. dhi(1) .and. slo(2) .lt. dlo(2) .and. slo(3) .lt. dlo(3)) then
       if (bchi(1) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dlo(2)-1,dlo(3)-1) = sigma(dhi(1),dlo(2)-1,dlo(3)-1)
       else if (bclo(2) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dlo(2)-1,dlo(3)-1) = sigma(dhi(1)+1,dlo(2),dlo(3)-1)
       else if (bclo(3) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dlo(2)-1,dlo(3)-1) = sigma(dhi(1)+1,dlo(2)-1,dlo(3))
       end if
    end if

    if (slo(1) .lt. dlo(1) .and. shi(2) .gt. dhi(2) .and. slo(3) .lt. dlo(3)) then
       if (bclo(1) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dhi(2)+1,dlo(3)-1) = sigma(dlo(1),dhi(2)+1,dlo(3)-1)
       else if (bchi(2) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dhi(2)+1,dlo(3)-1) = sigma(dlo(1)-1,dhi(2),dlo(3)-1)
       else if (bclo(3) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dhi(2)+1,dlo(3)-1) = sigma(dlo(1)-1,dhi(2)+1,dlo(3))
       end if
    end if

    if (shi(1) .gt. dhi(1) .and. shi(2) .gt. dhi(2) .and. slo(3) .lt. dlo(3)) then
       if (bchi(1) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dhi(2)+1,dlo(3)-1) = sigma(dhi(1),dhi(2)+1,dlo(3)-1)
       else if (bchi(2) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dhi(2)+1,dlo(3)-1) = sigma(dhi(1)+1,dhi(2),dlo(3)-1)
       else if (bclo(3) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dhi(2)+1,dlo(3)-1) = sigma(dhi(1)+1,dhi(2)+1,dlo(3))
       end if
    end if

    if (slo(1) .lt. dlo(1) .and. slo(2) .lt. dlo(2) .and. shi(3) .gt. dhi(3)) then
       if (bclo(1) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dlo(2)-1,dhi(3)+1) = sigma(dlo(1),dlo(2)-1,dhi(3)+1)
       else if (bclo(2) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dlo(2)-1,dhi(3)+1) = sigma(dlo(1)-1,dlo(2),dhi(3)+1)
       else if (bchi(3) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dlo(2)-1,dhi(3)+1) = sigma(dlo(1)-1,dlo(2)-1,dhi(3))
       end if
    end if

    if (shi(1) .gt. dhi(1) .and. slo(2) .lt. dlo(2) .and. shi(3) .gt. dhi(3)) then
       if (bchi(1) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dlo(2)-1,dhi(3)+1) = sigma(dhi(1),dlo(2)-1,dhi(3)+1)
       else if (bclo(2) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dlo(2)-1,dhi(3)+1) = sigma(dhi(1)+1,dlo(2),dhi(3)+1)
       else if (bchi(3) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dlo(2)-1,dhi(3)+1) = sigma(dhi(1)+1,dlo(2)-1,dhi(3))
       end if
    end if

    if (slo(1) .lt. dlo(1) .and. shi(2) .gt. dhi(2) .and. shi(3) .gt. dhi(3)) then
       if (bclo(1) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dhi(2)+1,dhi(3)+1) = sigma(dlo(1),dhi(2)+1,dhi(3)+1)
       else if (bchi(2) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dhi(2)+1,dhi(3)+1) = sigma(dlo(1)-1,dhi(2),dhi(3)+1)
       else if (bchi(3) .ne. amrex_lo_periodic) then
          sigma(dlo(1)-1,dhi(2)+1,dhi(3)+1) = sigma(dlo(1)-1,dhi(2)+1,dhi(3))
       end if
    end if

    if (shi(1) .gt. dhi(1) .and. shi(2) .gt. dhi(2) .and. shi(3) .gt. dhi(3)) then
       if (bchi(1) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dhi(2)+1,dhi(3)+1) = sigma(dhi(1),dhi(2)+1,dhi(3)+1)
       else if (bchi(2) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dhi(2)+1,dhi(3)+1) = sigma(dhi(1)+1,dhi(2),dhi(3)+1)
       else if (bchi(3) .ne. amrex_lo_periodic) then
          sigma(dhi(1)+1,dhi(2)+1,dhi(3)+1) = sigma(dhi(1)+1,dhi(2)+1,dhi(3))
       end if
    end if

  end subroutine amrex_mlndlap_fillbc_cc


  subroutine amrex_mlndlap_applybc (phi, hlo, hhi, dlo, dhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_applybc')
    integer, dimension(3) :: hlo, hhi, dlo, dhi, bclo, bchi
    real(amrex_real), intent(inout) :: phi(hlo(1):hhi(1),hlo(2):hhi(2),hlo(3):hhi(3))

    integer :: ilo, ihi, jlo, jhi, klo, khi

    if (any(bclo.eq.amrex_lo_inflow) .or. any(bchi.eq.amrex_lo_inflow)) then
       call amrex_error("amrex_mlndlap_applybc: inflow not supported yet")
    end if

    ilo = max(dlo(1), hlo(1))
    ihi = min(dhi(1), hhi(1))
    jlo = max(dlo(2), hlo(2))
    jhi = min(dhi(2), hhi(2))
    klo = max(dlo(3), hlo(3))
    khi = min(dhi(3), hhi(3))

    ! neumann

    if (bclo(1) .eq. amrex_lo_neumann .and. hlo(1) .lt. dlo(1)) then
       phi(dlo(1)-1,jlo:jhi,klo:khi) = phi(dlo(1)+1,jlo:jhi,klo:khi)
    end if

    if (bchi(1) .eq. amrex_lo_neumann .and. hhi(1) .gt. dhi(1)) then
       phi(dhi(1)+1,jlo:jhi,klo:khi) = phi(dhi(1)-1,jlo:jhi,klo:khi)
    end if

    if (bclo(2) .eq. amrex_lo_neumann .and. hlo(2) .lt. dlo(2)) then
       phi(ilo:ihi,dlo(2)-1,klo:khi) = phi(ilo:ihi,dlo(2)+1,klo:khi)
    end if

    if (bchi(2) .eq. amrex_lo_neumann .and. hhi(2) .gt. dhi(2)) then
       phi(ilo:ihi,dhi(2)+1,klo:khi) = phi(ilo:ihi,dhi(2)-1,klo:khi)
    end if

    if (bclo(3) .eq. amrex_lo_neumann .and. hlo(3) .lt. dlo(3)) then
       phi(ilo:ihi,jlo:jhi,dlo(3)-1) = phi(ilo:ihi,jlo:jhi,dlo(3)+1)
    end if
    
    if (bchi(3) .eq. amrex_lo_neumann .and. hhi(3) .gt. dhi(3)) then
       phi(ilo:ihi,jlo:jhi,dhi(3)+1) = phi(ilo:ihi,jlo:jhi,dhi(3)-1)
    end if

    !edges

    if (hlo(1) .lt. dlo(1) .and. hlo(2) .lt. dlo(2)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dlo(2)-1,klo:khi) = phi(dlo(1)+1,dlo(2)-1,klo:khi)
       else if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dlo(2)-1,klo:khi) = phi(dlo(1)-1,dlo(2)+1,klo:khi)
       end if
    end if

    if (hhi(1) .gt. dhi(1) .and. hlo(2) .lt. dlo(2)) then
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dlo(2)-1,klo:khi) = phi(dhi(1)-1,dlo(2)-1,klo:khi)
       else if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dlo(2)-1,klo:khi) = phi(dhi(1)+1,dlo(2)+1,klo:khi)
       end if
    end if

    if (hlo(1) .lt. dlo(1) .and. hhi(2) .gt. dhi(2)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dhi(2)+1,klo:khi) = phi(dlo(1)+1,dhi(2)+1,klo:khi)
       else if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dhi(2)+1,klo:khi) = phi(dlo(1)-1,dhi(2)-1,klo:khi)
       end if
    end if

    if (hhi(1) .gt. dhi(1) .and. hhi(2) .gt. dhi(2)) then
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dhi(2)+1,klo:khi) = phi(dhi(1)-1,dhi(2)+1,klo:khi)
       else  if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dhi(2)+1,klo:khi) = phi(dhi(1)+1,dhi(2)-1,klo:khi)
       end if
    end if

    if (hlo(1) .lt. dlo(1) .and. hlo(3) .lt. dlo(3)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,jlo:jhi,dlo(3)-1) = phi(dlo(1)+1,jlo:jhi,dlo(3)-1)
       else if (bclo(3) .ne. amrex_lo_neumann .or. bclo(3) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,jlo:jhi,dlo(3)-1) = phi(dlo(1)-1,jlo:jhi,dlo(3)+1)
       end if
    end if
    
    if (hhi(1) .gt. dhi(1) .and. hlo(3) .lt. dlo(3)) then
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,jlo:jhi,dlo(3)-1) = phi(dhi(1)-1,jlo:jhi,dlo(3)-1)
       else if (bclo(3) .eq. amrex_lo_neumann .or. bclo(3) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,jlo:jhi,dlo(3)-1) = phi(dhi(1)+1,jlo:jhi,dlo(3)+1)
       end if
    end if
    
    if (hlo(1) .lt. dlo(1) .and. hhi(3) .gt. dhi(3)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,jlo:jhi,dhi(3)+1) = phi(dlo(1)+1,jlo:jhi,dhi(3)+1)
       else if (bchi(3) .eq. amrex_lo_neumann .or. bchi(3) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,jlo:jhi,dhi(3)+1) = phi(dlo(1)-1,jlo:jhi,dhi(3)-1)
       end if
    end if

    if (hhi(1) .gt. dhi(1) .and. hhi(3) .gt. dhi(3)) then
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,jlo:jhi,dhi(3)+1) = phi(dhi(1)-1,jlo:jhi,dhi(3)+1)
       else if (bchi(3) .eq. amrex_lo_neumann .or. bchi(3) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,jlo:jhi,dhi(3)+1) = phi(dhi(1)+1,jlo:jhi,dhi(3)-1)
       end if
    end if

    if (hlo(2) .lt. dlo(2) .and. hlo(3) .lt. dlo(3)) then
       if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          phi(ilo:ihi,dlo(2)-1,dlo(3)-1) = phi(ilo:ihi,dlo(2)+1,dlo(3)-1)
       else if (bclo(3) .eq. amrex_lo_neumann .or. bclo(3) .eq. amrex_lo_inflow) then
          phi(ilo:ihi,dlo(2)-1,dlo(3)-1) = phi(ilo:ihi,dlo(2)-1,dlo(3)+1)
       end if
    end if

    if (hhi(2) .gt. dhi(2) .and. hlo(3) .lt. dlo(3)) then
       if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          phi(ilo:ihi,dhi(2)+1,dlo(3)-1) = phi(ilo:ihi,dhi(2)-1,dlo(3)-1)
       else if (bclo(3) .eq. amrex_lo_neumann .or. bclo(3) .eq. amrex_lo_inflow) then
          phi(ilo:ihi,dhi(2)+1,dlo(3)-1) = phi(ilo:ihi,dhi(2)+1,dlo(3)+1)
       end if
    end if

    if (hlo(2) .lt. dlo(2) .and. hhi(3) .gt. dhi(3)) then
       if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          phi(ilo:ihi,dlo(2)-1,dhi(3)+1) = phi(ilo:ihi,dlo(2)+1,dhi(3)+1)
       else if (bchi(3) .eq. amrex_lo_neumann .or. bchi(3) .eq. amrex_lo_inflow) then
          phi(ilo:ihi,dlo(2)-1,dhi(3)+1) = phi(ilo:ihi,dlo(2)-1,dhi(3)-1)
       end if
    end if

    if (hhi(2) .gt. dhi(2) .and. hhi(3) .gt. dhi(3)) then
       if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          phi(ilo:ihi,dhi(2)+1,dhi(3)+1) = phi(ilo:ihi,dhi(2)-1,dhi(3)+1)
       else if (bchi(3) .eq. amrex_lo_neumann .or. bchi(3) .eq. amrex_lo_inflow) then
          phi(ilo:ihi,dhi(2)+1,dhi(3)+1) = phi(ilo:ihi,dhi(2)+1,dhi(3)-1)
       end if
    end if

    ! corners
    if (hlo(1) .lt. dlo(1) .and. hlo(2) .lt. dlo(2) .and. hlo(3) .lt. dlo(3)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dlo(2)-1,dlo(3)-1) = phi(dlo(1)+1,dlo(2)-1,dlo(3)-1)
       else if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dlo(2)-1,dlo(3)-1) = phi(dlo(1)-1,dlo(2)+1,dlo(3)-1)
       else if (bclo(3) .eq. amrex_lo_neumann .or. bclo(3) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dlo(2)-1,dlo(3)-1) = phi(dlo(1)-1,dlo(2)-1,dlo(3)+1)
       end if
    end if

    if (hhi(1) .gt. dhi(1) .and. hlo(2) .lt. dlo(2) .and. hlo(3) .lt. dlo(3)) then
       if (bchi(1) .ne. amrex_lo_neumann .or. bchi(1) .ne. amrex_lo_inflow) then
          phi(dhi(1)+1,dlo(2)-1,dlo(3)-1) = phi(dhi(1)-1,dlo(2)-1,dlo(3)-1)
       else if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dlo(2)-1,dlo(3)-1) = phi(dhi(1)+1,dlo(2)+1,dlo(3)-1)
       else if (bclo(3) .eq. amrex_lo_neumann .or. bclo(3) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dlo(2)-1,dlo(3)-1) = phi(dhi(1)+1,dlo(2)-1,dlo(3)+1)
       end if
    end if

    if (hlo(1) .lt. dlo(1) .and. hhi(2) .gt. dhi(2) .and. hlo(3) .lt. dlo(3)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dhi(2)+1,dlo(3)-1) = phi(dlo(1)+1,dhi(2)+1,dlo(3)-1)
       else if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dhi(2)+1,dlo(3)-1) = phi(dlo(1)-1,dhi(2)-1,dlo(3)-1)
       else if (bclo(3) .eq. amrex_lo_neumann .or. bclo(3) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dhi(2)+1,dlo(3)-1) = phi(dlo(1)-1,dhi(2)+1,dlo(3)+1)
       end if
    end if

    if (hhi(1) .gt. dhi(1) .and. hhi(2) .gt. dhi(2) .and. hlo(3) .lt. dlo(3)) then
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dhi(2)+1,dlo(3)-1) = phi(dhi(1)-1,dhi(2)+1,dlo(3)-1)
       else if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dhi(2)+1,dlo(3)-1) = phi(dhi(1)+1,dhi(2)-1,dlo(3)-1)
       else if (bclo(3) .eq. amrex_lo_neumann .or. bclo(3) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dhi(2)+1,dlo(3)-1) = phi(dhi(1)+1,dhi(2)+1,dlo(3)+1)
       end if
    end if

    if (hlo(1) .lt. dlo(1) .and. hlo(2) .lt. dlo(2) .and. hhi(3) .gt. dhi(3)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dlo(2)-1,dhi(3)+1) = phi(dlo(1)+1,dlo(2)-1,dhi(3)+1)
       else if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dlo(2)-1,dhi(3)+1) = phi(dlo(1)-1,dlo(2)+1,dhi(3)+1)
       else if (bchi(3) .eq. amrex_lo_neumann .or. bchi(3) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dlo(2)-1,dhi(3)+1) = phi(dlo(1)-1,dlo(2)-1,dhi(3)-1)
       end if
    end if

    if (hhi(1) .gt. dhi(1) .and. hlo(2) .lt. dlo(2) .and. hhi(3) .gt. dhi(3)) then
       if (bchi(1) .ne. amrex_lo_neumann .or. bchi(1) .ne. amrex_lo_inflow) then
          phi(dhi(1)+1,dlo(2)-1,dhi(3)+1) = phi(dhi(1)-1,dlo(2)-1,dhi(3)+1)
       else if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dlo(2)-1,dhi(3)+1) = phi(dhi(1)+1,dlo(2)+1,dhi(3)+1)
       else if (bchi(3) .eq. amrex_lo_neumann .or. bchi(3) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dlo(2)-1,dhi(3)+1) = phi(dhi(1)+1,dlo(2)-1,dhi(3)-1)
       end if
    end if

    if (hlo(1) .lt. dlo(1) .and. hhi(2) .gt. dhi(2) .and. hhi(3) .gt. dhi(3)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dhi(2)+1,dhi(3)+1) = phi(dlo(1)+1,dhi(2)+1,dhi(3)+1)
       else if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dhi(2)+1,dhi(3)+1) = phi(dlo(1)-1,dhi(2)-1,dhi(3)+1)
       else if (bchi(3) .eq. amrex_lo_neumann .or. bchi(3) .eq. amrex_lo_inflow) then
          phi(dlo(1)-1,dhi(2)+1,dhi(3)+1) = phi(dlo(1)-1,dhi(2)+1,dhi(3)-1)
       end if
    end if

    if (hhi(1) .gt. dhi(1) .and. hhi(2) .gt. dhi(2) .and. hhi(3) .gt. dhi(3)) then
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dhi(2)+1,dhi(3)+1) = phi(dhi(1)-1,dhi(2)+1,dhi(3)+1)
       else if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dhi(2)+1,dhi(3)+1) = phi(dhi(1)+1,dhi(2)-1,dhi(3)+1)
       else if (bchi(3) .eq. amrex_lo_neumann .or. bchi(3) .eq. amrex_lo_inflow) then
          phi(dhi(1)+1,dhi(2)+1,dhi(3)+1) = phi(dhi(1)+1,dhi(2)+1,dhi(3)-1)
       end if
    end if

  end subroutine amrex_mlndlap_applybc


  subroutine amrex_mlndlap_adotx_ha (lo, hi, y, ylo, yhi, x, xlo, xhi, &
       sx, sxlo, sxhi, sy, sylo, syhi, sz, szlo, szhi, msk, mlo, mhi, &
       dxinv, domlo, domhi, bclo, bchi) bind(c,name='amrex_mlndlap_adotx_ha')
    integer, dimension(3), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, sxlo, sxhi, &
         sylo, syhi, szlo, szhi, mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) ::  y( ylo(1): yhi(1), ylo(2): yhi(2), ylo(3): yhi(3))
    real(amrex_real), intent(in   ) ::  x( xlo(1): xhi(1), xlo(2): xhi(2), xlo(3): xhi(3))
    real(amrex_real), intent(in   ) :: sx(sxlo(1):sxhi(1),sxlo(2):sxhi(2),sxlo(3):sxhi(3))
    real(amrex_real), intent(in   ) :: sy(sylo(1):syhi(1),sylo(2):syhi(2),sylo(3):syhi(3))
    real(amrex_real), intent(in   ) :: sz(szlo(1):szhi(1),szlo(2):szhi(2),szlo(3):szhi(3))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k
    real(amrex_real) :: facx, facy, facz
    
    facx = (1.d0/36.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/36.d0)*dxinv(2)*dxinv(2)
    facz = (1.d0/36.d0)*dxinv(3)*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                y(i,j,k) = x(i,j,k)*(-4.d0)*(facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1) &
                     &                            +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  )) &
                     &                      +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1) &
                     &                            +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  )) &
                     &                      +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1) &
                     &                            +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )))
                y(i,j,k) = y(i,j,k) &
                     &   + x(i-1,j-1,k-1)*(facx*sx(i-1,j-1,k-1) &
                     &                    +facy*sy(i-1,j-1,k-1) &
                     &                    +facz*sz(i-1,j-1,k-1)) &
                     &   + x(i+1,j-1,k-1)*(facx*sx(i  ,j-1,k-1) &
                     &                    +facy*sy(i  ,j-1,k-1) &
                     &                    +facz*sz(i  ,j-1,k-1)) &
                     &   + x(i-1,j+1,k-1)*(facx*sx(i-1,j  ,k-1) &
                     &                    +facy*sy(i-1,j  ,k-1) &
                     &                    +facz*sz(i-1,j  ,k-1)) &
                     &   + x(i+1,j+1,k-1)*(facx*sx(i  ,j  ,k-1) &
                     &                    +facy*sy(i  ,j  ,k-1) &
                     &                    +facz*sz(i  ,j  ,k-1)) &
                     &   + x(i-1,j-1,k+1)*(facx*sx(i-1,j-1,k  ) &
                     &                    +facy*sy(i-1,j-1,k  ) &
                     &                    +facz*sz(i-1,j-1,k  )) &
                     &   + x(i+1,j-1,k+1)*(facx*sx(i  ,j-1,k  ) &
                     &                    +facy*sy(i  ,j-1,k  ) &
                     &                    +facz*sz(i  ,j-1,k  )) &
                     &   + x(i-1,j+1,k+1)*(facx*sx(i-1,j  ,k  ) &
                     &                    +facy*sy(i-1,j  ,k  ) &
                     &                    +facz*sz(i-1,j  ,k  )) &
                     &   + x(i+1,j+1,k+1)*(facx*sx(i  ,j  ,k  ) &
                     &                    +facy*sy(i  ,j  ,k  ) &
                     &                    +facz*sz(i  ,j  ,k  ))
                y(i,j,k) = y(i,j,k) &
                     + x(i  ,j-1,k-1)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)) &
                     &                +2.d0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)) &
                     &                +2.d0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1))) &
                     + x(i  ,j+1,k-1)*(    -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1)) &
                     &                +2.d0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1)) &
                     &                +2.d0*facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1))) &
                     + x(i  ,j-1,k+1)*(    -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  )) &
                     &                +2.d0*facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  )) &
                     &                +2.d0*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  ))) &
                     + x(i  ,j+1,k+1)*(    -facx*(sx(i-1,j  ,k  )+sx(i,j  ,k  )) &
                     &                +2.d0*facy*(sy(i-1,j  ,k  )+sy(i,j  ,k  )) &
                     &                +2.d0*facz*(sz(i-1,j  ,k  )+sz(i,j  ,k  ))) &
                     !
                     + x(i-1,j  ,k-1)*(2.d0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1)) &
                     &                     -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1)) &
                     &                +2.d0*facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1))) &
                     + x(i+1,j  ,k-1)*(2.d0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1)) &
                     &                     -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1)) &
                     &                +2.d0*facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1))) &
                     + x(i-1,j  ,k+1)*(2.d0*facx*(sx(i-1,j-1,k  )+sx(i-1,j,k  )) &
                     &                     -facy*(sy(i-1,j-1,k  )+sy(i-1,j,k  )) &
                     &                +2.d0*facz*(sz(i-1,j-1,k  )+sz(i-1,j,k  ))) &
                     + x(i+1,j  ,k+1)*(2.d0*facx*(sx(i  ,j-1,k  )+sx(i  ,j,k  )) &
                     &                     -facy*(sy(i  ,j-1,k  )+sy(i  ,j,k  )) &
                     &                +2.d0*facz*(sz(i  ,j-1,k  )+sz(i  ,j,k  ))) &
                     !
                     + x(i-1,j-1,k  )*(2.d0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j-1,k)) &
                     &                +2.d0*facy*(sy(i-1,j-1,k-1)+sy(i-1,j-1,k)) &
                     &                     -facz*(sz(i-1,j-1,k-1)+sz(i-1,j-1,k))) &
                     + x(i+1,j-1,k  )*(2.d0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j-1,k)) &
                     &                +2.d0*facy*(sy(i  ,j-1,k-1)+sy(i  ,j-1,k)) &
                     &                     -facz*(sz(i  ,j-1,k-1)+sz(i  ,j-1,k))) &
                     + x(i-1,j+1,k  )*(2.d0*facx*(sx(i-1,j  ,k-1)+sx(i-1,j  ,k)) &
                     &                +2.d0*facy*(sy(i-1,j  ,k-1)+sy(i-1,j  ,k)) &
                     &                     -facz*(sz(i-1,j  ,k-1)+sz(i-1,j  ,k))) &
                     + x(i+1,j+1,k  )*(2.d0*facx*(sx(i  ,j  ,k-1)+sx(i  ,j  ,k)) &
                     &                +2.d0*facy*(sy(i  ,j  ,k-1)+sy(i  ,j  ,k)) &
                     &                     -facz*(sz(i  ,j  ,k-1)+sz(i  ,j  ,k)))
                y(i,j,k) = y(i,j,k) &
                     + 2.d0*x(i-1,j,k)*(2.d0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1)+sx(i-1,j-1,k)+sx(i-1,j,k)) &
                     &                      -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1)+sy(i-1,j-1,k)+sy(i-1,j,k)) &
                     &                      -facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)+sz(i-1,j-1,k)+sz(i-1,j,k))) &
                     + 2.d0*x(i  ,j,k)*(2.d0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1)+sx(i  ,j-1,k)+sx(i  ,j,k)) &
                     &                      -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1)+sy(i  ,j-1,k)+sy(i  ,j,k)) &
                     &                      -facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)+sz(i  ,j-1,k)+sz(i  ,j,k))) &
                     + 2.d0*x(i,j-1,k)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j-1,k)+sx(i,j-1,k)) &
                     &                 +2.d0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j-1,k)+sy(i,j-1,k)) &
                     &                      -facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j-1,k)+sz(i,j-1,k))) &
                     + 2.d0*x(i,j  ,k)*(    -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1)+sx(i-1,j  ,k)+sx(i,j  ,k)) &
                     &                 +2.d0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1)+sy(i-1,j  ,k)+sy(i,j  ,k)) &
                     &                      -facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)+sz(i-1,j  ,k)+sz(i,j  ,k))) &
                     + 2.d0*x(i,j,k-1)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)) &
                     &                      -facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)) &
                     &                 +2.d0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1))) &
                     + 2.d0*x(i,j,k  )*(    -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  )) &
                     &                      -facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  )) &
                     &                 +2.d0*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )))
             else
                y(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_adotx_ha


  subroutine amrex_mlndlap_adotx_aa (lo, hi, y, ylo, yhi, x, xlo, xhi, &
       sig, slo, shi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_adotx_aa')
    integer, dimension(3), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, slo, shi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) ::   y(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    real(amrex_real), intent(in   ) ::   x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k
    real(amrex_real) :: facx, facy, facz, fxyz, fmx2y2z, f2xmy2z, f2x2ymz
    real(amrex_real) :: f4xm2ym2z, fm2x4ym2z, fm2xm2y4z
    
    facx = (1.d0/36.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/36.d0)*dxinv(2)*dxinv(2)
    facz = (1.d0/36.d0)*dxinv(3)*dxinv(3)
    fxyz = facx + facy + facz
    fmx2y2z = -facx + 2.d0*facy + 2.d0*facz
    f2xmy2z = 2.d0*facx - facy + 2.d0*facz
    f2x2ymz = 2.d0*facx + 2.d0*facy - facz
    f4xm2ym2z = 4.d0*facx - 2.d0*facy - 2.d0*facz
    fm2x4ym2z = -2.d0*facx + 4.d0*facy - 2.d0*facz
    fm2xm2y4z = -2.d0*facx - 2.d0*facy + 4.d0*facz

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                y(i,j,k) = x(i,j,k)*(-4.d0)*fxyz* &
                     (sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1) &
                     +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )) &
                     !
                     + fxyz*(x(i-1,j-1,k-1)*sig(i-1,j-1,k-1) &
                     &     + x(i+1,j-1,k-1)*sig(i  ,j-1,k-1) &
                     &     + x(i-1,j+1,k-1)*sig(i-1,j  ,k-1) &
                     &     + x(i+1,j+1,k-1)*sig(i  ,j  ,k-1) &
                     &     + x(i-1,j-1,k+1)*sig(i-1,j-1,k  ) &
                     &     + x(i+1,j-1,k+1)*sig(i  ,j-1,k  ) &
                     &     + x(i-1,j+1,k+1)*sig(i-1,j  ,k  ) &
                     &     + x(i+1,j+1,k+1)*sig(i  ,j  ,k  )) &
                     !
                     + fmx2y2z*(x(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)) &
                     &        + x(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)) &
                     &        + x(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )) &
                     &        + x(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  ))) &
                     !
                     + f2xmy2z*(x(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)) &
                     &        + x(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)) &
                     &        + x(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  )) &
                     &        + x(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  ))) &
                     !
                     + f2x2ymz*(x(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k)) &
                     &        + x(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k)) &
                     &        + x(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k)) &
                     &        + x(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k))) &
                     !
                     + f4xm2ym2z*(x(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k)) &
                     &          + x(i  ,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k))) &
                     + fm2x4ym2z*(x(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k)) &
                     &          + x(i,j  ,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k))) &
                     + fm2xm2y4z*(x(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)) &
                     &          + x(i,j,k  )*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )))
             else
                y(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_adotx_aa


  subroutine amrex_mlndlap_jacobi_ha (lo, hi, sol, slo, shi, Ax, alo, ahi, rhs, rlo, rhi, &
       sx, sxlo, sxhi, sy, sylo, syhi, sz, szlo, szhi, msk, mlo, mhi, &
       dxinv, domlo, domhi, bclo, bchi) bind(c,name='amrex_mlndlap_jacobi_ha')
    integer, dimension(3),intent(in) :: lo,hi,slo,shi,alo,ahi,rlo,rhi,sxlo,sxhi,sylo,syhi, &
         szlo, szhi, mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    real(amrex_real), intent(in   ) :: Ax ( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) :: sx (sxlo(1):sxhi(1),sxlo(2):sxhi(2),sxlo(3):sxhi(3))
    real(amrex_real), intent(in   ) :: sy (sylo(1):syhi(1),sylo(2):syhi(2),sylo(3):syhi(3))
    real(amrex_real), intent(in   ) :: sz (szlo(1):szhi(1),szlo(2):szhi(2),szlo(3):szhi(3))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k
    real(amrex_real) :: facx, facy, facz
    real(amrex_real), parameter :: omega = 2.d0/3.d0

    facx = -4.d0 * (1.d0/36.d0)*dxinv(1)*dxinv(1)
    facy = -4.d0 * (1.d0/36.d0)*dxinv(2)*dxinv(2)
    facz = -4.d0 * (1.d0/36.d0)*dxinv(3)*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                sol(i,j,k) = sol(i,j,k) + omega * (rhs(i,j,k) - Ax(i,j,k)) &
                     / (facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1) &
                     &       +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  )) &
                     & +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1) &
                     &       +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  )) &
                     & +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1) &
                     &       +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )))
             else
                sol(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_jacobi_ha


  subroutine amrex_mlndlap_jacobi_aa (lo, hi, sol, slo, shi, Ax, alo, ahi, rhs, rlo, rhi, &
       sig, sglo, sghi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_jacobi_aa')
    integer, dimension(3),intent(in) :: lo,hi,slo,shi,alo,ahi,rlo,rhi,sglo,sghi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    real(amrex_real), intent(in   ) :: Ax ( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) :: sig(sglo(1):sghi(1),sglo(2):sghi(2),sglo(3):sghi(3))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k
    real(amrex_real) :: fxyz
    real(amrex_real), parameter :: omega = 2.d0/3.d0

    fxyz = -4.d0 * (1.d0/36.d0)*(dxinv(1)*dxinv(1) + dxinv(2)*dxinv(2) + dxinv(3)*dxinv(3))

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                sol(i,j,k) = sol(i,j,k) + omega * (rhs(i,j,k) - Ax(i,j,k)) &
                     / (fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1) &
                     &       +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )))

             else
                sol(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_jacobi_aa


  subroutine amrex_mlndlap_gauss_seidel_ha (lo, hi, sol, slo, shi, rhs, rlo, rhi, &
       sx, sxlo, sxhi, sy, sylo, syhi, sz, szlo, szhi, msk, mlo, mhi, &
       dxinv, domlo, domhi, bclo, bchi) bind(c,name='amrex_mlndlap_gauss_seidel_ha')
    integer, dimension(3),intent(in) :: lo,hi,slo,shi,rlo,rhi,sxlo,sxhi,sylo,syhi, &
         szlo, szhi, mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) :: sx (sxlo(1):sxhi(1),sxlo(2):sxhi(2),sxlo(3):sxhi(3))
    real(amrex_real), intent(in   ) :: sy (sylo(1):syhi(1),sylo(2):syhi(2),sylo(3):syhi(3))
    real(amrex_real), intent(in   ) :: sz (szlo(1):szhi(1),szlo(2):szhi(2),szlo(3):szhi(3))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k
    real(amrex_real) :: facx, facy, facz, Ax, s0
    
    facx = (1.d0/36.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/36.d0)*dxinv(2)*dxinv(2)
    facz = (1.d0/36.d0)*dxinv(3)*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                s0 = (-4.d0)*(facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1) &
                     &             +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  )) &
                     &       +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1) &
                     &             +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  )) &
                     &       +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1) &
                     &             +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )))
                Ax = sol(i,j,k)*s0 &
                     !
                     +  sol(i-1,j-1,k-1)*(facx*sx(i-1,j-1,k-1) &
                                        +facy*sy(i-1,j-1,k-1) &
                                        +facz*sz(i-1,j-1,k-1)) &
                     + sol(i+1,j-1,k-1)*(facx*sx(i  ,j-1,k-1) &
                                        +facy*sy(i  ,j-1,k-1) &
                                        +facz*sz(i  ,j-1,k-1)) &
                     + sol(i-1,j+1,k-1)*(facx*sx(i-1,j  ,k-1) &
                                        +facy*sy(i-1,j  ,k-1) &
                                        +facz*sz(i-1,j  ,k-1)) &
                     + sol(i+1,j+1,k-1)*(facx*sx(i  ,j  ,k-1) &
                                        +facy*sy(i  ,j  ,k-1) &
                                        +facz*sz(i  ,j  ,k-1)) &
                     + sol(i-1,j-1,k+1)*(facx*sx(i-1,j-1,k  ) &
                                        +facy*sy(i-1,j-1,k  ) &
                                        +facz*sz(i-1,j-1,k  )) &
                     + sol(i+1,j-1,k+1)*(facx*sx(i  ,j-1,k  ) &
                                        +facy*sy(i  ,j-1,k  ) &
                                        +facz*sz(i  ,j-1,k  )) &
                     + sol(i-1,j+1,k+1)*(facx*sx(i-1,j  ,k  ) &
                                        +facy*sy(i-1,j  ,k  ) &
                                        +facz*sz(i-1,j  ,k  )) &
                     + sol(i+1,j+1,k+1)*(facx*sx(i  ,j  ,k  ) &
                                        +facy*sy(i  ,j  ,k  ) &
                                        +facz*sz(i  ,j  ,k  ))
                Ax = Ax  &
                     +sol(i  ,j-1,k-1)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)) &
                     &                 +2.d0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)) &
                     &                 +2.d0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1))) &
                     +sol(i  ,j+1,k-1)*(    -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1)) &
                     &                 +2.d0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1)) &
                     &                 +2.d0*facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1))) &
                     +sol(i  ,j-1,k+1)*(    -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  )) &
                     &                 +2.d0*facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  )) &
                     &                 +2.d0*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  ))) &
                     +sol(i  ,j+1,k+1)*(    -facx*(sx(i-1,j  ,k  )+sx(i,j  ,k  )) &
                     &                 +2.d0*facy*(sy(i-1,j  ,k  )+sy(i,j  ,k  )) &
                     &                 +2.d0*facz*(sz(i-1,j  ,k  )+sz(i,j  ,k  ))) &
                     !
                     +sol(i-1,j  ,k-1)*(2.d0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1)) &
                     &                      -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1)) &
                     &                 +2.d0*facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1))) &
                     +sol(i+1,j  ,k-1)*(2.d0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1)) &
                     &                      -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1)) &
                     &                 +2.d0*facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1))) &
                     +sol(i-1,j  ,k+1)*(2.d0*facx*(sx(i-1,j-1,k  )+sx(i-1,j,k  )) &
                     &                      -facy*(sy(i-1,j-1,k  )+sy(i-1,j,k  )) &
                     &                 +2.d0*facz*(sz(i-1,j-1,k  )+sz(i-1,j,k  ))) &
                     +sol(i+1,j  ,k+1)*(2.d0*facx*(sx(i  ,j-1,k  )+sx(i  ,j,k  )) &
                     &                      -facy*(sy(i  ,j-1,k  )+sy(i  ,j,k  )) &
                     &                 +2.d0*facz*(sz(i  ,j-1,k  )+sz(i  ,j,k  ))) &
                     !
                     +sol(i-1,j-1,k  )*(2.d0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j-1,k)) &
                     &                 +2.d0*facy*(sy(i-1,j-1,k-1)+sy(i-1,j-1,k)) &
                     &                      -facz*(sz(i-1,j-1,k-1)+sz(i-1,j-1,k))) &
                     +sol(i+1,j-1,k  )*(2.d0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j-1,k)) &
                     &                 +2.d0*facy*(sy(i  ,j-1,k-1)+sy(i  ,j-1,k)) &
                     &                      -facz*(sz(i  ,j-1,k-1)+sz(i  ,j-1,k))) &
                     +sol(i-1,j+1,k  )*(2.d0*facx*(sx(i-1,j  ,k-1)+sx(i-1,j  ,k)) &
                     &                 +2.d0*facy*(sy(i-1,j  ,k-1)+sy(i-1,j  ,k)) &
                     &                      -facz*(sz(i-1,j  ,k-1)+sz(i-1,j  ,k))) &
                     +sol(i+1,j+1,k  )*(2.d0*facx*(sx(i  ,j  ,k-1)+sx(i  ,j  ,k)) &
                     &                 +2.d0*facy*(sy(i  ,j  ,k-1)+sy(i  ,j  ,k)) &
                     &                      -facz*(sz(i  ,j  ,k-1)+sz(i  ,j  ,k)))
                Ax = Ax &
                     + 2.d0*sol(i-1,j,k)*(2.d0*facx*(sx(i-1,j-1,k-1)+sx(i-1,j,k-1)+sx(i-1,j-1,k)+sx(i-1,j,k)) &
                     &                        -facy*(sy(i-1,j-1,k-1)+sy(i-1,j,k-1)+sy(i-1,j-1,k)+sy(i-1,j,k)) &
                     &                        -facz*(sz(i-1,j-1,k-1)+sz(i-1,j,k-1)+sz(i-1,j-1,k)+sz(i-1,j,k))) &
                     + 2.d0*sol(i  ,j,k)*(2.d0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1)+sx(i  ,j-1,k)+sx(i  ,j,k)) &
                     &                        -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1)+sy(i  ,j-1,k)+sy(i  ,j,k)) &
                     &                        -facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)+sz(i  ,j-1,k)+sz(i  ,j,k))) &
                     + 2.d0*sol(i,j-1,k)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j-1,k)+sx(i,j-1,k)) &
                     &                   +2.d0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j-1,k)+sy(i,j-1,k)) &
                     &                        -facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j-1,k)+sz(i,j-1,k))) &
                     + 2.d0*sol(i,j  ,k)*(    -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1)+sx(i-1,j  ,k)+sx(i,j  ,k)) &
                     &                   +2.d0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1)+sy(i-1,j  ,k)+sy(i,j  ,k)) &
                     &                        -facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)+sz(i-1,j  ,k)+sz(i,j  ,k))) &
                     + 2.d0*sol(i,j,k-1)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)) &
                     &                        -facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)) &
                     &                   +2.d0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1))) &
                     + 2.d0*sol (i,j,k  )*(    -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  )) &
                     &                        -facy*(sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  )) &
                     &                   +2.d0*facz*(sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  )))
                
                sol(i,j,k) = sol(i,j,k) + (rhs(i,j,k) - Ax) / s0
             else
                sol(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_gauss_seidel_ha


  subroutine amrex_mlndlap_gauss_seidel_aa (lo, hi, sol, slo, shi, rhs, rlo, rhi, &
       sig, sglo, sghi, msk, mlo, mhi, dxinv, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_gauss_seidel_aa')
    integer, dimension(3),intent(in) :: lo,hi,slo,shi,rlo,rhi,sglo,sghi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) :: sig(sglo(1):sghi(1),sglo(2):sghi(2),sglo(3):sghi(3))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k
    real(amrex_real) :: Ax, s0, facx, facy, facz, fxyz, fmx2y2z, f2xmy2z, f2x2ymz
    real(amrex_real) :: f4xm2ym2z, fm2x4ym2z, fm2xm2y4z
    
    facx = (1.d0/36.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/36.d0)*dxinv(2)*dxinv(2)
    facz = (1.d0/36.d0)*dxinv(3)*dxinv(3)
    fxyz = facx + facy + facz
    fmx2y2z = -facx + 2.d0*facy + 2.d0*facz
    f2xmy2z = 2.d0*facx - facy + 2.d0*facz
    f2x2ymz = 2.d0*facx + 2.d0*facy - facz
    f4xm2ym2z = 4.d0*facx - 2.d0*facy - 2.d0*facz
    fm2x4ym2z = -2.d0*facx + 4.d0*facy - 2.d0*facz
    fm2xm2y4z = -2.d0*facx - 2.d0*facy + 4.d0*facz

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                s0 = (-4.d0)*fxyz*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1) &
                     &            +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  ))
                Ax = sol(i,j,k)*s0 &
                     + fxyz*(sol(i-1,j-1,k-1)*sig(i-1,j-1,k-1) &
                     &     + sol(i+1,j-1,k-1)*sig(i  ,j-1,k-1) &
                     &     + sol(i-1,j+1,k-1)*sig(i-1,j  ,k-1) &
                     &     + sol(i+1,j+1,k-1)*sig(i  ,j  ,k-1) &
                     &     + sol(i-1,j-1,k+1)*sig(i-1,j-1,k  ) &
                     &     + sol(i+1,j-1,k+1)*sig(i  ,j-1,k  ) &
                     &     + sol(i-1,j+1,k+1)*sig(i-1,j  ,k  ) &
                     &     + sol(i+1,j+1,k+1)*sig(i  ,j  ,k  )) &
                     !
                     + fmx2y2z*(sol(i  ,j-1,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)) &
                     &        + sol(i  ,j+1,k-1)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)) &
                     &        + sol(i  ,j-1,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )) &
                     &        + sol(i  ,j+1,k+1)*(sig(i-1,j  ,k  )+sig(i,j  ,k  ))) &
                     !
                     + f2xmy2z*(sol(i-1,j  ,k-1)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)) &
                     &        + sol(i+1,j  ,k-1)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)) &
                     &        + sol(i-1,j  ,k+1)*(sig(i-1,j-1,k  )+sig(i-1,j,k  )) &
                     &        + sol(i+1,j  ,k+1)*(sig(i  ,j-1,k  )+sig(i  ,j,k  ))) &
                     !
                     + f2x2ymz*(sol(i-1,j-1,k  )*(sig(i-1,j-1,k-1)+sig(i-1,j-1,k)) &
                     &        + sol(i+1,j-1,k  )*(sig(i  ,j-1,k-1)+sig(i  ,j-1,k)) &
                     &        + sol(i-1,j+1,k  )*(sig(i-1,j  ,k-1)+sig(i-1,j  ,k)) &
                     &        + sol(i+1,j+1,k  )*(sig(i  ,j  ,k-1)+sig(i  ,j  ,k))) &
                     !
                     + f4xm2ym2z*(sol(i-1,j,k)*(sig(i-1,j-1,k-1)+sig(i-1,j,k-1)+sig(i-1,j-1,k)+sig(i-1,j,k)) &
                     &          + sol(i  ,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k))) &
                     + fm2x4ym2z*(sol(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k)) &
                     &          + sol(i,j  ,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k))) &
                     + fm2xm2y4z*(sol(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)) &
                     &          + sol(i,j,k  )*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )))

                sol(i,j,k) = sol(i,j,k) + (rhs(i,j,k) - Ax) / s0
             else
                sol(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_gauss_seidel_aa


  subroutine amrex_mlndlap_restriction (lo, hi, crse, clo, chi, fine, flo, fhi, msk, mlo, mhi, &
       domlo, domhi, bclo, bchi) bind(c,name='amrex_mlndlap_restriction')
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, flo, fhi, mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i, j, k, ii, jj, kk
    real(amrex_real), parameter :: fac1 = 1.d0/64.d0
    real(amrex_real), parameter :: fac2 = 1.d0/32.d0
    real(amrex_real), parameter :: fac3 = 1.d0/16.d0
    real(amrex_real), parameter :: fac4 = 1.d0/8.d0

    do k = lo(3), hi(3)
       kk = 2*k
       do j = lo(2), hi(2)
          jj = 2*j
          do i = lo(1), hi(1)
             ii = 2*i
             if (msk(ii,jj,kk) .ne. dirichlet) then
                crse(i,j,k) = fac1*(fine(ii-1,jj-1,kk-1)+fine(ii+1,jj-1,kk-1) &
                     &             +fine(ii-1,jj+1,kk-1)+fine(ii+1,jj+1,kk-1) &
                     &             +fine(ii-1,jj-1,kk+1)+fine(ii+1,jj-1,kk+1) &
                     &             +fine(ii-1,jj+1,kk+1)+fine(ii+1,jj+1,kk+1)) &
                     !
                     &      + fac2*(fine(ii  ,jj-1,kk-1)+fine(ii  ,jj+1,kk-1) &
                     &             +fine(ii  ,jj-1,kk+1)+fine(ii  ,jj+1,kk+1) &
                     &             +fine(ii-1,jj  ,kk-1)+fine(ii+1,jj  ,kk-1) &
                     &             +fine(ii-1,jj  ,kk+1)+fine(ii+1,jj  ,kk+1) &
                     &             +fine(ii-1,jj-1,kk  )+fine(ii+1,jj-1,kk  ) &
                     &             +fine(ii-1,jj+1,kk  )+fine(ii+1,jj+1,kk  )) &
                     !
                     &      + fac3*(fine(ii-1,jj,kk)+fine(ii+1,jj,kk) &
                     &             +fine(ii,jj-1,kk)+fine(ii,jj+1,kk) &
                     &             +fine(ii,jj,kk-1)+fine(ii,jj,kk+1)) &
                     &      + fac4*fine(ii,jj,kk)
             else
                crse(i,j,k) = 0.d0
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_restriction


  subroutine amrex_mlndlap_interpolation_ha (clo, chi, fine, fflo, ffhi, crse, cflo, cfhi, &
       sigx, sxlo, sxhi, sigy, sylo, syhi, sigz, szlo, szhi, msk, mlo, mhi, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_interpolation_ha')
    integer, dimension(3), intent(in) :: clo,chi,fflo,ffhi,cflo,cfhi,sxlo,sxhi,sylo,syhi, &
         szlo, szhi, mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in   ) :: crse(cflo(1):cfhi(1),cflo(2):cfhi(2),cflo(3):cfhi(3))
    real(amrex_real), intent(inout) :: fine(fflo(1):ffhi(1),fflo(2):ffhi(2),fflo(3):ffhi(3))
    real(amrex_real), intent(in   ) :: sigx(sxlo(1):sxhi(1),sxlo(2):sxhi(2),sxlo(3):sxhi(3))
    real(amrex_real), intent(in   ) :: sigy(sylo(1):syhi(1),sylo(2):syhi(2),sylo(3):syhi(3))
    real(amrex_real), intent(in   ) :: sigz(szlo(1):szhi(1),szlo(2):szhi(2),szlo(3):szhi(3))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: flo(3), fhi(3), i, j, k, ii, jj, kk
    real(amrex_real) :: w1, w2, w3, w4, w5, w6

    flo = 2*clo
    fhi = 2*chi

    do k = clo(3), chi(3)
       kk = 2*k
       do j = clo(2), chi(2)
          jj = 2*j
          do i = clo(1), chi(1)
             ii = 2*i

             if (msk(ii,jj,kk) .ne. dirichlet) then
                fine(ii,jj,kk) = crse(i,j,k)
             else
                fine(ii,jj,kk) = 0.d0
             end if

             if (ii+1 .le. flo(1)) then
                if (msk(ii+1,jj,kk) .ne. dirichlet) then
                   w1 = sum(sigx(ii  ,jj-1:jj,kk-1:kk))
                   w2 = sum(sigx(ii+1,jj-1:jj,kk-1:kk))
                   fine(ii+1,jj,kk) = (w1*crse(i,j,k)+w2*crse(i+1,j,k))/(w1+w2)
                else
                   fine(ii+1,jj,kk) = 0.d0
                end if
             end if

             if (jj+1 .lt. flo(2)) then
                if (msk(ii,jj+1,kk) .ne. dirichlet) then
                   w1 = sum(sigy(ii-1:ii,jj  ,kk-1:kk))
                   w2 = sum(sigy(ii-1:ii,jj+1,kk-1:kk))
                   fine(ii,jj+1,kk) = (w1*crse(i,j,k)+w2*crse(i,j+1,k))/(w1+w2)
                else
                   fine(ii,jj+1,kk) = 0.d0
                end if
             end if

             if (kk+1 .lt. flo(3)) then
                if (msk(ii,jj,kk+1) .ne. dirichlet) then
                   w1 = sum(sigz(ii-1:ii,jj-1:jj,kk  ))
                   w2 = sum(sigz(ii-1:ii,jj-1:jj,kk+1))
                   fine(ii,jj,kk+1) = (w1*crse(i,j,k)+w2*crse(i,j,k+1))/(w1+w2)
                else
                   fine(ii,jj,kk+1) = 0.d0
                end if
             end if
          end do
       end do
    end do

    do k = clo(3), chi(3)
       kk = 2*k
       do j = clo(2), chi(2)
          jj = 2*j
          do i = clo(1), chi(1)
             ii = 2*i

             if (ii+1 .le. flo(1) .and. jj+1 .le. flo(2)) then
                if (msk(ii+1,jj+1,kk) .ne. dirichlet) then
                   w1 = sum(sigx(ii     ,jj:jj+1,kk-1:kk))
                   w2 = sum(sigx(ii+1   ,jj:jj+1,kk-1:kk))
                   w3 = sum(sigy(ii:ii+1,jj     ,kk-1:kk))
                   w4 = sum(sigy(ii:ii+1,jj+1   ,kk-1:kk))
                   fine(ii+1,jj+1,kk) = (w1*fine(ii,jj+1,kk) + w2*fine(ii+2,jj+1,kk) &
                        + w3*fine(ii+1,jj,kk) + w4*fine(ii+1,jj+2,kk)) / (w1+w2+w3+w4)
                else
                   fine(ii+1,jj+1,kk) = 0.d0
                end if
             end if

             if (ii+1 .le. flo(1) .and. kk+1 .le. flo(3)) then
                if (msk(ii+1,jj,kk+1) .ne. dirichlet) then
                   w1 = sum(sigx(ii     ,jj-1:jj,kk:kk+1))
                   w2 = sum(sigx(ii+1   ,jj-1:jj,kk:kk+1))
                   w3 = sum(sigz(ii:ii+1,jj-1:jj,kk     ))
                   w4 = sum(sigz(ii:ii+1,jj-1:jj,kk+1   ))
                   fine(ii+1,jj,kk+1) = (w1*fine(ii,jj,kk+1) + w2*fine(ii+2,jj,kk+1) &
                        + w3*fine(ii+1,jj,kk) + w4*fine(ii+1,jj,kk+2)) / (w1+w2+w3+w4)
                else
                   fine(ii+1,jj,kk+1) = 0.d0
                end if
             end if

             if (jj+1 .le. flo(2) .and. kk+1 .lt. flo(3)) then
                if (msk(ii,jj+1,kk+1) .ne. dirichlet) then
                   w1 = sum(sigy(ii-1:ii,jj     ,kk:kk+1))
                   w2 = sum(sigy(ii-1:ii,jj+1   ,kk:kk+1))
                   w3 = sum(sigz(ii-1:ii,jj:jj+1,kk     ))
                   w4 = sum(sigz(ii-1:ii,jj:jj+1,kk+1   ))
                   fine(ii,jj+1,kk+1) = (w1*fine(ii,jj,kk+1) + w2*fine(ii,jj+2,kk+1) &
                        + w3*fine(ii,jj+1,kk) + w4*fine(ii,jj+1,kk+2)) / (w1+w2+w3+w4)
                else
                   fine(ii,jj+1,kk+1) = 0.d0
                end if
             end if
          end do
       end do
    end do
    
    do       kk = flo(3)+1, fhi(3)-1, 2
       do    jj = flo(2)+1, fhi(2)-1, 2
          do ii = flo(1)+1, fhi(1)-1, 2
             w1 = sum(sigx(ii-1,jj-1:jj,kk-1:kk))
             w2 = sum(sigx(ii  ,jj-1:jj,kk-1:kk))
             w3 = sum(sigy(ii-1:ii,jj-1,kk-1:kk))
             w4 = sum(sigy(ii-1:ii,jj  ,kk-1:kk))
             w5 = sum(sigz(ii-1:ii,jj-1:jj,kk-1))
             w5 = sum(sigz(ii-1:ii,jj-1:jj,kk  ))
             fine(ii,jj,kk) = (w1*fine(ii-1,jj,kk) + w2*fine(ii+1,jj,kk) &
                  + w3*fine(ii,jj-1,kk) + w4*fine(ii,jj+1,kk) &
                  + w5*fine(ii,jj,kk-1) + w6*fine(ii,jj,kk+1)) &
                  / (w1+w2+w3+w4+w5+w6)
          end do
       end do
    end do

  end subroutine amrex_mlndlap_interpolation_ha


  subroutine amrex_mlndlap_interpolation_aa (clo, chi, fine, fflo, ffhi, crse, cflo, cfhi, &
       sig, sglo, sghi, msk, mlo, mhi, domlo, domhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_interpolation_aa')
    integer, dimension(3), intent(in) :: clo,chi,fflo,ffhi,cflo,cfhi,sglo,sghi, &
         mlo, mhi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(in   ) :: crse(cflo(1):cfhi(1),cflo(2):cfhi(2),cflo(3):cfhi(3))
    real(amrex_real), intent(inout) :: fine(fflo(1):ffhi(1),fflo(2):ffhi(2),fflo(3):ffhi(3))
    real(amrex_real), intent(in   ) :: sig (sglo(1):sghi(1),sglo(2):sghi(2),sglo(3):sghi(3))
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: flo(3), fhi(3), i, j, k, ii, jj, kk
    real(amrex_real) :: w1, w2, w3, w4, w5, w6

    flo = 2*clo
    fhi = 2*chi

    do k = clo(3), chi(3)
       kk = 2*k
       do j = clo(2), chi(2)
          jj = 2*j
          do i = clo(1), chi(1)
             ii = 2*i

             if (msk(ii,jj,kk) .ne. dirichlet) then
                fine(ii,jj,kk) = crse(i,j,k)
             else
                fine(ii,jj,kk) = 0.d0
             end if

             if (ii+1 .le. flo(1)) then
                if (msk(ii+1,jj,kk) .ne. dirichlet) then
                   w1 = sum(sig(ii  ,jj-1:jj,kk-1:kk))
                   w2 = sum(sig(ii+1,jj-1:jj,kk-1:kk))
                   fine(ii+1,jj,kk) = (w1*crse(i,j,k)+w2*crse(i+1,j,k))/(w1+w2)
                else
                   fine(ii+1,jj,kk) = 0.d0
                end if
             end if

             if (jj+1 .lt. flo(2)) then
                if (msk(ii,jj+1,kk) .ne. dirichlet) then
                   w1 = sum(sig(ii-1:ii,jj  ,kk-1:kk))
                   w2 = sum(sig(ii-1:ii,jj+1,kk-1:kk))
                   fine(ii,jj+1,kk) = (w1*crse(i,j,k)+w2*crse(i,j+1,k))/(w1+w2)
                else
                   fine(ii,jj+1,kk) = 0.d0
                end if
             end if

             if (kk+1 .lt. flo(3)) then
                if (msk(ii,jj,kk+1) .ne. dirichlet) then
                   w1 = sum(sig(ii-1:ii,jj-1:jj,kk  ))
                   w2 = sum(sig(ii-1:ii,jj-1:jj,kk+1))
                   fine(ii,jj,kk+1) = (w1*crse(i,j,k)+w2*crse(i,j,k+1))/(w1+w2)
                else
                   fine(ii,jj,kk+1) = 0.d0
                end if
             end if
          end do
       end do
    end do

    do k = clo(3), chi(3)
       kk = 2*k
       do j = clo(2), chi(2)
          jj = 2*j
          do i = clo(1), chi(1)
             ii = 2*i

             if (ii+1 .le. flo(1) .and. jj+1 .le. flo(2)) then
                if (msk(ii+1,jj+1,kk) .ne. dirichlet) then
                   w1 = sum(sig(ii     ,jj:jj+1,kk-1:kk))
                   w2 = sum(sig(ii+1   ,jj:jj+1,kk-1:kk))
                   w3 = sum(sig(ii:ii+1,jj     ,kk-1:kk))
                   w4 = sum(sig(ii:ii+1,jj+1   ,kk-1:kk))
                   fine(ii+1,jj+1,kk) = (w1*fine(ii,jj+1,kk) + w2*fine(ii+2,jj+1,kk) &
                        + w3*fine(ii+1,jj,kk) + w4*fine(ii+1,jj+2,kk)) / (w1+w2+w3+w4)
                else
                   fine(ii+1,jj+1,kk) = 0.d0
                end if
             end if

             if (ii+1 .le. flo(1) .and. kk+1 .le. flo(3)) then
                if (msk(ii+1,jj,kk+1) .ne. dirichlet) then
                   w1 = sum(sig(ii     ,jj-1:jj,kk:kk+1))
                   w2 = sum(sig(ii+1   ,jj-1:jj,kk:kk+1))
                   w3 = sum(sig(ii:ii+1,jj-1:jj,kk     ))
                   w4 = sum(sig(ii:ii+1,jj-1:jj,kk+1   ))
                   fine(ii+1,jj,kk+1) = (w1*fine(ii,jj,kk+1) + w2*fine(ii+2,jj,kk+1) &
                        + w3*fine(ii+1,jj,kk) + w4*fine(ii+1,jj,kk+2)) / (w1+w2+w3+w4)
                else
                   fine(ii+1,jj,kk+1) = 0.d0
                end if
             end if

             if (jj+1 .le. flo(2) .and. kk+1 .lt. flo(3)) then
                if (msk(ii,jj+1,kk+1) .ne. dirichlet) then
                   w1 = sum(sig(ii-1:ii,jj     ,kk:kk+1))
                   w2 = sum(sig(ii-1:ii,jj+1   ,kk:kk+1))
                   w3 = sum(sig(ii-1:ii,jj:jj+1,kk     ))
                   w4 = sum(sig(ii-1:ii,jj:jj+1,kk+1   ))
                   fine(ii,jj+1,kk+1) = (w1*fine(ii,jj,kk+1) + w2*fine(ii,jj+2,kk+1) &
                        + w3*fine(ii,jj+1,kk) + w4*fine(ii,jj+1,kk+2)) / (w1+w2+w3+w4)
                else
                   fine(ii,jj+1,kk+1) = 0.d0
                end if
             end if
          end do
       end do
    end do
    
    do       kk = flo(3)+1, fhi(3)-1, 2
       do    jj = flo(2)+1, fhi(2)-1, 2
          do ii = flo(1)+1, fhi(1)-1, 2
             w1 = sum(sig(ii-1,jj-1:jj,kk-1:kk))
             w2 = sum(sig(ii  ,jj-1:jj,kk-1:kk))
             w3 = sum(sig(ii-1:ii,jj-1,kk-1:kk))
             w4 = sum(sig(ii-1:ii,jj  ,kk-1:kk))
             w5 = sum(sig(ii-1:ii,jj-1:jj,kk-1))
             w5 = sum(sig(ii-1:ii,jj-1:jj,kk  ))
             fine(ii,jj,kk) = (w1*fine(ii-1,jj,kk) + w2*fine(ii+1,jj,kk) &
                  + w3*fine(ii,jj-1,kk) + w4*fine(ii,jj+1,kk) &
                  + w5*fine(ii,jj,kk-1) + w6*fine(ii,jj,kk+1)) &
                  / (w1+w2+w3+w4+w5+w6)
          end do
       end do
    end do

  end subroutine amrex_mlndlap_interpolation_aa


  subroutine amrex_mlndlap_divu (lo, hi, rhs, rlo, rhi, vel, vlo, vhi, msk, mlo, mhi, &
       dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_divu')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, vlo, vhi, mlo, mhi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k
    real(amrex_real) :: facx, facy, facz

    facx = 0.25d0*dxinv(1)
    facy = 0.25d0*dxinv(2)
    facz = 0.25d0*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                rhs(i,j,k) = facx*(-vel(i-1,j-1,k-1,1)+vel(i,j-1,k-1,1) &
                     &             -vel(i-1,j  ,k-1,1)+vel(i,j  ,k-1,1) &
                     &             -vel(i-1,j-1,k  ,1)+vel(i,j-1,k  ,1) &
                     &             -vel(i-1,j  ,k  ,1)+vel(i,j  ,k  ,1)) &
                     &     + facy*(-vel(i-1,j-1,k-1,2)-vel(i,j-1,k-1,2) &
                     &             +vel(i-1,j  ,k-1,2)+vel(i,j  ,k-1,2) &
                     &             -vel(i-1,j-1,k  ,2)-vel(i,j-1,k  ,2) &
                     &             +vel(i-1,j  ,k  ,2)+vel(i,j  ,k  ,2)) &
                     &     + facy*(-vel(i-1,j-1,k-1,3)-vel(i,j-1,k-1,3) &
                     &             -vel(i-1,j  ,k-1,3)-vel(i,j  ,k-1,3) &
                     &             +vel(i-1,j-1,k  ,3)+vel(i,j-1,k  ,3) &
                     &             +vel(i-1,j  ,k  ,3)+vel(i,j  ,k  ,3))
             else
                rhs(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

    if (lo(1) .eq. ndlo(1)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then 
          rhs(lo(1),lo(2):hi(2),lo(3):hi(3)) = 2.d0*rhs(lo(1),lo(2):hi(2),lo(3):hi(3))
       end if
    end if

    if (hi(1) .eq. ndhi(1)) then
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          rhs(hi(1),lo(2):hi(2),lo(3):hi(3)) = 2.d0*rhs(hi(1),lo(2):hi(2),lo(3):hi(3))
       end if
    end if

    if (lo(2) .eq. ndlo(2)) then
       if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          rhs(lo(1):hi(1),lo(2),lo(3):hi(3)) = 2.d0*rhs(lo(1):hi(1),lo(2),lo(3):hi(3))
       end if
    end if

    if (hi(2) .eq. ndhi(2)) then
       if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          rhs(lo(1):hi(1),hi(2),lo(3):hi(3)) = 2.d0*rhs(lo(1):hi(1),hi(2),lo(3):hi(3))
       end if
    end if

    if (lo(3) .eq. ndlo(3)) then
       if (bclo(3) .eq. amrex_lo_neumann .or. bclo(3) .eq. amrex_lo_inflow) then
          rhs(lo(1):hi(1),lo(2):hi(2),lo(3)) = 2.d0*rhs(lo(1):hi(1),lo(2):hi(2),lo(3))
       end if
    end if

    if (hi(3) .eq. ndhi(3)) then
       if (bchi(3) .eq. amrex_lo_neumann .or. bchi(3) .eq. amrex_lo_inflow) then
          rhs(lo(1):hi(1),lo(2):hi(2),hi(3)) = 2.d0*rhs(lo(1):hi(1),lo(2):hi(2),hi(3))
       end if
    end if

  end subroutine amrex_mlndlap_divu


  subroutine amrex_mlndlap_add_rhcc (lo, hi, rhs, rlo, rhi, rhcc, clo, chi, msk, mlo, mhi) &
       bind(c,name='amrex_mlndlap_add_rhcc')
    integer, dimension(3) :: lo, hi, rlo, rhi, clo, chi, mlo, mhi
    real(amrex_real), intent(inout) :: rhs (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(in   ) :: rhcc(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    integer,          intent(in   ) :: msk (mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                rhs(i,j,k) = rhs(i,j,k) + 0.125d0* &
                     (rhcc(i-1,j-1,k-1)+rhcc(i,j-1,k-1)+rhcc(i-1,j,k-1)+rhcc(i,j,k-1) &
                     +rhcc(i-1,j-1,k  )+rhcc(i,j-1,k  )+rhcc(i-1,j,k  )+rhcc(i,j,k  ))
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_add_rhcc


  subroutine amrex_mlndlap_mknewu (lo, hi, u, ulo, uhi, p, plo, phi, sig, slo, shi, dxinv) &
       bind(c,name='amrex_mlndlap_mknewu')
    integer, dimension(3), intent(in) :: lo, hi, ulo, uhi, plo, phi, slo, shi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) ::   u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)
    real(amrex_real), intent(in   ) ::   p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))

    integer :: i, j, k
    real(amrex_real) :: facx, facy, facz

    facx = 0.25d0*dxinv(1)
    facy = 0.25d0*dxinv(2)
    facz = 0.25d0*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             u(i,j,k,1) = u(i,j,k,1) - sig(i,j,k)*facx &
                  * (-p(i,j,k  )+p(i+1,j,k  )-p(i,j+1,k  )+p(i+1,j+1,k  ) &
                  &  -p(i,j,k+1)+p(i+1,j,k+1)-p(i,j+1,k+1)+p(i+1,j+1,k+1))
             u(i,j,k,2) = u(i,j,k,2) - sig(i,j,k)*facy &
                  * (-p(i,j,k  )-p(i+1,j,k  )+p(i,j+1,k  )+p(i+1,j+1,k  ) &
                  &  -p(i,j,k+1)-p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1))
             u(i,j,k,3) = u(i,j,k,3) - sig(i,j,k)*facy &
                  * (-p(i,j,k  )-p(i+1,j,k  )-p(i,j+1,k  )-p(i+1,j+1,k  ) &
                  &  +p(i,j,k+1)+p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1))
          end do
       end do
    end do

  end subroutine amrex_mlndlap_mknewu


  subroutine amrex_mlndlap_divu_fine_contrib (clo, chi, cglo, cghi, rhs, rlo, rhi, &
       vel, vlo, vhi, frh, flo, fhi, msk, mlo, mhi, dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_divu_fine_contrib')
    integer, dimension(3), intent(in) :: clo, chi, cglo, cghi, rlo, rhi, vlo, vhi, &
         flo, fhi, mlo, mhi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)
    real(amrex_real), intent(inout) :: frh(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer, dimension(3) :: lo, hi, glo, ghi
    integer :: i, j, k, ii, jj, kk, step
    real(amrex_real) :: facx, facy, facz
    real(amrex_real), parameter :: rfd = 0.125d0
    real(amrex_real), parameter :: chip = 0.5d0
    real(amrex_real), parameter :: chip2 = 0.25d0
    real(amrex_real), parameter :: chip3 = 0.125d0

    ! note that dxinv is fine dxinv
    facx = 0.25d0*dxinv(1)
    facy = 0.25d0*dxinv(2)
    facz = 0.25d0*dxinv(3)

    lo = 2*clo
    hi = 2*chi
    glo = 2*cglo
    ghi = 2*cghi

    do    kk = lo(3), hi(3)
       do jj = lo(2), hi(2)
          if (jj.eq.glo(2) .or. jj.eq.ghi(2) .or. kk.eq.glo(3) .or. kk.eq.ghi(3)) then
             step = 1
          else
             step = hi(1)-lo(1)
          end if
          do ii = lo(1), hi(1), step
             if (ii.eq.glo(1) .or. ii.eq.ghi(1) .or. step.eq.1) then
                frh(ii,jj,kk) = facx*(-vel(ii-1,jj-1,kk-1,1)+vel(ii,jj-1,kk-1,1) &
                     &                -vel(ii-1,jj  ,kk-1,1)+vel(ii,jj  ,kk-1,1) &
                     &                -vel(ii-1,jj-1,kk  ,1)+vel(ii,jj-1,kk  ,1) &
                     &                -vel(ii-1,jj  ,kk  ,1)+vel(ii,jj  ,kk  ,1)) &
                     &        + facy*(-vel(ii-1,jj-1,kk-1,2)-vel(ii,jj-1,kk-1,2) &
                     &                +vel(ii-1,jj  ,kk-1,2)+vel(ii,jj  ,kk-1,2) &
                     &                -vel(ii-1,jj-1,kk  ,2)-vel(ii,jj-1,kk  ,2) &
                     &                +vel(ii-1,jj  ,kk  ,2)+vel(ii,jj  ,kk  ,2)) &
                     &        + facy*(-vel(ii-1,jj-1,kk-1,3)-vel(ii,jj-1,kk-1,3) &
                     &                -vel(ii-1,jj  ,kk-1,3)-vel(ii,jj  ,kk-1,3) &
                     &                +vel(ii-1,jj-1,kk  ,3)+vel(ii,jj-1,kk  ,3) &
                     &                +vel(ii-1,jj  ,kk  ,3)+vel(ii,jj  ,kk  ,3))
             end if
          end do
       end do
    end do

    do k = clo(3), chi(3)
       kk = 2*k
       do j = clo(2), chi(2)
          jj = 2*j
          if (jj.eq.glo(2) .or. jj.eq.ghi(2) .or. kk.eq.glo(3) .or. kk.eq.ghi(3)) then
             step = 1
          else
             step = chi(1)-clo(1)
          end if
          do i = clo(1), chi(1), step
             ii = 2*i
             if (msk(ii,jj,kk) .eq. dirichlet) then
                rhs(i,j,k) = rhs(i,j,k) + rfd*(frh(ii,jj,kk) &
                     + chip*(frh(ii,jj,kk-1)+frh(ii,jj,kk+1) &
                     &      +frh(ii,jj-1,kk)+frh(ii,jj+1,kk) &
                     &      +frh(ii-1,jj,kk)+frh(ii+1,jj,kk)) &
                     + chip2*(frh(ii,jj-1,kk-1)+frh(ii,jj+1,kk-1)+frh(ii,jj-1,kk+1)+frh(ii,jj+1,kk+1) &
                     &       +frh(ii-1,jj,kk-1)+frh(ii+1,jj,kk-1)+frh(ii-1,jj,kk+1)+frh(ii+1,jj,kk+1) &
                     &       +frh(ii-1,jj-1,kk)+frh(ii+1,jj-1,kk)+frh(ii-1,jj+1,kk)+frh(ii+1,jj+1,kk)) &
                     + chip3*(frh(ii-1,jj-1,kk-1)+frh(ii+1,jj-1,kk-1) &
                     &       +frh(ii-1,jj+1,kk-1)+frh(ii+1,jj+1,kk-1) &
                     &       +frh(ii-1,jj-1,kk+1)+frh(ii+1,jj-1,kk+1) &
                     &       +frh(ii-1,jj+1,kk+1)+frh(ii+1,jj+1,kk+1)))
             end if
          end do
       end do
    end do

    ! xxxxx what do we do at physical boundaries?

  end subroutine amrex_mlndlap_divu_fine_contrib


  subroutine amrex_mlndlap_divu_cf_contrib (lo, hi,  rhs, rlo, rhi, vel, vlo, vhi, dmsk, mlo, mhi, &
       ndmsk, nmlo, nmhi, ccmsk, cmlo, cmhi, fc, clo, chi, dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_divu_cf_contrib')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, vlo, vhi, mlo, mhi, &
         nmlo, nmhi, cmlo, cmhi, clo, chi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: rhs  ( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) :: vel  ( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3),3)
    real(amrex_real), intent(in   ) :: fc   ( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))
    integer         , intent(in   ) :: dmsk ( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3))
    integer         , intent(in   ) :: ndmsk(nmlo(1):nmhi(1),nmlo(2):nmhi(2),nmlo(3):nmhi(3))
    integer         , intent(in   ) :: ccmsk(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3))

    integer :: i,j,k
    real(amrex_real) :: facx, facy, facz

    facx = 0.25d0*dxinv(1)
    facy = 0.25d0*dxinv(2)
    facz = 0.25d0*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (dmsk(i,j,k) .ne. dirichlet) then
                if (ndmsk(i,j,k) .eq. crse_fine_node) then
                   rhs(i,j,k) = fc(i,j,k)
                   if (ccmsk(i-1,j-1,k-1) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) - facx*vel(i-1,j-1,k-1,1) &
                           &                  - facy*vel(i-1,j-1,k-1,2) &
                           &                  - facz*vel(i-1,j-1,k-1,3)
                   end if
                   if (ccmsk(i,j-1,k-1) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) + facx*vel(i  ,j-1,k-1,1) &
                           &                  - facy*vel(i  ,j-1,k-1,2) &
                           &                  - facz*vel(i  ,j-1,k-1,3)
                   end if
                   if (ccmsk(i-1,j,k-1) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) - facx*vel(i-1,j  ,k-1,1) &
                           &                  + facy*vel(i-1,j  ,k-1,2) &
                           &                  - facz*vel(i-1,j  ,k-1,3)
                   end if
                   if (ccmsk(i,j,k-1) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) + facx*vel(i  ,j  ,k-1,1) &
                           &                  + facy*vel(i  ,j  ,k-1,2) &
                           &                  - facz*vel(i  ,j  ,k-1,3)
                   end if
                   if (ccmsk(i-1,j-1,k) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) - facx*vel(i-1,j-1,k  ,1) &
                           &                  - facy*vel(i-1,j-1,k  ,2) &
                        &                     + facz*vel(i-1,j-1,k  ,3)
                   end if
                   if (ccmsk(i,j-1,k) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) + facx*vel(i  ,j-1,k  ,1) &
                           &                  - facy*vel(i  ,j-1,k  ,2) &
                           &                  + facz*vel(i  ,j-1,k  ,3)
                   end if
                   if (ccmsk(i-1,j,k) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) - facx*vel(i-1,j  ,k  ,1) &
                           &                  + facy*vel(i-1,j  ,k  ,2) &
                           &                  + facz*vel(i-1,j  ,k  ,3)
                   end if
                   if (ccmsk(i,j,k) .eq. crse_cell) then
                      rhs(i,j,k) = rhs(i,j,k) + facx*vel(i  ,j  ,k  ,1) &
                           &                  + facy*vel(i  ,j  ,k  ,2) &
                           &                  + facz*vel(i  ,j  ,k  ,3)
                   end if
                end if
             end if
          end do
       end do
    end do

    ! xxxxx what do we do at physical boundaries?

  end subroutine amrex_mlndlap_divu_cf_contrib


  subroutine amrex_mlndlap_crse_resid (lo, hi, resid, rslo, rshi, rhs, rhlo, rhhi, msk, mlo, mhi) &
       bind(c, name='amrex_mlndlap_crse_resid')
    integer, dimension(3), intent(in) :: lo, hi, rslo, rshi, rhlo, rhhi, mlo, mhi
    real(amrex_real), intent(inout) :: resid(rslo(1):rshi(1),rslo(2):rshi(2),rslo(3):rshi(3))
    real(amrex_real), intent(in   ) :: rhs  (rhlo(1):rhhi(1),rhlo(2):rhhi(2),rhlo(3):rhhi(3))
    integer         , intent(in   ) :: msk  ( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (any(msk(i-1:i,j-1:j,k-1:k).eq.0) .and. any(msk(i-1:i,j-1:j,k-1:k).eq.1)) then
                resid(i,j,k) = rhs(i,j,k) - resid(i,j,k)
             else
                resid(i,j,k) = 0.d0
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_crse_resid


  subroutine amrex_mlndlap_res_fine_contrib (clo, chi, cglo, cghi, f, flo, fhi, &
       x, xlo, xhi, sig, slo, shi, Ax, alo, ahi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_res_fine_contrib')
    integer, dimension(3), intent(in) :: clo, chi, cglo, cghi, flo, fhi, xlo, xhi, &
         slo, shi, alo, ahi, mlo, mhi
    real(amrex_real), intent(inout) :: f  (flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(amrex_real), intent(in   ) :: x  (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(amrex_real), intent(inout) :: Ax (alo(1):ahi(1),alo(2):ahi(2),alo(3):ahi(3))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    real(amrex_real), intent(in) :: dxinv(3)

    integer, dimension(3) :: lo, hi, glo, ghi
    integer :: i, j, k, ii, jj, kk, step
    real(amrex_real) :: facx, facy, facz, fxyz, fmx2y2z, f2xmy2z, f2x2ymz
    real(amrex_real) :: f4xm2ym2z, fm2x4ym2z, fm2xm2y4z
    real(amrex_real), parameter :: rfd = 0.125d0
    real(amrex_real), parameter :: chip = 0.5d0
    real(amrex_real), parameter :: chip2 = 0.25d0
    real(amrex_real), parameter :: chip3 = 0.125d0

    facx = (1.d0/36.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/36.d0)*dxinv(2)*dxinv(2)
    facz = (1.d0/36.d0)*dxinv(3)*dxinv(3)
    fxyz = facx + facy + facz
    fmx2y2z = -facx + 2.d0*facy + 2.d0*facz
    f2xmy2z = 2.d0*facx - facy + 2.d0*facz
    f2x2ymz = 2.d0*facx + 2.d0*facy - facz
    f4xm2ym2z = 4.d0*facx - 2.d0*facy - 2.d0*facz
    fm2x4ym2z = -2.d0*facx + 4.d0*facy - 2.d0*facz
    fm2xm2y4z = -2.d0*facx - 2.d0*facy + 4.d0*facz

    lo = 2*clo
    hi = 2*chi
    glo = 2*cglo
    ghi = 2*cghi

    do    kk = lo(3), hi(3)
       do jj = lo(2), hi(2)
          if (jj.eq.glo(2) .or. jj.eq.ghi(2) .or. kk.eq.glo(3) .or. kk.eq.ghi(3)) then
             step = 1
          else
             step = hi(1)-lo(1)
          end if
          do ii = lo(1), hi(1), step
             if (ii.eq.glo(1) .or. ii.eq.ghi(1) .or. step.eq.1) then
                Ax(ii,jj,kk) = x(ii,jj,kk)*(-4.d0)*fxyz* &
                     (sig(ii-1,jj-1,kk-1)+sig(ii,jj-1,kk-1)+sig(ii-1,jj,kk-1)+sig(ii,jj,kk-1) &
                     +sig(ii-1,jj-1,kk  )+sig(ii,jj-1,kk  )+sig(ii-1,jj,kk  )+sig(ii,jj,kk  )) &
                     !
                     + fxyz*(x(ii-1,jj-1,kk-1)*sig(ii-1,jj-1,kk-1) &
                     &     + x(ii+1,jj-1,kk-1)*sig(ii  ,jj-1,kk-1) &
                     &     + x(ii-1,jj+1,kk-1)*sig(ii-1,jj  ,kk-1) &
                     &     + x(ii+1,jj+1,kk-1)*sig(ii  ,jj  ,kk-1) &
                     &     + x(ii-1,jj-1,kk+1)*sig(ii-1,jj-1,kk  ) &
                     &     + x(ii+1,jj-1,kk+1)*sig(ii  ,jj-1,kk  ) &
                     &     + x(ii-1,jj+1,kk+1)*sig(ii-1,jj  ,kk  ) &
                     &     + x(ii+1,jj+1,kk+1)*sig(ii  ,jj  ,kk  )) &
                     !
                     + fmx2y2z*(x(ii  ,jj-1,kk-1)*(sig(ii-1,jj-1,kk-1)+sig(ii,jj-1,kk-1)) &
                     &        + x(ii  ,jj+1,kk-1)*(sig(ii-1,jj  ,kk-1)+sig(ii,jj  ,kk-1)) &
                     &        + x(ii  ,jj-1,kk+1)*(sig(ii-1,jj-1,kk  )+sig(ii,jj-1,kk  )) &
                     &        + x(ii  ,jj+1,kk+1)*(sig(ii-1,jj  ,kk  )+sig(ii,jj  ,kk  ))) &
                     !
                     + f2xmy2z*(x(ii-1,jj  ,kk-1)*(sig(ii-1,jj-1,kk-1)+sig(ii-1,jj,kk-1)) &
                     &        + x(ii+1,jj  ,kk-1)*(sig(ii  ,jj-1,kk-1)+sig(ii  ,jj,kk-1)) &
                     &        + x(ii-1,jj  ,kk+1)*(sig(ii-1,jj-1,kk  )+sig(ii-1,jj,kk  )) &
                     &        + x(ii+1,jj  ,kk+1)*(sig(ii  ,jj-1,kk  )+sig(ii  ,jj,kk  ))) &
                     !
                     + f2x2ymz*(x(ii-1,jj-1,kk  )*(sig(ii-1,jj-1,kk-1)+sig(ii-1,jj-1,kk)) &
                     &        + x(ii+1,jj-1,kk  )*(sig(ii  ,jj-1,kk-1)+sig(ii  ,jj-1,kk)) &
                     &        + x(ii-1,jj+1,kk  )*(sig(ii-1,jj  ,kk-1)+sig(ii-1,jj  ,kk)) &
                     &        + x(ii+1,jj+1,kk  )*(sig(ii  ,jj  ,kk-1)+sig(ii  ,jj  ,kk))) &
                     !
                     + f4xm2ym2z*(x(ii-1,jj,kk)*(sig(ii-1,jj-1,kk-1)+sig(ii-1,jj,kk-1)+sig(ii-1,jj-1,kk)+sig(ii-1,jj,kk)) &
                     &          + x(ii  ,jj,kk)*(sig(ii  ,jj-1,kk-1)+sig(ii  ,jj,kk-1)+sig(ii  ,jj-1,kk)+sig(ii  ,jj,kk))) &
                     + fm2x4ym2z*(x(ii,jj-1,kk)*(sig(ii-1,jj-1,kk-1)+sig(ii,jj-1,kk-1)+sig(ii-1,jj-1,kk)+sig(ii,jj-1,kk)) &
                     &          + x(ii,jj  ,kk)*(sig(ii-1,jj  ,kk-1)+sig(ii,jj  ,kk-1)+sig(ii-1,jj  ,kk)+sig(ii,jj  ,kk))) &
                     + fm2xm2y4z*(x(ii,jj,kk-1)*(sig(ii-1,jj-1,kk-1)+sig(ii,jj-1,kk-1)+sig(ii-1,jj,kk-1)+sig(ii,jj,kk-1)) &
                     &          + x(ii,jj,kk  )*(sig(ii-1,jj-1,kk  )+sig(ii,jj-1,kk  )+sig(ii-1,jj,kk  )+sig(ii,jj,kk  )))

             end if
          end do
       end do
    end do

    do k = clo(3), chi(3)
       kk = 2*k
       do j = clo(2), chi(2)
          jj = 2*j
          if (jj.eq.glo(2) .or. jj.eq.ghi(2) .or. kk.eq.glo(3) .or. kk.eq.ghi(3)) then
             step = 1
          else
             step = chi(1)-clo(1)
          end if
          do i = clo(1), chi(1), step
             ii = 2*i
             if (msk(ii,jj,kk) .eq. dirichlet) then
                f(i,j,k) = f(i,j,k) + rfd*(Ax(ii,jj,kk) &
                     + chip*(Ax(ii,jj,kk-1)+Ax(ii,jj,kk+1) &
                     &      +Ax(ii,jj-1,kk)+Ax(ii,jj+1,kk) &
                     &      +Ax(ii-1,jj,kk)+Ax(ii+1,jj,kk)) &
                     + chip2*(Ax(ii,jj-1,kk-1)+Ax(ii,jj+1,kk-1)+Ax(ii,jj-1,kk+1)+Ax(ii,jj+1,kk+1) &
                     &       +Ax(ii-1,jj,kk-1)+Ax(ii+1,jj,kk-1)+Ax(ii-1,jj,kk+1)+Ax(ii+1,jj,kk+1) &
                     &       +Ax(ii-1,jj-1,kk)+Ax(ii+1,jj-1,kk)+Ax(ii-1,jj+1,kk)+Ax(ii+1,jj+1,kk)) &
                     + chip3*(Ax(ii-1,jj-1,kk-1)+Ax(ii+1,jj-1,kk-1) &
                     &       +Ax(ii-1,jj+1,kk-1)+Ax(ii+1,jj+1,kk-1) &
                     &       +Ax(ii-1,jj-1,kk+1)+Ax(ii+1,jj-1,kk+1) &
                     &       +Ax(ii-1,jj+1,kk+1)+Ax(ii+1,jj+1,kk+1)))
             end if
          end do
       end do
    end do

    ! xxxxx what do we do at physical boundaries?

  end subroutine amrex_mlndlap_res_fine_contrib


  subroutine amrex_mlndlap_res_cf_contrib (lo, hi, res, rlo, rhi, phi, phlo, phhi, &
       rhs, rhlo, rhhi, sig, slo, shi, dmsk, mlo, mhi, ndmsk, nmlo, nmhi, ccmsk, cmlo, cmhi, &
       fc, clo, chi, dxinv) &
       bind(c,name='amrex_mlndlap_res_cf_contrib')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, phlo, phhi, rhlo, rhhi, slo, shi, &
         mlo, mhi, nmlo, nmhi, cmlo, cmhi, clo, chi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: res( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) :: phi(phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))
    real(amrex_real), intent(in   ) :: rhs(rhlo(1):rhhi(1),rhlo(2):rhhi(2),rhlo(3):rhhi(3))
    real(amrex_real), intent(in   ) :: sig( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    real(amrex_real), intent(inout) :: fc ( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))
    integer, intent(in) ::  dmsk( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3))
    integer, intent(in) :: ndmsk(nmlo(1):nmhi(1),nmlo(2):nmhi(2),nmlo(3):nmhi(3))
    integer, intent(in) :: ccmsk(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3))

    integer :: i,j,k
    real(amrex_real) :: Ax, facx, facy, facz

    facx = (1.d0/36.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/36.d0)*dxinv(2)*dxinv(2)
    facz = (1.d0/36.d0)*dxinv(3)*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (dmsk(i,j,k) .ne. dirichlet) then
                if (ndmsk(i,j,k) .eq. crse_fine_node) then
                   Ax = 0.d0
                   if (1.d0-ccmsk(i-1,j-1,k-1) .eq. crse_cell) then
                      Ax = Ax + sig(i-1,j-1,k-1)*(facx*(4.d0*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                           +2.d0*(phi(i-1,j-1,k  )-phi(i  ,j-1,k  )) &
                           &                           +2.d0*(phi(i-1,j  ,k-1)-phi(i  ,j  ,k-1)) &
                           &                           +     (phi(i-1,j-1,k-1)-phi(i  ,j-1,k-1))) &
                           &                    + facy*(4.d0*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  )) &
                           &                           +2.d0*(phi(i-1,j-1,k  )-phi(i-1,j  ,k  )) &
                           &                           +2.d0*(phi(i  ,j-1,k-1)-phi(i  ,j  ,k-1)) &
                           &                           +     (phi(i-1,j-1,k-1)-phi(i-1,j  ,k-1))) &
                           &                    + facz*(4.d0*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  )) &
                           &                           +2.d0*(phi(i-1,j  ,k-1)-phi(i-1,j  ,k  )) &
                           &                           +2.d0*(phi(i  ,j-1,k-1)-phi(i  ,j-1,k  )) &
                           &                           +     (phi(i-1,j-1,k-1)-phi(i-1,j-1,k  ))))
                   end if
                   if (1.d0-ccmsk(i,j-1,k-1) .eq. crse_cell) then
                      Ax = Ax + sig(i,j-1,k-1)*(facx*(4.d0*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i+1,j-1,k  )-phi(i  ,j-1,k  )) &
                           &                         +2.d0*(phi(i+1,j  ,k-1)-phi(i  ,j  ,k-1)) &
                           &                         +     (phi(i+1,j-1,k-1)-phi(i  ,j-1,k-1))) &
                           &                  + facy*(4.d0*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i+1,j-1,k  )-phi(i+1,j  ,k  )) &
                           &                         +2.d0*(phi(i  ,j-1,k-1)-phi(i  ,j  ,k-1)) &
                           &                         +     (phi(i+1,j-1,k-1)-phi(i+1,j  ,k-1))) &
                           &                  + facz*(4.d0*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i+1,j  ,k-1)-phi(i+1,j  ,k  )) &
                           &                         +2.d0*(phi(i  ,j-1,k-1)-phi(i  ,j-1,k  )) &
                           &                         +     (phi(i+1,j-1,k-1)-phi(i+1,j-1,k  ))))
                   end if
                   if (1.d0-ccmsk(i-1,j,k-1) .eq. crse_cell) then
                      Ax = Ax + sig(i-1,j,k-1)*(facx*(4.d0*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i-1,j+1,k  )-phi(i  ,j+1,k  )) &
                           &                         +2.d0*(phi(i-1,j  ,k-1)-phi(i  ,j  ,k-1)) &
                           &                         +     (phi(i-1,j+1,k-1)-phi(i  ,j+1,k-1))) &
                           &                  + facy*(4.d0*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i-1,j+1,k  )-phi(i-1,j  ,k  )) &
                           &                         +2.d0*(phi(i  ,j+1,k-1)-phi(i  ,j  ,k-1)) &
                           &                         +     (phi(i-1,j+1,k-1)-phi(i-1,j  ,k-1))) &
                           &                  + facz*(4.d0*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i-1,j  ,k-1)-phi(i-1,j  ,k  )) &
                           &                         +2.d0*(phi(i  ,j+1,k-1)-phi(i  ,j+1,k  )) &
                           &                         +     (phi(i-1,j+1,k-1)-phi(i-1,j+1,k  ))))
                   end if
                   if (1.d0-ccmsk(i,j,k-1) .eq. crse_cell) then
                      Ax = Ax + sig(i,j,k-1)*(facx*(4.d0*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i+1,j+1,k  )-phi(i  ,j+1,k  )) &
                           &                       +2.d0*(phi(i+1,j  ,k-1)-phi(i  ,j  ,k-1)) &
                           &                       +     (phi(i+1,j+1,k-1)-phi(i  ,j+1,k-1))) &
                           &                + facy*(4.d0*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i+1,j+1,k  )-phi(i+1,j  ,k  )) &
                           &                       +2.d0*(phi(i  ,j+1,k-1)-phi(i  ,j  ,k-1)) &
                           &                       +     (phi(i+1,j+1,k-1)-phi(i+1,j  ,k-1))) &
                           &                + facz*(4.d0*(phi(i  ,j  ,k-1)-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i+1,j  ,k-1)-phi(i+1,j  ,k  )) &
                           &                       +2.d0*(phi(i  ,j+1,k-1)-phi(i  ,j+1,k  )) &
                           &                       +     (phi(i+1,j+1,k-1)-phi(i+1,j+1,k  ))))
                   end if
                   if (1.d0-ccmsk(i-1,j-1,k) .eq. crse_cell) then
                      Ax = Ax + sig(i-1,j-1,k)*(facx*(4.d0*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i-1,j-1,k  )-phi(i  ,j-1,k  )) &
                           &                         +2.d0*(phi(i-1,j  ,k+1)-phi(i  ,j  ,k+1)) &
                           &                         +     (phi(i-1,j-1,k+1)-phi(i  ,j-1,k+1))) &
                           &                  + facy*(4.d0*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i-1,j-1,k  )-phi(i-1,j  ,k  )) &
                           &                         +2.d0*(phi(i  ,j-1,k+1)-phi(i  ,j  ,k+1)) &
                           &                         +     (phi(i-1,j-1,k+1)-phi(i-1,j  ,k+1))) &
                           &                  + facz*(4.d0*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  )) &
                           &                         +2.d0*(phi(i-1,j  ,k+1)-phi(i-1,j  ,k  )) &
                           &                         +2.d0*(phi(i  ,j-1,k+1)-phi(i  ,j-1,k  )) &
                           &                         +     (phi(i-1,j-1,k+1)-phi(i-1,j-1,k  ))))
                   end if
                   if (1.d0-ccmsk(i,j-1,k) .eq. crse_cell) then
                      Ax = Ax + sig(i,j-1,k)*(facx*(4.d0*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i+1,j-1,k  )-phi(i  ,j-1,k  )) &
                           &                       +2.d0*(phi(i+1,j  ,k+1)-phi(i  ,j  ,k+1)) &
                           &                       +     (phi(i+1,j-1,k+1)-phi(i  ,j-1,k+1))) &
                           &                + facy*(4.d0*(phi(i  ,j-1,k  )-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i+1,j-1,k  )-phi(i+1,j  ,k  )) &
                           &                       +2.d0*(phi(i  ,j-1,k+1)-phi(i  ,j  ,k+1)) &
                           &                       +     (phi(i+1,j-1,k+1)-phi(i+1,j  ,k+1))) &
                           &                + facz*(4.d0*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i+1,j  ,k+1)-phi(i+1,j  ,k  )) &
                           &                       +2.d0*(phi(i  ,j-1,k+1)-phi(i  ,j-1,k  )) &
                           &                       +     (phi(i+1,j-1,k+1)-phi(i+1,j-1,k  ))))
                   end if
                   if (1.d0-ccmsk(i-1,j,k) .eq. crse_cell) then
                      Ax = Ax + sig(i-1,j,k)*(facx*(4.d0*(phi(i-1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i-1,j+1,k  )-phi(i  ,j+1,k  )) &
                           &                       +2.d0*(phi(i-1,j  ,k+1)-phi(i  ,j  ,k+1)) &
                           &                       +     (phi(i-1,j+1,k+1)-phi(i  ,j+1,k+1))) &
                           &                + facy*(4.d0*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i-1,j+1,k  )-phi(i-1,j  ,k  )) &
                           &                       +2.d0*(phi(i  ,j+1,k+1)-phi(i  ,j  ,k+1)) &
                           &                       +     (phi(i-1,j+1,k+1)-phi(i-1,j  ,k+1))) &
                           &                + facz*(4.d0*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  )) &
                           &                       +2.d0*(phi(i-1,j  ,k+1)-phi(i-1,j  ,k  )) &
                           &                       +2.d0*(phi(i  ,j+1,k+1)-phi(i  ,j+1,k  )) &
                           &                       +     (phi(i-1,j+1,k+1)-phi(i-1,j+1,k  ))))
                   end if
                   if (1.d0-ccmsk(i,j,k) .eq. crse_cell) then
                      Ax = Ax + sig(i,j,k)*(facx*(4.d0*(phi(i+1,j  ,k  )-phi(i  ,j  ,k  )) &
                           &                     +2.d0*(phi(i+1,j+1,k  )-phi(i  ,j+1,k  )) &
                           &                     +2.d0*(phi(i+1,j  ,k+1)-phi(i  ,j  ,k+1)) &
                           &                     +     (phi(i+1,j+1,k+1)-phi(i  ,j+1,k+1))) &
                           &              + facy*(4.d0*(phi(i  ,j+1,k  )-phi(i  ,j  ,k  )) &
                           &                     +2.d0*(phi(i+1,j+1,k  )-phi(i+1,j  ,k  )) &
                           &                     +2.d0*(phi(i  ,j+1,k+1)-phi(i  ,j  ,k+1)) &
                           &                     +     (phi(i+1,j+1,k+1)-phi(i+1,j  ,k+1))) &
                           &              + facz*(4.d0*(phi(i  ,j  ,k+1)-phi(i  ,j  ,k  )) &
                           &                     +2.d0*(phi(i+1,j  ,k+1)-phi(i+1,j  ,k  )) &
                           &                     +2.d0*(phi(i  ,j+1,k+1)-phi(i  ,j+1,k  )) &
                           &                     +     (phi(i+1,j+1,k+1)-phi(i+1,j+1,k  ))))
                   end if
                   res(i,j,k) = rhs(i,j,k) - Ax - fc(i,j,k)
                end if
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_res_cf_contrib


  subroutine amrex_mlndlap_zero_fine (lo, hi, phi, dlo, dhi, msk, mlo, mhi, fine_flag) &
       bind(c, name='amrex_mlndlap_zero_fine')
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, mlo, mhi
    real(amrex_real), intent(inout) :: phi(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    integer, intent(in) :: fine_flag

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ! Testing if the node is covered by a fine level in computing
             ! coarse sync residual
             if (all(msk(i-1:i,j-1:j,k-1:k).eq.fine_flag)) then
                phi(i,j,k) = 0.d0
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_zero_fine

end module amrex_mlnodelap_3d_module
