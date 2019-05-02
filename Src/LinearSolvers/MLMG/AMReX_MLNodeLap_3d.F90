
module amrex_mlnodelap_3d_module

  use amrex_error_module
  use amrex_fort_module, only : amrex_real
  use amrex_constants_module
  use amrex_lo_bctypes_module, only : amrex_lo_dirichlet, amrex_lo_neumann, amrex_lo_inflow, amrex_lo_periodic
  implicit none

  ! external dirichlet at physical boundary or internal dirichlet at crse/fine boundary
  integer, parameter :: dirichlet = 1

  integer, parameter :: crse_cell = 0
  integer, parameter :: fine_cell = 1
  integer, parameter :: crse_node = 0
  integer, parameter :: crse_fine_node = 1
  integer, parameter :: fine_node = 2

  real(amrex_real), private, parameter :: eps = 1.d-100

  integer, private, parameter :: ist_000 = 1
  integer, private, parameter :: ist_p00 = 2
  integer, private, parameter :: ist_0p0 = 3
  integer, private, parameter :: ist_00p = 4
  integer, private, parameter :: ist_pp0 = 5
  integer, private, parameter :: ist_p0p = 6
  integer, private, parameter :: ist_0pp = 7
  integer, private, parameter :: ist_ppp = 8
  integer, private, parameter :: ist_inv = 9
  integer, private, parameter :: n_sten = 9

#ifdef AMREX_USE_EB
  integer, private, parameter :: i_S_x     = 1
  integer, private, parameter :: i_S_y     = 2
  integer, private, parameter :: i_S_z     = 3
  integer, private, parameter :: i_S_x2    = 4
  integer, private, parameter :: i_S_y2    = 5
  integer, private, parameter :: i_S_z2    = 6
  integer, private, parameter :: i_S_x_y   = 7
  integer, private, parameter :: i_S_x_z   = 8
  integer, private, parameter :: i_S_y_z   = 9
  integer, private, parameter :: i_S_x2_y  = 10
  integer, private, parameter :: i_S_x2_z  = 11
  integer, private, parameter :: i_S_x_y2  = 12
  integer, private, parameter :: i_S_y2_z  = 13
  integer, private, parameter :: i_S_x_z2  = 14
  integer, private, parameter :: i_S_y_z2  = 15
  integer, private, parameter :: i_S_x2_y2 = 16
  integer, private, parameter :: i_S_x2_z2 = 17
  integer, private, parameter :: i_S_y2_z2 = 18
  integer, private, parameter :: n_Sintg = 18

  integer, private, parameter :: i_c_xmym = 1
  integer, private, parameter :: i_c_xmyb = 2
  integer, private, parameter :: i_c_xmyp = 3
  integer, private, parameter :: i_c_xbym = 4
  integer, private, parameter :: i_c_xbyb = 5
  integer, private, parameter :: i_c_xbyp = 6
  integer, private, parameter :: i_c_xpym = 7
  integer, private, parameter :: i_c_xpyb = 8
  integer, private, parameter :: i_c_xpyp = 9
  integer, private, parameter :: i_c_xmzm = 10
  integer, private, parameter :: i_c_xmzb = 11
  integer, private, parameter :: i_c_xmzp = 12
  integer, private, parameter :: i_c_xbzm = 13
  integer, private, parameter :: i_c_xbzb = 14
  integer, private, parameter :: i_c_xbzp = 15
  integer, private, parameter :: i_c_xpzm = 16
  integer, private, parameter :: i_c_xpzb = 17
  integer, private, parameter :: i_c_xpzp = 18
  integer, private, parameter :: i_c_ymzm = 19
  integer, private, parameter :: i_c_ymzb = 20
  integer, private, parameter :: i_c_ymzp = 21
  integer, private, parameter :: i_c_ybzm = 22
  integer, private, parameter :: i_c_ybzb = 23
  integer, private, parameter :: i_c_ybzp = 24
  integer, private, parameter :: i_c_ypzm = 25
  integer, private, parameter :: i_c_ypzb = 26
  integer, private, parameter :: i_c_ypzp = 27
  integer, private, parameter :: n_conn = 27
#endif

  private
  public :: &
       ! masks
       amrex_mlndlap_set_nodal_mask, amrex_mlndlap_set_dirichlet_mask, &
       amrex_mlndlap_fixup_res_mask, amrex_mlndlap_set_dot_mask, &
       amrex_mlndlap_any_fine_sync_cells, &
       ! coeffs
       amrex_mlndlap_avgdown_coeff, amrex_mlndlap_fillbc_cc, amrex_mlndlap_fillbc_cc_i, &
       ! bc
       amrex_mlndlap_applybc, amrex_mlndlap_impose_neumann_bc, &
       ! operator
       amrex_mlndlap_adotx_ha, amrex_mlndlap_adotx_aa, &
       amrex_mlndlap_normalize_ha, amrex_mlndlap_normalize_aa, &
       amrex_mlndlap_jacobi_ha, amrex_mlndlap_jacobi_aa, &
       amrex_mlndlap_gauss_seidel_ha, amrex_mlndlap_gauss_seidel_aa, &
       ! restriction
       amrex_mlndlap_restriction, &
       ! interpolation
       amrex_mlndlap_interpolation_ha, amrex_mlndlap_interpolation_aa, &
       ! rhs & u
       amrex_mlndlap_divu, amrex_mlndlap_rhcc, amrex_mlndlap_mknewu, &
       amrex_mlndlap_divu_fine_contrib, amrex_mlndlap_divu_cf_contrib, &
       amrex_mlndlap_rhcc_fine_contrib, amrex_mlndlap_rhcc_crse_contrib, &
       ! residual
       amrex_mlndlap_crse_resid, &
       amrex_mlndlap_res_fine_contrib, amrex_mlndlap_res_cf_contrib, &
       ! sync residual
       amrex_mlndlap_zero_fine

  ! RAP
  public:: amrex_mlndlap_set_stencil, amrex_mlndlap_set_stencil_s0, &
       amrex_mlndlap_adotx_sten, amrex_mlndlap_normalize_sten, &
       amrex_mlndlap_gauss_seidel_sten, amrex_mlndlap_jacobi_sten, &
       amrex_mlndlap_interpolation_rap, &
       amrex_mlndlap_restriction_rap, &
       amrex_mlndlap_stencil_rap

#ifdef AMREX_USE_EB
  public:: amrex_mlndlap_set_integral, amrex_mlndlap_set_integral_eb, &
       amrex_mlndlap_set_connection, amrex_mlndlap_set_stencil_eb, &
       amrex_mlndlap_divu_eb, amrex_mlndlap_mknewu_eb
#endif

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


  subroutine amrex_mlndlap_set_dot_mask (lo, hi, dmsk, dlo, dhi, omsk, olo, ohi, &
       domlo, domhi, bclo, bchi) bind(c,name='amrex_mlndlap_set_dot_mask')
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, olo, ohi, domlo, domhi, bclo, bchi
    real(amrex_real), intent(inout) :: dmsk(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    integer         , intent(in   ) :: omsk(olo(1):ohi(1),olo(2):ohi(2),olo(3):ohi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dmsk(i,j,k) = omsk(i,j,k)
          end do
       end do
    end do

    if (lo(1) .eq. domlo(1)) then
       if (bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) then
          dmsk(lo(1),lo(2):hi(2),lo(3):hi(3)) = 0.5d0*dmsk(lo(1),lo(2):hi(2),lo(3):hi(3))
       end if
    end if

    if (hi(1) .eq. domhi(1)) then
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
          dmsk(hi(1),lo(2):hi(2),lo(3):hi(3)) = 0.5d0*dmsk(hi(1),lo(2):hi(2),lo(3):hi(3))
       end if
    end if

    if (lo(2) .eq. domlo(2)) then
       if (bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) then
          dmsk(lo(1):hi(1),lo(2),lo(3):hi(3)) = 0.5d0*dmsk(lo(1):hi(1),lo(2),lo(3):hi(3))
       end if
    end if

    if (hi(2) .eq. domhi(2)) then
       if (bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) then
          dmsk(lo(1):hi(1),hi(2),lo(3):hi(3)) = 0.5d0*dmsk(lo(1):hi(1),hi(2),lo(3):hi(3))
       end if
    end if

    if (lo(3) .eq. domlo(3)) then
       if (bclo(3) .eq. amrex_lo_neumann .or. bclo(3) .eq. amrex_lo_inflow) then
          dmsk(lo(1):hi(1),lo(2):hi(2),lo(3)) = 0.5d0*dmsk(lo(1):hi(1),lo(2):hi(2),lo(3))
       end if
    end if

    if (hi(3) .eq. domhi(3)) then
       if (bchi(3) .eq. amrex_lo_neumann .or. bchi(3) .eq. amrex_lo_inflow) then
          dmsk(lo(1):hi(1),lo(2):hi(2),hi(3)) = 0.5d0*dmsk(lo(1):hi(1),lo(2):hi(2),hi(3))
       end if
    end if

  end subroutine amrex_mlndlap_set_dot_mask


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
                cr = 0.25d0*(fine(2*i,2*j+1,2*k  )+fine(2*i+1,2*j+1,2*k  ) &
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


  subroutine amrex_mlndlap_fillbc_cc_i (sigma, slo, shi, dlo, dhi, bclo, bchi) &
       bind(c, name='amrex_mlndlap_fillbc_cc_i')
    integer, dimension(3), intent(in) :: slo, shi, dlo, dhi, bclo, bchi
    integer, intent(inout) :: sigma(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    
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

  end subroutine amrex_mlndlap_fillbc_cc_i


  subroutine amrex_mlndlap_applybc (phi, hlo, hhi, dlo, dhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_applybc')
    integer, dimension(3) :: hlo, hhi, dlo, dhi, bclo, bchi
    real(amrex_real), intent(inout) :: phi(hlo(1):hhi(1),hlo(2):hhi(2),hlo(3):hhi(3))

    integer :: ilo, ihi, jlo, jhi, klo, khi

    ilo = max(dlo(1), hlo(1))
    ihi = min(dhi(1), hhi(1))
    jlo = max(dlo(2), hlo(2))
    jhi = min(dhi(2), hhi(2))
    klo = max(dlo(3), hlo(3))
    khi = min(dhi(3), hhi(3))

    ! neumann

    if ((bclo(1) .eq. amrex_lo_neumann .or. bclo(1) .eq. amrex_lo_inflow) &
         .and. hlo(1) .lt. dlo(1)) then
       phi(dlo(1)-1,jlo:jhi,klo:khi) = phi(dlo(1)+1,jlo:jhi,klo:khi)
    end if

    if ((bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) &
         .and. hhi(1) .gt. dhi(1)) then
       phi(dhi(1)+1,jlo:jhi,klo:khi) = phi(dhi(1)-1,jlo:jhi,klo:khi)
    end if

    if ((bclo(2) .eq. amrex_lo_neumann .or. bclo(2) .eq. amrex_lo_inflow) &
         .and. hlo(2) .lt. dlo(2)) then
       phi(ilo:ihi,dlo(2)-1,klo:khi) = phi(ilo:ihi,dlo(2)+1,klo:khi)
    end if

    if ((bchi(2) .eq. amrex_lo_neumann .or. bchi(2) .eq. amrex_lo_inflow) &
         .and. hhi(2) .gt. dhi(2)) then
       phi(ilo:ihi,dhi(2)+1,klo:khi) = phi(ilo:ihi,dhi(2)-1,klo:khi)
    end if

    if ((bclo(3) .eq. amrex_lo_neumann .or. bclo(3) .eq. amrex_lo_inflow) &
         .and. hlo(3) .lt. dlo(3)) then
       phi(ilo:ihi,jlo:jhi,dlo(3)-1) = phi(ilo:ihi,jlo:jhi,dlo(3)+1)
    end if
    
    if ((bchi(3) .eq. amrex_lo_neumann .or. bchi(3) .eq. amrex_lo_inflow) &
         .and. hhi(3) .gt. dhi(3)) then
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
       else if (bclo(3) .eq. amrex_lo_neumann .or. bclo(3) .eq. amrex_lo_inflow) then
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
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
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
       if (bchi(1) .eq. amrex_lo_neumann .or. bchi(1) .eq. amrex_lo_inflow) then
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


  subroutine amrex_mlndlap_impose_neumann_bc (lo, hi, rhs, rlo, rhi, ndlo, ndhi, bclo, bchi) &
       bind(c, name='amrex_mlndlap_impose_neumann_bc')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))

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

  end subroutine amrex_mlndlap_impose_neumann_bc


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
    integer         , intent(in   ) ::msk( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3))

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
                     + 2.d0*x(i+1,j,k)*(2.d0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1)+sx(i  ,j-1,k)+sx(i  ,j,k)) &
                     &                      -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1)+sy(i  ,j-1,k)+sy(i  ,j,k)) &
                     &                      -facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)+sz(i  ,j-1,k)+sz(i  ,j,k))) &
                     + 2.d0*x(i,j-1,k)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j-1,k)+sx(i,j-1,k)) &
                     &                 +2.d0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j-1,k)+sy(i,j-1,k)) &
                     &                      -facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j-1,k)+sz(i,j-1,k))) &
                     + 2.d0*x(i,j+1,k)*(    -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1)+sx(i-1,j  ,k)+sx(i,j  ,k)) &
                     &                 +2.d0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1)+sy(i-1,j  ,k)+sy(i,j  ,k)) &
                     &                      -facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)+sz(i-1,j  ,k)+sz(i,j  ,k))) &
                     + 2.d0*x(i,j,k-1)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)) &
                     &                      -facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)) &
                     &                 +2.d0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1))) &
                     + 2.d0*x(i,j,k+1)*(    -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  )) &
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
                     &          + x(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k))) &
                     + fm2x4ym2z*(x(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k)) &
                     &          + x(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k))) &
                     + fm2xm2y4z*(x(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)) &
                     &          + x(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )))
             else
                y(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_adotx_aa


  subroutine amrex_mlndlap_normalize_ha (lo, hi, x, xlo, xhi, &
       sx, sxlo, sxhi, sy, sylo, syhi, sz, szlo, szhi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_normalize_ha')
    integer, dimension(3), intent(in) :: lo, hi, xlo, xhi, sxlo, sxhi, &
         sylo, syhi, szlo, szhi, mlo, mhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) ::  x( xlo(1): xhi(1), xlo(2): xhi(2), xlo(3): xhi(3))
    real(amrex_real), intent(in   ) :: sx(sxlo(1):sxhi(1),sxlo(2):sxhi(2),sxlo(3):sxhi(3))
    real(amrex_real), intent(in   ) :: sy(sylo(1):syhi(1),sylo(2):syhi(2),sylo(3):syhi(3))
    real(amrex_real), intent(in   ) :: sz(szlo(1):szhi(1),szlo(2):szhi(2),szlo(3):szhi(3))
    integer         , intent(in   ) ::msk( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3))

    integer :: i,j,k
    real(amrex_real) :: facx, facy, facz

    facx = (1.d0/36.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/36.d0)*dxinv(2)*dxinv(2)
    facz = (1.d0/36.d0)*dxinv(3)*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                x(i,j,k) = x(i,j,k)/((-4.d0)*(facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1) &
                     &                             +sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  )) &
                     &                       +facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1) &
                     &                             +sy(i-1,j-1,k  )+sy(i,j-1,k  )+sy(i-1,j,k  )+sy(i,j,k  )) &
                     &                       +facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1) &
                     &                             +sz(i-1,j-1,k  )+sz(i,j-1,k  )+sz(i-1,j,k  )+sz(i,j,k  ))))

             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_normalize_ha


  subroutine amrex_mlndlap_normalize_aa (lo, hi, x, xlo, xhi, sig, slo, shi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_normalize_aa')
    integer, dimension(3), intent(in) :: lo, hi, xlo, xhi, slo, shi, mlo, mhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) ::   x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k
    real(amrex_real) :: facx, facy, facz, fxyz
    
    facx = (1.d0/36.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/36.d0)*dxinv(2)*dxinv(2)
    facz = (1.d0/36.d0)*dxinv(3)*dxinv(3)
    fxyz = facx + facy + facz

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                x(i,j,k) = x(i,j,k) / ((-4.d0)*fxyz* &
                     (sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1) &
                     +sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )))
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_normalize_aa


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
                     + sol(i-1,j-1,k-1)*(facx*sx(i-1,j-1,k-1) &
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
                     + 2.d0*sol(i+1,j,k)*(2.d0*facx*(sx(i  ,j-1,k-1)+sx(i  ,j,k-1)+sx(i  ,j-1,k)+sx(i  ,j,k)) &
                     &                        -facy*(sy(i  ,j-1,k-1)+sy(i  ,j,k-1)+sy(i  ,j-1,k)+sy(i  ,j,k)) &
                     &                        -facz*(sz(i  ,j-1,k-1)+sz(i  ,j,k-1)+sz(i  ,j-1,k)+sz(i  ,j,k))) &
                     + 2.d0*sol(i,j-1,k)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j-1,k)+sx(i,j-1,k)) &
                     &                   +2.d0*facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j-1,k)+sy(i,j-1,k)) &
                     &                        -facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j-1,k)+sz(i,j-1,k))) &
                     + 2.d0*sol(i,j+1,k)*(    -facx*(sx(i-1,j  ,k-1)+sx(i,j  ,k-1)+sx(i-1,j  ,k)+sx(i,j  ,k)) &
                     &                   +2.d0*facy*(sy(i-1,j  ,k-1)+sy(i,j  ,k-1)+sy(i-1,j  ,k)+sy(i,j  ,k)) &
                     &                        -facz*(sz(i-1,j  ,k-1)+sz(i,j  ,k-1)+sz(i-1,j  ,k)+sz(i,j  ,k))) &
                     + 2.d0*sol(i,j,k-1)*(    -facx*(sx(i-1,j-1,k-1)+sx(i,j-1,k-1)+sx(i-1,j,k-1)+sx(i,j,k-1)) &
                     &                        -facy*(sy(i-1,j-1,k-1)+sy(i,j-1,k-1)+sy(i-1,j,k-1)+sy(i,j,k-1)) &
                     &                   +2.d0*facz*(sz(i-1,j-1,k-1)+sz(i,j-1,k-1)+sz(i-1,j,k-1)+sz(i,j,k-1))) &
                     + 2.d0*sol(i,j,k+1)*(    -facx*(sx(i-1,j-1,k  )+sx(i,j-1,k  )+sx(i-1,j,k  )+sx(i,j,k  )) &
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
                     &          + sol(i+1,j,k)*(sig(i  ,j-1,k-1)+sig(i  ,j,k-1)+sig(i  ,j-1,k)+sig(i  ,j,k))) &
                     + fm2x4ym2z*(sol(i,j-1,k)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j-1,k)+sig(i,j-1,k)) &
                     &          + sol(i,j+1,k)*(sig(i-1,j  ,k-1)+sig(i,j  ,k-1)+sig(i-1,j  ,k)+sig(i,j  ,k))) &
                     + fm2xm2y4z*(sol(i,j,k-1)*(sig(i-1,j-1,k-1)+sig(i,j-1,k-1)+sig(i-1,j,k-1)+sig(i,j,k-1)) &
                     &          + sol(i,j,k+1)*(sig(i-1,j-1,k  )+sig(i,j-1,k  )+sig(i-1,j,k  )+sig(i,j,k  )))

                sol(i,j,k) = sol(i,j,k) + (rhs(i,j,k) - Ax) / s0
             else
                sol(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_gauss_seidel_aa


  subroutine amrex_mlndlap_restriction (lo, hi, crse, clo, chi, fine, flo, fhi, msk, mlo, mhi) &
       bind(c,name='amrex_mlndlap_restriction')
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, flo, fhi, mlo, mhi
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

             if (ii+1 .lt. fhi(1)) then
                if (msk(ii+1,jj,kk) .ne. dirichlet) then
                   w1 = sum(sigx(ii  ,jj-1:jj,kk-1:kk))
                   w2 = sum(sigx(ii+1,jj-1:jj,kk-1:kk))
                   fine(ii+1,jj,kk) = (w1*crse(i,j,k)+w2*crse(i+1,j,k))/(w1+w2)
                else
                   fine(ii+1,jj,kk) = 0.d0
                end if
             end if

             if (jj+1 .lt. fhi(2)) then
                if (msk(ii,jj+1,kk) .ne. dirichlet) then
                   w1 = sum(sigy(ii-1:ii,jj  ,kk-1:kk))
                   w2 = sum(sigy(ii-1:ii,jj+1,kk-1:kk))
                   fine(ii,jj+1,kk) = (w1*crse(i,j,k)+w2*crse(i,j+1,k))/(w1+w2)
                else
                   fine(ii,jj+1,kk) = 0.d0
                end if
             end if

             if (kk+1 .lt. fhi(3)) then
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

             if (ii+1 .lt. fhi(1) .and. jj+1 .lt. fhi(2)) then
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

             if (ii+1 .lt. fhi(1) .and. kk+1 .lt. fhi(3)) then
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

             if (jj+1 .lt. fhi(2) .and. kk+1 .lt. fhi(3)) then
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
             w6 = sum(sigz(ii-1:ii,jj-1:jj,kk  ))
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

             if (ii+1 .lt. fhi(1)) then
                if (msk(ii+1,jj,kk) .ne. dirichlet) then
                   w1 = sum(sig(ii  ,jj-1:jj,kk-1:kk))
                   w2 = sum(sig(ii+1,jj-1:jj,kk-1:kk))
                   fine(ii+1,jj,kk) = (w1*crse(i,j,k)+w2*crse(i+1,j,k))/(w1+w2)
                else
                   fine(ii+1,jj,kk) = 0.d0
                end if
             end if

             if (jj+1 .lt. fhi(2)) then
                if (msk(ii,jj+1,kk) .ne. dirichlet) then
                   w1 = sum(sig(ii-1:ii,jj  ,kk-1:kk))
                   w2 = sum(sig(ii-1:ii,jj+1,kk-1:kk))
                   fine(ii,jj+1,kk) = (w1*crse(i,j,k)+w2*crse(i,j+1,k))/(w1+w2)
                else
                   fine(ii,jj+1,kk) = 0.d0
                end if
             end if

             if (kk+1 .lt. fhi(3)) then
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

             if (ii+1 .lt. fhi(1) .and. jj+1 .lt. fhi(2)) then
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

             if (ii+1 .lt. fhi(1) .and. kk+1 .lt. fhi(3)) then
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

             if (jj+1 .lt. fhi(2) .and. kk+1 .lt. fhi(3)) then
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
             w6 = sum(sig(ii-1:ii,jj-1:jj,kk  ))
             fine(ii,jj,kk) = (w1*fine(ii-1,jj,kk) + w2*fine(ii+1,jj,kk) &
                  + w3*fine(ii,jj-1,kk) + w4*fine(ii,jj+1,kk) &
                  + w5*fine(ii,jj,kk-1) + w6*fine(ii,jj,kk+1)) &
                  / (w1+w2+w3+w4+w5+w6)
          end do
       end do
    end do

  end subroutine amrex_mlndlap_interpolation_aa


  subroutine amrex_mlndlap_divu (lo, hi, rhs, rlo, rhi, vel, vlo, vhi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_divu')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, vlo, vhi, mlo, mhi
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
                     &     + facz*(-vel(i-1,j-1,k-1,3)-vel(i,j-1,k-1,3) &
                     &             -vel(i-1,j  ,k-1,3)-vel(i,j  ,k-1,3) &
                     &             +vel(i-1,j-1,k  ,3)+vel(i,j-1,k  ,3) &
                     &             +vel(i-1,j  ,k  ,3)+vel(i,j  ,k  ,3))
             else
                rhs(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_divu


  subroutine amrex_mlndlap_rhcc (lo, hi, rhs, rlo, rhi, rhcc, clo, chi, msk, mlo, mhi) &
       bind(c,name='amrex_mlndlap_rhcc')
    integer, dimension(3) :: lo, hi, rlo, rhi, clo, chi, mlo, mhi
    real(amrex_real), intent(inout) :: rhs (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(in   ) :: rhcc(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    integer,          intent(in   ) :: msk (mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                rhs(i,j,k) = 0.125d0* &
                     (rhcc(i-1,j-1,k-1)+rhcc(i,j-1,k-1)+rhcc(i-1,j,k-1)+rhcc(i,j,k-1) &
                     +rhcc(i-1,j-1,k  )+rhcc(i,j-1,k  )+rhcc(i-1,j,k  )+rhcc(i,j,k  ))
             else
                rhs(i,j,k) = 0.d0
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_rhcc


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
             u(i,j,k,3) = u(i,j,k,3) - sig(i,j,k)*facz &
                  * (-p(i,j,k  )-p(i+1,j,k  )-p(i,j+1,k  )-p(i+1,j+1,k  ) &
                  &  +p(i,j,k+1)+p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1))
          end do
       end do
    end do

  end subroutine amrex_mlndlap_mknewu


  subroutine amrex_mlndlap_divu_fine_contrib (clo, chi, cglo, cghi, rhs, rlo, rhi, &
       vel, vlo, vhi, frh, flo, fhi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_divu_fine_contrib')
    integer, dimension(3), intent(in) :: clo, chi, cglo, cghi, rlo, rhi, vlo, vhi, &
         flo, fhi, mlo, mhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)
    real(amrex_real), intent(inout) :: frh(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer, dimension(3) :: lo, hi, glo, ghi, gtlo, gthi
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

    gtlo = max(lo-1,glo)
    gthi = min(hi+1,ghi)

    do    kk = gtlo(3), gthi(3)
       do jj = gtlo(2), gthi(2)
          if (jj.eq.glo(2) .or. jj.eq.ghi(2) .or. kk.eq.glo(3) .or. kk.eq.ghi(3)) then
             step = 1
          else
             step = gthi(1)-gtlo(1)
          end if
          do ii = gtlo(1), gthi(1), step
             if (ii.eq.glo(1) .or. ii.eq.ghi(1) .or. step.eq.1) then
                frh(ii,jj,kk) = facx*(-vel(ii-1,jj-1,kk-1,1)+vel(ii,jj-1,kk-1,1) &
                     &                -vel(ii-1,jj  ,kk-1,1)+vel(ii,jj  ,kk-1,1) &
                     &                -vel(ii-1,jj-1,kk  ,1)+vel(ii,jj-1,kk  ,1) &
                     &                -vel(ii-1,jj  ,kk  ,1)+vel(ii,jj  ,kk  ,1)) &
                     &        + facy*(-vel(ii-1,jj-1,kk-1,2)-vel(ii,jj-1,kk-1,2) &
                     &                +vel(ii-1,jj  ,kk-1,2)+vel(ii,jj  ,kk-1,2) &
                     &                -vel(ii-1,jj-1,kk  ,2)-vel(ii,jj-1,kk  ,2) &
                     &                +vel(ii-1,jj  ,kk  ,2)+vel(ii,jj  ,kk  ,2)) &
                     &        + facz*(-vel(ii-1,jj-1,kk-1,3)-vel(ii,jj-1,kk-1,3) &
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

                   if (i .eq. ndlo(1) .and. &
                        (    bclo(1) .eq. amrex_lo_neumann &
                        .or. bclo(1) .eq. amrex_lo_inflow)) then
                      rhs(i,j,k) = 2.d0*rhs(i,j,k)
                   else if (i.eq. ndhi(1) .and. &
                        (    bchi(1) .eq. amrex_lo_neumann &
                        .or. bchi(1) .eq. amrex_lo_inflow)) then
                      rhs(i,j,k) = 2.d0*rhs(i,j,k)
                   end if
                   
                   if (j .eq. ndlo(2) .and. &
                        (    bclo(2) .eq. amrex_lo_neumann &
                        .or. bclo(2) .eq. amrex_lo_inflow)) then
                      rhs(i,j,k) = 2.d0*rhs(i,j,k)                   
                   else if (j .eq. ndhi(2) .and. &
                        (    bchi(2) .eq. amrex_lo_neumann &
                        .or. bchi(2) .eq. amrex_lo_inflow)) then
                      rhs(i,j,k) = 2.d0*rhs(i,j,k)
                   end if

                   if (k .eq. ndlo(3) .and. &
                        (    bclo(3) .eq. amrex_lo_neumann &
                        .or. bclo(3) .eq. amrex_lo_inflow)) then
                      rhs(i,j,k) = 2.d0*rhs(i,j,k)                   
                   else if (k .eq. ndhi(3) .and. &
                        (    bchi(3) .eq. amrex_lo_neumann &
                        .or. bchi(3) .eq. amrex_lo_inflow)) then
                      rhs(i,j,k) = 2.d0*rhs(i,j,k)
                   end if

                end if
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_divu_cf_contrib


  subroutine amrex_mlndlap_rhcc_fine_contrib (clo, chi, cglo, cghi, rhs, rlo, rhi, &
       cc, cclo, cchi, msk, mlo, mhi) bind(c,name='amrex_mlndlap_rhcc_fine_contrib')
    integer, dimension(3), intent(in) :: clo, chi, cglo, cghi, rlo, rhi, cclo, cchi, mlo, mhi
    real(amrex_real), intent(inout) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) :: cc (cclo(1):cchi(1),cclo(2):cchi(2),cclo(3):cchi(3))
    integer         , intent(in   ) :: msk( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3))

    integer, dimension(3) :: lo, hi, glo, ghi
    integer :: i, j, k, ii, jj, kk, step, ioff, joff, koff
    real(amrex_real), parameter :: fac(-2:1) = [1.d0/8.d0, 3.d0/8.d0, 3.d0/8.d0, 1.d0/8.d0]

    lo = 2*clo
    hi = 2*chi
    glo = 2*cglo
    ghi = 2*cghi

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
                do koff = -2, 1
                   do joff = -2, 1
                      do ioff = -2, 1
                         rhs(i,j,k) = rhs(i,j,k) + cc(ii+ioff,jj+joff,kk+koff) &
                              * fac(ioff)*fac(joff)*fac(koff)
                      end do
                   end do
                end do
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_rhcc_fine_contrib


  subroutine amrex_mlndlap_rhcc_crse_contrib (lo, hi, crhs, rlo, rhi, rhcc, clo, chi, &
       dmsk, mlo, mhi, ndmsk, nmlo, nmhi, ccmsk, cmlo, cmhi) &
       bind(c,name='amrex_mlndlap_rhcc_crse_contrib')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, clo, chi, mlo, mhi, &
         nmlo, nmhi, cmlo, cmhi
    real(amrex_real), intent(inout) ::  crhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) ::  rhcc( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3))
    integer         , intent(in   ) ::  dmsk( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3))
    integer         , intent(in   ) :: ndmsk(nmlo(1):nmhi(1),nmlo(2):nmhi(2),nmlo(3):nmhi(3))
    integer         , intent(in   ) :: ccmsk(cmlo(1):cmhi(1),cmlo(2):cmhi(2),cmlo(3):cmhi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (dmsk(i,j,k) .ne. dirichlet) then
                if (ndmsk(i,j,k) .eq. crse_fine_node) then
                   crhs(i,j,k) = crhs(i,j,k) + 0.125d0 * &
                        ( (1.d0-ccmsk(i-1,j-1,k-1)) * rhcc(i-1,j-1,k-1) &
                        + (1.d0-ccmsk(i  ,j-1,k-1)) * rhcc(i  ,j-1,k-1) &
                        + (1.d0-ccmsk(i-1,j  ,k-1)) * rhcc(i-1,j  ,k-1) &
                        + (1.d0-ccmsk(i  ,j  ,k-1)) * rhcc(i  ,j  ,k-1) &
                        + (1.d0-ccmsk(i-1,j-1,k  )) * rhcc(i-1,j-1,k  ) &
                        + (1.d0-ccmsk(i  ,j-1,k  )) * rhcc(i  ,j-1,k  ) &
                        + (1.d0-ccmsk(i-1,j  ,k  )) * rhcc(i-1,j  ,k  ) &
                        + (1.d0-ccmsk(i  ,j  ,k  )) * rhcc(i  ,j  ,k  ) )
                end if
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_rhcc_crse_contrib


  subroutine amrex_mlndlap_crse_resid (lo, hi, resid, rslo, rshi, rhs, rhlo, rhhi, msk, mlo, mhi, &
       ndlo, ndhi, bclo, bchi) bind(c, name='amrex_mlndlap_crse_resid')
    integer, dimension(3), intent(in) :: lo, hi, rslo, rshi, rhlo, rhhi, mlo, mhi, ndlo, ndhi, bclo, bchi
    real(amrex_real), intent(inout) :: resid(rslo(1):rshi(1),rslo(2):rshi(2),rslo(3):rshi(3))
    real(amrex_real), intent(in   ) :: rhs  (rhlo(1):rhhi(1),rhlo(2):rhhi(2),rhlo(3):rhhi(3))
    integer         , intent(in   ) :: msk  ( mlo(1): mhi(1), mlo(2): mhi(2), mlo(3): mhi(3))

    integer :: i,j,k
    real(amrex_real) :: fac

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (any(msk(i-1:i,j-1:j,k-1:k).eq.0) .and. any(msk(i-1:i,j-1:j,k-1:k).eq.1)) then

                fac = 1.d0

                if (i .eq. ndlo(1) .and. &
                     (    bclo(1) .eq. amrex_lo_neumann &
                     .or. bclo(1) .eq. amrex_lo_inflow)) then
                   fac = 2.d0*fac
                else if (i.eq. ndhi(1) .and. &
                     (    bchi(1) .eq. amrex_lo_neumann &
                     .or. bchi(1) .eq. amrex_lo_inflow)) then
                   fac = 2.d0*fac
                end if
                
                if (j .eq. ndlo(2) .and. &
                     (    bclo(2) .eq. amrex_lo_neumann &
                     .or. bclo(2) .eq. amrex_lo_inflow)) then
                   fac = 2.d0*fac                   
                else if (j .eq. ndhi(2) .and. &
                     (    bchi(2) .eq. amrex_lo_neumann &
                     .or. bchi(2) .eq. amrex_lo_inflow)) then
                   fac = 2.d0*fac
                end if
                
                if (k .eq. ndlo(3) .and. &
                     (    bclo(3) .eq. amrex_lo_neumann &
                     .or. bclo(3) .eq. amrex_lo_inflow)) then
                   fac = 2.d0*fac                   
                else if (k .eq. ndhi(3) .and. &
                     (    bchi(3) .eq. amrex_lo_neumann &
                     .or. bchi(3) .eq. amrex_lo_inflow)) then
                   fac = 2.d0*fac
                end if

                resid(i,j,k) = (rhs(i,j,k) - resid(i,j,k)) * fac
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

    integer, dimension(3) :: lo, hi, glo, ghi, gtlo, gthi
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

    gtlo = max(lo-1,glo)
    gthi = min(hi+1,ghi)

    do    kk = gtlo(3), gthi(3)
       do jj = gtlo(2), gthi(2)
          if (jj.eq.glo(2) .or. jj.eq.ghi(2) .or. kk.eq.glo(3) .or. kk.eq.ghi(3)) then
             step = 1
          else
             step = gthi(1)-gtlo(1)
          end if
          do ii = gtlo(1), gthi(1), step
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
                     &          + x(ii+1,jj,kk)*(sig(ii  ,jj-1,kk-1)+sig(ii  ,jj,kk-1)+sig(ii  ,jj-1,kk)+sig(ii  ,jj,kk))) &
                     + fm2x4ym2z*(x(ii,jj-1,kk)*(sig(ii-1,jj-1,kk-1)+sig(ii,jj-1,kk-1)+sig(ii-1,jj-1,kk)+sig(ii,jj-1,kk)) &
                     &          + x(ii,jj+1,kk)*(sig(ii-1,jj  ,kk-1)+sig(ii,jj  ,kk-1)+sig(ii-1,jj  ,kk)+sig(ii,jj  ,kk))) &
                     + fm2xm2y4z*(x(ii,jj,kk-1)*(sig(ii-1,jj-1,kk-1)+sig(ii,jj-1,kk-1)+sig(ii-1,jj,kk-1)+sig(ii,jj,kk-1)) &
                     &          + x(ii,jj,kk+1)*(sig(ii-1,jj-1,kk  )+sig(ii,jj-1,kk  )+sig(ii-1,jj,kk  )+sig(ii,jj,kk  )))

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

  end subroutine amrex_mlndlap_res_fine_contrib


  subroutine amrex_mlndlap_res_cf_contrib (lo, hi, res, rlo, rhi, phi, phlo, phhi, &
       rhs, rhlo, rhhi, sig, slo, shi, dmsk, mlo, mhi, ndmsk, nmlo, nmhi, ccmsk, cmlo, cmhi, &
       fc, clo, chi, dxinv, ndlo, ndhi, bclo, bchi) &
       bind(c,name='amrex_mlndlap_res_cf_contrib')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, phlo, phhi, rhlo, rhhi, slo, shi, &
         mlo, mhi, nmlo, nmhi, cmlo, cmhi, clo, chi, ndlo, ndhi, bclo, bchi
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
    real(amrex_real) :: Ax, Axf, facx, facy, facz

    facx = (1.d0/36.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/36.d0)*dxinv(2)*dxinv(2)
    facz = (1.d0/36.d0)*dxinv(3)*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (dmsk(i,j,k) .ne. dirichlet) then
                if (ndmsk(i,j,k) .eq. crse_fine_node) then
                   Ax = 0.d0
                   if (ccmsk(i-1,j-1,k-1) .eq. crse_cell) then
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
                   if (ccmsk(i,j-1,k-1) .eq. crse_cell) then
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
                   if (ccmsk(i-1,j,k-1) .eq. crse_cell) then
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
                   if (ccmsk(i,j,k-1) .eq. crse_cell) then
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
                   if (ccmsk(i-1,j-1,k) .eq. crse_cell) then
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
                   if (ccmsk(i,j-1,k) .eq. crse_cell) then
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
                   if (ccmsk(i-1,j,k) .eq. crse_cell) then
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
                   if (ccmsk(i,j,k) .eq. crse_cell) then
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

                   Axf = fc(i,j,k)

                   if (i .eq. ndlo(1) .and. &
                        (    bclo(1) .eq. amrex_lo_neumann &
                        .or. bclo(1) .eq. amrex_lo_inflow)) then
                      Axf = 2.d0*Axf
                   else if (i.eq. ndhi(1) .and. &
                        (    bchi(1) .eq. amrex_lo_neumann &
                        .or. bchi(1) .eq. amrex_lo_inflow)) then
                      Axf = 2.d0*Axf
                   end if

                   if (j .eq. ndlo(2) .and. &
                        (    bclo(2) .eq. amrex_lo_neumann &
                        .or. bclo(2) .eq. amrex_lo_inflow)) then
                      Axf = 2.d0*Axf
                   else if (j .eq. ndhi(2) .and. &
                        (    bchi(2) .eq. amrex_lo_neumann &
                        .or. bchi(2) .eq. amrex_lo_inflow)) then
                      Axf = 2.d0*Axf
                   end if

                   if (k .eq. ndlo(3) .and. &
                        (    bclo(3) .eq. amrex_lo_neumann &
                        .or. bclo(3) .eq. amrex_lo_inflow)) then
                      Axf = 2.d0*Axf
                   else if (k .eq. ndhi(3) .and. &
                        (    bchi(3) .eq. amrex_lo_neumann &
                        .or. bchi(3) .eq. amrex_lo_inflow)) then
                      Axf = 2.d0*Axf
                   end if

                   res(i,j,k) = rhs(i,j,k) - (Ax + Axf)
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


  subroutine amrex_mlndlap_set_stencil (lo, hi, sten, tlo, thi, sig, glo, ghi, dxinv) &
       bind(c,name='amrex_mlndlap_set_stencil')
    integer, dimension(3), intent(in) :: lo, hi, tlo, thi, glo, ghi
    real(amrex_real), intent(inout) :: sten(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),n_sten)
    real(amrex_real), intent(in   ) ::  sig(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3))
    real(amrex_real), intent(in) :: dxinv(3)

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

             ! i+1,j,k
             sten(i,j,k,ist_p00) = f4xm2ym2z * (sig(i,j-1,k-1)+sig(i,j,k-1)+sig(i,j-1,k)+sig(i,j,k))
             ! i-1,j,k: sten(i-1,j,k,ist_p00)

             ! i,j+1,k
             sten(i,j,k,ist_0p0) = fm2x4ym2z * (sig(i-1,j,k-1)+sig(i,j,k-1)+sig(i-1,j,k)+sig(i,j,k))
             ! i,j-1,k: sten(i,j-1,k,ist_0p0)

             ! i,j,k+1
             sten(i,j,k,ist_00p) = fm2xm2y4z * (sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k))
             ! i,j,k-1: sten(i,j,k-1,ist_00p)

             ! i+1,j+1,k
             sten(i,j,k,ist_pp0) = f2x2ymz * (sig(i,j,k-1)+sig(i,j,k))
             ! i-1,j-1,k: sten(i-1,j-1,k,ist_pp0)
             ! i+1,j-1,k: sten(i  ,j-1,k,ist_pp0)
             ! i-1,j+1,k: sten(i-1,j  ,k,ist_pp0)
             
             ! i+1,j,k+1
             sten(i,j,k,ist_p0p) = f2xmy2z * (sig(i,j-1,k)+sig(i,j,k))
             ! i-1,j,k-1: sten(i-1,j,k-1,ist_p0p)
             ! i+1,j,k-1: sten(i  ,j,k-1,ist_p0p)
             ! i-1,j,k+1: sten(i-1,j,k  ,ist_p0p)

             ! i,j+1,k+1
             sten(i,j,k,ist_0pp) = fmx2y2z * (sig(i-1,j,k)+sig(i,j,k))
             ! i,j-1,k-1: sten(i,j-1,k-1,ist_0pp)
             ! i,j+1,k-1: sten(i,j  ,k-1,ist_0pp)
             ! i,j-1,k+1: sten(i,j-1,k  ,ist_0pp)

             ! i+1,j+1,k+1
             sten(i,j,k,ist_ppp) = fxyz * sig(i,j,k)
             ! i-1,j-1,k-1: sten(i-1,j-1,k-1,ist_ppp)
             ! i+1,j-1,k-1: sten(i  ,j-1,k-1,ist_ppp)
             ! i-1,j+1,k-1: sten(i-1,j  ,k-1,ist_ppp)
             ! i+1,j+1,k-1: sten(i  ,j  ,k-1,ist_ppp)
             ! i-1,j-1,k+1: sten(i-1,j-1,k  ,ist_ppp)
             ! i+1,j-1,k+1: sten(i  ,j-1,k  ,ist_ppp)
             ! i-1,j+1,k+1: sten(i-1,j  ,k  ,ist_ppp)
          end do
       end do
    end do
  end subroutine amrex_mlndlap_set_stencil


  subroutine amrex_mlndlap_set_stencil_s0 (lo, hi, sten, tlo, thi) &
       bind(c,name='amrex_mlndlap_set_stencil_s0')
    integer, dimension(3), intent(in) :: lo, hi, tlo, thi
    real(amrex_real), intent(inout) ::  sten(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),n_sten)

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             sten(i,j,k,ist_000) = -( &
                    sten(i-1,j,k,ist_p00) + sten(i,j,k,ist_p00) &
                  + sten(i,j-1,k,ist_0p0) + sten(i,j,k,ist_0p0) &
                  + sten(i,j,k-1,ist_00p) + sten(i,j,k,ist_00p) &
                  + sten(i-1,j-1,k,ist_pp0) + sten(i,j-1,k,ist_pp0) &
                  + sten(i-1,j,k,ist_pp0) + sten(i,j,k,ist_pp0) &
                  + sten(i-1,j,k-1,ist_p0p) + sten(i,j,k-1,ist_p0p) &
                  + sten(i-1,j,k,ist_p0p) + sten(i,j,k,ist_p0p) &
                  + sten(i,j-1,k-1,ist_0pp) + sten(i,j,k-1,ist_0pp) &
                  + sten(i,j-1,k,ist_0pp) + sten(i,j,k,ist_0pp) &
                  + sten(i-1,j-1,k-1,ist_ppp) + sten(i,j-1,k-1,ist_ppp) &
                  + sten(i-1,j,k-1,ist_ppp) + sten(i,j,k-1,ist_ppp) &
                  + sten(i-1,j-1,k,ist_ppp) + sten(i,j-1,k,ist_ppp) &
                  + sten(i-1,j,k,ist_ppp) + sten(i,j,k,ist_ppp))
             sten(i,j,k,ist_inv) = 1.d0 / (abs(sten(i-1,j,k,ist_p00)) + abs(sten(i,j,k,ist_p00)) &
                  + abs(sten(i,j-1,k,ist_0p0)) + abs(sten(i,j,k,ist_0p0)) &
                  + abs(sten(i,j,k-1,ist_00p)) + abs(sten(i,j,k,ist_00p)) &
                  + abs(sten(i-1,j-1,k,ist_pp0)) + abs(sten(i,j-1,k,ist_pp0)) &
                  + abs(sten(i-1,j,k,ist_pp0)) + abs(sten(i,j,k,ist_pp0)) &
                  + abs(sten(i-1,j,k-1,ist_p0p)) + abs(sten(i,j,k-1,ist_p0p)) &
                  + abs(sten(i-1,j,k,ist_p0p)) + abs(sten(i,j,k,ist_p0p)) &
                  + abs(sten(i,j-1,k-1,ist_0pp)) + abs(sten(i,j,k-1,ist_0pp)) &
                  + abs(sten(i,j-1,k,ist_0pp)) + abs(sten(i,j,k,ist_0pp)) &
                  + abs(sten(i-1,j-1,k-1,ist_ppp)) + abs(sten(i,j-1,k-1,ist_ppp)) &
                  + abs(sten(i-1,j,k-1,ist_ppp)) + abs(sten(i,j,k-1,ist_ppp)) &
                  + abs(sten(i-1,j-1,k,ist_ppp)) + abs(sten(i,j-1,k,ist_ppp)) &
                  + abs(sten(i-1,j,k,ist_ppp)) + abs(sten(i,j,k,ist_ppp)) + eps)

          end do
       end do
    end do

  end subroutine amrex_mlndlap_set_stencil_s0


  subroutine amrex_mlndlap_adotx_sten (lo, hi, y, ylo, yhi, x, xlo, xhi, &
       sten, slo, shi, msk, mlo, mhi) bind(c,name='amrex_mlndlap_adotx_sten')
    integer, dimension(3), intent(in) :: lo, hi, ylo, yhi, xlo, xhi, slo, shi, mlo, mhi
    real(amrex_real), intent(inout) ::   y(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3))
    real(amrex_real), intent(in   ) ::   x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(amrex_real), intent(in   ) ::sten(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),n_sten)
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                y(i,j,k) = x(i  ,j  ,k  ) * sten(i  ,j  ,k  ,ist_000) &
                     !
                     +     x(i-1,j  ,k  ) * sten(i-1,j  ,k  ,ist_p00) &
                     +     x(i+1,j  ,k  ) * sten(i  ,j  ,k  ,ist_p00) &
                     !
                     +     x(i  ,j-1,k  ) * sten(i  ,j-1,k  ,ist_0p0) &
                     +     x(i  ,j+1,k  ) * sten(i  ,j  ,k  ,ist_0p0) &
                     !
                     +     x(i  ,j  ,k-1) * sten(i  ,j  ,k-1,ist_00p) &
                     +     x(i  ,j  ,k+1) * sten(i  ,j  ,k  ,ist_00p) &
                     !
                     +     x(i-1,j-1,k  ) * sten(i-1,j-1,k  ,ist_pp0) &
                     +     x(i+1,j-1,k  ) * sten(i  ,j-1,k  ,ist_pp0) &
                     +     x(i-1,j+1,k  ) * sten(i-1,j  ,k  ,ist_pp0) &
                     +     x(i+1,j+1,k  ) * sten(i  ,j  ,k  ,ist_pp0) &
                     !
                     +     x(i-1,j  ,k-1) * sten(i-1,j  ,k-1,ist_p0p) &
                     +     x(i+1,j  ,k-1) * sten(i  ,j  ,k-1,ist_p0p) &
                     +     x(i-1,j  ,k+1) * sten(i-1,j  ,k  ,ist_p0p) &
                     +     x(i+1,j  ,k+1) * sten(i  ,j  ,k  ,ist_p0p) &
                     !
                     +     x(i  ,j-1,k-1) * sten(i  ,j-1,k-1,ist_0pp) &
                     +     x(i  ,j+1,k-1) * sten(i  ,j  ,k-1,ist_0pp) &
                     +     x(i  ,j-1,k+1) * sten(i  ,j-1,k  ,ist_0pp) &
                     +     x(i  ,j+1,k+1) * sten(i  ,j  ,k  ,ist_0pp) &
                     !
                     +     x(i-1,j-1,k-1) * sten(i-1,j-1,k-1,ist_ppp) &
                     +     x(i+1,j-1,k-1) * sten(i  ,j-1,k-1,ist_ppp) &
                     +     x(i-1,j+1,k-1) * sten(i-1,j  ,k-1,ist_ppp) &
                     +     x(i+1,j+1,k-1) * sten(i  ,j  ,k-1,ist_ppp) &
                     +     x(i-1,j-1,k+1) * sten(i-1,j-1,k  ,ist_ppp) &
                     +     x(i+1,j-1,k+1) * sten(i  ,j-1,k  ,ist_ppp) &
                     +     x(i-1,j+1,k+1) * sten(i-1,j  ,k  ,ist_ppp) &
                     +     x(i+1,j+1,k+1) * sten(i  ,j  ,k  ,ist_ppp)
             else
                y(i,j,k) = 0.d0
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_adotx_sten

  subroutine amrex_mlndlap_normalize_sten (lo, hi, x, xlo, xhi, &
       sten, slo, shi, msk, mlo, mhi) bind(c,name='amrex_mlndlap_normalize_sten')
    integer, dimension(3), intent(in) :: lo, hi, xlo, xhi, slo, shi, mlo, mhi
    real(amrex_real), intent(inout) ::   x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    real(amrex_real), intent(in   ) ::sten(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),n_sten)
    integer         , intent(in   ) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet .and. sten(i,j,k,ist_000) .ne. 0.d0) then
                x(i,j,k) = x(i,j,k) / sten(i,j,k,ist_000)
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_normalize_sten


  subroutine amrex_mlndlap_gauss_seidel_sten (lo, hi, sol, slo, shi, rhs, rlo, rhi, &
       sten, stlo, sthi, msk, mlo, mhi) &
       bind(c,name='amrex_mlndlap_gauss_seidel_sten')
    integer, dimension(3),intent(in) :: lo,hi,slo,shi,rlo,rhi,stlo,sthi,mlo,mhi
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) ::sten(stlo(1):sthi(1),stlo(2):sthi(2),stlo(3):sthi(3),n_sten)
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k
    real(amrex_real) :: Ax

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                if (sten(i,j,k,ist_000) .ne. 0.d0) then
                   Ax   = sol(i  ,j  ,k  ) * sten(i  ,j  ,k  ,ist_000) &
                        !
                        + sol(i-1,j  ,k  ) * sten(i-1,j  ,k  ,ist_p00) &
                        + sol(i+1,j  ,k  ) * sten(i  ,j  ,k  ,ist_p00) &
                        !
                        + sol(i  ,j-1,k  ) * sten(i  ,j-1,k  ,ist_0p0) &
                        + sol(i  ,j+1,k  ) * sten(i  ,j  ,k  ,ist_0p0) &
                        !
                        + sol(i  ,j  ,k-1) * sten(i  ,j  ,k-1,ist_00p) &
                        + sol(i  ,j  ,k+1) * sten(i  ,j  ,k  ,ist_00p) &
                        !
                        + sol(i-1,j-1,k  ) * sten(i-1,j-1,k  ,ist_pp0) &
                        + sol(i+1,j-1,k  ) * sten(i  ,j-1,k  ,ist_pp0) &
                        + sol(i-1,j+1,k  ) * sten(i-1,j  ,k  ,ist_pp0) &
                        + sol(i+1,j+1,k  ) * sten(i  ,j  ,k  ,ist_pp0) &
                        !
                        + sol(i-1,j  ,k-1) * sten(i-1,j  ,k-1,ist_p0p) &
                        + sol(i+1,j  ,k-1) * sten(i  ,j  ,k-1,ist_p0p) &
                        + sol(i-1,j  ,k+1) * sten(i-1,j  ,k  ,ist_p0p) &
                        + sol(i+1,j  ,k+1) * sten(i  ,j  ,k  ,ist_p0p) &
                        !
                        + sol(i  ,j-1,k-1) * sten(i  ,j-1,k-1,ist_0pp) &
                        + sol(i  ,j+1,k-1) * sten(i  ,j  ,k-1,ist_0pp) &
                        + sol(i  ,j-1,k+1) * sten(i  ,j-1,k  ,ist_0pp) &
                        + sol(i  ,j+1,k+1) * sten(i  ,j  ,k  ,ist_0pp) &
                        !
                        + sol(i-1,j-1,k-1) * sten(i-1,j-1,k-1,ist_ppp) &
                        + sol(i+1,j-1,k-1) * sten(i  ,j-1,k-1,ist_ppp) &
                        + sol(i-1,j+1,k-1) * sten(i-1,j  ,k-1,ist_ppp) &
                        + sol(i+1,j+1,k-1) * sten(i  ,j  ,k-1,ist_ppp) &
                        + sol(i-1,j-1,k+1) * sten(i-1,j-1,k  ,ist_ppp) &
                        + sol(i+1,j-1,k+1) * sten(i  ,j-1,k  ,ist_ppp) &
                        + sol(i-1,j+1,k+1) * sten(i-1,j  ,k  ,ist_ppp) &
                        + sol(i+1,j+1,k+1) * sten(i  ,j  ,k  ,ist_ppp)

                   sol(i,j,k) = sol(i,j,k) + (rhs(i,j,k) - Ax) / sten(i,j,k,ist_000)
                end if
             else
                sol(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_gauss_seidel_sten

  subroutine amrex_mlndlap_jacobi_sten (lo, hi, sol, slo, shi, Ax, alo, ahi, &
       rhs, rlo, rhi, sten, stlo, sthi, msk, mlo, mhi) &
       bind(c,name='amrex_mlndlap_jacobi_sten')
    integer, dimension(3),intent(in) :: lo,hi,slo,shi,alo,ahi,rlo,rhi,stlo,sthi,mlo,mhi
    real(amrex_real), intent(inout) :: sol( slo(1): shi(1), slo(2): shi(2), slo(3): shi(3))
    real(amrex_real), intent(in   ) :: Ax ( alo(1): ahi(1), alo(2): ahi(2), alo(3): ahi(3))
    real(amrex_real), intent(in   ) :: rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(in   ) ::sten(stlo(1):sthi(1),stlo(2):sthi(2),stlo(3):sthi(3),1)
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: i,j,k
    real(amrex_real), parameter :: omega = 2.d0/3.d0

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (msk(i,j,k) .ne. dirichlet) then
                if (sten(i,j,k,ist_000) .ne. 0.d0) then
                   sol(i,j,k) = sol(i,j,k) + omega * (rhs(i,j,k) - Ax(i,j,k)) / sten(i,j,k,ist_000)
                end if
             else
                sol(i,j,k) = 0.d0
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_jacobi_sten


  subroutine amrex_mlndlap_interpolation_rap (clo, chi, fine, fflo, ffhi, crse, cflo, cfhi, &
       sten, stlo, sthi, msk, mlo, mhi) bind(c,name='amrex_mlndlap_interpolation_rap')
    integer, dimension(3), intent(in) :: clo,chi,fflo,ffhi,cflo,cfhi,stlo,sthi, mlo, mhi
    real(amrex_real), intent(in   ) :: crse(cflo(1):cfhi(1),cflo(2):cfhi(2),cflo(3):cfhi(3))
    real(amrex_real), intent(inout) :: fine(fflo(1):ffhi(1),fflo(2):ffhi(2),fflo(3):ffhi(3))
    real(amrex_real), intent(in   ) :: sten(stlo(1):sthi(1),stlo(2):sthi(2),stlo(3):sthi(3),n_sten)
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))

    integer :: flo(3), fhi(3), i,j,k, ic,jc,kc
    logical :: ieven, jeven, keven
    real(amrex_real) :: w1, w2, w1m, w1p, w2m, w2p, wmm, wpm, wmp, wpp, wtmp
    real(amrex_real) :: wmmm, wpmm, wmpm, wppm, wmmp, wpmp, wmpp, wppp

    flo = 2*clo
    fhi = 2*chi

    do k = flo(3), fhi(3)
       kc = (k-flo(3))/2 + clo(3)
       keven = kc*2 .eq. k
       do j = flo(2), fhi(2)
          jc = (j-flo(2))/2 + clo(2)
          jeven = jc*2 .eq. j
          do i = flo(1), fhi(1)
             if (msk(i,j,k) .ne. dirichlet .and. sten(i,j,k,ist_000) .ne. 0.d0) then
                ic = (i-flo(1))/2 + clo(1)
                ieven = ic*2 .eq. i
                if (ieven .and. jeven .and. keven) then
                   fine(i,j,k) = crse(ic,jc,kc)
                else if (ieven .and. jeven) then
                   w1 = abs(sten(i,j,k-1,ist_00p))
                   w2 = abs(sten(i,j,k  ,ist_00p))
                   if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
                      fine(i,j,k) = 0.5d0*(crse(ic,jc,kc)+crse(ic,jc,kc+1))
                   else
                      fine(i,j,k) = (w1*crse(ic,jc,kc) + w2*crse(ic,jc,kc+1)) / (w1+w2)
                   end if
                else if (ieven .and. keven) then
                   w1 = abs(sten(i,j-1,k,ist_0p0))
                   w2 = abs(sten(i,j  ,k,ist_0p0))
                   if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
                      fine(i,j,k) = 0.5d0*(crse(ic,jc,kc)+crse(ic,jc+1,kc))
                   else
                      fine(i,j,k) = (w1*crse(ic,jc,kc) + w2*crse(ic,jc+1,kc)) / (w1+w2)
                   end if
                else if (jeven .and. keven) then
                   w1 = abs(sten(i-1,j,k,ist_p00))
                   w2 = abs(sten(i  ,j,k,ist_p00))
                   if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
                      fine(i,j,k) = 0.5d0*(crse(ic,jc,kc)+crse(ic+1,jc,kc))
                   else
                      fine(i,j,k) = (w1*crse(ic,jc,kc) + w2*crse(ic+1,jc,kc)) / (w1+w2)
                   end if
                else if (ieven) then
                   w1m = abs(sten(i,j-1,k,ist_0p0)) / (abs(sten(i,j-1,k-1,ist_0pp)) &
                        &                             +abs(sten(i,j-1,k  ,ist_0pp)) + eps)
                   w1p = abs(sten(i,j  ,k,ist_0p0)) / (abs(sten(i,j  ,k-1,ist_0pp)) &
                        &                             +abs(sten(i,j  ,k  ,ist_0pp)) + eps)
                   w2m = abs(sten(i,j,k-1,ist_00p)) / (abs(sten(i,j-1,k-1,ist_0pp)) &
                        &                             +abs(sten(i,j  ,k-1,ist_0pp)) + eps)
                   w2p = abs(sten(i,j,k  ,ist_00p)) / (abs(sten(i,j-1,k  ,ist_0pp)) &
                        &                             +abs(sten(i,j  ,k  ,ist_0pp)) + eps)
                   wmm = abs(sten(i,j-1,k-1,ist_0pp)) * (1.d0 + w1m + w2m)
                   wpm = abs(sten(i,j  ,k-1,ist_0pp)) * (1.d0 + w1p + w2m)
                   wmp = abs(sten(i,j-1,k  ,ist_0pp)) * (1.d0 + w1m + w2p)
                   wpp = abs(sten(i,j  ,k  ,ist_0pp)) * (1.d0 + w1p + w2p)
                   fine(i,j,k) = (wmm*crse(ic,jc,kc) + wpm*crse(ic,jc+1,kc) &
                        + wmp*crse(ic,jc,kc+1) + wpp*crse(ic,jc+1,kc+1)) &
                        / (wmm+wpm+wmp+wpp+eps)
                else if (jeven) then
                   w1m = abs(sten(i-1,j,k,ist_p00)) / (abs(sten(i-1,j,k-1,ist_p0p)) &
                        &                             +abs(sten(i-1,j,k  ,ist_p0p)) + eps)
                   w1p = abs(sten(i  ,j,k,ist_p00)) / (abs(sten(i  ,j,k-1,ist_p0p)) &
                        &                             +abs(sten(i  ,j,k  ,ist_p0p)) + eps)
                   w2m = abs(sten(i,j,k-1,ist_00p)) / (abs(sten(i-1,j,k-1,ist_p0p)) &
                        &                             +abs(sten(i  ,j,k-1,ist_p0p)) + eps)
                   w2p = abs(sten(i,j,k  ,ist_00p)) / (abs(sten(i-1,j,k  ,ist_p0p)) &
                        &                             +abs(sten(i  ,j,k  ,ist_p0p)) + eps)
                   wmm = abs(sten(i-1,j,k-1,ist_p0p)) * (1.d0 + w1m + w2m)
                   wpm = abs(sten(i  ,j,k-1,ist_p0p)) * (1.d0 + w1p + w2m)
                   wmp = abs(sten(i-1,j,k  ,ist_p0p)) * (1.d0 + w1m + w2p)
                   wpp = abs(sten(i  ,j,k  ,ist_p0p)) * (1.d0 + w1p + w2p)
                   fine(i,j,k) = (wmm*crse(ic,jc,kc) + wpm*crse(ic+1,jc,kc) &
                        + wmp*crse(ic,jc,kc+1) + wpp*crse(ic+1,jc,kc+1)) &
                        / (wmm+wpm+wmp+wpp+eps)
                else if (keven) then
                   w1m = abs(sten(i-1,j,k,ist_p00)) / (abs(sten(i-1,j-1,k,ist_pp0)) &
                        &                             +abs(sten(i-1,j  ,k,ist_pp0)) + eps)
                   w1p = abs(sten(i  ,j,k,ist_p00)) / (abs(sten(i  ,j-1,k,ist_pp0)) &
                        &                             +abs(sten(i  ,j  ,k,ist_pp0)) + eps)
                   w2m = abs(sten(i,j-1,k,ist_0p0)) / (abs(sten(i-1,j-1,k,ist_pp0)) &
                        &                             +abs(sten(i  ,j-1,k,ist_pp0)) + eps)
                   w2p = abs(sten(i,j  ,k,ist_0p0)) / (abs(sten(i-1,j  ,k,ist_pp0)) &
                        &                             +abs(sten(i  ,j  ,k,ist_pp0)) + eps)
                   wmm = abs(sten(i-1,j-1,k,ist_pp0)) * (1.d0 + w1m + w2m)
                   wpm = abs(sten(i  ,j-1,k,ist_pp0)) * (1.d0 + w1p + w2m)
                   wmp = abs(sten(i-1,j  ,k,ist_pp0)) * (1.d0 + w1m + w2p)
                   wpp = abs(sten(i  ,j  ,k,ist_pp0)) * (1.d0 + w1p + w2p)
                   fine(i,j,k) = (wmm*crse(ic,jc,kc) + wpm*crse(ic+1,jc,kc) &
                        + wmp*crse(ic,jc+1,kc) + wpp*crse(ic+1,jc+1,kc)) &
                        / (wmm+wpm+wmp+wpp+eps)
                else
                   wmmm = 1.d0
                   wpmm = 1.d0
                   wmpm = 1.d0
                   wppm = 1.d0
                   wmmp = 1.d0
                   wpmp = 1.d0
                   wmpp = 1.d0
                   wppp = 1.d0

                   wtmp = abs(sten(i-1,j,k,ist_p00)) / &
                        &     ( abs(sten(i-1,j-1,k-1,ist_ppp)) &
                        &     + abs(sten(i-1,j  ,k-1,ist_ppp)) &
                        &     + abs(sten(i-1,j-1,k  ,ist_ppp)) &
                        &     + abs(sten(i-1,j  ,k  ,ist_ppp)) + eps)
                   wmmm = wmmm + wtmp
                   wmpm = wmpm + wtmp
                   wmmp = wmmp + wtmp
                   wmpp = wmpp + wtmp

                   wtmp = abs(sten(i,j,k,ist_p00)) / &
                        &     ( abs(sten(i,j-1,k-1,ist_ppp)) &
                        &     + abs(sten(i,j  ,k-1,ist_ppp)) &
                        &     + abs(sten(i,j-1,k  ,ist_ppp)) &
                        &     + abs(sten(i,j  ,k  ,ist_ppp)) + eps)
                   wpmm = wpmm + wtmp
                   wppm = wppm + wtmp
                   wpmp = wpmp + wtmp
                   wppp = wppp + wtmp

                   wtmp = abs(sten(i,j-1,k,ist_0p0)) / &
                        &     ( abs(sten(i-1,j-1,k-1,ist_ppp)) &
                        &     + abs(sten(i  ,j-1,k-1,ist_ppp)) &
                        &     + abs(sten(i-1,j-1,k  ,ist_ppp)) &
                        &     + abs(sten(i  ,j-1,k  ,ist_ppp)) + eps)
                   wmmm = wmmm + wtmp
                   wpmm = wpmm + wtmp
                   wmmp = wmmp + wtmp
                   wpmp = wpmp + wtmp

                   wtmp = abs(sten(i,j,k,ist_0p0)) / &
                        &     ( abs(sten(i-1,j,k-1,ist_ppp)) &
                        &     + abs(sten(i  ,j,k-1,ist_ppp)) &
                        &     + abs(sten(i-1,j,k  ,ist_ppp)) &
                        &     + abs(sten(i  ,j,k  ,ist_ppp)) + eps)
                   wmpm = wmpm + wtmp
                   wppm = wppm + wtmp
                   wmpp = wmpp + wtmp
                   wppp = wppp + wtmp

                   wtmp = abs(sten(i,j,k-1,ist_00p)) / &
                        &     ( abs(sten(i-1,j-1,k-1,ist_ppp)) &
                        &     + abs(sten(i  ,j-1,k-1,ist_ppp)) &
                        &     + abs(sten(i-1,j  ,k-1,ist_ppp)) &
                        &     + abs(sten(i  ,j  ,k-1,ist_ppp)) + eps)
                   wmmm = wmmm + wtmp
                   wpmm = wpmm + wtmp
                   wmpm = wmpm + wtmp
                   wppm = wppm + wtmp

                   wtmp = abs(sten(i,j,k,ist_00p)) / &
                        &     ( abs(sten(i-1,j-1,k,ist_ppp)) &
                        &     + abs(sten(i  ,j-1,k,ist_ppp)) &
                        &     + abs(sten(i-1,j  ,k,ist_ppp)) &
                        &     + abs(sten(i  ,j  ,k,ist_ppp)) + eps)
                   wmmp = wmmp + wtmp
                   wpmp = wpmp + wtmp
                   wmpp = wmpp + wtmp
                   wppp = wppp + wtmp

                   wtmp = abs(sten(i-1,j-1,k,ist_pp0)) / &
                        &     ( abs(sten(i-1,j-1,k-1,ist_ppp)) &
                        &     + abs(sten(i-1,j-1,k  ,ist_ppp)) + eps)
                   wmmm = wmmm + wtmp
                   wmmp = wmmp + wtmp

                   wtmp = abs(sten(i,j-1,k,ist_pp0)) / &
                        &     ( abs(sten(i,j-1,k-1,ist_ppp)) &
                        &     + abs(sten(i,j-1,k  ,ist_ppp)) + eps)
                   wpmm = wpmm + wtmp
                   wpmp = wpmp + wtmp

                   wtmp = abs(sten(i-1,j,k,ist_pp0)) / &
                        &     ( abs(sten(i-1,j,k-1,ist_ppp)) &
                        &     + abs(sten(i-1,j,k  ,ist_ppp)) + eps)
                   wmpm = wmpm + wtmp
                   wmpp = wmpp + wtmp

                   wtmp = abs(sten(i,j,k,ist_pp0)) / &
                        &     ( abs(sten(i,j,k-1,ist_ppp)) &
                        &     + abs(sten(i,j,k  ,ist_ppp)) + eps)
                   wppm = wppm + wtmp
                   wppp = wppp + wtmp

                   wtmp = abs(sten(i-1,j,k-1,ist_p0p)) / &
                        &     ( abs(sten(i-1,j-1,k-1,ist_ppp)) &
                        &     + abs(sten(i-1,j  ,k-1,ist_ppp)) + eps)
                   wmmm = wmmm + wtmp
                   wmpm = wmpm + wtmp

                   wtmp = abs(sten(i,j,k-1,ist_p0p)) / &
                        &     ( abs(sten(i,j-1,k-1,ist_ppp)) &
                        &     + abs(sten(i,j  ,k-1,ist_ppp)) + eps)
                   wpmm = wpmm + wtmp
                   wppm = wppm + wtmp

                   wtmp = abs(sten(i-1,j,k,ist_p0p)) / &
                        &     ( abs(sten(i-1,j-1,k,ist_ppp)) &
                        &     + abs(sten(i-1,j  ,k,ist_ppp)) + eps)
                   wmmp = wmmp + wtmp
                   wmpp = wmpp + wtmp

                   wtmp = abs(sten(i,j,k,ist_p0p)) / &
                        &     ( abs(sten(i,j-1,k,ist_ppp)) &
                        &     + abs(sten(i,j  ,k,ist_ppp)) + eps)
                   wpmp = wpmp + wtmp
                   wppp = wppp + wtmp

                   wtmp = abs(sten(i,j-1,k-1,ist_0pp)) / &
                        &     ( abs(sten(i-1,j-1,k-1,ist_ppp)) &
                        &     + abs(sten(i  ,j-1,k-1,ist_ppp)) + eps)
                   wmmm = wmmm + wtmp
                   wpmm = wpmm + wtmp

                   wtmp = abs(sten(i,j,k-1,ist_0pp)) / &
                        &     ( abs(sten(i-1,j,k-1,ist_ppp)) &
                        &     + abs(sten(i  ,j,k-1,ist_ppp)) + eps)
                   wmpm = wmpm + wtmp
                   wppm = wppm + wtmp

                   wtmp = abs(sten(i,j-1,k,ist_0pp)) / &
                        &     ( abs(sten(i-1,j-1,k,ist_ppp)) &
                        &     + abs(sten(i  ,j-1,k,ist_ppp)) + eps)
                   wmmp = wmmp + wtmp
                   wpmp = wpmp + wtmp

                   wtmp = abs(sten(i,j,k,ist_0pp)) / &
                        &     ( abs(sten(i-1,j,k,ist_ppp)) &
                        &     + abs(sten(i  ,j,k,ist_ppp)) + eps)
                   wmpp = wmpp + wtmp
                   wppp = wppp + wtmp

                   wmmm = wmmm * abs(sten(i-1,j-1,k-1,ist_ppp))
                   wpmm = wpmm * abs(sten(i  ,j-1,k-1,ist_ppp))
                   wmpm = wmpm * abs(sten(i-1,j  ,k-1,ist_ppp))
                   wppm = wppm * abs(sten(i  ,j  ,k-1,ist_ppp))
                   wmmp = wmmp * abs(sten(i-1,j-1,k  ,ist_ppp))
                   wpmp = wpmp * abs(sten(i  ,j-1,k  ,ist_ppp))
                   wmpp = wmpp * abs(sten(i-1,j  ,k  ,ist_ppp))
                   wppp = wppp * abs(sten(i  ,j  ,k  ,ist_ppp))
                   fine(i,j,k) = (wmmm*crse(ic,jc  ,kc  ) + wpmm*crse(ic+1,jc  ,kc  ) &
                        &       + wmpm*crse(ic,jc+1,kc  ) + wppm*crse(ic+1,jc+1,kc  ) &
                        &       + wmmp*crse(ic,jc  ,kc+1) + wpmp*crse(ic+1,jc  ,kc+1) &
                        &       + wmpp*crse(ic,jc+1,kc+1) + wppp*crse(ic+1,jc+1,kc+1)) &
                        / (wmmm + wpmm + wmpm + wppm + wmmp + wpmp + wmpp + wppp + eps)
                end if
             else
                fine(i,j,k) = 0.d0
             end if
          end do
       end do
    end do
  end subroutine amrex_mlndlap_interpolation_rap


  subroutine amrex_mlndlap_restriction_rap (lo, hi, crse, clo, chi, fine, flo, fhi, &
       sten, slo, shi, msk, mlo, mhi) bind(c,name='amrex_mlndlap_restriction_rap')
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, flo, fhi, slo, shi, mlo, mhi
    real(amrex_real), intent(inout) :: crse(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(amrex_real), intent(in   ) :: fine(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(amrex_real), intent(in   ) :: sten(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),n_sten)
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))


    integer :: i,j,k,ii,jj,kk
    real(amrex_real) :: w1m, w1p, w2m ,w2p, wmm, wpm, wmp, wpp
    real(amrex_real) :: wmmm, wpmm, wmpm, wppm, wmmp, wpmp, wmpp, wppp
    real(amrex_real) :: sten_lo, sten_hi

    do k = lo(3), hi(3)
       kk = 2*k
       do j = lo(2), hi(2)
          jj = 2*j
          do i = lo(1), hi(1)
             ii = 2*i

             if (msk(ii,jj,kk) .ne. dirichlet) then

                crse(i,j,k) = fine(ii,jj,kk)

                ! ************************************
                ! Adding fine(ii-1,jj,kk)
                ! ************************************

                sten_lo = abs(sten(ii-2,jj,kk,ist_p00))
                sten_hi = abs(sten(ii-1,jj,kk,ist_p00))

                if (sten_lo .eq. 0.d0 .and. sten_hi .eq. 0.d0) then
                   crse(i,j,k) = crse(i,j,k) + 0.5d0*fine(ii-1,jj,kk)
                else
                   crse(i,j,k) = crse(i,j,k) + fine(ii-1,jj,kk) * sten_hi / (sten_lo + sten_hi)
                end if

                ! ************************************
                ! Adding fine(ii+1,jj,kk)
                ! ************************************

                sten_lo = abs(sten(ii  ,jj,kk,ist_p00))
                sten_hi = abs(sten(ii+1,jj,kk,ist_p00))

                if (sten_lo .eq. 0.d0 .and. sten_hi .eq. 0.d0) then
                   crse(i,j,k) = crse(i,j,k) + 0.5d0*fine(ii+1,jj,kk)
                else
                   crse(i,j,k) = crse(i,j,k) + fine(ii+1,jj,kk) * sten_lo / (sten_lo + sten_hi)
                end if

                ! ************************************
                ! Adding fine(ii,jj-1,kk)
                ! ************************************

                sten_lo = abs(sten(ii,jj-2,kk,ist_0p0))
                sten_hi = abs(sten(ii,jj-1,kk,ist_0p0))

                if (sten_lo .eq. 0.d0 .and. sten_hi .eq. 0.d0) then
                   crse(i,j,k) = crse(i,j,k) + 0.5d0*fine(ii,jj-1,kk)
                else
                   crse(i,j,k) = crse(i,j,k) + fine(ii,jj-1,kk) * sten_hi / (sten_lo + sten_hi)
                end if

                ! ************************************
                ! Adding fine(ii,jj+1,kk)
                ! ************************************

                sten_lo = abs(sten(ii,jj  ,kk,ist_0p0))
                sten_hi = abs(sten(ii,jj+1,kk,ist_0p0))

                if (sten_lo .eq. 0.d0 .and. sten_hi .eq. 0.d0) then
                   crse(i,j,k) = crse(i,j,k) + 0.5d0*fine(ii,jj+1,kk)
                else
                   crse(i,j,k) = crse(i,j,k) + fine(ii,jj+1,kk) * sten_lo / (sten_lo + sten_hi)
                end if

                ! ************************************
                ! Adding fine(ii,jj,kk-1)
                ! ************************************

                sten_lo = abs(sten(ii,jj,kk-2,ist_00p))
                sten_hi = abs(sten(ii,jj,kk-1,ist_00p))

                if (sten_lo .eq. 0.d0 .and. sten_hi .eq. 0.d0) then
                   crse(i,j,k) = crse(i,j,k) + 0.5d0*fine(ii,jj,kk-1)
                else
                   crse(i,j,k) = crse(i,j,k) + fine(ii,jj,kk-1)*sten_hi / (sten_lo + sten_hi)
                end if

                ! ************************************
                ! Adding fine(ii,jj,kk+1)
                ! ************************************

                sten_lo = abs(sten(ii,jj,kk  ,ist_00p))
                sten_hi = abs(sten(ii,jj,kk+1,ist_00p))

                if (sten_lo .eq. 0.d0 .and. sten_hi .eq. 0.d0) then
                   crse(i,j,k) = crse(i,j,k) + 0.5d0*fine(ii,jj,kk+1)
                else
                   crse(i,j,k) = crse(i,j,k) + fine(ii,jj,kk+1)*sten_lo  / (sten_lo + sten_hi)
                end if

                ! ************************************
                ! Adding fine(ii-1,jj-1,kk)
                ! ************************************

                ! keven
                w1m = abs(sten(ii-2,jj-1,kk,ist_p00)) / (abs(sten(ii-2,jj-2,kk,ist_pp0)) &
                     &                                  +abs(sten(ii-2,jj-1,kk,ist_pp0)) + eps)
                w1p = abs(sten(ii-1,jj-1,kk,ist_p00)) / (abs(sten(ii-1,jj-2,kk,ist_pp0)) &
                     &                                  +abs(sten(ii-1,jj-1,kk,ist_pp0)) + eps)
                w2m = abs(sten(ii-1,jj-2,kk,ist_0p0)) / (abs(sten(ii-2,jj-2,kk,ist_pp0)) &
                     &                                  +abs(sten(ii-1,jj-2,kk,ist_pp0)) + eps)
                w2p = abs(sten(ii-1,jj-1,kk,ist_0p0)) / (abs(sten(ii-2,jj-1,kk,ist_pp0)) &
                     &                                  +abs(sten(ii-1,jj-1,kk,ist_pp0)) + eps)
                wmm = abs(sten(ii-2,jj-2,kk,ist_pp0)) * (1.d0 + w1m + w2m)
                wpm = abs(sten(ii-1,jj-2,kk,ist_pp0)) * (1.d0 + w1p + w2m)
                wmp = abs(sten(ii-2,jj-1,kk,ist_pp0)) * (1.d0 + w1m + w2p)
                wpp = abs(sten(ii-1,jj-1,kk,ist_pp0)) * (1.d0 + w1p + w2p)
                crse(i,j,k) = crse(i,j,k) + fine(ii-1,jj-1,kk)*wpp/(wmm+wpm+wmp+wpp+eps)

                ! ************************************
                ! Adding fine(ii+1,jj-1,kk)
                ! ************************************

                w1m = abs(sten(ii  ,jj-1,kk,ist_p00)) / (abs(sten(ii  ,jj-2,kk,ist_pp0)) &
                     &                                  +abs(sten(ii  ,jj-1,kk,ist_pp0)) + eps)
                w1p = abs(sten(ii+1,jj-1,kk,ist_p00)) / (abs(sten(ii+1,jj-2,kk,ist_pp0)) &
                     &                                  +abs(sten(ii+1,jj-1,kk,ist_pp0)) + eps)
                w2m = abs(sten(ii+1,jj-2,kk,ist_0p0)) / (abs(sten(ii  ,jj-2,kk,ist_pp0)) &
                     &                                  +abs(sten(ii+1,jj-2,kk,ist_pp0)) + eps) 
                w2p = abs(sten(ii+1,jj-1,kk,ist_0p0)) / (abs(sten(ii  ,jj-1,kk,ist_pp0)) &
                     &                                  +abs(sten(ii+1,jj-1,kk,ist_pp0)) + eps)
                wmm = abs(sten(ii  ,jj-2,kk,ist_pp0)) * (1.d0 + w1m + w2m)
                wpm = abs(sten(ii+1,jj-2,kk,ist_pp0)) * (1.d0 + w1p + w2m)
                wmp = abs(sten(ii  ,jj-1,kk,ist_pp0)) * (1.d0 + w1m + w2p)
                wpp = abs(sten(ii+1,jj-1,kk,ist_pp0)) * (1.d0 + w1p + w2p)
                crse(i,j,k) = crse(i,j,k) + fine(ii+1,jj-1,kk)*wmp/(wmm+wpm+wmp+wpp+eps)

                ! ************************************
                ! Adding fine(ii-1,jj+1,kk)
                ! ************************************

                w1m = abs(sten(ii-2,jj+1,kk,ist_p00)) / (abs(sten(ii-2,jj  ,kk,ist_pp0)) &
                     &                                  +abs(sten(ii-2,jj+1,kk,ist_pp0)) + eps)
                w1p = abs(sten(ii-1,jj+1,kk,ist_p00)) / (abs(sten(ii-1,jj  ,kk,ist_pp0)) &
                     &                                  +abs(sten(ii-1,jj+1,kk,ist_pp0)) + eps)
                w2m = abs(sten(ii-1,jj  ,kk,ist_0p0)) / (abs(sten(ii-2,jj  ,kk,ist_pp0)) &
                     &                                  +abs(sten(ii-1,jj  ,kk,ist_pp0)) + eps)
                w2p = abs(sten(ii-1,jj+1,kk,ist_0p0)) / (abs(sten(ii-2,jj+1,kk,ist_pp0)) &
                     &                                  +abs(sten(ii-1,jj+1,kk,ist_pp0)) + eps)
                wmm = abs(sten(ii-2,jj  ,kk,ist_pp0)) * (1.d0 + w1m + w2m)
                wpm = abs(sten(ii-1,jj  ,kk,ist_pp0)) * (1.d0 + w1p + w2m)
                wmp = abs(sten(ii-2,jj+1,kk,ist_pp0)) * (1.d0 + w1m + w2p)
                wpp = abs(sten(ii-1,jj+1,kk,ist_pp0)) * (1.d0 + w1p + w2p)
                crse(i,j,k) = crse(i,j,k) + fine(ii-1,jj+1,kk)*wpm/(wmm+wpm+wmp+wpp+eps)

                ! ************************************
                ! Adding fine(ii+1,jj+1,kk)
                ! ************************************

                w1m = abs(sten(ii  ,jj+1,kk,ist_p00)) / (abs(sten(ii  ,jj+1,kk,ist_pp0)) &
                     &                                  +abs(sten(ii  ,jj  ,kk,ist_pp0)) + eps)
                w1p = abs(sten(ii+1,jj+1,kk,ist_p00)) / (abs(sten(ii+1,jj+1,kk,ist_pp0)) &
                     &                                  +abs(sten(ii+1,jj  ,kk,ist_pp0)) + eps)
                w2m = abs(sten(ii+1,jj  ,kk,ist_0p0)) / (abs(sten(ii+1,jj  ,kk,ist_pp0)) &
                     &                                  +abs(sten(ii  ,jj  ,kk,ist_pp0)) + eps)
                w2p = abs(sten(ii+1,jj+1,kk,ist_0p0)) / (abs(sten(ii+1,jj+1,kk,ist_pp0)) &
                     &                                  +abs(sten(ii  ,jj+1,kk,ist_pp0)) + eps)
                wmm = abs(sten(ii  ,jj  ,kk,ist_pp0)) * (1.d0 + w1m + w2m)
                wpm = abs(sten(ii+1,jj  ,kk,ist_pp0)) * (1.d0 + w1p + w2m)
                wmp = abs(sten(ii  ,jj+1,kk,ist_pp0)) * (1.d0 + w1m + w2p)
                wpp = abs(sten(ii+1,jj+1,kk,ist_pp0)) * (1.d0 + w1p + w2p)
                crse(i,j,k) = crse(i,j,k) + fine(ii+1,jj+1,kk)*wmm/(wmm+wpm+wmp+wpp+eps)

                ! ************************************
                ! Adding fine(ii-1,jj,kk-1)
                ! ************************************

                ! jeven
                w1m = abs(sten(ii-2,jj,kk-1,ist_p00)) / (abs(sten(ii-2,jj,kk-2,ist_p0p)) &
                     &                                  +abs(sten(ii-2,jj,kk-1,ist_p0p)) + eps)
                w1p = abs(sten(ii-1,jj,kk-1,ist_p00)) / (abs(sten(ii-1,jj,kk-2,ist_p0p)) &
                     &                                  +abs(sten(ii-1,jj,kk-1,ist_p0p)) + eps)
                w2m = abs(sten(ii-1,jj,kk-2,ist_00p)) / (abs(sten(ii-2,jj,kk-2,ist_p0p)) &
                     &                                  +abs(sten(ii-1,jj,kk-2,ist_p0p)) + eps)
                w2p = abs(sten(ii-1,jj,kk-1,ist_00p)) / (abs(sten(ii-2,jj,kk-1,ist_p0p)) &
                     &                                  +abs(sten(ii-1,jj,kk-1,ist_p0p)) + eps)
                wmm = abs(sten(ii-2,jj,kk-2,ist_p0p)) * (1.d0 + w1m + w2m)
                wpm = abs(sten(ii-1,jj,kk-2,ist_p0p)) * (1.d0 + w1p + w2m)
                wmp = abs(sten(ii-2,jj,kk-1,ist_p0p)) * (1.d0 + w1m + w2p)
                wpp = abs(sten(ii-1,jj,kk-1,ist_p0p)) * (1.d0 + w1p + w2p)
                crse(i,j,k) = crse(i,j,k) + fine(ii-1,jj,kk-1)*wpp/(wmm+wpm+wmp+wpp+eps)

                ! ************************************
                ! Adding fine(ii+1,jj,kk-1)
                ! ************************************

                w1m = abs(sten(ii  ,jj,kk-1,ist_p00)) / (abs(sten(ii  ,jj,kk-2,ist_p0p)) &
                     &                                  +abs(sten(ii  ,jj,kk-1,ist_p0p)) + eps)
                w1p = abs(sten(ii+1,jj,kk-1,ist_p00)) / (abs(sten(ii+1,jj,kk-2,ist_p0p)) &
                     &                                  +abs(sten(ii+1,jj,kk-1,ist_p0p)) + eps)
                w2m = abs(sten(ii+1,jj,kk-2,ist_00p)) / (abs(sten(ii+1,jj,kk-2,ist_p0p)) &
                     &                                  +abs(sten(ii  ,jj,kk-2,ist_p0p)) + eps)
                w2p = abs(sten(ii+1,jj,kk-1,ist_00p)) / (abs(sten(ii+1,jj,kk-1,ist_p0p)) &
                     &                                  +abs(sten(ii  ,jj,kk-1,ist_p0p)) + eps)
                wmm = abs(sten(ii  ,jj,kk-2,ist_p0p)) * (1.d0 + w1m + w2m)
                wpm = abs(sten(ii+1,jj,kk-2,ist_p0p)) * (1.d0 + w1p + w2m) 
                wmp = abs(sten(ii  ,jj,kk-1,ist_p0p)) * (1.d0 + w1m + w2p)
                wpp = abs(sten(ii+1,jj,kk-1,ist_p0p)) * (1.d0 + w1p + w2p)
                crse(i,j,k) = crse(i,j,k) + fine(ii+1,jj,kk-1)*wmp/(wmm+wpm+wmp+wpp+eps)

                ! ************************************
                ! Adding fine(ii-1,jj,kk+1)
                ! ************************************

                w1m = abs(sten(ii-2,jj,kk+1,ist_p00)) / (abs(sten(ii-2,jj,kk+1,ist_p0p)) &
                     &                                  +abs(sten(ii-2,jj,kk  ,ist_p0p)) + eps)
                w1p = abs(sten(ii-1,jj,kk+1,ist_p00)) / (abs(sten(ii-1,jj,kk+1,ist_p0p)) &
                     &                                  +abs(sten(ii-1,jj,kk  ,ist_p0p)) + eps)
                w2m = abs(sten(ii-1,jj,kk  ,ist_00p)) / (abs(sten(ii-2,jj,kk  ,ist_p0p)) &
                     &                                  +abs(sten(ii-1,jj,kk  ,ist_p0p)) + eps)
                w2p = abs(sten(ii-1,jj,kk+1,ist_00p)) / (abs(sten(ii-2,jj,kk+1,ist_p0p)) &
                     &                                  +abs(sten(ii-1,jj,kk+1,ist_p0p)) + eps)
                wmm = abs(sten(ii-2,jj,kk  ,ist_p0p)) * (1.d0 + w1m + w2m)
                wpm = abs(sten(ii-1,jj,kk  ,ist_p0p)) * (1.d0 + w1p + w2m)
                wmp = abs(sten(ii-2,jj,kk+1,ist_p0p)) * (1.d0 + w1m + w2p)
                wpp = abs(sten(ii-1,jj,kk+1,ist_p0p)) * (1.d0 + w1p + w2p)
                crse(i,j,k) = crse(i,j,k) + fine(ii-1,jj,kk+1)*wpm/(wmm+wpm+wmp+wpp+eps)

                ! ************************************
                ! Adding fine(ii+1,jj,kk+1)
                ! ************************************

                w1m = abs(sten(ii  ,jj,kk+1,ist_p00)) / (abs(sten(ii  ,jj,kk+1,ist_p0p)) &
                     &                                  +abs(sten(ii  ,jj,kk  ,ist_p0p)) + eps)
                w1p = abs(sten(ii+1,jj,kk+1,ist_p00)) / (abs(sten(ii+1,jj,kk+1,ist_p0p)) &
                     &                                  +abs(sten(ii+1,jj,kk  ,ist_p0p)) + eps)
                w2m = abs(sten(ii+1,jj,kk  ,ist_00p)) / (abs(sten(ii+1,jj,kk  ,ist_p0p)) &
                     &                                  +abs(sten(ii  ,jj,kk  ,ist_p0p)) + eps)
                w2p = abs(sten(ii+1,jj,kk+1,ist_00p)) / (abs(sten(ii+1,jj,kk+1,ist_p0p)) &
                     &                                  +abs(sten(ii  ,jj,kk+1,ist_p0p)) + eps)
                wmm = abs(sten(ii  ,jj,kk  ,ist_p0p)) * (1.d0 + w1m + w2m)
                wpm = abs(sten(ii+1,jj,kk  ,ist_p0p)) * (1.d0 + w1p + w2m)
                wmp = abs(sten(ii  ,jj,kk+1,ist_p0p)) * (1.d0 + w1m + w2p)
                wpp = abs(sten(ii+1,jj,kk+1,ist_p0p)) * (1.d0 + w1p + w2p)
                crse(i,j,k) = crse(i,j,k) + fine(ii+1,jj,kk+1)*wmm/(wmm+wpm+wmp+wpp+eps)

                ! ************************************
                ! Adding fine(ii,jj-1,kk-1)
                ! ************************************

                ! ieven
                w1m = abs(sten(ii,jj-2,kk-1,ist_0p0)) / (abs(sten(ii,jj-2,kk-2,ist_0pp)) &
                     &                                  +abs(sten(ii,jj-2,kk-1,ist_0pp)) + eps)
                w2m = abs(sten(ii,jj-1,kk-2,ist_00p)) / (abs(sten(ii,jj-2,kk-2,ist_0pp)) &
                     &                                  +abs(sten(ii,jj-1,kk-2,ist_0pp)) + eps)
                w1p = abs(sten(ii,jj-1,kk-1,ist_0p0)) / (abs(sten(ii,jj-1,kk-2,ist_0pp)) &
                     &                                  +abs(sten(ii,jj-1,kk-1,ist_0pp)) + eps)
                w2p = abs(sten(ii,jj-1,kk-1,ist_00p)) / (abs(sten(ii,jj-2,kk-1,ist_0pp)) &
                     &                                  +abs(sten(ii,jj-1,kk-1,ist_0pp)) + eps)
                wmm = abs(sten(ii,jj-2,kk-2,ist_0pp)) * (1.d0 + w1m + w2m)
                wpm = abs(sten(ii,jj-1,kk-2,ist_0pp)) * (1.d0 + w1p + w2m)
                wmp = abs(sten(ii,jj-2,kk-1,ist_0pp)) * (1.d0 + w1m + w2p)
                wpp = abs(sten(ii,jj-1,kk-1,ist_0pp)) * (1.d0 + w1p + w2p)
                crse(i,j,k) = crse(i,j,k) + fine(ii,jj-1,kk-1)*wpp/(wmm+wpm+wmp+wpp+eps)

                ! ************************************
                ! Adding fine(ii,jj+1,kk-1)
                ! ************************************

                w1m = abs(sten(ii,jj  ,kk-1,ist_0p0)) / (abs(sten(ii,jj  ,kk-2,ist_0pp)) &
                     &                                  +abs(sten(ii,jj  ,kk-1,ist_0pp)) + eps)
                w1p = abs(sten(ii,jj+1,kk-1,ist_0p0)) / (abs(sten(ii,jj+1,kk-2,ist_0pp)) &
                     &                                  +abs(sten(ii,jj+1,kk-1,ist_0pp)) + eps)
                w2m = abs(sten(ii,jj+1,kk-2,ist_00p)) / (abs(sten(ii,jj+1,kk-2,ist_0pp)) &
                     &                                  +abs(sten(ii,jj  ,kk-2,ist_0pp)) + eps)
                w2p = abs(sten(ii,jj+1,kk-1,ist_00p)) / (abs(sten(ii,jj+1,kk-1,ist_0pp)) &
                     &                                  +abs(sten(ii,jj  ,kk-1,ist_0pp)) + eps)
                wmm = abs(sten(ii,jj  ,kk-2,ist_0pp)) * (1.d0 + w1m + w2m)
                wpm = abs(sten(ii,jj+1,kk-2,ist_0pp)) * (1.d0 + w1p + w2m)
                wmp = abs(sten(ii,jj  ,kk-1,ist_0pp)) * (1.d0 + w1m + w2p)
                wpp = abs(sten(ii,jj+1,kk-1,ist_0pp)) * (1.d0 + w1p + w2p)

                crse(i,j,k) = crse(i,j,k) + fine(ii,jj+1,kk-1)*wmp/(wmm+wpm+wmp+wpp+eps)

                ! ************************************
                ! Adding fine(ii,jj-1,kk+1)
                ! ************************************

                w1m = abs(sten(ii,jj-2,kk+1,ist_0p0)) / (abs(sten(ii,jj-2,kk+1,ist_0pp)) &
                     &                                  +abs(sten(ii,jj-2,kk  ,ist_0pp)) + eps)
                w1p = abs(sten(ii,jj-1,kk+1,ist_0p0)) / (abs(sten(ii,jj-1,kk+1,ist_0pp)) &
                     &                                  +abs(sten(ii,jj-1,kk  ,ist_0pp)) + eps)
                w2m = abs(sten(ii,jj-1,kk  ,ist_00p)) / (abs(sten(ii,jj-2,kk  ,ist_0pp)) &
                     &                                  +abs(sten(ii,jj-1,kk  ,ist_0pp)) + eps)
                w2p = abs(sten(ii,jj-1,kk+1,ist_00p)) / (abs(sten(ii,jj-2,kk+1,ist_0pp)) &
                     &                                  +abs(sten(ii,jj-1,kk+1,ist_0pp)) + eps)
                wmm = abs(sten(ii,jj-2,kk  ,ist_0pp)) * (1.d0 + w1m + w2m)
                wpm = abs(sten(ii,jj-1,kk  ,ist_0pp)) * (1.d0 + w1p + w2m)
                wmp = abs(sten(ii,jj-2,kk+1,ist_0pp)) * (1.d0 + w1m + w2p)
                wpp = abs(sten(ii,jj-1,kk+1,ist_0pp)) * (1.d0 + w1p + w2p)
                crse(i,j,k) = crse(i,j,k) + fine(ii,jj-1,kk+1)*wpm/(wmm+wpm+wmp+wpp+eps)

                ! ************************************
                ! Adding fine(ii,jj+1,kk+1)
                ! ************************************

                w1m = abs(sten(ii,jj  ,kk+1,ist_0p0)) / (abs(sten(ii,jj  ,kk+1,ist_0pp)) &
                     &                                  +abs(sten(ii,jj  ,kk  ,ist_0pp)) + eps)
                w1p = abs(sten(ii,jj+1,kk+1,ist_0p0)) / (abs(sten(ii,jj+1,kk+1,ist_0pp)) &
                     &                                  +abs(sten(ii,jj+1,kk  ,ist_0pp)) + eps)
                w2m = abs(sten(ii,jj+1,kk  ,ist_00p)) / (abs(sten(ii,jj+1,kk  ,ist_0pp)) &
                     &                                  +abs(sten(ii,jj  ,kk  ,ist_0pp)) + eps)
                w2p = abs(sten(ii,jj+1,kk+1,ist_00p)) / (abs(sten(ii,jj+1,kk+1,ist_0pp)) &
                     &                                  +abs(sten(ii,jj  ,kk+1,ist_0pp)) + eps)
                wmm = abs(sten(ii,jj  ,kk  ,ist_0pp)) * (1.d0 + w1m + w2m)
                wpm = abs(sten(ii,jj+1,kk  ,ist_0pp)) * (1.d0 + w1p + w2m)
                wmp = abs(sten(ii,jj  ,kk+1,ist_0pp)) * (1.d0 + w1m + w2p)
                wpp = abs(sten(ii,jj+1,kk+1,ist_0pp)) * (1.d0 + w1p + w2p)
                crse(i,j,k) = crse(i,j,k) + fine(ii,jj+1,kk+1)*wmm/(wmm+wpm+wmp+wpp+eps)

                ! ************************************
                ! Adding fine at corners
                ! ************************************

                wmmm = 1.d0 &
                     +   abs(sten(ii  ,jj+1,kk+1,ist_p00)) / &
                     & ( abs(sten(ii  ,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii  ,jj+1,kk  ,ist_ppp)) &
                     & + abs(sten(ii  ,jj  ,kk+1,ist_ppp)) &
                     & + abs(sten(ii  ,jj+1,kk+1,ist_ppp)) + eps) &
                     +   abs(sten(ii+1,jj  ,kk+1,ist_0p0)) / &
                     & ( abs(sten(ii  ,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii+1,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii  ,jj  ,kk+1,ist_ppp)) &
                     & + abs(sten(ii+1,jj  ,kk+1,ist_ppp)) + eps) &
                     +   abs(sten(ii+1,jj+1,kk  ,ist_00p)) / &
                     & ( abs(sten(ii  ,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii+1,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii  ,jj+1,kk  ,ist_ppp)) &
                     & + abs(sten(ii+1,jj+1,kk  ,ist_ppp)) + eps) &
                     +   abs(sten(ii  ,jj  ,kk+1,ist_pp0)) / &
                     & ( abs(sten(ii  ,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii  ,jj  ,kk+1,ist_ppp)) + eps) &
                     +   abs(sten(ii  ,jj+1,kk  ,ist_p0p)) / &
                     & ( abs(sten(ii  ,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii  ,jj+1,kk  ,ist_ppp)) + eps) &
                      +  abs(sten(ii+1,jj  ,kk  ,ist_0pp)) / &
                     & ( abs(sten(ii  ,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii+1,jj  ,kk  ,ist_ppp)) + eps)
                wmmm = wmmm * abs(sten(ii,jj,kk,ist_ppp))
                crse(i,j,k) = crse(i,j,k) + wmmm*fine(ii+1,jj+1,kk+1)*sten(ii+1,jj+1,kk+1,ist_inv)

                wpmm = 1.d0 &
                     +   abs(sten(ii-1,jj+1,kk+1,ist_p00)) / &
                     & ( abs(sten(ii-1,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj+1,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj  ,kk+1,ist_ppp)) &
                     & + abs(sten(ii-1,jj+1,kk+1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj  ,kk+1,ist_0p0)) / &
                     & ( abs(sten(ii-2,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii-2,jj  ,kk+1,ist_ppp)) &
                     & + abs(sten(ii-1,jj  ,kk+1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj+1,kk  ,ist_00p)) / &
                     & ( abs(sten(ii-2,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii-2,jj+1,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj+1,kk  ,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj  ,kk+1,ist_pp0)) / &
                     & ( abs(sten(ii-1,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj  ,kk+1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj+1,kk  ,ist_p0p)) / &
                     & ( abs(sten(ii-1,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj+1,kk  ,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj  ,kk  ,ist_0pp)) / &
                     & ( abs(sten(ii-2,jj  ,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj  ,kk  ,ist_ppp)) + eps)
                wpmm = wpmm * abs(sten(ii-1,jj,kk,ist_ppp))
                crse(i,j,k) = crse(i,j,k) + wpmm*fine(ii-1,jj+1,kk+1)*sten(ii-1,jj+1,kk+1,ist_inv)

                wmpm = 1.d0 &
                     +   abs(sten(ii  ,jj-1,kk+1,ist_p00)) / &
                     & ( abs(sten(ii  ,jj-2,kk  ,ist_ppp)) &
                     & + abs(sten(ii  ,jj-1,kk  ,ist_ppp)) &
                     & + abs(sten(ii  ,jj-2,kk+1,ist_ppp)) &
                     & + abs(sten(ii  ,jj-1,kk+1,ist_ppp)) + eps) &
                     +   abs(sten(ii+1,jj-1,kk+1,ist_0p0)) / &
                     & ( abs(sten(ii  ,jj-1,kk  ,ist_ppp)) &
                     & + abs(sten(ii+1,jj-1,kk  ,ist_ppp)) &
                     & + abs(sten(ii  ,jj-1,kk+1,ist_ppp)) &
                     & + abs(sten(ii+1,jj-1,kk+1,ist_ppp)) + eps) &
                     +   abs(sten(ii+1,jj-1,kk  ,ist_00p)) / &
                     & ( abs(sten(ii  ,jj-2,kk  ,ist_ppp)) &
                     & + abs(sten(ii+1,jj-2,kk  ,ist_ppp)) &
                     & + abs(sten(ii  ,jj-1,kk  ,ist_ppp)) &
                     & + abs(sten(ii+1,jj-1,kk  ,ist_ppp)) + eps) &
                     +   abs(sten(ii  ,jj-1,kk+1,ist_pp0)) / &
                     & ( abs(sten(ii  ,jj-1,kk  ,ist_ppp)) &
                     & + abs(sten(ii  ,jj-1,kk+1,ist_ppp)) + eps) &
                     +   abs(sten(ii  ,jj-1,kk  ,ist_p0p)) / &
                     & ( abs(sten(ii  ,jj-2,kk  ,ist_ppp)) &
                     & + abs(sten(ii  ,jj-1,kk  ,ist_ppp)) + eps) &
                     +   abs(sten(ii+1,jj-1,kk  ,ist_0pp)) / &
                     & ( abs(sten(ii  ,jj-1,kk  ,ist_ppp)) &
                     & + abs(sten(ii+1,jj-1,kk  ,ist_ppp)) + eps)
                   wmpm = wmpm * abs(sten(ii  ,jj-1,kk  ,ist_ppp))
                crse(i,j,k) = crse(i,j,k) + wmpm*fine(ii+1,jj-1,kk+1)*sten(ii+1,jj-1,kk+1,ist_inv)

                wppm = 1.d0 &
                     +   abs(sten(ii-1,jj-1,kk+1,ist_p00)) / &
                     & ( abs(sten(ii-1,jj-2,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj-2,kk+1,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj-1,kk+1,ist_0p0)) / &
                     & ( abs(sten(ii-2,jj-1,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk  ,ist_ppp)) &
                     & + abs(sten(ii-2,jj-1,kk+1,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj-1,kk  ,ist_00p)) / &
                     & ( abs(sten(ii-2,jj-2,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj-2,kk  ,ist_ppp)) &
                     & + abs(sten(ii-2,jj-1,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj-1,kk+1,ist_pp0)) / &
                     & ( abs(sten(ii-1,jj-1,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk+1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj-1,kk  ,ist_p0p)) / &
                     & ( abs(sten(ii-1,jj-2,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj-1,kk  ,ist_0pp)) / &
                     & ( abs(sten(ii-2,jj-1,kk  ,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk  ,ist_ppp)) + eps)
                wppm = wppm * abs(sten(ii-1,jj-1,kk  ,ist_ppp))
                crse(i,j,k) = crse(i,j,k) + wppm*fine(ii-1,jj-1,kk+1)*sten(ii-1,jj-1,kk+1,ist_inv)

                wmmp = 1.d0 &
                     +   abs(sten(ii  ,jj+1,kk-1,ist_p00)) / &
                     & ( abs(sten(ii  ,jj  ,kk-2,ist_ppp)) &
                     & + abs(sten(ii  ,jj+1,kk-2,ist_ppp)) &
                     & + abs(sten(ii  ,jj  ,kk-1,ist_ppp)) &
                     & + abs(sten(ii  ,jj+1,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii+1,jj  ,kk-1,ist_0p0)) / &
                     & ( abs(sten(ii  ,jj  ,kk-2,ist_ppp)) &
                     & + abs(sten(ii+1,jj  ,kk-2,ist_ppp)) &
                     & + abs(sten(ii  ,jj  ,kk-1,ist_ppp)) &
                     & + abs(sten(ii+1,jj  ,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii+1,jj+1,kk-1,ist_00p)) / &
                     & ( abs(sten(ii  ,jj  ,kk-1,ist_ppp)) &
                     & + abs(sten(ii+1,jj  ,kk-1,ist_ppp)) &
                     & + abs(sten(ii  ,jj+1,kk-1,ist_ppp)) &
                     & + abs(sten(ii+1,jj+1,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii  ,jj  ,kk-1,ist_pp0)) / &
                     & ( abs(sten(ii  ,jj  ,kk-2,ist_ppp)) &
                     & + abs(sten(ii  ,jj  ,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii  ,jj+1,kk-1,ist_p0p)) / &
                     & ( abs(sten(ii  ,jj  ,kk-1,ist_ppp)) &
                     & + abs(sten(ii  ,jj+1,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii+1,jj  ,kk-1,ist_0pp)) / &
                     & ( abs(sten(ii  ,jj  ,kk-1,ist_ppp)) &
                     & + abs(sten(ii+1,jj  ,kk-1,ist_ppp)) + eps)
                wmmp = wmmp * abs(sten(ii  ,jj  ,kk-1,ist_ppp))
                crse(i,j,k) = crse(i,j,k) + wmmp*fine(ii+1,jj+1,kk-1)*sten(ii+1,jj+1,kk-1,ist_inv)

                wpmp = 1.d0 &
                     +   abs(sten(ii-1,jj+1,kk-1,ist_p00)) / &
                     & ( abs(sten(ii-1,jj  ,kk-2,ist_ppp)) &
                     & + abs(sten(ii-1,jj+1,kk-2,ist_ppp)) &
                     & + abs(sten(ii-1,jj  ,kk-1,ist_ppp)) &
                     & + abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj  ,kk-1,ist_0p0)) / &
                     & ( abs(sten(ii-2,jj  ,kk-2,ist_ppp)) &
                     & + abs(sten(ii-1,jj  ,kk-2,ist_ppp)) &
                     & + abs(sten(ii-2,jj  ,kk-1,ist_ppp)) &
                     & + abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj+1,kk-1,ist_00p)) / &
                     & ( abs(sten(ii-2,jj  ,kk-1,ist_ppp)) &
                     & + abs(sten(ii-1,jj  ,kk-1,ist_ppp)) &
                     & + abs(sten(ii-2,jj+1,kk-1,ist_ppp)) &
                     & + abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj  ,kk-1,ist_pp0)) / &
                     & ( abs(sten(ii-1,jj  ,kk-2,ist_ppp)) &
                     & + abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj+1,kk-1,ist_p0p)) / &
                     & ( abs(sten(ii-1,jj  ,kk-1,ist_ppp)) &
                     & + abs(sten(ii-1,jj+1,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj  ,kk-1,ist_0pp)) / &
                     & ( abs(sten(ii-2,jj  ,kk-1,ist_ppp)) &
                     & + abs(sten(ii-1,jj  ,kk-1,ist_ppp)) + eps)
                wpmp = wpmp * abs(sten(ii-1,jj  ,kk-1,ist_ppp))
                crse(i,j,k) = crse(i,j,k) + wpmp*fine(ii-1,jj+1,kk-1)*sten(ii-1,jj+1,kk-1,ist_inv)

                wmpp = 1.d0 &
                     +       abs(sten(ii  ,jj-1,kk-1,ist_p00)) / &
                     &     ( abs(sten(ii  ,jj-2,kk-2,ist_ppp)) &
                     &     + abs(sten(ii  ,jj-1,kk-2,ist_ppp)) &
                     &     + abs(sten(ii  ,jj-2,kk-1,ist_ppp)) &
                     &     + abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps) &
                     +       abs(sten(ii+1,jj-1,kk-1,ist_0p0)) / &
                     &     ( abs(sten(ii  ,jj-1,kk-2,ist_ppp)) &
                     &     + abs(sten(ii+1,jj-1,kk-2,ist_ppp)) &
                     &     + abs(sten(ii  ,jj-1,kk-1,ist_ppp)) &
                     &     + abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps) &
                     +       abs(sten(ii+1,jj-1,kk-1,ist_00p)) / &
                     &     ( abs(sten(ii  ,jj-2,kk-1,ist_ppp)) &
                     &     + abs(sten(ii+1,jj-2,kk-1,ist_ppp)) &
                     &     + abs(sten(ii  ,jj-1,kk-1,ist_ppp)) &
                     &     + abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps) &
                     +       abs(sten(ii  ,jj-1,kk-1,ist_pp0)) / &
                     &     ( abs(sten(ii  ,jj-1,kk-2,ist_ppp)) &
                     &     + abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps) &
                     +       abs(sten(ii  ,jj-1,kk-1,ist_p0p)) / &
                     &     ( abs(sten(ii  ,jj-2,kk-1,ist_ppp)) &
                     &     + abs(sten(ii  ,jj-1,kk-1,ist_ppp)) + eps) &
                     +       abs(sten(ii+1,jj-1,kk-1,ist_0pp)) / &
                     &     ( abs(sten(ii  ,jj-1,kk-1,ist_ppp)) &
                     &     + abs(sten(ii+1,jj-1,kk-1,ist_ppp)) + eps)
                wmpp = wmpp * abs(sten(ii  ,jj-1,kk-1,ist_ppp))
                crse(i,j,k) = crse(i,j,k) + wmpp*fine(ii+1,jj-1,kk-1)*sten(ii+1,jj-1,kk-1,ist_inv)

                wppp = 1.d0 &
                     +   abs(sten(ii-1,jj-1,kk-1,ist_p00)) / &
                     & ( abs(sten(ii-1,jj-2,kk-2,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk-2,ist_ppp)) &
                     & + abs(sten(ii-1,jj-2,kk-1,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj-1,kk-1,ist_0p0)) / &
                     & ( abs(sten(ii-2,jj-1,kk-2,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk-2,ist_ppp)) &
                     & + abs(sten(ii-2,jj-1,kk-1,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj-1,kk-1,ist_00p)) / &
                     & ( abs(sten(ii-2,jj-2,kk-1,ist_ppp)) &
                     & + abs(sten(ii-1,jj-2,kk-1,ist_ppp)) &
                     & + abs(sten(ii-2,jj-1,kk-1,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj-1,kk-1,ist_pp0)) / &
                     & ( abs(sten(ii-1,jj-1,kk-2,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj-1,kk-1,ist_p0p)) / &
                     & ( abs(sten(ii-1,jj-2,kk-1,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps) &
                     +   abs(sten(ii-1,jj-1,kk-1,ist_0pp)) / &
                     & ( abs(sten(ii-2,jj-1,kk-1,ist_ppp)) &
                     & + abs(sten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
                wppp = wppp * abs(sten(ii-1,jj-1,kk-1,ist_ppp))
                crse(i,j,k) = crse(i,j,k) + wppp*fine(ii-1,jj-1,kk-1)*sten(ii-1,jj-1,kk-1,ist_inv)

                crse(i,j,k) = crse(i,j,k) * 0.125d0

             else
                crse(i,j,k) = 0.d0
             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_restriction_rap


  subroutine amrex_mlndlap_stencil_rap (lo, hi, csten, clo, chi, fsten, flo, fhi) &
       bind(c,name='amrex_mlndlap_stencil_rap')
    integer, dimension(3), intent(in) :: lo, hi, clo, chi, flo, fhi
    real(amrex_real), intent(inout) :: csten(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),n_sten)
    real(amrex_real), intent(in   ) :: fsten(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3),n_sten)

    integer :: i,j,k,ii,jj,kk,iii,jjj,kkk
    real(amrex_real) :: ap(-1:1,-1:1,-1:1), p(-1:1,-1:1,-1:1)
    real(amrex_real) :: cs1, cs2, cs3, cs4

    do k = lo(3), hi(3)
       kk = 2*k
       do j = lo(2), hi(2)
          jj = 2*j
          do i = lo(1), hi(1)
             ii = 2*i

             ! csten(i,j,k,ist_p00)
             iii = ii
             jjj = jj
             kkk = kk
             p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj-1,kkk-1)
             p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj-1,kkk-1)
             p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj  ,kkk-1)
             p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj  ,kkk-1)
             p(-1,+1,-1) = interp_from_pmp_to(iii+1,jjj+1,kkk-1)
             p( 0,+1,-1) = interp_from_0mp_to(iii+2,jjj+1,kkk-1)
             p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj-1,kkk  )
             p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj-1,kkk  )
             p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj  ,kkk  )
             p( 0, 0, 0) = 1.d0
             p(-1,+1, 0) = interp_from_pm0_to(iii+1,jjj+1,kkk  )
             p( 0,+1, 0) = interp_from_0m0_to(iii+2,jjj+1,kkk  )
             p(-1,-1,+1) = interp_from_ppm_to(iii+1,jjj-1,kkk+1)
             p( 0,-1,+1) = interp_from_0pm_to(iii+2,jjj-1,kkk+1)
             p(-1, 0,+1) = interp_from_p0m_to(iii+1,jjj  ,kkk+1)
             p( 0, 0,+1) = interp_from_00m_to(iii+2,jjj  ,kkk+1)
             p(-1,+1,+1) = interp_from_pmm_to(iii+1,jjj+1,kkk+1)
             p( 0,+1,+1) = interp_from_0mm_to(iii+2,jjj+1,kkk+1)
             ap(0,-1,-1) = &
               &              Ap00(iii,jjj-1,kkk-1) * p(-1,-1,-1) &
               +              App0(iii,jjj-1,kkk-1) * p(-1, 0,-1) &
               +              Ap0p(iii,jjj-1,kkk-1) * p(-1,-1, 0) &
               +              Appp(iii,jjj-1,kkk-1) * p(-1, 0, 0)
             ap(1,-1,-1) = &
               &              A000(iii+1,jjj-1,kkk-1) * p(-1,-1,-1) &
               +              Ap00(iii+1,jjj-1,kkk-1) * p( 0,-1,-1) &
               +              A0p0(iii+1,jjj-1,kkk-1) * p(-1, 0,-1) &
               +              App0(iii+1,jjj-1,kkk-1) * p( 0, 0,-1) &
               +              A00p(iii+1,jjj-1,kkk-1) * p(-1,-1, 0) &
               +              Ap0p(iii+1,jjj-1,kkk-1) * p( 0,-1, 0) &
               +              A0pp(iii+1,jjj-1,kkk-1) * p(-1, 0, 0) &
               +              Appp(iii+1,jjj-1,kkk-1) * p( 0, 0, 0)
             ap(0,0,-1) = &
               &              Apm0(iii,jjj,kkk-1) * p(-1,-1,-1) &
               +              Ap00(iii,jjj,kkk-1) * p(-1, 0,-1) &
               +              App0(iii,jjj,kkk-1) * p(-1,+1,-1) &
               +              Apmp(iii,jjj,kkk-1) * p(-1,-1, 0) &
               +              Ap0p(iii,jjj,kkk-1) * p(-1, 0, 0) &
               +              Appp(iii,jjj,kkk-1) * p(-1,+1, 0)
             ap(1,0,-1) = &
               &              A0m0(iii+1,jjj,kkk-1) * p(-1,-1,-1) &
               +              Apm0(iii+1,jjj,kkk-1) * p( 0,-1,-1) &
               +              A000(iii+1,jjj,kkk-1) * p(-1, 0,-1) &
               +              Ap00(iii+1,jjj,kkk-1) * p( 0, 0,-1) &
               +              A0p0(iii+1,jjj,kkk-1) * p(-1,+1,-1) &
               +              App0(iii+1,jjj,kkk-1) * p( 0,+1,-1) &
               +              A0mp(iii+1,jjj,kkk-1) * p(-1,-1, 0) &
               +              Apmp(iii+1,jjj,kkk-1) * p( 0,-1, 0) &
               +              A00p(iii+1,jjj,kkk-1) * p(-1, 0, 0) &
               +              Ap0p(iii+1,jjj,kkk-1) * p( 0, 0, 0) &
               +              A0pp(iii+1,jjj,kkk-1) * p(-1,+1, 0) &
               +              Appp(iii+1,jjj,kkk-1) * p( 0,+1, 0)
             ap(0,1,-1) = &
               &              Apm0(iii,jjj+1,kkk-1) * p(-1, 0,-1) &
               +              Ap00(iii,jjj+1,kkk-1) * p(-1,+1,-1) &
               +              Apmp(iii,jjj+1,kkk-1) * p(-1, 0, 0) &
               +              Ap0p(iii,jjj+1,kkk-1) * p(-1,+1, 0)
             ap(1,1,-1) = &
               &              A0m0(iii+1,jjj+1,kkk-1) * p(-1, 0,-1) &
               +              Apm0(iii+1,jjj+1,kkk-1) * p( 0, 0,-1) &
               +              A000(iii+1,jjj+1,kkk-1) * p(-1,+1,-1) &
               +              Ap00(iii+1,jjj+1,kkk-1) * p( 0,+1,-1) &
               +              A0mp(iii+1,jjj+1,kkk-1) * p(-1, 0, 0) &
               +              Apmp(iii+1,jjj+1,kkk-1) * p( 0, 0, 0) &
               +              A00p(iii+1,jjj+1,kkk-1) * p(-1,+1, 0) &
               +              Ap0p(iii+1,jjj+1,kkk-1) * p( 0,+1, 0)
             ap(0,-1,0) = &
               &              Ap0m(iii,jjj-1,kkk) * p(-1,-1,-1) &
               +              Appm(iii,jjj-1,kkk) * p(-1, 0,-1) &
               +              Ap00(iii,jjj-1,kkk) * p(-1,-1, 0) &
               +              App0(iii,jjj-1,kkk) * p(-1, 0, 0) &
               +              Ap0p(iii,jjj-1,kkk) * p(-1,-1,+1) &
               +              Appp(iii,jjj-1,kkk) * p(-1, 0,+1)
             ap(1,-1,0) = &
               &              A00m(iii+1,jjj-1,kkk) * p(-1,-1,-1) &
               +              Ap0m(iii+1,jjj-1,kkk) * p( 0,-1,-1) &
               +              A0pm(iii+1,jjj-1,kkk) * p(-1, 0,-1) &
               +              Appm(iii+1,jjj-1,kkk) * p( 0, 0,-1) &
               +              A000(iii+1,jjj-1,kkk) * p(-1,-1, 0) &
               +              Ap00(iii+1,jjj-1,kkk) * p( 0,-1, 0) &
               +              A0p0(iii+1,jjj-1,kkk) * p(-1, 0, 0) &
               +              App0(iii+1,jjj-1,kkk) * p( 0, 0, 0) &
               +              A00p(iii+1,jjj-1,kkk) * p(-1,-1,+1) &
               +              Ap0p(iii+1,jjj-1,kkk) * p( 0,-1,+1) &
               +              A0pp(iii+1,jjj-1,kkk) * p(-1, 0,+1) &
               +              Appp(iii+1,jjj-1,kkk) * p( 0, 0,+1)
             ap(0,0,0) = &
               &              Apmm(iii,jjj,kkk) * p(-1,-1,-1) &
               +              Ap0m(iii,jjj,kkk) * p(-1, 0,-1) &
               +              Appm(iii,jjj,kkk) * p(-1,+1,-1) &
               +              Apm0(iii,jjj,kkk) * p(-1,-1, 0) &
               +              Ap00(iii,jjj,kkk) * p(-1, 0, 0) &
               +              App0(iii,jjj,kkk) * p(-1,+1, 0) &
               +              Apmp(iii,jjj,kkk) * p(-1,-1,+1) &
               +              Ap0p(iii,jjj,kkk) * p(-1, 0,+1) &
               +              Appp(iii,jjj,kkk) * p(-1,+1,+1)
             ap(1,0,0) = &
               &              A0mm(iii+1,jjj,kkk) * p(-1,-1,-1) &
               +              Apmm(iii+1,jjj,kkk) * p( 0,-1,-1) &
               +              A00m(iii+1,jjj,kkk) * p(-1, 0,-1) &
               +              Ap0m(iii+1,jjj,kkk) * p( 0, 0,-1) &
               +              A0pm(iii+1,jjj,kkk) * p(-1,+1,-1) &
               +              Appm(iii+1,jjj,kkk) * p( 0,+1,-1) &
               +              A0m0(iii+1,jjj,kkk) * p(-1,-1, 0) &
               +              Apm0(iii+1,jjj,kkk) * p( 0,-1, 0) &
               +              A000(iii+1,jjj,kkk) * p(-1, 0, 0) &
               +              Ap00(iii+1,jjj,kkk) * p( 0, 0, 0) &
               +              A0p0(iii+1,jjj,kkk) * p(-1,+1, 0) &
               +              App0(iii+1,jjj,kkk) * p( 0,+1, 0) &
               +              A0mp(iii+1,jjj,kkk) * p(-1,-1,+1) &
               +              Apmp(iii+1,jjj,kkk) * p( 0,-1,+1) &
               +              A00p(iii+1,jjj,kkk) * p(-1, 0,+1) &
               +              Ap0p(iii+1,jjj,kkk) * p( 0, 0,+1) &
               +              A0pp(iii+1,jjj,kkk) * p(-1,+1,+1) &
               +              Appp(iii+1,jjj,kkk) * p( 0,+1,+1)
             ap(0,1,0) = &
               &              Apmm(iii,jjj+1,kkk) * p(-1, 0,-1) &
               +              Ap0m(iii,jjj+1,kkk) * p(-1,+1,-1) &
               +              Apm0(iii,jjj+1,kkk) * p(-1, 0, 0) &
               +              Ap00(iii,jjj+1,kkk) * p(-1,+1, 0) &
               +              Apmp(iii,jjj+1,kkk) * p(-1, 0,+1) &
               +              Ap0p(iii,jjj+1,kkk) * p(-1,+1,+1)
             ap(1,1,0) = &
               &              A0mm(iii+1,jjj+1,kkk) * p(-1, 0,-1) &
               +              Apmm(iii+1,jjj+1,kkk) * p( 0, 0,-1) &
               +              A00m(iii+1,jjj+1,kkk) * p(-1,+1,-1) &
               +              Ap0m(iii+1,jjj+1,kkk) * p( 0,+1,-1) &
               +              A0m0(iii+1,jjj+1,kkk) * p(-1, 0, 0) &
               +              Apm0(iii+1,jjj+1,kkk) * p( 0, 0, 0) &
               +              A000(iii+1,jjj+1,kkk) * p(-1,+1, 0) &
               +              Ap00(iii+1,jjj+1,kkk) * p( 0,+1, 0) &
               +              A0mp(iii+1,jjj+1,kkk) * p(-1, 0,+1) &
               +              Apmp(iii+1,jjj+1,kkk) * p( 0, 0,+1) &
               +              A00p(iii+1,jjj+1,kkk) * p(-1,+1,+1) &
               +              Ap0p(iii+1,jjj+1,kkk) * p( 0,+1,+1)
             ap(0,-1,1) = &
               &              Ap0m(iii,jjj-1,kkk+1) * p(-1,-1, 0) &
               +              Appm(iii,jjj-1,kkk+1) * p(-1, 0, 0) &
               +              Ap00(iii,jjj-1,kkk+1) * p(-1,-1,+1) &
               +              App0(iii,jjj-1,kkk+1) * p(-1, 0,+1)
             ap(1,-1,1) = &
               &              A00m(iii+1,jjj-1,kkk+1) * p(-1,-1, 0) &
               +              Ap0m(iii+1,jjj-1,kkk+1) * p( 0,-1, 0) &
               +              A0pm(iii+1,jjj-1,kkk+1) * p(-1, 0, 0) &
               +              Appm(iii+1,jjj-1,kkk+1) * p( 0, 0, 0) &
               +              A000(iii+1,jjj-1,kkk+1) * p(-1,-1,+1) &
               +              Ap00(iii+1,jjj-1,kkk+1) * p( 0,-1,+1) &
               +              A0p0(iii+1,jjj-1,kkk+1) * p(-1, 0,+1) &
               +              App0(iii+1,jjj-1,kkk+1) * p( 0, 0,+1)
             ap(0,0,1) = &
               &              Apmm(iii,jjj,kkk+1) * p(-1,-1, 0) &
               +              Ap0m(iii,jjj,kkk+1) * p(-1, 0, 0) &
               +              Appm(iii,jjj,kkk+1) * p(-1,+1, 0) &
               +              Apm0(iii,jjj,kkk+1) * p(-1,-1,+1) &
               +              Ap00(iii,jjj,kkk+1) * p(-1, 0,+1) &
               +              App0(iii,jjj,kkk+1) * p(-1,+1,+1)
             ap(1,0,1) = &
               &              A0mm(iii+1,jjj,kkk+1) * p(-1,-1, 0) &
               +              Apmm(iii+1,jjj,kkk+1) * p( 0,-1, 0) &
               +              A00m(iii+1,jjj,kkk+1) * p(-1, 0, 0) &
               +              Ap0m(iii+1,jjj,kkk+1) * p( 0, 0, 0) &
               +              A0pm(iii+1,jjj,kkk+1) * p(-1,+1, 0) &
               +              Appm(iii+1,jjj,kkk+1) * p( 0,+1, 0) &
               +              A0m0(iii+1,jjj,kkk+1) * p(-1,-1,+1) &
               +              Apm0(iii+1,jjj,kkk+1) * p( 0,-1,+1) &
               +              A000(iii+1,jjj,kkk+1) * p(-1, 0,+1) &
               +              Ap00(iii+1,jjj,kkk+1) * p( 0, 0,+1) &
               +              A0p0(iii+1,jjj,kkk+1) * p(-1,+1,+1) &
               +              App0(iii+1,jjj,kkk+1) * p( 0,+1,+1)
             ap(0,1,1) = &
               &              Apmm(iii,jjj+1,kkk+1) * p(-1, 0, 0) &
               +              Ap0m(iii,jjj+1,kkk+1) * p(-1,+1, 0) &
               +              Apm0(iii,jjj+1,kkk+1) * p(-1, 0,+1) &
               +              Ap00(iii,jjj+1,kkk+1) * p(-1,+1,+1)
             ap(1,1,1) = &
               &              A0mm(iii+1,jjj+1,kkk+1) * p(-1, 0, 0) &
               +              Apmm(iii+1,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              A00m(iii+1,jjj+1,kkk+1) * p(-1,+1, 0) &
               +              Ap0m(iii+1,jjj+1,kkk+1) * p( 0,+1, 0) &
               +              A0m0(iii+1,jjj+1,kkk+1) * p(-1, 0,+1) &
               +              Apm0(iii+1,jjj+1,kkk+1) * p( 0, 0,+1) &
               +              A000(iii+1,jjj+1,kkk+1) * p(-1,+1,+1) &
               +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,+1,+1)
             csten(i,j,k,ist_p00) = 0.125d0 * &
               ( restrict_from_0mm_to(iii,jjj,kkk) * ap( 0,-1,-1) &
               + restrict_from_pmm_to(iii,jjj,kkk) * ap(+1,-1,-1) &
               + restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1) &
               + restrict_from_p0m_to(iii,jjj,kkk) * ap(+1, 0,-1) &
               + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1) &
               + restrict_from_ppm_to(iii,jjj,kkk) * ap(+1,+1,-1) &
               + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0) &
               + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0) &
               + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0) &
               + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0) &
               + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0) &
               + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0) &
               + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1) &
               + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1) &
               + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1) &
               + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1) &
               + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1) &
               + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1))

             ! csten(i,j,k,ist_0p0)
             iii = ii
             jjj = jj
             kkk = kk
             p(-1,-1,-1) = interp_from_ppp_to(iii-1,jjj+1,kkk-1)
             p( 0,-1,-1) = interp_from_0pp_to(iii  ,jjj+1,kkk-1)
             p(+1,-1,-1) = interp_from_mpp_to(iii+1,jjj+1,kkk-1)
             p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj+2,kkk-1)
             p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj+2,kkk-1)
             p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj+2,kkk-1)
             p(-1,-1, 0) = interp_from_pp0_to(iii-1,jjj+1,kkk  )
             p( 0,-1, 0) = interp_from_0p0_to(iii  ,jjj+1,kkk  )
             p(+1,-1, 0) = interp_from_mp0_to(iii+1,jjj+1,kkk  )
             p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj+2,kkk  )
             p( 0, 0, 0) = 1.d0
             p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj+2,kkk  )
             p(-1,-1,+1) = interp_from_ppm_to(iii-1,jjj+1,kkk+1)
             p( 0,-1,+1) = interp_from_0pm_to(iii  ,jjj+1,kkk+1)
             p(+1,-1,+1) = interp_from_mpm_to(iii+1,jjj+1,kkk+1)
             p(-1, 0,+1) = interp_from_p0m_to(iii-1,jjj+2,kkk+1)
             p( 0, 0,+1) = interp_from_00m_to(iii  ,jjj+2,kkk+1)
             p(+1, 0,+1) = interp_from_m0m_to(iii+1,jjj+2,kkk+1)
             ap(-1,0,-1) = &
               &              A0p0(iii-1,jjj,kkk-1) * p(-1,-1,-1) &
               +              App0(iii-1,jjj,kkk-1) * p( 0,-1,-1) &
               +              A0pp(iii-1,jjj,kkk-1) * p(-1,-1, 0) &
               +              Appp(iii-1,jjj,kkk-1) * p( 0,-1, 0)
             ap(0,0,-1) = &
               &              Amp0(iii,jjj,kkk-1) * p(-1,-1,-1) &
               +              A0p0(iii,jjj,kkk-1) * p( 0,-1,-1) &
               +              App0(iii,jjj,kkk-1) * p(+1,-1,-1) &
               +              Ampp(iii,jjj,kkk-1) * p(-1,-1, 0) &
               +              A0pp(iii,jjj,kkk-1) * p( 0,-1, 0) &
               +              Appp(iii,jjj,kkk-1) * p(+1,-1, 0)
             ap(1,0,-1) = &
               &              Amp0(iii+1,jjj,kkk-1) * p( 0,-1,-1) &
               +              A0p0(iii+1,jjj,kkk-1) * p(+1,-1,-1) &
               +              Ampp(iii+1,jjj,kkk-1) * p( 0,-1, 0) &
               +              A0pp(iii+1,jjj,kkk-1) * p(+1,-1, 0)
             ap(-1,1,-1) = &
               &              A000(iii-1,jjj+1,kkk-1) * p(-1,-1,-1) &
               +              Ap00(iii-1,jjj+1,kkk-1) * p( 0,-1,-1) &
               +              A0p0(iii-1,jjj+1,kkk-1) * p(-1, 0,-1) &
               +              App0(iii-1,jjj+1,kkk-1) * p( 0, 0,-1) &
               +              A00p(iii-1,jjj+1,kkk-1) * p(-1,-1, 0) &
               +              Ap0p(iii-1,jjj+1,kkk-1) * p( 0,-1, 0) &
               +              A0pp(iii-1,jjj+1,kkk-1) * p(-1, 0, 0) &
               +              Appp(iii-1,jjj+1,kkk-1) * p( 0, 0, 0)
             ap(0,1,-1) = &
               &              Am00(iii,jjj+1,kkk-1) * p(-1,-1,-1) &
               +              A000(iii,jjj+1,kkk-1) * p( 0,-1,-1) &
               +              Ap00(iii,jjj+1,kkk-1) * p(+1,-1,-1) &
               +              Amp0(iii,jjj+1,kkk-1) * p(-1, 0,-1) &
               +              A0p0(iii,jjj+1,kkk-1) * p( 0, 0,-1) &
               +              App0(iii,jjj+1,kkk-1) * p(+1, 0,-1) &
               +              Am0p(iii,jjj+1,kkk-1) * p(-1,-1, 0) &
               +              A00p(iii,jjj+1,kkk-1) * p( 0,-1, 0) &
               +              Ap0p(iii,jjj+1,kkk-1) * p(+1,-1, 0) &
               +              Ampp(iii,jjj+1,kkk-1) * p(-1, 0, 0) &
               +              A0pp(iii,jjj+1,kkk-1) * p( 0, 0, 0) &
               +              Appp(iii,jjj+1,kkk-1) * p(+1, 0, 0)
             ap(1,1,-1) = &
               &              Am00(iii+1,jjj+1,kkk-1) * p( 0,-1,-1) &
               +              A000(iii+1,jjj+1,kkk-1) * p(+1,-1,-1) &
               +              Amp0(iii+1,jjj+1,kkk-1) * p( 0, 0,-1) &
               +              A0p0(iii+1,jjj+1,kkk-1) * p(+1, 0,-1) &
               +              Am0p(iii+1,jjj+1,kkk-1) * p( 0,-1, 0) &
               +              A00p(iii+1,jjj+1,kkk-1) * p(+1,-1, 0) &
               +              Ampp(iii+1,jjj+1,kkk-1) * p( 0, 0, 0) &
               +              A0pp(iii+1,jjj+1,kkk-1) * p(+1, 0, 0)
             ap(-1,0,0) = &
               &              A0pm(iii-1,jjj,kkk) * p(-1,-1,-1) &
               +              Appm(iii-1,jjj,kkk) * p( 0,-1,-1) &
               +              A0p0(iii-1,jjj,kkk) * p(-1,-1, 0) &
               +              App0(iii-1,jjj,kkk) * p( 0,-1, 0) &
               +              A0pp(iii-1,jjj,kkk) * p(-1,-1,+1) &
               +              Appp(iii-1,jjj,kkk) * p( 0,-1,+1)
             ap(0,0,0) = &
               &              Ampm(iii,jjj,kkk) * p(-1,-1,-1) &
               +              A0pm(iii,jjj,kkk) * p( 0,-1,-1) &
               +              Appm(iii,jjj,kkk) * p(+1,-1,-1) &
               +              Amp0(iii,jjj,kkk) * p(-1,-1, 0) &
               +              A0p0(iii,jjj,kkk) * p( 0,-1, 0) &
               +              App0(iii,jjj,kkk) * p(+1,-1, 0) &
               +              Ampp(iii,jjj,kkk) * p(-1,-1,+1) &
               +              A0pp(iii,jjj,kkk) * p( 0,-1,+1) &
               +              Appp(iii,jjj,kkk) * p(+1,-1,+1)
             ap(1,0,0) = &
               &              Ampm(iii+1,jjj,kkk) * p( 0,-1,-1) &
               +              A0pm(iii+1,jjj,kkk) * p(+1,-1,-1) &
               +              Amp0(iii+1,jjj,kkk) * p( 0,-1, 0) &
               +              A0p0(iii+1,jjj,kkk) * p(+1,-1, 0) &
               +              Ampp(iii+1,jjj,kkk) * p( 0,-1,+1) &
               +              A0pp(iii+1,jjj,kkk) * p(+1,-1,+1)
             ap(-1,1,0) = &
               &              A00m(iii-1,jjj+1,kkk) * p(-1,-1,-1) &
               +              Ap0m(iii-1,jjj+1,kkk) * p( 0,-1,-1) &
               +              A0pm(iii-1,jjj+1,kkk) * p(-1, 0,-1) &
               +              Appm(iii-1,jjj+1,kkk) * p( 0, 0,-1) &
               +              A000(iii-1,jjj+1,kkk) * p(-1,-1, 0) &
               +              Ap00(iii-1,jjj+1,kkk) * p( 0,-1, 0) &
               +              A0p0(iii-1,jjj+1,kkk) * p(-1, 0, 0) &
               +              App0(iii-1,jjj+1,kkk) * p( 0, 0, 0) &
               +              A00p(iii-1,jjj+1,kkk) * p(-1,-1,+1) &
               +              Ap0p(iii-1,jjj+1,kkk) * p( 0,-1,+1) &
               +              A0pp(iii-1,jjj+1,kkk) * p(-1, 0,+1) &
               +              Appp(iii-1,jjj+1,kkk) * p( 0, 0,+1)
             ap(0,1,0) = &
               &              Am0m(iii,jjj+1,kkk) * p(-1,-1,-1) &
               +              A00m(iii,jjj+1,kkk) * p( 0,-1,-1) &
               +              Ap0m(iii,jjj+1,kkk) * p(+1,-1,-1) &
               +              Ampm(iii,jjj+1,kkk) * p(-1, 0,-1) &
               +              A0pm(iii,jjj+1,kkk) * p( 0, 0,-1) &
               +              Appm(iii,jjj+1,kkk) * p(+1, 0,-1) &
               +              Am00(iii,jjj+1,kkk) * p(-1,-1, 0) &
               +              A000(iii,jjj+1,kkk) * p( 0,-1, 0) &
               +              Ap00(iii,jjj+1,kkk) * p(+1,-1, 0) &
               +              Amp0(iii,jjj+1,kkk) * p(-1, 0, 0) &
               +              A0p0(iii,jjj+1,kkk) * p( 0, 0, 0) &
               +              App0(iii,jjj+1,kkk) * p(+1, 0, 0) &
               +              Am0p(iii,jjj+1,kkk) * p(-1,-1,+1) &
               +              A00p(iii,jjj+1,kkk) * p( 0,-1,+1) &
               +              Ap0p(iii,jjj+1,kkk) * p(+1,-1,+1) &
               +              Ampp(iii,jjj+1,kkk) * p(-1, 0,+1) &
               +              A0pp(iii,jjj+1,kkk) * p( 0, 0,+1) &
               +              Appp(iii,jjj+1,kkk) * p(+1, 0,+1)
             ap(1,1,0) = &
               &              Am0m(iii+1,jjj+1,kkk) * p( 0,-1,-1) &
               +              A00m(iii+1,jjj+1,kkk) * p(+1,-1,-1) &
               +              Ampm(iii+1,jjj+1,kkk) * p( 0, 0,-1) &
               +              A0pm(iii+1,jjj+1,kkk) * p(+1, 0,-1) &
               +              Am00(iii+1,jjj+1,kkk) * p( 0,-1, 0) &
               +              A000(iii+1,jjj+1,kkk) * p(+1,-1, 0) &
               +              Amp0(iii+1,jjj+1,kkk) * p( 0, 0, 0) &
               +              A0p0(iii+1,jjj+1,kkk) * p(+1, 0, 0) &
               +              Am0p(iii+1,jjj+1,kkk) * p( 0,-1,+1) &
               +              A00p(iii+1,jjj+1,kkk) * p(+1,-1,+1) &
               +              Ampp(iii+1,jjj+1,kkk) * p( 0, 0,+1) &
               +              A0pp(iii+1,jjj+1,kkk) * p(+1, 0,+1)
             ap(-1,0,1) = &
               &              A0pm(iii-1,jjj,kkk+1) * p(-1,-1, 0) &
               +              Appm(iii-1,jjj,kkk+1) * p( 0,-1, 0) &
               +              A0p0(iii-1,jjj,kkk+1) * p(-1,-1,+1) &
               +              App0(iii-1,jjj,kkk+1) * p( 0,-1,+1)
             ap(0,0,1) = &
               &              Ampm(iii,jjj,kkk+1) * p(-1,-1, 0) &
               +              A0pm(iii,jjj,kkk+1) * p( 0,-1, 0) &
               +              Appm(iii,jjj,kkk+1) * p(+1,-1, 0) &
               +              Amp0(iii,jjj,kkk+1) * p(-1,-1,+1) &
               +              A0p0(iii,jjj,kkk+1) * p( 0,-1,+1) &
               +              App0(iii,jjj,kkk+1) * p(+1,-1,+1)
             ap(1,0,1) = &
               &              Ampm(iii+1,jjj,kkk+1) * p( 0,-1, 0) &
               +              A0pm(iii+1,jjj,kkk+1) * p(+1,-1, 0) &
               +              Amp0(iii+1,jjj,kkk+1) * p( 0,-1,+1) &
               +              A0p0(iii+1,jjj,kkk+1) * p(+1,-1,+1)
             ap(-1,1,1) = &
               &              A00m(iii-1,jjj+1,kkk+1) * p(-1,-1, 0) &
               +              Ap0m(iii-1,jjj+1,kkk+1) * p( 0,-1, 0) &
               +              A0pm(iii-1,jjj+1,kkk+1) * p(-1, 0, 0) &
               +              Appm(iii-1,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              A000(iii-1,jjj+1,kkk+1) * p(-1,-1,+1) &
               +              Ap00(iii-1,jjj+1,kkk+1) * p( 0,-1,+1) &
               +              A0p0(iii-1,jjj+1,kkk+1) * p(-1, 0,+1) &
               +              App0(iii-1,jjj+1,kkk+1) * p( 0, 0,+1)
             ap(0,1,1) = &
               &              Am0m(iii,jjj+1,kkk+1) * p(-1,-1, 0) &
               +              A00m(iii,jjj+1,kkk+1) * p( 0,-1, 0) &
               +              Ap0m(iii,jjj+1,kkk+1) * p(+1,-1, 0) &
               +              Ampm(iii,jjj+1,kkk+1) * p(-1, 0, 0) &
               +              A0pm(iii,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              Appm(iii,jjj+1,kkk+1) * p(+1, 0, 0) &
               +              Am00(iii,jjj+1,kkk+1) * p(-1,-1,+1) &
               +              A000(iii,jjj+1,kkk+1) * p( 0,-1,+1) &
               +              Ap00(iii,jjj+1,kkk+1) * p(+1,-1,+1) &
               +              Amp0(iii,jjj+1,kkk+1) * p(-1, 0,+1) &
               +              A0p0(iii,jjj+1,kkk+1) * p( 0, 0,+1) &
               +              App0(iii,jjj+1,kkk+1) * p(+1, 0,+1)
             ap(1,1,1) = &
               &              Am0m(iii+1,jjj+1,kkk+1) * p( 0,-1, 0) &
               +              A00m(iii+1,jjj+1,kkk+1) * p(+1,-1, 0) &
               +              Ampm(iii+1,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              A0pm(iii+1,jjj+1,kkk+1) * p(+1, 0, 0) &
               +              Am00(iii+1,jjj+1,kkk+1) * p( 0,-1,+1) &
               +              A000(iii+1,jjj+1,kkk+1) * p(+1,-1,+1) &
               +              Amp0(iii+1,jjj+1,kkk+1) * p( 0, 0,+1) &
               +              A0p0(iii+1,jjj+1,kkk+1) * p(+1, 0,+1)
             csten(i,j,k,ist_0p0) = 0.125d0 * &
               ( restrict_from_m0m_to(iii,jjj,kkk) * ap(-1, 0,-1) &
               + restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1) &
               + restrict_from_p0m_to(iii,jjj,kkk) * ap(+1, 0,-1) &
               + restrict_from_mpm_to(iii,jjj,kkk) * ap(-1,+1,-1) &
               + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1) &
               + restrict_from_ppm_to(iii,jjj,kkk) * ap(+1,+1,-1) &
               + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0) &
               + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0) &
               + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0) &
               + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0) &
               + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0) &
               + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0) &
               + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1) &
               + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1) &
               + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1) &
               + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1) &
               + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1) &
               + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1))

             ! csten(i,j,k,ist_00p)
             iii = ii
             jjj = jj
             kkk = kk
             p(-1,-1,-1) = interp_from_ppp_to(iii-1,jjj-1,kkk+1)
             p( 0,-1,-1) = interp_from_0pp_to(iii  ,jjj-1,kkk+1)
             p(+1,-1,-1) = interp_from_mpp_to(iii+1,jjj-1,kkk+1)
             p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj  ,kkk+1)
             p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj  ,kkk+1)
             p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj  ,kkk+1)
             p(-1,+1,-1) = interp_from_pmp_to(iii-1,jjj+1,kkk+1)
             p( 0,+1,-1) = interp_from_0mp_to(iii  ,jjj+1,kkk+1)
             p(+1,+1,-1) = interp_from_mmp_to(iii+1,jjj+1,kkk+1)
             p(-1,-1, 0) = interp_from_pp0_to(iii-1,jjj-1,kkk+2)
             p( 0,-1, 0) = interp_from_0p0_to(iii  ,jjj-1,kkk+2)
             p(+1,-1, 0) = interp_from_mp0_to(iii+1,jjj-1,kkk+2)
             p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj  ,kkk+2)
             p( 0, 0, 0) = 1.d0
             p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj  ,kkk+2)
             p(-1,+1, 0) = interp_from_pm0_to(iii-1,jjj+1,kkk+2)
             p( 0,+1, 0) = interp_from_0m0_to(iii  ,jjj+1,kkk+2)
             p(+1,+1, 0) = interp_from_mm0_to(iii+1,jjj+1,kkk+2)
             ap(-1,-1,0) = &
               &              A00p(iii-1,jjj-1,kkk) * p(-1,-1,-1) &
               +              Ap0p(iii-1,jjj-1,kkk) * p( 0,-1,-1) &
               +              A0pp(iii-1,jjj-1,kkk) * p(-1, 0,-1) &
               +              Appp(iii-1,jjj-1,kkk) * p( 0, 0,-1)
             ap(0,-1,0) = &
               &              Am0p(iii,jjj-1,kkk) * p(-1,-1,-1) &
               +              A00p(iii,jjj-1,kkk) * p( 0,-1,-1) &
               +              Ap0p(iii,jjj-1,kkk) * p(+1,-1,-1) &
               +              Ampp(iii,jjj-1,kkk) * p(-1, 0,-1) &
               +              A0pp(iii,jjj-1,kkk) * p( 0, 0,-1) &
               +              Appp(iii,jjj-1,kkk) * p(+1, 0,-1)
             ap(1,-1,0) = &
               &              Am0p(iii+1,jjj-1,kkk) * p( 0,-1,-1) &
               +              A00p(iii+1,jjj-1,kkk) * p(+1,-1,-1) &
               +              Ampp(iii+1,jjj-1,kkk) * p( 0, 0,-1) &
               +              A0pp(iii+1,jjj-1,kkk) * p(+1, 0,-1)
             ap(-1,0,0) = &
               &              A0mp(iii-1,jjj,kkk) * p(-1,-1,-1) &
               +              Apmp(iii-1,jjj,kkk) * p( 0,-1,-1) &
               +              A00p(iii-1,jjj,kkk) * p(-1, 0,-1) &
               +              Ap0p(iii-1,jjj,kkk) * p( 0, 0,-1) &
               +              A0pp(iii-1,jjj,kkk) * p(-1,+1,-1) &
               +              Appp(iii-1,jjj,kkk) * p( 0,+1,-1)
             ap(0,0,0) = &
               &              Ammp(iii,jjj,kkk) * p(-1,-1,-1) &
               +              A0mp(iii,jjj,kkk) * p( 0,-1,-1) &
               +              Apmp(iii,jjj,kkk) * p(+1,-1,-1) &
               +              Am0p(iii,jjj,kkk) * p(-1, 0,-1) &
               +              A00p(iii,jjj,kkk) * p( 0, 0,-1) &
               +              Ap0p(iii,jjj,kkk) * p(+1, 0,-1) &
               +              Ampp(iii,jjj,kkk) * p(-1,+1,-1) &
               +              A0pp(iii,jjj,kkk) * p( 0,+1,-1) &
               +              Appp(iii,jjj,kkk) * p(+1,+1,-1)
             ap(1,0,0) = &
               &              Ammp(iii+1,jjj,kkk) * p( 0,-1,-1) &
               +              A0mp(iii+1,jjj,kkk) * p(+1,-1,-1) &
               +              Am0p(iii+1,jjj,kkk) * p( 0, 0,-1) &
               +              A00p(iii+1,jjj,kkk) * p(+1, 0,-1) &
               +              Ampp(iii+1,jjj,kkk) * p( 0,+1,-1) &
               +              A0pp(iii+1,jjj,kkk) * p(+1,+1,-1)
             ap(-1,1,0) = &
               &              A0mp(iii-1,jjj+1,kkk) * p(-1, 0,-1) &
               +              Apmp(iii-1,jjj+1,kkk) * p( 0, 0,-1) &
               +              A00p(iii-1,jjj+1,kkk) * p(-1,+1,-1) &
               +              Ap0p(iii-1,jjj+1,kkk) * p( 0,+1,-1)
             ap(0,1,0) = &
               &              Ammp(iii,jjj+1,kkk) * p(-1, 0,-1) &
               +              A0mp(iii,jjj+1,kkk) * p( 0, 0,-1) &
               +              Apmp(iii,jjj+1,kkk) * p(+1, 0,-1) &
               +              Am0p(iii,jjj+1,kkk) * p(-1,+1,-1) &
               +              A00p(iii,jjj+1,kkk) * p( 0,+1,-1) &
               +              Ap0p(iii,jjj+1,kkk) * p(+1,+1,-1)
             ap(1,1,0) = &
               &              Ammp(iii+1,jjj+1,kkk) * p( 0, 0,-1) &
               +              A0mp(iii+1,jjj+1,kkk) * p(+1, 0,-1) &
               +              Am0p(iii+1,jjj+1,kkk) * p( 0,+1,-1) &
               +              A00p(iii+1,jjj+1,kkk) * p(+1,+1,-1)
             ap(-1,-1,1) = &
               &              A000(iii-1,jjj-1,kkk+1) * p(-1,-1,-1) &
               +              Ap00(iii-1,jjj-1,kkk+1) * p( 0,-1,-1) &
               +              A0p0(iii-1,jjj-1,kkk+1) * p(-1, 0,-1) &
               +              App0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1) &
               +              A00p(iii-1,jjj-1,kkk+1) * p(-1,-1, 0) &
               +              Ap0p(iii-1,jjj-1,kkk+1) * p( 0,-1, 0) &
               +              A0pp(iii-1,jjj-1,kkk+1) * p(-1, 0, 0) &
               +              Appp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0)
             ap(0,-1,1) = &
               &              Am00(iii,jjj-1,kkk+1) * p(-1,-1,-1) &
               +              A000(iii,jjj-1,kkk+1) * p( 0,-1,-1) &
               +              Ap00(iii,jjj-1,kkk+1) * p(+1,-1,-1) &
               +              Amp0(iii,jjj-1,kkk+1) * p(-1, 0,-1) &
               +              A0p0(iii,jjj-1,kkk+1) * p( 0, 0,-1) &
               +              App0(iii,jjj-1,kkk+1) * p(+1, 0,-1) &
               +              Am0p(iii,jjj-1,kkk+1) * p(-1,-1, 0) &
               +              A00p(iii,jjj-1,kkk+1) * p( 0,-1, 0) &
               +              Ap0p(iii,jjj-1,kkk+1) * p(+1,-1, 0) &
               +              Ampp(iii,jjj-1,kkk+1) * p(-1, 0, 0) &
               +              A0pp(iii,jjj-1,kkk+1) * p( 0, 0, 0) &
               +              Appp(iii,jjj-1,kkk+1) * p(+1, 0, 0)
             ap(1,-1,1) = &
               &              Am00(iii+1,jjj-1,kkk+1) * p( 0,-1,-1) &
               +              A000(iii+1,jjj-1,kkk+1) * p(+1,-1,-1) &
               +              Amp0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1) &
               +              A0p0(iii+1,jjj-1,kkk+1) * p(+1, 0,-1) &
               +              Am0p(iii+1,jjj-1,kkk+1) * p( 0,-1, 0) &
               +              A00p(iii+1,jjj-1,kkk+1) * p(+1,-1, 0) &
               +              Ampp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0) &
               +              A0pp(iii+1,jjj-1,kkk+1) * p(+1, 0, 0)
             ap(-1,0,1) = &
               &              A0m0(iii-1,jjj,kkk+1) * p(-1,-1,-1) &
               +              Apm0(iii-1,jjj,kkk+1) * p( 0,-1,-1) &
               +              A000(iii-1,jjj,kkk+1) * p(-1, 0,-1) &
               +              Ap00(iii-1,jjj,kkk+1) * p( 0, 0,-1) &
               +              A0p0(iii-1,jjj,kkk+1) * p(-1,+1,-1) &
               +              App0(iii-1,jjj,kkk+1) * p( 0,+1,-1) &
               +              A0mp(iii-1,jjj,kkk+1) * p(-1,-1, 0) &
               +              Apmp(iii-1,jjj,kkk+1) * p( 0,-1, 0) &
               +              A00p(iii-1,jjj,kkk+1) * p(-1, 0, 0) &
               +              Ap0p(iii-1,jjj,kkk+1) * p( 0, 0, 0) &
               +              A0pp(iii-1,jjj,kkk+1) * p(-1,+1, 0) &
               +              Appp(iii-1,jjj,kkk+1) * p( 0,+1, 0)
             ap(0,0,1) = &
               &              Amm0(iii,jjj,kkk+1) * p(-1,-1,-1) &
               +              A0m0(iii,jjj,kkk+1) * p( 0,-1,-1) &
               +              Apm0(iii,jjj,kkk+1) * p(+1,-1,-1) &
               +              Am00(iii,jjj,kkk+1) * p(-1, 0,-1) &
               +              A000(iii,jjj,kkk+1) * p( 0, 0,-1) &
               +              Ap00(iii,jjj,kkk+1) * p(+1, 0,-1) &
               +              Amp0(iii,jjj,kkk+1) * p(-1,+1,-1) &
               +              A0p0(iii,jjj,kkk+1) * p( 0,+1,-1) &
               +              App0(iii,jjj,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii,jjj,kkk+1) * p(-1,-1, 0) &
               +              A0mp(iii,jjj,kkk+1) * p( 0,-1, 0) &
               +              Apmp(iii,jjj,kkk+1) * p(+1,-1, 0) &
               +              Am0p(iii,jjj,kkk+1) * p(-1, 0, 0) &
               +              A00p(iii,jjj,kkk+1) * p( 0, 0, 0) &
               +              Ap0p(iii,jjj,kkk+1) * p(+1, 0, 0) &
               +              Ampp(iii,jjj,kkk+1) * p(-1,+1, 0) &
               +              A0pp(iii,jjj,kkk+1) * p( 0,+1, 0) &
               +              Appp(iii,jjj,kkk+1) * p(+1,+1, 0)
             ap(1,0,1) = &
               &              Amm0(iii+1,jjj,kkk+1) * p( 0,-1,-1) &
               +              A0m0(iii+1,jjj,kkk+1) * p(+1,-1,-1) &
               +              Am00(iii+1,jjj,kkk+1) * p( 0, 0,-1) &
               +              A000(iii+1,jjj,kkk+1) * p(+1, 0,-1) &
               +              Amp0(iii+1,jjj,kkk+1) * p( 0,+1,-1) &
               +              A0p0(iii+1,jjj,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii+1,jjj,kkk+1) * p( 0,-1, 0) &
               +              A0mp(iii+1,jjj,kkk+1) * p(+1,-1, 0) &
               +              Am0p(iii+1,jjj,kkk+1) * p( 0, 0, 0) &
               +              A00p(iii+1,jjj,kkk+1) * p(+1, 0, 0) &
               +              Ampp(iii+1,jjj,kkk+1) * p( 0,+1, 0) &
               +              A0pp(iii+1,jjj,kkk+1) * p(+1,+1, 0)
             ap(-1,1,1) = &
               &              A0m0(iii-1,jjj+1,kkk+1) * p(-1, 0,-1) &
               +              Apm0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1) &
               +              A000(iii-1,jjj+1,kkk+1) * p(-1,+1,-1) &
               +              Ap00(iii-1,jjj+1,kkk+1) * p( 0,+1,-1) &
               +              A0mp(iii-1,jjj+1,kkk+1) * p(-1, 0, 0) &
               +              Apmp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              A00p(iii-1,jjj+1,kkk+1) * p(-1,+1, 0) &
               +              Ap0p(iii-1,jjj+1,kkk+1) * p( 0,+1, 0)
             ap(0,1,1) = &
               &              Amm0(iii,jjj+1,kkk+1) * p(-1, 0,-1) &
               +              A0m0(iii,jjj+1,kkk+1) * p( 0, 0,-1) &
               +              Apm0(iii,jjj+1,kkk+1) * p(+1, 0,-1) &
               +              Am00(iii,jjj+1,kkk+1) * p(-1,+1,-1) &
               +              A000(iii,jjj+1,kkk+1) * p( 0,+1,-1) &
               +              Ap00(iii,jjj+1,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii,jjj+1,kkk+1) * p(-1, 0, 0) &
               +              A0mp(iii,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              Apmp(iii,jjj+1,kkk+1) * p(+1, 0, 0) &
               +              Am0p(iii,jjj+1,kkk+1) * p(-1,+1, 0) &
               +              A00p(iii,jjj+1,kkk+1) * p( 0,+1, 0) &
               +              Ap0p(iii,jjj+1,kkk+1) * p(+1,+1, 0)
             ap(1,1,1) = &
               &              Amm0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1) &
               +              A0m0(iii+1,jjj+1,kkk+1) * p(+1, 0,-1) &
               +              Am00(iii+1,jjj+1,kkk+1) * p( 0,+1,-1) &
               +              A000(iii+1,jjj+1,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              A0mp(iii+1,jjj+1,kkk+1) * p(+1, 0, 0) &
               +              Am0p(iii+1,jjj+1,kkk+1) * p( 0,+1, 0) &
               +              A00p(iii+1,jjj+1,kkk+1) * p(+1,+1, 0)
             csten(i,j,k,ist_00p) = 0.125d0 * &
               ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0) &
               + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0) &
               + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0) &
               + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0) &
               + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0) &
               + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0) &
               + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0) &
               + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0) &
               + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0) &
               + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1) &
               + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1) &
               + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1) &
               + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1) &
               + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1) &
               + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1) &
               + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1) &
               + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1) &
               + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1))

             ! csten(i,j,k,ist_pp0)
             iii = ii
             jjj = jj
             kkk = kk
             p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj+1,kkk-1)
             p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj+1,kkk-1)
             p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj+2,kkk-1)
             p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj+2,kkk-1)
             p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj+1,kkk  )
             p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj+1,kkk  )
             p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj+2,kkk  )
             p( 0, 0, 0) = 1.d0
             p(-1,-1,+1) = interp_from_ppm_to(iii+1,jjj+1,kkk+1)
             p( 0,-1,+1) = interp_from_0pm_to(iii+2,jjj+1,kkk+1)
             p(-1, 0,+1) = interp_from_p0m_to(iii+1,jjj+2,kkk+1)
             p( 0, 0,+1) = interp_from_00m_to(iii+2,jjj+2,kkk+1)
             ap(0,0,-1) = &
               &              App0(iii,jjj,kkk-1) * p(-1,-1,-1) &
               +              Appp(iii,jjj,kkk-1) * p(-1,-1, 0)
             ap(1,0,-1) = &
               &              A0p0(iii+1,jjj,kkk-1) * p(-1,-1,-1) &
               +              App0(iii+1,jjj,kkk-1) * p( 0,-1,-1) &
               +              A0pp(iii+1,jjj,kkk-1) * p(-1,-1, 0) &
               +              Appp(iii+1,jjj,kkk-1) * p( 0,-1, 0)
             ap(0,1,-1) = &
               &              Ap00(iii,jjj+1,kkk-1) * p(-1,-1,-1) &
               +              App0(iii,jjj+1,kkk-1) * p(-1, 0,-1) &
               +              Ap0p(iii,jjj+1,kkk-1) * p(-1,-1, 0) &
               +              Appp(iii,jjj+1,kkk-1) * p(-1, 0, 0)
             ap(1,1,-1) = &
               &              A000(iii+1,jjj+1,kkk-1) * p(-1,-1,-1) &
               +              Ap00(iii+1,jjj+1,kkk-1) * p( 0,-1,-1) &
               +              A0p0(iii+1,jjj+1,kkk-1) * p(-1, 0,-1) &
               +              App0(iii+1,jjj+1,kkk-1) * p( 0, 0,-1) &
               +              A00p(iii+1,jjj+1,kkk-1) * p(-1,-1, 0) &
               +              Ap0p(iii+1,jjj+1,kkk-1) * p( 0,-1, 0) &
               +              A0pp(iii+1,jjj+1,kkk-1) * p(-1, 0, 0) &
               +              Appp(iii+1,jjj+1,kkk-1) * p( 0, 0, 0)
             ap(0,0,0) = &
               &              Appm(iii,jjj,kkk) * p(-1,-1,-1) &
               +              App0(iii,jjj,kkk) * p(-1,-1, 0) &
               +              Appp(iii,jjj,kkk) * p(-1,-1,+1)
             ap(1,0,0) = &
               &              A0pm(iii+1,jjj,kkk) * p(-1,-1,-1) &
               +              Appm(iii+1,jjj,kkk) * p( 0,-1,-1) &
               +              A0p0(iii+1,jjj,kkk) * p(-1,-1, 0) &
               +              App0(iii+1,jjj,kkk) * p( 0,-1, 0) &
               +              A0pp(iii+1,jjj,kkk) * p(-1,-1,+1) &
               +              Appp(iii+1,jjj,kkk) * p( 0,-1,+1)
             ap(0,1,0) = &
               &              Ap0m(iii,jjj+1,kkk) * p(-1,-1,-1) &
               +              Appm(iii,jjj+1,kkk) * p(-1, 0,-1) &
               +              Ap00(iii,jjj+1,kkk) * p(-1,-1, 0) &
               +              App0(iii,jjj+1,kkk) * p(-1, 0, 0) &
               +              Ap0p(iii,jjj+1,kkk) * p(-1,-1,+1) &
               +              Appp(iii,jjj+1,kkk) * p(-1, 0,+1)
             ap(1,1,0) = &
               &              A00m(iii+1,jjj+1,kkk) * p(-1,-1,-1) &
               +              Ap0m(iii+1,jjj+1,kkk) * p( 0,-1,-1) &
               +              A0pm(iii+1,jjj+1,kkk) * p(-1, 0,-1) &
               +              Appm(iii+1,jjj+1,kkk) * p( 0, 0,-1) &
               +              A000(iii+1,jjj+1,kkk) * p(-1,-1, 0) &
               +              Ap00(iii+1,jjj+1,kkk) * p( 0,-1, 0) &
               +              A0p0(iii+1,jjj+1,kkk) * p(-1, 0, 0) &
               +              App0(iii+1,jjj+1,kkk) * p( 0, 0, 0) &
               +              A00p(iii+1,jjj+1,kkk) * p(-1,-1,+1) &
               +              Ap0p(iii+1,jjj+1,kkk) * p( 0,-1,+1) &
               +              A0pp(iii+1,jjj+1,kkk) * p(-1, 0,+1) &
               +              Appp(iii+1,jjj+1,kkk) * p( 0, 0,+1)
             ap(0,0,1) = &
               &              Appm(iii,jjj,kkk+1) * p(-1,-1, 0) &
               +              App0(iii,jjj,kkk+1) * p(-1,-1,+1)
             ap(1,0,1) = &
               &              A0pm(iii+1,jjj,kkk+1) * p(-1,-1, 0) &
               +              Appm(iii+1,jjj,kkk+1) * p( 0,-1, 0) &
               +              A0p0(iii+1,jjj,kkk+1) * p(-1,-1,+1) &
               +              App0(iii+1,jjj,kkk+1) * p( 0,-1,+1)
             ap(0,1,1) = &
               &              Ap0m(iii,jjj+1,kkk+1) * p(-1,-1, 0) &
               +              Appm(iii,jjj+1,kkk+1) * p(-1, 0, 0) &
               +              Ap00(iii,jjj+1,kkk+1) * p(-1,-1,+1) &
               +              App0(iii,jjj+1,kkk+1) * p(-1, 0,+1)
             ap(1,1,1) = &
               &              A00m(iii+1,jjj+1,kkk+1) * p(-1,-1, 0) &
               +              Ap0m(iii+1,jjj+1,kkk+1) * p( 0,-1, 0) &
               +              A0pm(iii+1,jjj+1,kkk+1) * p(-1, 0, 0) &
               +              Appm(iii+1,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              A000(iii+1,jjj+1,kkk+1) * p(-1,-1,+1) &
               +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,-1,+1) &
               +              A0p0(iii+1,jjj+1,kkk+1) * p(-1, 0,+1) &
               +              App0(iii+1,jjj+1,kkk+1) * p( 0, 0,+1)
             cs1 = 0.125d0 * &
               ( restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1) &
               + restrict_from_p0m_to(iii,jjj,kkk) * ap(+1, 0,-1) &
               + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1) &
               + restrict_from_ppm_to(iii,jjj,kkk) * ap(+1,+1,-1) &
               + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0) &
               + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0) &
               + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0) &
               + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0) &
               + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1) &
               + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1) &
               + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1) &
               + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1))

             ! alternative: csten(i+1,j,k,ist_mp0)
             iii = ii+2
             jjj = jj
             kkk = kk
             p( 0,-1,-1) = interp_from_0pp_to(iii-2,jjj+1,kkk-1)
             p(+1,-1,-1) = interp_from_mpp_to(iii-1,jjj+1,kkk-1)
             p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj+2,kkk-1)
             p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj+2,kkk-1)
             p( 0,-1, 0) = interp_from_0p0_to(iii-2,jjj+1,kkk  )
             p(+1,-1, 0) = interp_from_mp0_to(iii-1,jjj+1,kkk  )
             p( 0, 0, 0) = 1.d0
             p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj+2,kkk  )
             p( 0,-1,+1) = interp_from_0pm_to(iii-2,jjj+1,kkk+1)
             p(+1,-1,+1) = interp_from_mpm_to(iii-1,jjj+1,kkk+1)
             p( 0, 0,+1) = interp_from_00m_to(iii-2,jjj+2,kkk+1)
             p(+1, 0,+1) = interp_from_m0m_to(iii-1,jjj+2,kkk+1)
             ap(-1,0,-1) = &
               &              Amp0(iii-1,jjj,kkk-1) * p( 0,-1,-1) &
               +              A0p0(iii-1,jjj,kkk-1) * p(+1,-1,-1) &
               +              Ampp(iii-1,jjj,kkk-1) * p( 0,-1, 0) &
               +              A0pp(iii-1,jjj,kkk-1) * p(+1,-1, 0)
             ap(0,0,-1) = &
               &              Amp0(iii,jjj,kkk-1) * p(+1,-1,-1) &
               +              Ampp(iii,jjj,kkk-1) * p(+1,-1, 0)
             ap(-1,1,-1) = &
               &              Am00(iii-1,jjj+1,kkk-1) * p( 0,-1,-1) &
               +              A000(iii-1,jjj+1,kkk-1) * p(+1,-1,-1) &
               +              Amp0(iii-1,jjj+1,kkk-1) * p( 0, 0,-1) &
               +              A0p0(iii-1,jjj+1,kkk-1) * p(+1, 0,-1) &
               +              Am0p(iii-1,jjj+1,kkk-1) * p( 0,-1, 0) &
               +              A00p(iii-1,jjj+1,kkk-1) * p(+1,-1, 0) &
               +              Ampp(iii-1,jjj+1,kkk-1) * p( 0, 0, 0) &
               +              A0pp(iii-1,jjj+1,kkk-1) * p(+1, 0, 0)
             ap(0,1,-1) = &
               &              Am00(iii,jjj+1,kkk-1) * p(+1,-1,-1) &
               +              Amp0(iii,jjj+1,kkk-1) * p(+1, 0,-1) &
               +              Am0p(iii,jjj+1,kkk-1) * p(+1,-1, 0) &
               +              Ampp(iii,jjj+1,kkk-1) * p(+1, 0, 0)
             ap(-1,0,0) = &
               &              Ampm(iii-1,jjj,kkk) * p( 0,-1,-1) &
               +              A0pm(iii-1,jjj,kkk) * p(+1,-1,-1) &
               +              Amp0(iii-1,jjj,kkk) * p( 0,-1, 0) &
               +              A0p0(iii-1,jjj,kkk) * p(+1,-1, 0) &
               +              Ampp(iii-1,jjj,kkk) * p( 0,-1,+1) &
               +              A0pp(iii-1,jjj,kkk) * p(+1,-1,+1)
             ap(0,0,0) = &
               &              Ampm(iii,jjj,kkk) * p(+1,-1,-1) &
               +              Amp0(iii,jjj,kkk) * p(+1,-1, 0) &
               +              Ampp(iii,jjj,kkk) * p(+1,-1,+1)
             ap(-1,1,0) = &
               &              Am0m(iii-1,jjj+1,kkk) * p( 0,-1,-1) &
               +              A00m(iii-1,jjj+1,kkk) * p(+1,-1,-1) &
               +              Ampm(iii-1,jjj+1,kkk) * p( 0, 0,-1) &
               +              A0pm(iii-1,jjj+1,kkk) * p(+1, 0,-1) &
               +              Am00(iii-1,jjj+1,kkk) * p( 0,-1, 0) &
               +              A000(iii-1,jjj+1,kkk) * p(+1,-1, 0) &
               +              Amp0(iii-1,jjj+1,kkk) * p( 0, 0, 0) &
               +              A0p0(iii-1,jjj+1,kkk) * p(+1, 0, 0) &
               +              Am0p(iii-1,jjj+1,kkk) * p( 0,-1,+1) &
               +              A00p(iii-1,jjj+1,kkk) * p(+1,-1,+1) &
               +              Ampp(iii-1,jjj+1,kkk) * p( 0, 0,+1) &
               +              A0pp(iii-1,jjj+1,kkk) * p(+1, 0,+1)
             ap(0,1,0) = &
               &              Am0m(iii,jjj+1,kkk) * p(+1,-1,-1) &
               +              Ampm(iii,jjj+1,kkk) * p(+1, 0,-1) &
               +              Am00(iii,jjj+1,kkk) * p(+1,-1, 0) &
               +              Amp0(iii,jjj+1,kkk) * p(+1, 0, 0) &
               +              Am0p(iii,jjj+1,kkk) * p(+1,-1,+1) &
               +              Ampp(iii,jjj+1,kkk) * p(+1, 0,+1)
             ap(-1,0,1) = &
               &              Ampm(iii-1,jjj,kkk+1) * p( 0,-1, 0) &
               +              A0pm(iii-1,jjj,kkk+1) * p(+1,-1, 0) &
               +              Amp0(iii-1,jjj,kkk+1) * p( 0,-1,+1) &
               +              A0p0(iii-1,jjj,kkk+1) * p(+1,-1,+1)
             ap(0,0,1) = &
               &              Ampm(iii,jjj,kkk+1) * p(+1,-1, 0) &
               +              Amp0(iii,jjj,kkk+1) * p(+1,-1,+1)
             ap(-1,1,1) = &
               &              Am0m(iii-1,jjj+1,kkk+1) * p( 0,-1, 0) &
               +              A00m(iii-1,jjj+1,kkk+1) * p(+1,-1, 0) &
               +              Ampm(iii-1,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              A0pm(iii-1,jjj+1,kkk+1) * p(+1, 0, 0) &
               +              Am00(iii-1,jjj+1,kkk+1) * p( 0,-1,+1) &
               +              A000(iii-1,jjj+1,kkk+1) * p(+1,-1,+1) &
               +              Amp0(iii-1,jjj+1,kkk+1) * p( 0, 0,+1) &
               +              A0p0(iii-1,jjj+1,kkk+1) * p(+1, 0,+1)
             ap(0,1,1) = &
               &              Am0m(iii,jjj+1,kkk+1) * p(+1,-1, 0) &
               +              Ampm(iii,jjj+1,kkk+1) * p(+1, 0, 0) &
               +              Am00(iii,jjj+1,kkk+1) * p(+1,-1,+1) &
               +              Amp0(iii,jjj+1,kkk+1) * p(+1, 0,+1)
             cs2 = 0.125d0 * &
               ( restrict_from_m0m_to(iii,jjj,kkk) * ap(-1, 0,-1) &
               + restrict_from_00m_to(iii,jjj,kkk) * ap( 0, 0,-1) &
               + restrict_from_mpm_to(iii,jjj,kkk) * ap(-1,+1,-1) &
               + restrict_from_0pm_to(iii,jjj,kkk) * ap( 0,+1,-1) &
               + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0) &
               + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0) &
               + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0) &
               + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0) &
               + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1) &
               + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1) &
               + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1) &
               + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1))

             csten(i,j,k,ist_pp0) = 0.5d0*(cs1 + cs2)

             ! csten(i,j,k,ist_p0p)
             iii = ii
             jjj = jj
             kkk = kk
             p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj-1,kkk+1)
             p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj-1,kkk+1)
             p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj  ,kkk+1)
             p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj  ,kkk+1)
             p(-1,+1,-1) = interp_from_pmp_to(iii+1,jjj+1,kkk+1)
             p( 0,+1,-1) = interp_from_0mp_to(iii+2,jjj+1,kkk+1)
             p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj-1,kkk+2)
             p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj-1,kkk+2)
             p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj  ,kkk+2)
             p( 0, 0, 0) = 1.d0
             p(-1,+1, 0) = interp_from_pm0_to(iii+1,jjj+1,kkk+2)
             p( 0,+1, 0) = interp_from_0m0_to(iii+2,jjj+1,kkk+2)
             ap(0,-1,0) = &
               &              Ap0p(iii,jjj-1,kkk) * p(-1,-1,-1) &
               +              Appp(iii,jjj-1,kkk) * p(-1, 0,-1)
             ap(1,-1,0) = &
               &              A00p(iii+1,jjj-1,kkk) * p(-1,-1,-1) &
               +              Ap0p(iii+1,jjj-1,kkk) * p( 0,-1,-1) &
               +              A0pp(iii+1,jjj-1,kkk) * p(-1, 0,-1) &
               +              Appp(iii+1,jjj-1,kkk) * p( 0, 0,-1)
             ap(0,0,0) = &
               &              Apmp(iii,jjj,kkk) * p(-1,-1,-1) &
               +              Ap0p(iii,jjj,kkk) * p(-1, 0,-1) &
               +              Appp(iii,jjj,kkk) * p(-1,+1,-1)
             ap(1,0,0) = &
               &              A0mp(iii+1,jjj,kkk) * p(-1,-1,-1) &
               +              Apmp(iii+1,jjj,kkk) * p( 0,-1,-1) &
               +              A00p(iii+1,jjj,kkk) * p(-1, 0,-1) &
               +              Ap0p(iii+1,jjj,kkk) * p( 0, 0,-1) &
               +              A0pp(iii+1,jjj,kkk) * p(-1,+1,-1) &
               +              Appp(iii+1,jjj,kkk) * p( 0,+1,-1)
             ap(0,1,0) = &
               &              Apmp(iii,jjj+1,kkk) * p(-1, 0,-1) &
               +              Ap0p(iii,jjj+1,kkk) * p(-1,+1,-1)
             ap(1,1,0) = &
               &              A0mp(iii+1,jjj+1,kkk) * p(-1, 0,-1) &
               +              Apmp(iii+1,jjj+1,kkk) * p( 0, 0,-1) &
               +              A00p(iii+1,jjj+1,kkk) * p(-1,+1,-1) &
               +              Ap0p(iii+1,jjj+1,kkk) * p( 0,+1,-1)
             ap(0,-1,1) = &
               &              Ap00(iii,jjj-1,kkk+1) * p(-1,-1,-1) &
               +              App0(iii,jjj-1,kkk+1) * p(-1, 0,-1) &
               +              Ap0p(iii,jjj-1,kkk+1) * p(-1,-1, 0) &
               +              Appp(iii,jjj-1,kkk+1) * p(-1, 0, 0)
             ap(1,-1,1) = &
               &              A000(iii+1,jjj-1,kkk+1) * p(-1,-1,-1) &
               +              Ap00(iii+1,jjj-1,kkk+1) * p( 0,-1,-1) &
               +              A0p0(iii+1,jjj-1,kkk+1) * p(-1, 0,-1) &
               +              App0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1) &
               +              A00p(iii+1,jjj-1,kkk+1) * p(-1,-1, 0) &
               +              Ap0p(iii+1,jjj-1,kkk+1) * p( 0,-1, 0) &
               +              A0pp(iii+1,jjj-1,kkk+1) * p(-1, 0, 0) &
               +              Appp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0)
             ap(0,0,1) = &
               &              Apm0(iii,jjj,kkk+1) * p(-1,-1,-1) &
               +              Ap00(iii,jjj,kkk+1) * p(-1, 0,-1) &
               +              App0(iii,jjj,kkk+1) * p(-1,+1,-1) &
               +              Apmp(iii,jjj,kkk+1) * p(-1,-1, 0) &
               +              Ap0p(iii,jjj,kkk+1) * p(-1, 0, 0) &
               +              Appp(iii,jjj,kkk+1) * p(-1,+1, 0)
             ap(1,0,1) = &
               &              A0m0(iii+1,jjj,kkk+1) * p(-1,-1,-1) &
               +              Apm0(iii+1,jjj,kkk+1) * p( 0,-1,-1) &
               +              A000(iii+1,jjj,kkk+1) * p(-1, 0,-1) &
               +              Ap00(iii+1,jjj,kkk+1) * p( 0, 0,-1) &
               +              A0p0(iii+1,jjj,kkk+1) * p(-1,+1,-1) &
               +              App0(iii+1,jjj,kkk+1) * p( 0,+1,-1) &
               +              A0mp(iii+1,jjj,kkk+1) * p(-1,-1, 0) &
               +              Apmp(iii+1,jjj,kkk+1) * p( 0,-1, 0) &
               +              A00p(iii+1,jjj,kkk+1) * p(-1, 0, 0) &
               +              Ap0p(iii+1,jjj,kkk+1) * p( 0, 0, 0) &
               +              A0pp(iii+1,jjj,kkk+1) * p(-1,+1, 0) &
               +              Appp(iii+1,jjj,kkk+1) * p( 0,+1, 0)
             ap(0,1,1) = &
               &              Apm0(iii,jjj+1,kkk+1) * p(-1, 0,-1) &
               +              Ap00(iii,jjj+1,kkk+1) * p(-1,+1,-1) &
               +              Apmp(iii,jjj+1,kkk+1) * p(-1, 0, 0) &
               +              Ap0p(iii,jjj+1,kkk+1) * p(-1,+1, 0)
             ap(1,1,1) = &
               &              A0m0(iii+1,jjj+1,kkk+1) * p(-1, 0,-1) &
               +              Apm0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1) &
               +              A000(iii+1,jjj+1,kkk+1) * p(-1,+1,-1) &
               +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,+1,-1) &
               +              A0mp(iii+1,jjj+1,kkk+1) * p(-1, 0, 0) &
               +              Apmp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              A00p(iii+1,jjj+1,kkk+1) * p(-1,+1, 0) &
               +              Ap0p(iii+1,jjj+1,kkk+1) * p( 0,+1, 0)
             cs1 = 0.125d0 * &
               ( restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0) &
               + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0) &
               + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0) &
               + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0) &
               + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0) &
               + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0) &
               + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1) &
               + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1) &
               + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1) &
               + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1) &
               + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1) &
               + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1))

             ! alternative: csten(i+1,j,k,ist_m0p)
             iii = ii+2
             jjj = jj
             kkk = kk
             p( 0,-1,-1) = interp_from_0pp_to(iii-2,jjj-1,kkk+1)
             p(+1,-1,-1) = interp_from_mpp_to(iii-1,jjj-1,kkk+1)
             p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj  ,kkk+1)
             p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj  ,kkk+1)
             p( 0,+1,-1) = interp_from_0mp_to(iii-2,jjj+1,kkk+1)
             p(+1,+1,-1) = interp_from_mmp_to(iii-1,jjj+1,kkk+1)
             p( 0,-1, 0) = interp_from_0p0_to(iii-2,jjj-1,kkk+2)
             p(+1,-1, 0) = interp_from_mp0_to(iii-1,jjj-1,kkk+2)
             p( 0, 0, 0) = 1.d0
             p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj  ,kkk+2)
             p( 0,+1, 0) = interp_from_0m0_to(iii-2,jjj+1,kkk+2)
             p(+1,+1, 0) = interp_from_mm0_to(iii-1,jjj+1,kkk+2)

             ap(-1,-1,0) = &
               &              Am0p(iii-1,jjj-1,kkk) * p( 0,-1,-1) &
               +              A00p(iii-1,jjj-1,kkk) * p(+1,-1,-1) &
               +              Ampp(iii-1,jjj-1,kkk) * p( 0, 0,-1) &
               +              A0pp(iii-1,jjj-1,kkk) * p(+1, 0,-1)
             ap(0,-1,0) = &
               &              Am0p(iii,jjj-1,kkk) * p(+1,-1,-1) &
               +              Ampp(iii,jjj-1,kkk) * p(+1, 0,-1)
             ap(-1,0,0) = &
               &              Ammp(iii-1,jjj,kkk) * p( 0,-1,-1) &
               +              A0mp(iii-1,jjj,kkk) * p(+1,-1,-1) &
               +              Am0p(iii-1,jjj,kkk) * p( 0, 0,-1) &
               +              A00p(iii-1,jjj,kkk) * p(+1, 0,-1) &
               +              Ampp(iii-1,jjj,kkk) * p( 0,+1,-1) &
               +              A0pp(iii-1,jjj,kkk) * p(+1,+1,-1)
             ap(0,0,0) = &
               &              Ammp(iii,jjj,kkk) * p(+1,-1,-1) &
               +              Am0p(iii,jjj,kkk) * p(+1, 0,-1) &
               +              Ampp(iii,jjj,kkk) * p(+1,+1,-1)
             ap(-1,1,0) = &
               &              Ammp(iii-1,jjj+1,kkk) * p( 0, 0,-1) &
               +              A0mp(iii-1,jjj+1,kkk) * p(+1, 0,-1) &
               +              Am0p(iii-1,jjj+1,kkk) * p( 0,+1,-1) &
               +              A00p(iii-1,jjj+1,kkk) * p(+1,+1,-1)
             ap(0,1,0) = &
               &              Ammp(iii,jjj+1,kkk) * p(+1, 0,-1) &
               +              Am0p(iii,jjj+1,kkk) * p(+1,+1,-1)
             ap(-1,-1,1) = &
               &              Am00(iii-1,jjj-1,kkk+1) * p( 0,-1,-1) &
               +              A000(iii-1,jjj-1,kkk+1) * p(+1,-1,-1) &
               +              Amp0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1) &
               +              A0p0(iii-1,jjj-1,kkk+1) * p(+1, 0,-1) &
               +              Am0p(iii-1,jjj-1,kkk+1) * p( 0,-1, 0) &
               +              A00p(iii-1,jjj-1,kkk+1) * p(+1,-1, 0) &
               +              Ampp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0) &
               +              A0pp(iii-1,jjj-1,kkk+1) * p(+1, 0, 0)
             ap(0,-1,1) = &
               &              Am00(iii,jjj-1,kkk+1) * p(+1,-1,-1) &
               +              Amp0(iii,jjj-1,kkk+1) * p(+1, 0,-1) &
               +              Am0p(iii,jjj-1,kkk+1) * p(+1,-1, 0) &
               +              Ampp(iii,jjj-1,kkk+1) * p(+1, 0, 0)
             ap(-1,0,1) = &
               &              Amm0(iii-1,jjj,kkk+1) * p( 0,-1,-1) &
               +              A0m0(iii-1,jjj,kkk+1) * p(+1,-1,-1) &
               +              Am00(iii-1,jjj,kkk+1) * p( 0, 0,-1) &
               +              A000(iii-1,jjj,kkk+1) * p(+1, 0,-1) &
               +              Amp0(iii-1,jjj,kkk+1) * p( 0,+1,-1) &
               +              A0p0(iii-1,jjj,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii-1,jjj,kkk+1) * p( 0,-1, 0) &
               +              A0mp(iii-1,jjj,kkk+1) * p(+1,-1, 0) &
               +              Am0p(iii-1,jjj,kkk+1) * p( 0, 0, 0) &
               +              A00p(iii-1,jjj,kkk+1) * p(+1, 0, 0) &
               +              Ampp(iii-1,jjj,kkk+1) * p( 0,+1, 0) &
               +              A0pp(iii-1,jjj,kkk+1) * p(+1,+1, 0)
             ap(0,0,1) = &
               &              Amm0(iii,jjj,kkk+1) * p(+1,-1,-1) &
               +              Am00(iii,jjj,kkk+1) * p(+1, 0,-1) &
               +              Amp0(iii,jjj,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii,jjj,kkk+1) * p(+1,-1, 0) &
               +              Am0p(iii,jjj,kkk+1) * p(+1, 0, 0) &
               +              Ampp(iii,jjj,kkk+1) * p(+1,+1, 0)
             ap(-1,1,1) = &
               &              Amm0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1) &
               +              A0m0(iii-1,jjj+1,kkk+1) * p(+1, 0,-1) &
               +              Am00(iii-1,jjj+1,kkk+1) * p( 0,+1,-1) &
               +              A000(iii-1,jjj+1,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              A0mp(iii-1,jjj+1,kkk+1) * p(+1, 0, 0) &
               +              Am0p(iii-1,jjj+1,kkk+1) * p( 0,+1, 0) &
               +              A00p(iii-1,jjj+1,kkk+1) * p(+1,+1, 0)
             ap(0,1,1) = &
               &              Amm0(iii,jjj+1,kkk+1) * p(+1, 0,-1) &
               +              Am00(iii,jjj+1,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii,jjj+1,kkk+1) * p(+1, 0, 0) &
               +              Am0p(iii,jjj+1,kkk+1) * p(+1,+1, 0)
             cs2 = 0.125d0 * &
               ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0) &
               + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0) &
               + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0) &
               + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0) &
               + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0) &
               + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0) &
               + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1) &
               + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1) &
               + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1) &
               + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1) &
               + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1) &
               + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1))

             csten(i,j,k,ist_p0p) = 0.5d0*(cs1+cs2)

             ! csten(i,j,k,ist_0pp)
             iii = ii
             jjj = jj
             kkk = kk
             p(-1,-1,-1) = interp_from_ppp_to(iii-1,jjj+1,kkk+1)
             p( 0,-1,-1) = interp_from_0pp_to(iii  ,jjj+1,kkk+1)
             p(+1,-1,-1) = interp_from_mpp_to(iii+1,jjj+1,kkk+1)
             p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj+2,kkk+1)
             p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj+2,kkk+1)
             p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj+2,kkk+1)
             p(-1,-1, 0) = interp_from_pp0_to(iii-1,jjj+1,kkk+2)
             p( 0,-1, 0) = interp_from_0p0_to(iii  ,jjj+1,kkk+2)
             p(+1,-1, 0) = interp_from_mp0_to(iii+1,jjj+1,kkk+2)
             p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj+2,kkk+2)
             p( 0, 0, 0) = 1.d0
             p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj+2,kkk+2)
             ap(-1,0,0) = &
               &              A0pp(iii-1,jjj,kkk) * p(-1,-1,-1) &
               +              Appp(iii-1,jjj,kkk) * p( 0,-1,-1)
             ap(0,0,0) = &
               &              Ampp(iii,jjj,kkk) * p(-1,-1,-1) &
               +              A0pp(iii,jjj,kkk) * p( 0,-1,-1) &
               +              Appp(iii,jjj,kkk) * p(+1,-1,-1)
             ap(1,0,0) = &
               &              Ampp(iii+1,jjj,kkk) * p( 0,-1,-1) &
               +              A0pp(iii+1,jjj,kkk) * p(+1,-1,-1)
             ap(-1,1,0) = &
               &              A00p(iii-1,jjj+1,kkk) * p(-1,-1,-1) &
               +              Ap0p(iii-1,jjj+1,kkk) * p( 0,-1,-1) &
               +              A0pp(iii-1,jjj+1,kkk) * p(-1, 0,-1) &
               +              Appp(iii-1,jjj+1,kkk) * p( 0, 0,-1)
             ap(0,1,0) = &
               &              Am0p(iii,jjj+1,kkk) * p(-1,-1,-1) &
               +              A00p(iii,jjj+1,kkk) * p( 0,-1,-1) &
               +              Ap0p(iii,jjj+1,kkk) * p(+1,-1,-1) &
               +              Ampp(iii,jjj+1,kkk) * p(-1, 0,-1) &
               +              A0pp(iii,jjj+1,kkk) * p( 0, 0,-1) &
               +              Appp(iii,jjj+1,kkk) * p(+1, 0,-1)
             ap(1,1,0) = &
               &              Am0p(iii+1,jjj+1,kkk) * p( 0,-1,-1) &
               +              A00p(iii+1,jjj+1,kkk) * p(+1,-1,-1) &
               +              Ampp(iii+1,jjj+1,kkk) * p( 0, 0,-1) &
               +              A0pp(iii+1,jjj+1,kkk) * p(+1, 0,-1)
             ap(-1,0,1) = &
               &              A0p0(iii-1,jjj,kkk+1) * p(-1,-1,-1) &
               +              App0(iii-1,jjj,kkk+1) * p( 0,-1,-1) &
               +              A0pp(iii-1,jjj,kkk+1) * p(-1,-1, 0) &
               +              Appp(iii-1,jjj,kkk+1) * p( 0,-1, 0)
             ap(0,0,1) = &
               &              Amp0(iii,jjj,kkk+1) * p(-1,-1,-1) &
               +              A0p0(iii,jjj,kkk+1) * p( 0,-1,-1) &
               +              App0(iii,jjj,kkk+1) * p(+1,-1,-1) &
               +              Ampp(iii,jjj,kkk+1) * p(-1,-1, 0) &
               +              A0pp(iii,jjj,kkk+1) * p( 0,-1, 0) &
               +              Appp(iii,jjj,kkk+1) * p(+1,-1, 0)
             ap(1,0,1) = &
               &              Amp0(iii+1,jjj,kkk+1) * p( 0,-1,-1) &
               +              A0p0(iii+1,jjj,kkk+1) * p(+1,-1,-1) &
               +              Ampp(iii+1,jjj,kkk+1) * p( 0,-1, 0) &
               +              A0pp(iii+1,jjj,kkk+1) * p(+1,-1, 0)
             ap(-1,1,1) = &
               &              A000(iii-1,jjj+1,kkk+1) * p(-1,-1,-1) &
               +              Ap00(iii-1,jjj+1,kkk+1) * p( 0,-1,-1) &
               +              A0p0(iii-1,jjj+1,kkk+1) * p(-1, 0,-1) &
               +              App0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1) &
               +              A00p(iii-1,jjj+1,kkk+1) * p(-1,-1, 0) &
               +              Ap0p(iii-1,jjj+1,kkk+1) * p( 0,-1, 0) &
               +              A0pp(iii-1,jjj+1,kkk+1) * p(-1, 0, 0) &
               +              Appp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0)
             ap(0,1,1) = &
               &              Am00(iii,jjj+1,kkk+1) * p(-1,-1,-1) &
               +              A000(iii,jjj+1,kkk+1) * p( 0,-1,-1) &
               +              Ap00(iii,jjj+1,kkk+1) * p(+1,-1,-1) &
               +              Amp0(iii,jjj+1,kkk+1) * p(-1, 0,-1) &
               +              A0p0(iii,jjj+1,kkk+1) * p( 0, 0,-1) &
               +              App0(iii,jjj+1,kkk+1) * p(+1, 0,-1) &
               +              Am0p(iii,jjj+1,kkk+1) * p(-1,-1, 0) &
               +              A00p(iii,jjj+1,kkk+1) * p( 0,-1, 0) &
               +              Ap0p(iii,jjj+1,kkk+1) * p(+1,-1, 0) &
               +              Ampp(iii,jjj+1,kkk+1) * p(-1, 0, 0) &
               +              A0pp(iii,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              Appp(iii,jjj+1,kkk+1) * p(+1, 0, 0)
             ap(1,1,1) = &
               &              Am00(iii+1,jjj+1,kkk+1) * p( 0,-1,-1) &
               +              A000(iii+1,jjj+1,kkk+1) * p(+1,-1,-1) &
               +              Amp0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1) &
               +              A0p0(iii+1,jjj+1,kkk+1) * p(+1, 0,-1) &
               +              Am0p(iii+1,jjj+1,kkk+1) * p( 0,-1, 0) &
               +              A00p(iii+1,jjj+1,kkk+1) * p(+1,-1, 0) &
               +              Ampp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              A0pp(iii+1,jjj+1,kkk+1) * p(+1, 0, 0)
             cs1 = 0.125d0 * &
               ( restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0) &
               + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0) &
               + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0) &
               + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0) &
               + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0) &
               + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0) &
               + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1) &
               + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1) &
               + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1) &
               + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1) &
               + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1) &
               + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1))

             ! alternative: csten(i,j+1,k,ist_0mp)
             iii = ii
             jjj = jj+2
             kkk = kk
             p(-1, 0,-1) = interp_from_p0p_to(iii-1,jjj-2,kkk+1)
             p( 0, 0,-1) = interp_from_00p_to(iii  ,jjj-2,kkk+1)
             p(+1, 0,-1) = interp_from_m0p_to(iii+1,jjj-2,kkk+1)
             p(-1,+1,-1) = interp_from_pmp_to(iii-1,jjj-1,kkk+1)
             p( 0,+1,-1) = interp_from_0mp_to(iii  ,jjj-1,kkk+1)
             p(+1,+1,-1) = interp_from_mmp_to(iii+1,jjj-1,kkk+1)
             p(-1, 0, 0) = interp_from_p00_to(iii-1,jjj-2,kkk+2)
             p( 0, 0, 0) = 1.d0
             p(+1, 0, 0) = interp_from_m00_to(iii+1,jjj-2,kkk+2)
             p(-1,+1, 0) = interp_from_pm0_to(iii-1,jjj-1,kkk+2)
             p( 0,+1, 0) = interp_from_0m0_to(iii  ,jjj-1,kkk+2)
             p(+1,+1, 0) = interp_from_mm0_to(iii+1,jjj-1,kkk+2)
             ap(-1,-1,0) = &
               &              A0mp(iii-1,jjj-1,kkk) * p(-1, 0,-1) &
               +              Apmp(iii-1,jjj-1,kkk) * p( 0, 0,-1) &
               +              A00p(iii-1,jjj-1,kkk) * p(-1,+1,-1) &
               +              Ap0p(iii-1,jjj-1,kkk) * p( 0,+1,-1)
             ap(0,-1,0) = &
               &              Ammp(iii,jjj-1,kkk) * p(-1, 0,-1) &
               +              A0mp(iii,jjj-1,kkk) * p( 0, 0,-1) &
               +              Apmp(iii,jjj-1,kkk) * p(+1, 0,-1) &
               +              Am0p(iii,jjj-1,kkk) * p(-1,+1,-1) &
               +              A00p(iii,jjj-1,kkk) * p( 0,+1,-1) &
               +              Ap0p(iii,jjj-1,kkk) * p(+1,+1,-1)
             ap(1,-1,0) = &
               &              Ammp(iii+1,jjj-1,kkk) * p( 0, 0,-1) &
               +              A0mp(iii+1,jjj-1,kkk) * p(+1, 0,-1) &
               +              Am0p(iii+1,jjj-1,kkk) * p( 0,+1,-1) &
               +              A00p(iii+1,jjj-1,kkk) * p(+1,+1,-1)
             ap(-1,0,0) = &
               &              A0mp(iii-1,jjj,kkk) * p(-1,+1,-1) &
               +              Apmp(iii-1,jjj,kkk) * p( 0,+1,-1)
             ap(0,0,0) = &
               &              Ammp(iii,jjj,kkk) * p(-1,+1,-1) &
               +              A0mp(iii,jjj,kkk) * p( 0,+1,-1) &
               +              Apmp(iii,jjj,kkk) * p(+1,+1,-1)
             ap(1,0,0) = &
               &              Ammp(iii+1,jjj,kkk) * p( 0,+1,-1) &
               +              A0mp(iii+1,jjj,kkk) * p(+1,+1,-1)
             ap(-1,-1,1) = &
               &              A0m0(iii-1,jjj-1,kkk+1) * p(-1, 0,-1) &
               +              Apm0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1) &
               +              A000(iii-1,jjj-1,kkk+1) * p(-1,+1,-1) &
               +              Ap00(iii-1,jjj-1,kkk+1) * p( 0,+1,-1) &
               +              A0mp(iii-1,jjj-1,kkk+1) * p(-1, 0, 0) &
               +              Apmp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0) &
               +              A00p(iii-1,jjj-1,kkk+1) * p(-1,+1, 0) &
               +              Ap0p(iii-1,jjj-1,kkk+1) * p( 0,+1, 0)
             ap(0,-1,1) = &
               &              Amm0(iii,jjj-1,kkk+1) * p(-1, 0,-1) &
               +              A0m0(iii,jjj-1,kkk+1) * p( 0, 0,-1) &
               +              Apm0(iii,jjj-1,kkk+1) * p(+1, 0,-1) &
               +              Am00(iii,jjj-1,kkk+1) * p(-1,+1,-1) &
               +              A000(iii,jjj-1,kkk+1) * p( 0,+1,-1) &
               +              Ap00(iii,jjj-1,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii,jjj-1,kkk+1) * p(-1, 0, 0) &
               +              A0mp(iii,jjj-1,kkk+1) * p( 0, 0, 0) &
               +              Apmp(iii,jjj-1,kkk+1) * p(+1, 0, 0) &
               +              Am0p(iii,jjj-1,kkk+1) * p(-1,+1, 0) &
               +              A00p(iii,jjj-1,kkk+1) * p( 0,+1, 0) &
               +              Ap0p(iii,jjj-1,kkk+1) * p(+1,+1, 0)
             ap(1,-1,1) = &
               &              Amm0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1) &
               +              A0m0(iii+1,jjj-1,kkk+1) * p(+1, 0,-1) &
               +              Am00(iii+1,jjj-1,kkk+1) * p( 0,+1,-1) &
               +              A000(iii+1,jjj-1,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0) &
               +              A0mp(iii+1,jjj-1,kkk+1) * p(+1, 0, 0) &
               +              Am0p(iii+1,jjj-1,kkk+1) * p( 0,+1, 0) &
               +              A00p(iii+1,jjj-1,kkk+1) * p(+1,+1, 0)
             ap(-1,0,1) = &
               &              A0m0(iii-1,jjj,kkk+1) * p(-1,+1,-1) &
               +              Apm0(iii-1,jjj,kkk+1) * p( 0,+1,-1) &
               +              A0mp(iii-1,jjj,kkk+1) * p(-1,+1, 0) &
               +              Apmp(iii-1,jjj,kkk+1) * p( 0,+1, 0)
             ap(0,0,1) = &
               &              Amm0(iii,jjj,kkk+1) * p(-1,+1,-1) &
               +              A0m0(iii,jjj,kkk+1) * p( 0,+1,-1) &
               +              Apm0(iii,jjj,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii,jjj,kkk+1) * p(-1,+1, 0) &
               +              A0mp(iii,jjj,kkk+1) * p( 0,+1, 0) &
               +              Apmp(iii,jjj,kkk+1) * p(+1,+1, 0)
             ap(1,0,1) = &
               &              Amm0(iii+1,jjj,kkk+1) * p( 0,+1,-1) &
               +              A0m0(iii+1,jjj,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii+1,jjj,kkk+1) * p( 0,+1, 0) &
               +              A0mp(iii+1,jjj,kkk+1) * p(+1,+1, 0)
             cs2 = 0.125d0 * &
               ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0) &
               + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0) &
               + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0) &
               + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0) &
               + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0) &
               + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0) &
               + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1) &
               + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1) &
               + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1) &
               + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1) &
               + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1) &
               + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1))

             csten(i,j,k,ist_0pp) = 0.5d0*(cs1+cs2)

             ! csten(i,j,k,ist_ppp)
             iii = ii
             jjj = jj
             kkk = kk
             p(-1,-1,-1) = interp_from_ppp_to(iii+1,jjj+1,kkk+1)
             p( 0,-1,-1) = interp_from_0pp_to(iii+2,jjj+1,kkk+1)
             p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj+2,kkk+1)
             p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj+2,kkk+1)
             p(-1,-1, 0) = interp_from_pp0_to(iii+1,jjj+1,kkk+2)
             p( 0,-1, 0) = interp_from_0p0_to(iii+2,jjj+1,kkk+2)
             p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj+2,kkk+2)
             p( 0, 0, 0) = 1.d0
             ap(0,0,0) = &
               &              Appp(iii,jjj,kkk) * p(-1,-1,-1)
             ap(1,0,0) = &
               &              A0pp(iii+1,jjj,kkk) * p(-1,-1,-1) &
               +              Appp(iii+1,jjj,kkk) * p( 0,-1,-1)
             ap(0,1,0) = &
               &              Ap0p(iii,jjj+1,kkk) * p(-1,-1,-1) &
               +              Appp(iii,jjj+1,kkk) * p(-1, 0,-1)
             ap(1,1,0) = &
               &              A00p(iii+1,jjj+1,kkk) * p(-1,-1,-1) &
               +              Ap0p(iii+1,jjj+1,kkk) * p( 0,-1,-1) &
               +              A0pp(iii+1,jjj+1,kkk) * p(-1, 0,-1) &
               +              Appp(iii+1,jjj+1,kkk) * p( 0, 0,-1)
             ap(0,0,1) = &
               &              App0(iii,jjj,kkk+1) * p(-1,-1,-1) &
               +              Appp(iii,jjj,kkk+1) * p(-1,-1, 0)
             ap(1,0,1) = &
               &              A0p0(iii+1,jjj,kkk+1) * p(-1,-1,-1) &
               +              App0(iii+1,jjj,kkk+1) * p( 0,-1,-1) &
               +              A0pp(iii+1,jjj,kkk+1) * p(-1,-1, 0) &
               +              Appp(iii+1,jjj,kkk+1) * p( 0,-1, 0)
             ap(0,1,1) = &
               &              Ap00(iii,jjj+1,kkk+1) * p(-1,-1,-1) &
               +              App0(iii,jjj+1,kkk+1) * p(-1, 0,-1) &
               +              Ap0p(iii,jjj+1,kkk+1) * p(-1,-1, 0) &
               +              Appp(iii,jjj+1,kkk+1) * p(-1, 0, 0)
             ap(1,1,1) = &
               &              A000(iii+1,jjj+1,kkk+1) * p(-1,-1,-1) &
               +              Ap00(iii+1,jjj+1,kkk+1) * p( 0,-1,-1) &
               +              A0p0(iii+1,jjj+1,kkk+1) * p(-1, 0,-1) &
               +              App0(iii+1,jjj+1,kkk+1) * p( 0, 0,-1) &
               +              A00p(iii+1,jjj+1,kkk+1) * p(-1,-1, 0) &
               +              Ap0p(iii+1,jjj+1,kkk+1) * p( 0,-1, 0) &
               +              A0pp(iii+1,jjj+1,kkk+1) * p(-1, 0, 0) &
               +              Appp(iii+1,jjj+1,kkk+1) * p( 0, 0, 0)
             cs1 = 0.125d0 * &
               ( restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0) &
               + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0) &
               + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0) &
               + restrict_from_pp0_to(iii,jjj,kkk) * ap(+1,+1, 0) &
               + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1) &
               + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1) &
               + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1) &
               + restrict_from_ppp_to(iii,jjj,kkk) * ap(+1,+1,+1))

             ! alternative: csten(i+1,j,k,ist_mpp)
             iii = ii+2
             jjj = jj
             kkk = kk
             p( 0,-1,-1) = interp_from_0pp_to(iii-2,jjj+1,kkk+1)
             p(+1,-1,-1) = interp_from_mpp_to(iii-1,jjj+1,kkk+1)
             p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj+2,kkk+1)
             p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj+2,kkk+1)
             p( 0,-1, 0) = interp_from_0p0_to(iii-2,jjj+1,kkk+2)
             p(+1,-1, 0) = interp_from_mp0_to(iii-1,jjj+1,kkk+2)
             p( 0, 0, 0) = 1.d0
             p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj+2,kkk+2)
             ap(-1,0,0) = &
               &              Ampp(iii-1,jjj,kkk) * p( 0,-1,-1) &
               +              A0pp(iii-1,jjj,kkk) * p(+1,-1,-1)
             ap(0,0,0) = &
               &              Ampp(iii,jjj,kkk) * p(+1,-1,-1)
             ap(-1,1,0) = &
               &              Am0p(iii-1,jjj+1,kkk) * p( 0,-1,-1) &
               +              A00p(iii-1,jjj+1,kkk) * p(+1,-1,-1) &
               +              Ampp(iii-1,jjj+1,kkk) * p( 0, 0,-1) &
               +              A0pp(iii-1,jjj+1,kkk) * p(+1, 0,-1)
             ap(0,1,0) = &
               &              Am0p(iii,jjj+1,kkk) * p(+1,-1,-1) &
               +              Ampp(iii,jjj+1,kkk) * p(+1, 0,-1)
             ap(-1,0,1) = &
               &              Amp0(iii-1,jjj,kkk+1) * p( 0,-1,-1) &
               +              A0p0(iii-1,jjj,kkk+1) * p(+1,-1,-1) &
               +              Ampp(iii-1,jjj,kkk+1) * p( 0,-1, 0) &
               +              A0pp(iii-1,jjj,kkk+1) * p(+1,-1, 0)
             ap(0,0,1) = &
               &              Amp0(iii,jjj,kkk+1) * p(+1,-1,-1) &
               +              Ampp(iii,jjj,kkk+1) * p(+1,-1, 0)
             ap(-1,1,1) = &
               &              Am00(iii-1,jjj+1,kkk+1) * p( 0,-1,-1) &
               +              A000(iii-1,jjj+1,kkk+1) * p(+1,-1,-1) &
               +              Amp0(iii-1,jjj+1,kkk+1) * p( 0, 0,-1) &
               +              A0p0(iii-1,jjj+1,kkk+1) * p(+1, 0,-1) &
               +              Am0p(iii-1,jjj+1,kkk+1) * p( 0,-1, 0) &
               +              A00p(iii-1,jjj+1,kkk+1) * p(+1,-1, 0) &
               +              Ampp(iii-1,jjj+1,kkk+1) * p( 0, 0, 0) &
               +              A0pp(iii-1,jjj+1,kkk+1) * p(+1, 0, 0)
             ap(0,1,1) = &
               &              Am00(iii,jjj+1,kkk+1) * p(+1,-1,-1) &
               +              Amp0(iii,jjj+1,kkk+1) * p(+1, 0,-1) &
               +              Am0p(iii,jjj+1,kkk+1) * p(+1,-1, 0) &
               +              Ampp(iii,jjj+1,kkk+1) * p(+1, 0, 0)
             cs2 = 0.125d0 * &
               ( restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0) &
               + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0) &
               + restrict_from_mp0_to(iii,jjj,kkk) * ap(-1,+1, 0) &
               + restrict_from_0p0_to(iii,jjj,kkk) * ap( 0,+1, 0) &
               + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1) &
               + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1) &
               + restrict_from_mpp_to(iii,jjj,kkk) * ap(-1,+1,+1) &
               + restrict_from_0pp_to(iii,jjj,kkk) * ap( 0,+1,+1))

             ! alternative: csten(i,j+1,k,ist_pmp)
             iii = ii
             jjj = jj+2
             kkk = kk
             p(-1, 0,-1) = interp_from_p0p_to(iii+1,jjj-2,kkk+1)
             p( 0, 0,-1) = interp_from_00p_to(iii+2,jjj-2,kkk+1)
             p(-1,+1,-1) = interp_from_pmp_to(iii+1,jjj-1,kkk+1)
             p( 0,+1,-1) = interp_from_0mp_to(iii+2,jjj-1,kkk+1)
             p(-1, 0, 0) = interp_from_p00_to(iii+1,jjj-2,kkk+2)
             p( 0, 0, 0) = 1.d0
             p(-1,+1, 0) = interp_from_pm0_to(iii+1,jjj-1,kkk+2)
             p( 0,+1, 0) = interp_from_0m0_to(iii+2,jjj-1,kkk+2)
             ap(0,-1,0) = &
               &              Apmp(iii,jjj-1,kkk) * p(-1, 0,-1) &
               +              Ap0p(iii,jjj-1,kkk) * p(-1,+1,-1)
             ap(1,-1,0) = &
               &              A0mp(iii+1,jjj-1,kkk) * p(-1, 0,-1) &
               +              Apmp(iii+1,jjj-1,kkk) * p( 0, 0,-1) &
               +              A00p(iii+1,jjj-1,kkk) * p(-1,+1,-1) &
               +              Ap0p(iii+1,jjj-1,kkk) * p( 0,+1,-1)
             ap(0,0,0) = &
               &              Apmp(iii,jjj,kkk) * p(-1,+1,-1)
             ap(1,0,0) = &
               &              A0mp(iii+1,jjj,kkk) * p(-1,+1,-1) &
               +              Apmp(iii+1,jjj,kkk) * p( 0,+1,-1)
             ap(0,-1,1) = &
               &              Apm0(iii,jjj-1,kkk+1) * p(-1, 0,-1) &
               +              Ap00(iii,jjj-1,kkk+1) * p(-1,+1,-1) &
               +              Apmp(iii,jjj-1,kkk+1) * p(-1, 0, 0) &
               +              Ap0p(iii,jjj-1,kkk+1) * p(-1,+1, 0)
             ap(1,-1,1) = &
               &              A0m0(iii+1,jjj-1,kkk+1) * p(-1, 0,-1) &
               +              Apm0(iii+1,jjj-1,kkk+1) * p( 0, 0,-1) &
               +              A000(iii+1,jjj-1,kkk+1) * p(-1,+1,-1) &
               +              Ap00(iii+1,jjj-1,kkk+1) * p( 0,+1,-1) &
               +              A0mp(iii+1,jjj-1,kkk+1) * p(-1, 0, 0) &
               +              Apmp(iii+1,jjj-1,kkk+1) * p( 0, 0, 0) &
               +              A00p(iii+1,jjj-1,kkk+1) * p(-1,+1, 0) &
               +              Ap0p(iii+1,jjj-1,kkk+1) * p( 0,+1, 0)
             ap(0,0,1) = &
               &              Apm0(iii,jjj,kkk+1) * p(-1,+1,-1) &
               +              Apmp(iii,jjj,kkk+1) * p(-1,+1, 0)
             ap(1,0,1) = &
               &              A0m0(iii+1,jjj,kkk+1) * p(-1,+1,-1) &
               +              Apm0(iii+1,jjj,kkk+1) * p( 0,+1,-1) &
               +              A0mp(iii+1,jjj,kkk+1) * p(-1,+1, 0) &
               +              Apmp(iii+1,jjj,kkk+1) * p( 0,+1, 0)
             cs3 = 0.125d0 * &
               ( restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0) &
               + restrict_from_pm0_to(iii,jjj,kkk) * ap(+1,-1, 0) &
               + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0) &
               + restrict_from_p00_to(iii,jjj,kkk) * ap(+1, 0, 0) &
               + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1) &
               + restrict_from_pmp_to(iii,jjj,kkk) * ap(+1,-1,+1) &
               + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1) &
               + restrict_from_p0p_to(iii,jjj,kkk) * ap(+1, 0,+1))

             ! alternative: csten(i+1,j+1,k,ist_mmp)
             iii = ii+2
             jjj = jj+2
             kkk = kk
             p( 0, 0,-1) = interp_from_00p_to(iii-2,jjj-2,kkk+1)
             p(+1, 0,-1) = interp_from_m0p_to(iii-1,jjj-2,kkk+1)
             p( 0,+1,-1) = interp_from_0mp_to(iii-2,jjj-1,kkk+1)
             p(+1,+1,-1) = interp_from_mmp_to(iii-1,jjj-1,kkk+1)
             p( 0, 0, 0) = 1.d0
             p(+1, 0, 0) = interp_from_m00_to(iii-1,jjj-2,kkk+2)
             p( 0,+1, 0) = interp_from_0m0_to(iii-2,jjj-1,kkk+2)
             p(+1,+1, 0) = interp_from_mm0_to(iii-1,jjj-1,kkk+2)
             ap(-1,-1,0) = &
               &              Ammp(iii-1,jjj-1,kkk) * p( 0, 0,-1) &
               +              A0mp(iii-1,jjj-1,kkk) * p(+1, 0,-1) &
               +              Am0p(iii-1,jjj-1,kkk) * p( 0,+1,-1) &
               +              A00p(iii-1,jjj-1,kkk) * p(+1,+1,-1)
             ap(0,-1,0) = &
               &              Ammp(iii,jjj-1,kkk) * p(+1, 0,-1) &
               +              Am0p(iii,jjj-1,kkk) * p(+1,+1,-1)
             ap(-1,0,0) = &
               &              Ammp(iii-1,jjj,kkk) * p( 0,+1,-1) &
               +              A0mp(iii-1,jjj,kkk) * p(+1,+1,-1)
             ap(0,0,0) = &
               &              Ammp(iii,jjj,kkk) * p(+1,+1,-1)
             ap(-1,-1,1) = &
               &              Amm0(iii-1,jjj-1,kkk+1) * p( 0, 0,-1) &
               +              A0m0(iii-1,jjj-1,kkk+1) * p(+1, 0,-1) &
               +              Am00(iii-1,jjj-1,kkk+1) * p( 0,+1,-1) &
               +              A000(iii-1,jjj-1,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii-1,jjj-1,kkk+1) * p( 0, 0, 0) &
               +              A0mp(iii-1,jjj-1,kkk+1) * p(+1, 0, 0) &
               +              Am0p(iii-1,jjj-1,kkk+1) * p( 0,+1, 0) &
               +              A00p(iii-1,jjj-1,kkk+1) * p(+1,+1, 0)
             ap(0,-1,1) = &
               &              Amm0(iii,jjj-1,kkk+1) * p(+1, 0,-1) &
               +              Am00(iii,jjj-1,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii,jjj-1,kkk+1) * p(+1, 0, 0) &
               +              Am0p(iii,jjj-1,kkk+1) * p(+1,+1, 0)
             ap(-1,0,1) = &
               &              Amm0(iii-1,jjj,kkk+1) * p( 0,+1,-1) &
               +              A0m0(iii-1,jjj,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii-1,jjj,kkk+1) * p( 0,+1, 0) &
               +              A0mp(iii-1,jjj,kkk+1) * p(+1,+1, 0)
             ap(0,0,1) = &
               &              Amm0(iii,jjj,kkk+1) * p(+1,+1,-1) &
               +              Ammp(iii,jjj,kkk+1) * p(+1,+1, 0)
             cs4 = 0.125d0 * &
               ( restrict_from_mm0_to(iii,jjj,kkk) * ap(-1,-1, 0) &
               + restrict_from_0m0_to(iii,jjj,kkk) * ap( 0,-1, 0) &
               + restrict_from_m00_to(iii,jjj,kkk) * ap(-1, 0, 0) &
               + restrict_from_000_to(iii,jjj,kkk) * ap( 0, 0, 0) &
               + restrict_from_mmp_to(iii,jjj,kkk) * ap(-1,-1,+1) &
               + restrict_from_0mp_to(iii,jjj,kkk) * ap( 0,-1,+1) &
               + restrict_from_m0p_to(iii,jjj,kkk) * ap(-1, 0,+1) &
               + restrict_from_00p_to(iii,jjj,kkk) * ap( 0, 0,+1))

             csten(i,j,k,ist_ppp) = 0.25d0*(cs1+cs2+cs3+cs4)

          end do
       end do
    end do

  contains

    elemental function interp_from_mmm_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p
      p = 1.d0
      p = p      + abs(fsten(i-1,j  ,k  ,ist_p00)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k  ,ist_ppp)) + eps)
      p = p      + abs(fsten(i  ,j-1,k  ,ist_0p0)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k  ,ist_ppp)) + eps)
      p = p      + abs(fsten(i  ,j  ,k-1,ist_00p)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k-1,ist_ppp)) + eps)
      p = p      + abs(fsten(i-1,j-1,k  ,ist_pp0)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j-1,k  ,ist_ppp)) + eps)
      p = p      + abs(fsten(i-1,j  ,k-1,ist_p0p)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k-1,ist_ppp)) + eps)
      p = p      + abs(fsten(i  ,j-1,k-1,ist_0pp)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k-1,ist_ppp)) + eps)
      p = p * abs(fsten(i-1,j-1,k-1,ist_ppp)) * fsten(i,j,k,ist_inv)
    end function interp_from_mmm_to

    elemental function interp_from_pmm_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p
      p = 1.d0
      p = p + abs(fsten(i,j,k,ist_p00)) / &
           &     ( abs(fsten(i,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i,j  ,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i,j-1,k,ist_0p0)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i,j,k-1,ist_00p)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k-1,ist_ppp)) + eps)
      p = p + abs(fsten(i,j-1,k,ist_pp0)) / &
           &     ( abs(fsten(i,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i,j-1,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i,j,k-1,ist_p0p)) / &
           &     ( abs(fsten(i,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i,j  ,k-1,ist_ppp)) + eps)
      p = p + abs(fsten(i,j-1,k-1,ist_0pp)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k-1,ist_ppp)) + eps)
      p = p * abs(fsten(i  ,j-1,k-1,ist_ppp)) * fsten(i,j,k,ist_inv)
    end function interp_from_pmm_to

    elemental function interp_from_mpm_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p
      p = 1.d0
      p = p + abs(fsten(i-1,j,k,ist_p00)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i,j,k,ist_0p0)) / &
           &     ( abs(fsten(i-1,j,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i,j,k-1,ist_00p)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k-1,ist_ppp)) + eps)
      p = p + abs(fsten(i-1,j,k,ist_pp0)) / &
           &     ( abs(fsten(i-1,j,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i-1,j,k-1,ist_p0p)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k-1,ist_ppp)) + eps)
      p = p + abs(fsten(i,j,k-1,ist_0pp)) / &
           &     ( abs(fsten(i-1,j,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j,k-1,ist_ppp)) + eps)
      p = p * abs(fsten(i-1,j  ,k-1,ist_ppp)) * fsten(i,j,k,ist_inv)
    end function interp_from_mpm_to

    elemental function interp_from_ppm_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p
      p = 1.d0
      p = p + abs(fsten(i,j,k,ist_p00)) / &
           &     ( abs(fsten(i,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i,j  ,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i,j,k,ist_0p0)) / &
           &     ( abs(fsten(i-1,j,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i,j,k-1,ist_00p)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k-1,ist_ppp)) + eps)
      p = p + abs(fsten(i,j,k,ist_pp0)) / &
           &     ( abs(fsten(i,j,k-1,ist_ppp)) &
           &     + abs(fsten(i,j,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i,j,k-1,ist_p0p)) / &
           &     ( abs(fsten(i,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i,j  ,k-1,ist_ppp)) + eps)
      p = p + abs(fsten(i,j,k-1,ist_0pp)) / &
           &     ( abs(fsten(i-1,j,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j,k-1,ist_ppp)) + eps)
      p = p * abs(fsten(i  ,j  ,k-1,ist_ppp)) * fsten(i,j,k,ist_inv)
    end function interp_from_ppm_to

    elemental function interp_from_mmp_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p
      p = 1.d0
      p = p + abs(fsten(i-1,j,k,ist_p00)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i,j-1,k,ist_0p0)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i,j,k,ist_00p)) / &
           &     ( abs(fsten(i-1,j-1,k,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k,ist_ppp)) + eps)
      p = p + abs(fsten(i-1,j-1,k,ist_pp0)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j-1,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i-1,j,k,ist_p0p)) / &
           &     ( abs(fsten(i-1,j-1,k,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k,ist_ppp)) + eps)
      p = p + abs(fsten(i,j-1,k,ist_0pp)) / &
           &     ( abs(fsten(i-1,j-1,k,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k,ist_ppp)) + eps)
      p = p * abs(fsten(i-1,j-1,k  ,ist_ppp)) * fsten(i,j,k,ist_inv)
    end function interp_from_mmp_to

    elemental function interp_from_pmp_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p
      p = 1.d0
      p = p + abs(fsten(i,j,k,ist_p00)) / &
           &     ( abs(fsten(i,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i,j  ,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i,j-1,k,ist_0p0)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i,j,k,ist_00p)) / &
           &     ( abs(fsten(i-1,j-1,k,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k,ist_ppp)) + eps)
      p = p + abs(fsten(i,j-1,k,ist_pp0)) / &
           &     ( abs(fsten(i,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i,j-1,k  ,ist_ppp)) + eps)
      p = p + abs(fsten(i,j,k,ist_p0p)) / &
           &     ( abs(fsten(i,j-1,k,ist_ppp)) &
           &     + abs(fsten(i,j  ,k,ist_ppp)) + eps)
      p = p + abs(fsten(i,j-1,k,ist_0pp)) / &
           &     ( abs(fsten(i-1,j-1,k,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k,ist_ppp)) + eps)
      p = p * abs(fsten(i  ,j-1,k  ,ist_ppp)) * fsten(i,j,k,ist_inv)
    end function interp_from_pmp_to

    elemental function interp_from_mpp_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p
      p = 1.d0
      p = p      + abs(fsten(i-1,j  ,k  ,ist_p00)) / &
           &     ( abs(fsten(i-1,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k  ,ist_ppp)) + eps)
      p = p      + abs(fsten(i  ,j  ,k  ,ist_0p0)) / &
           &     ( abs(fsten(i-1,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k  ,ist_ppp)) + eps)
      p = p      + abs(fsten(i  ,j  ,k  ,ist_00p)) / &
           &     ( abs(fsten(i-1,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k  ,ist_ppp)) + eps)
      p = p      + abs(fsten(i-1,j  ,k  ,ist_pp0)) / &
           &     ( abs(fsten(i-1,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k  ,ist_ppp)) + eps)
      p = p      + abs(fsten(i-1,j  ,k  ,ist_p0p)) / &
           &     ( abs(fsten(i-1,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k  ,ist_ppp)) + eps)
      p = p      + abs(fsten(i  ,j  ,k  ,ist_0pp)) / &
           &     ( abs(fsten(i-1,j  ,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k  ,ist_ppp)) + eps)
      p = p * abs(fsten(i-1,j  ,k  ,ist_ppp)) * fsten(i,j,k,ist_inv)
    end function interp_from_mpp_to

    elemental function interp_from_ppp_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p
      p = 1.d0
      p = p      + abs(fsten(i  ,j  ,k  ,ist_p00)) / &
           &     ( abs(fsten(i  ,j-1,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k  ,ist_ppp)) + eps)
      p = p      + abs(fsten(i  ,j  ,k  ,ist_0p0)) / &
           &     ( abs(fsten(i-1,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k  ,ist_ppp)) + eps)
      p = p      + abs(fsten(i  ,j  ,k  ,ist_00p)) / &
           &     ( abs(fsten(i-1,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i-1,j  ,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k  ,ist_ppp)) + eps)
      p = p      + abs(fsten(i  ,j  ,k  ,ist_pp0)) / &
           &     ( abs(fsten(i  ,j  ,k-1,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k  ,ist_ppp)) + eps)
      p = p      + abs(fsten(i  ,j  ,k  ,ist_p0p)) / &
           &     ( abs(fsten(i  ,j-1,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k  ,ist_ppp)) + eps)
      p = p      + abs(fsten(i  ,j  ,k  ,ist_0pp)) / &
           &     ( abs(fsten(i-1,j  ,k  ,ist_ppp)) &
           &     + abs(fsten(i  ,j  ,k  ,ist_ppp)) + eps)
      p = p * abs(fsten(i,j,k,ist_ppp)) * fsten(i,j,k,ist_inv)
    end function interp_from_ppp_to

    elemental function interp_from_0mm_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1m, w2m, w1p, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(i,j-1,k,ist_0p0)) / (abs(fsten(i,j-1,k-1,ist_0pp)) &
           &                              +abs(fsten(i,j-1,k  ,ist_0pp)) + eps)
      w1p = abs(fsten(i,j  ,k,ist_0p0)) / (abs(fsten(i,j  ,k-1,ist_0pp)) &
           &                              +abs(fsten(i,j  ,k  ,ist_0pp)) + eps)
      w2m = abs(fsten(i,j,k-1,ist_00p)) / (abs(fsten(i,j-1,k-1,ist_0pp)) &
           &                              +abs(fsten(i,j  ,k-1,ist_0pp)) + eps)
      w2p = abs(fsten(i,j,k  ,ist_00p)) / (abs(fsten(i,j-1,k  ,ist_0pp)) &
           &                              +abs(fsten(i,j  ,k  ,ist_0pp)) + eps)
      wmm = abs(fsten(i,j-1,k-1,ist_0pp)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(i,j  ,k-1,ist_0pp)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(i,j-1,k  ,ist_0pp)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(i,j  ,k  ,ist_0pp)) * (1.d0 + w1p + w2p)
      p = wmm / (wmm+wpm+wmp+wpp+eps)
    end function interp_from_0mm_to

    elemental function interp_from_0mp_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1m, w2m, w1p, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(i,j-1,k,ist_0p0)) / (abs(fsten(i,j-1,k-1,ist_0pp)) &
           &                              +abs(fsten(i,j-1,k  ,ist_0pp)) + eps)
      w1p = abs(fsten(i,j  ,k,ist_0p0)) / (abs(fsten(i,j  ,k-1,ist_0pp)) &
           &                              +abs(fsten(i,j  ,k  ,ist_0pp)) + eps)
      w2m = abs(fsten(i,j,k-1,ist_00p)) / (abs(fsten(i,j-1,k-1,ist_0pp)) &
           &                              +abs(fsten(i,j  ,k-1,ist_0pp)) + eps)
      w2p = abs(fsten(i,j,k  ,ist_00p)) / (abs(fsten(i,j-1,k  ,ist_0pp)) &
           &                              +abs(fsten(i,j  ,k  ,ist_0pp)) + eps)
      wmm = abs(fsten(i,j-1,k-1,ist_0pp)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(i,j  ,k-1,ist_0pp)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(i,j-1,k  ,ist_0pp)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(i,j  ,k  ,ist_0pp)) * (1.d0 + w1p + w2p)
      p = wmp / (wmm+wpm+wmp+wpp+eps)
    end function interp_from_0mp_to

    elemental function interp_from_0pm_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1m, w2m, w1p, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(i,j-1,k,ist_0p0)) / (abs(fsten(i,j-1,k-1,ist_0pp)) &
           &                              +abs(fsten(i,j-1,k  ,ist_0pp)) + eps)
      w1p = abs(fsten(i,j  ,k,ist_0p0)) / (abs(fsten(i,j  ,k-1,ist_0pp)) &
           &                              +abs(fsten(i,j  ,k  ,ist_0pp)) + eps)
      w2m = abs(fsten(i,j,k-1,ist_00p)) / (abs(fsten(i,j-1,k-1,ist_0pp)) &
           &                              +abs(fsten(i,j  ,k-1,ist_0pp)) + eps)
      w2p = abs(fsten(i,j,k  ,ist_00p)) / (abs(fsten(i,j-1,k  ,ist_0pp)) &
           &                              +abs(fsten(i,j  ,k  ,ist_0pp)) + eps)
      wmm = abs(fsten(i,j-1,k-1,ist_0pp)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(i,j  ,k-1,ist_0pp)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(i,j-1,k  ,ist_0pp)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(i,j  ,k  ,ist_0pp)) * (1.d0 + w1p + w2p)
      p = wpm / (wmm+wpm+wmp+wpp+eps)
    end function interp_from_0pm_to

    elemental function interp_from_0pp_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1m, w2m, w1p, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(i,j-1,k,ist_0p0)) / (abs(fsten(i,j-1,k-1,ist_0pp)) &
           &                              +abs(fsten(i,j-1,k  ,ist_0pp)) + eps)
      w1p = abs(fsten(i,j  ,k,ist_0p0)) / (abs(fsten(i,j  ,k-1,ist_0pp)) &
           &                              +abs(fsten(i,j  ,k  ,ist_0pp)) + eps)
      w2m = abs(fsten(i,j,k-1,ist_00p)) / (abs(fsten(i,j-1,k-1,ist_0pp)) &
           &                              +abs(fsten(i,j  ,k-1,ist_0pp)) + eps)
      w2p = abs(fsten(i,j,k  ,ist_00p)) / (abs(fsten(i,j-1,k  ,ist_0pp)) &
           &                              +abs(fsten(i,j  ,k  ,ist_0pp)) + eps)
      wmm = abs(fsten(i,j-1,k-1,ist_0pp)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(i,j  ,k-1,ist_0pp)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(i,j-1,k  ,ist_0pp)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(i,j  ,k  ,ist_0pp)) * (1.d0 + w1p + w2p)
      p = wpp / (wmm+wpm+wmp+wpp+eps)
    end function interp_from_0pp_to

    elemental function interp_from_m0m_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1m, w2m, w1p, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(i-1,j,k,ist_p00)) / (abs(fsten(i-1,j,k-1,ist_p0p)) &
           &                              +abs(fsten(i-1,j,k  ,ist_p0p)) + eps)
      w1p = abs(fsten(i  ,j,k,ist_p00)) / (abs(fsten(i  ,j,k-1,ist_p0p)) &
           &                              +abs(fsten(i  ,j,k  ,ist_p0p)) + eps)
      w2m = abs(fsten(i,j,k-1,ist_00p)) / (abs(fsten(i-1,j,k-1,ist_p0p)) &
           &                              +abs(fsten(i  ,j,k-1,ist_p0p)) + eps)
      w2p = abs(fsten(i,j,k  ,ist_00p)) / (abs(fsten(i-1,j,k  ,ist_p0p)) &
           &                              +abs(fsten(i  ,j,k  ,ist_p0p)) + eps)
      wmm = abs(fsten(i-1,j,k-1,ist_p0p)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(i  ,j,k-1,ist_p0p)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(i-1,j,k  ,ist_p0p)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(i  ,j,k  ,ist_p0p)) * (1.d0 + w1p + w2p)
      p = wmm / (wmm+wpm+wmp+wpp+eps)
    end function interp_from_m0m_to

    elemental function interp_from_p0m_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1m, w2m, w1p, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(i-1,j,k,ist_p00)) / (abs(fsten(i-1,j,k-1,ist_p0p)) &
           &                              +abs(fsten(i-1,j,k  ,ist_p0p)) + eps)
      w1p = abs(fsten(i  ,j,k,ist_p00)) / (abs(fsten(i  ,j,k-1,ist_p0p)) &
           &                              +abs(fsten(i  ,j,k  ,ist_p0p)) + eps)
      w2m = abs(fsten(i,j,k-1,ist_00p)) / (abs(fsten(i-1,j,k-1,ist_p0p)) &
           &                              +abs(fsten(i  ,j,k-1,ist_p0p)) + eps)
      w2p = abs(fsten(i,j,k  ,ist_00p)) / (abs(fsten(i-1,j,k  ,ist_p0p)) &
           &                              +abs(fsten(i  ,j,k  ,ist_p0p)) + eps)
      wmm = abs(fsten(i-1,j,k-1,ist_p0p)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(i  ,j,k-1,ist_p0p)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(i-1,j,k  ,ist_p0p)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(i  ,j,k  ,ist_p0p)) * (1.d0 + w1p + w2p)
      p = wpm / (wmm+wpm+wmp+wpp+eps)
    end function interp_from_p0m_to

    elemental function interp_from_m0p_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1m, w2m, w1p, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(i-1,j,k,ist_p00)) / (abs(fsten(i-1,j,k-1,ist_p0p)) &
           &                              +abs(fsten(i-1,j,k  ,ist_p0p)) + eps)
      w1p = abs(fsten(i  ,j,k,ist_p00)) / (abs(fsten(i  ,j,k-1,ist_p0p)) &
           &                              +abs(fsten(i  ,j,k  ,ist_p0p)) + eps)
      w2m = abs(fsten(i,j,k-1,ist_00p)) / (abs(fsten(i-1,j,k-1,ist_p0p)) &
           &                              +abs(fsten(i  ,j,k-1,ist_p0p)) + eps)
      w2p = abs(fsten(i,j,k  ,ist_00p)) / (abs(fsten(i-1,j,k  ,ist_p0p)) &
           &                              +abs(fsten(i  ,j,k  ,ist_p0p)) + eps)
      wmm = abs(fsten(i-1,j,k-1,ist_p0p)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(i  ,j,k-1,ist_p0p)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(i-1,j,k  ,ist_p0p)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(i  ,j,k  ,ist_p0p)) * (1.d0 + w1p + w2p)
      p = wmp / (wmm+wpm+wmp+wpp+eps)
    end function interp_from_m0p_to

    elemental function interp_from_p0p_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1m, w2m, w1p, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(i-1,j,k,ist_p00)) / (abs(fsten(i-1,j,k-1,ist_p0p)) &
           &                              +abs(fsten(i-1,j,k  ,ist_p0p)) + eps)
      w1p = abs(fsten(i  ,j,k,ist_p00)) / (abs(fsten(i  ,j,k-1,ist_p0p)) &
           &                              +abs(fsten(i  ,j,k  ,ist_p0p)) + eps)
      w2m = abs(fsten(i,j,k-1,ist_00p)) / (abs(fsten(i-1,j,k-1,ist_p0p)) &
           &                              +abs(fsten(i  ,j,k-1,ist_p0p)) + eps)
      w2p = abs(fsten(i,j,k  ,ist_00p)) / (abs(fsten(i-1,j,k  ,ist_p0p)) &
           &                              +abs(fsten(i  ,j,k  ,ist_p0p)) + eps)
      wmm = abs(fsten(i-1,j,k-1,ist_p0p)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(i  ,j,k-1,ist_p0p)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(i-1,j,k  ,ist_p0p)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(i  ,j,k  ,ist_p0p)) * (1.d0 + w1p + w2p)
      p = wpp / (wmm+wpm+wmp+wpp+eps)
    end function interp_from_p0p_to

    elemental function interp_from_mm0_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1m, w2m, w1p, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(i-1,j,k,ist_p00)) / (abs(fsten(i-1,j-1,k,ist_pp0)) &
           &                              +abs(fsten(i-1,j  ,k,ist_pp0)) + eps)
      w1p = abs(fsten(i  ,j,k,ist_p00)) / (abs(fsten(i  ,j-1,k,ist_pp0)) &
           &                              +abs(fsten(i  ,j  ,k,ist_pp0)) + eps)
      w2m = abs(fsten(i,j-1,k,ist_0p0)) / (abs(fsten(i-1,j-1,k,ist_pp0)) &
           &                              +abs(fsten(i  ,j-1,k,ist_pp0)) + eps)
      w2p = abs(fsten(i,j  ,k,ist_0p0)) / (abs(fsten(i-1,j  ,k,ist_pp0)) &
           &                              +abs(fsten(i  ,j  ,k,ist_pp0)) + eps)
      wmm = abs(fsten(i-1,j-1,k,ist_pp0)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(i  ,j-1,k,ist_pp0)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(i-1,j  ,k,ist_pp0)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(i  ,j  ,k,ist_pp0)) * (1.d0 + w1p + w2p)
      p = wmm / (wmm+wpm+wmp+wpp+eps)
    end function interp_from_mm0_to

    elemental function interp_from_mp0_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1m, w2m, w1p, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(i-1,j,k,ist_p00)) / (abs(fsten(i-1,j-1,k,ist_pp0)) &
           &                              +abs(fsten(i-1,j  ,k,ist_pp0)) + eps)
      w1p = abs(fsten(i  ,j,k,ist_p00)) / (abs(fsten(i  ,j-1,k,ist_pp0)) &
           &                              +abs(fsten(i  ,j  ,k,ist_pp0)) + eps)
      w2m = abs(fsten(i,j-1,k,ist_0p0)) / (abs(fsten(i-1,j-1,k,ist_pp0)) &
           &                              +abs(fsten(i  ,j-1,k,ist_pp0)) + eps)
      w2p = abs(fsten(i,j  ,k,ist_0p0)) / (abs(fsten(i-1,j  ,k,ist_pp0)) &
           &                              +abs(fsten(i  ,j  ,k,ist_pp0)) + eps)
      wmm = abs(fsten(i-1,j-1,k,ist_pp0)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(i  ,j-1,k,ist_pp0)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(i-1,j  ,k,ist_pp0)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(i  ,j  ,k,ist_pp0)) * (1.d0 + w1p + w2p)
      p = wmp / (wmm+wpm+wmp+wpp+eps)
    end function interp_from_mp0_to

    elemental function interp_from_pm0_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1m, w2m, w1p, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(i-1,j,k,ist_p00)) / (abs(fsten(i-1,j-1,k,ist_pp0)) &
           &                              +abs(fsten(i-1,j  ,k,ist_pp0)) + eps)
      w1p = abs(fsten(i  ,j,k,ist_p00)) / (abs(fsten(i  ,j-1,k,ist_pp0)) &
           &                              +abs(fsten(i  ,j  ,k,ist_pp0)) + eps)
      w2m = abs(fsten(i,j-1,k,ist_0p0)) / (abs(fsten(i-1,j-1,k,ist_pp0)) &
           &                              +abs(fsten(i  ,j-1,k,ist_pp0)) + eps)
      w2p = abs(fsten(i,j  ,k,ist_0p0)) / (abs(fsten(i-1,j  ,k,ist_pp0)) &
           &                              +abs(fsten(i  ,j  ,k,ist_pp0)) + eps)
      wmm = abs(fsten(i-1,j-1,k,ist_pp0)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(i  ,j-1,k,ist_pp0)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(i-1,j  ,k,ist_pp0)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(i  ,j  ,k,ist_pp0)) * (1.d0 + w1p + w2p)
      p = wpm / (wmm+wpm+wmp+wpp+eps)
    end function interp_from_pm0_to

    elemental function interp_from_pp0_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1m, w2m, w1p, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(i-1,j,k,ist_p00)) / (abs(fsten(i-1,j-1,k,ist_pp0)) &
           &                              +abs(fsten(i-1,j  ,k,ist_pp0)) + eps)
      w1p = abs(fsten(i  ,j,k,ist_p00)) / (abs(fsten(i  ,j-1,k,ist_pp0)) &
           &                              +abs(fsten(i  ,j  ,k,ist_pp0)) + eps)
      w2m = abs(fsten(i,j-1,k,ist_0p0)) / (abs(fsten(i-1,j-1,k,ist_pp0)) &
           &                              +abs(fsten(i  ,j-1,k,ist_pp0)) + eps)
      w2p = abs(fsten(i,j  ,k,ist_0p0)) / (abs(fsten(i-1,j  ,k,ist_pp0)) &
           &                              +abs(fsten(i  ,j  ,k,ist_pp0)) + eps)
      wmm = abs(fsten(i-1,j-1,k,ist_pp0)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(i  ,j-1,k,ist_pp0)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(i-1,j  ,k,ist_pp0)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(i  ,j  ,k,ist_pp0)) * (1.d0 + w1p + w2p)
      p = wpp / (wmm+wpm+wmp+wpp+eps)
    end function interp_from_pp0_to

    elemental function interp_from_00m_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1, w2
      w1 = abs(fsten(i,j,k-1,ist_00p))
      w2 = abs(fsten(i,j,k  ,ist_00p))
      if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
         p = 0.5d0
      else
         p = w1 / (w1+w2)
      end if
    end function interp_from_00m_to

    elemental function interp_from_00p_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1, w2
      w1 = abs(fsten(i,j,k-1,ist_00p))
      w2 = abs(fsten(i,j,k  ,ist_00p))
      if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
         p = 0.5d0
      else
         p = w2 / (w1+w2)
      end if
    end function interp_from_00p_to

    elemental function interp_from_0m0_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1, w2
      w1 = abs(fsten(i,j-1,k,ist_0p0))
      w2 = abs(fsten(i,j  ,k,ist_0p0))
      if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
         p = 0.5d0
      else
         p = w1 / (w1+w2)
      end if
    end function interp_from_0m0_to

    elemental function interp_from_0p0_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1, w2
      w1 = abs(fsten(i,j-1,k,ist_0p0))
      w2 = abs(fsten(i,j  ,k,ist_0p0))
      if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
         p = 0.5d0
      else
         p = w2 / (w1+w2)
      end if
    end function interp_from_0p0_to

    elemental function interp_from_m00_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1, w2
      w1 = abs(fsten(i-1,j,k,ist_p00))
      w2 = abs(fsten(i  ,j,k,ist_p00))
      if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
         p = 0.5d0
      else
         p = w1 / (w1+w2)
      end if
    end function interp_from_m00_to

    elemental function interp_from_p00_to (i,j,k) result(p)
      integer, intent(in) :: i,j,k
      real(amrex_real) :: p, w1, w2
      w1 = abs(fsten(i-1,j,k,ist_p00))
      w2 = abs(fsten(i  ,j,k,ist_p00))
      if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
         p = 0.5d0
      else
         p = w2 / (w1+w2)
      end if
    end function interp_from_p00_to

    elemental real(amrex_real) function Ammm (i,j,k)
      integer, intent(in) :: i,j,k
      Ammm = fsten(i-1,j-1,k-1,ist_ppp)
    end function Ammm

    elemental real(amrex_real) function A0mm (i,j,k)
      integer, intent(in) :: i,j,k
      A0mm = fsten(i  ,j-1,k-1,ist_0pp)
    end function A0mm

    elemental real(amrex_real) function Apmm (i,j,k)
      integer, intent(in) :: i,j,k
      Apmm = fsten(i  ,j-1,k-1,ist_ppp)
    end function Apmm

    elemental real(amrex_real) function Am0m (i,j,k)
      integer, intent(in) :: i,j,k
      Am0m = fsten(i-1,j  ,k-1,ist_p0p)
    end function Am0m

    elemental real(amrex_real) function A00m (i,j,k)
      integer, intent(in) :: i,j,k
      A00m = fsten(i  ,j  ,k-1,ist_00p)
    end function A00m

    elemental real(amrex_real) function Ap0m (i,j,k)
      integer, intent(in) :: i,j,k
      Ap0m = fsten(i  ,j  ,k-1,ist_p0p)
    end function Ap0m

    elemental real(amrex_real) function Ampm (i,j,k)
      integer, intent(in) :: i,j,k
      Ampm = fsten(i-1,j  ,k-1,ist_ppp)
    end function Ampm

    elemental real(amrex_real) function A0pm (i,j,k)
      integer, intent(in) :: i,j,k
      A0pm = fsten(i  ,j  ,k-1,ist_0pp)
    end function A0pm

    elemental real(amrex_real) function Appm (i,j,k)
      integer, intent(in) :: i,j,k
      Appm = fsten(i  ,j  ,k-1,ist_ppp)
    end function Appm

    elemental real(amrex_real) function Amm0 (i,j,k)
      integer, intent(in) :: i,j,k
      Amm0 = fsten(i-1,j-1,k  ,ist_pp0)
    end function Amm0

    elemental real(amrex_real) function A0m0 (i,j,k)
      integer, intent(in) :: i,j,k
      A0m0 = fsten(i  ,j-1,k  ,ist_0p0)
    end function A0m0

    elemental real(amrex_real) function Apm0 (i,j,k)
      integer, intent(in) :: i,j,k
      Apm0 = fsten(i  ,j-1,k  ,ist_pp0)
    end function Apm0

    elemental real(amrex_real) function Am00 (i,j,k)
      integer, intent(in) :: i,j,k
      Am00 = fsten(i-1,j  ,k  ,ist_p00)
    end function Am00

    elemental real(amrex_real) function A000 (i,j,k)
      integer, intent(in) :: i,j,k
      A000 = fsten(i  ,j  ,k  ,ist_000)
    end function A000

    elemental real(amrex_real) function Ap00 (i,j,k)
      integer, intent(in) :: i,j,k
      Ap00 = fsten(i  ,j  ,k  ,ist_p00)
    end function Ap00

    elemental real(amrex_real) function Amp0 (i,j,k)
      integer, intent(in) :: i,j,k
      Amp0 = fsten(i-1,j  ,k  ,ist_pp0)
    end function Amp0

    elemental real(amrex_real) function A0p0 (i,j,k)
      integer, intent(in) :: i,j,k
      A0p0 = fsten(i  ,j  ,k  ,ist_0p0)
    end function A0p0

    elemental real(amrex_real) function App0 (i,j,k)
      integer, intent(in) :: i,j,k
      App0 = fsten(i  ,j  ,k  ,ist_pp0)
    end function App0

    elemental real(amrex_real) function Ammp (i,j,k)
      integer, intent(in) :: i,j,k
      Ammp = fsten(i-1,j-1,k  ,ist_ppp)
    end function Ammp

    elemental real(amrex_real) function A0mp (i,j,k)
      integer, intent(in) :: i,j,k
      A0mp = fsten(i  ,j-1,k  ,ist_0pp)
    end function A0mp

    elemental real(amrex_real) function Apmp (i,j,k)
      integer, intent(in) :: i,j,k
      Apmp = fsten(i  ,j-1,k  ,ist_ppp)
    end function Apmp

    elemental real(amrex_real) function Am0p (i,j,k)
      integer, intent(in) :: i,j,k
      Am0p = fsten(i-1,j  ,k  ,ist_p0p)
    end function Am0p

    elemental real(amrex_real) function A00p (i,j,k)
      integer, intent(in) :: i,j,k
      A00p = fsten(i  ,j  ,k  ,ist_00p)
    end function A00p

    elemental real(amrex_real) function Ap0p (i,j,k)
      integer, intent(in) :: i,j,k
      Ap0p = fsten(i  ,j  ,k  ,ist_p0p)
    end function Ap0p

    elemental real(amrex_real) function Ampp (i,j,k)
      integer, intent(in) :: i,j,k
      Ampp = fsten(i-1,j  ,k  ,ist_ppp)
    end function Ampp

    elemental real(amrex_real) function A0pp (i,j,k)
      integer, intent(in) :: i,j,k
      A0pp = fsten(i  ,j  ,k  ,ist_0pp)
    end function A0pp

    elemental real(amrex_real) function Appp (i,j,k)
      integer, intent(in) :: i,j,k
      Appp = fsten(i  ,j  ,k  ,ist_ppp)
    end function Appp

    elemental function restrict_from_mmm_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r
      r = 1.d0
      r = r      + abs(fsten(ii-1,jj-1,kk-1,ist_p00)) / &
           &     ( abs(fsten(ii-1,jj-2,kk-2,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk-2,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-2,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj-1,kk-1,ist_0p0)) / &
           &     ( abs(fsten(ii-2,jj-1,kk-2,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk-2,ist_ppp)) &
           &     + abs(fsten(ii-2,jj-1,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj-1,kk-1,ist_00p)) / &
           &     ( abs(fsten(ii-2,jj-2,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-2,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-2,jj-1,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj-1,kk-1,ist_pp0)) / &
           &     ( abs(fsten(ii-1,jj-1,kk-2,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj-1,kk-1,ist_p0p)) / &
           &     ( abs(fsten(ii-1,jj-2,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj-1,kk-1,ist_0pp)) / &
           &     ( abs(fsten(ii-2,jj-1,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk-1,ist_ppp)) + eps)
      r = r * abs(fsten(ii-1,jj-1,kk-1,ist_ppp)) * fsten(ii-1,jj-1,kk-1,ist_inv)
    end function restrict_from_mmm_to

    elemental function restrict_from_0mm_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1m, w1p, w2m, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(ii,jj-2,kk-1,ist_0p0)) / (abs(fsten(ii,jj-2,kk-2,ist_0pp)) &
           &                                   +abs(fsten(ii,jj-2,kk-1,ist_0pp)) + eps)
      w1p = abs(fsten(ii,jj-1,kk-1,ist_0p0)) / (abs(fsten(ii,jj-1,kk-2,ist_0pp)) &
           &                                   +abs(fsten(ii,jj-1,kk-1,ist_0pp)) + eps)
      w2m = abs(fsten(ii,jj-1,kk-2,ist_00p)) / (abs(fsten(ii,jj-2,kk-2,ist_0pp)) &
           &                                   +abs(fsten(ii,jj-1,kk-2,ist_0pp)) + eps)
      w2p = abs(fsten(ii,jj-1,kk-1,ist_00p)) / (abs(fsten(ii,jj-2,kk-1,ist_0pp)) &
           &                                   +abs(fsten(ii,jj-1,kk-1,ist_0pp)) + eps)
      wmm = abs(fsten(ii,jj-2,kk-2,ist_0pp)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(ii,jj-1,kk-2,ist_0pp)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(ii,jj-2,kk-1,ist_0pp)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(ii,jj-1,kk-1,ist_0pp)) * (1.d0 + w1p + w2p)
      r = wpp / (wmm+wpm+wmp+wpp+eps)
    end function restrict_from_0mm_to

    elemental function restrict_from_pmm_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r
      r = 1.d0
      r = r      + abs(fsten(ii  ,jj-1,kk-1,ist_p00)) / &
           &     ( abs(fsten(ii  ,jj-2,kk-2,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-1,kk-2,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-2,kk-1,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii+1,jj-1,kk-1,ist_0p0)) / &
           &     ( abs(fsten(ii  ,jj-1,kk-2,ist_ppp)) &
           &     + abs(fsten(ii+1,jj-1,kk-2,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-1,kk-1,ist_ppp)) &
           &     + abs(fsten(ii+1,jj-1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii+1,jj-1,kk-1,ist_00p)) / &
           &     ( abs(fsten(ii  ,jj-2,kk-1,ist_ppp)) &
           &     + abs(fsten(ii+1,jj-2,kk-1,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-1,kk-1,ist_ppp)) &
           &     + abs(fsten(ii+1,jj-1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii  ,jj-1,kk-1,ist_pp0)) / &
           &     ( abs(fsten(ii  ,jj-1,kk-2,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii  ,jj-1,kk-1,ist_p0p)) / &
           &     ( abs(fsten(ii  ,jj-2,kk-1,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii+1,jj-1,kk-1,ist_0pp)) / &
           &     ( abs(fsten(ii  ,jj-1,kk-1,ist_ppp)) &
           &     + abs(fsten(ii+1,jj-1,kk-1,ist_ppp)) + eps)
      r = r * abs(fsten(ii  ,jj-1,kk-1,ist_ppp)) * fsten(ii+1,jj-1,kk-1,ist_inv)
    end function restrict_from_pmm_to

    elemental function restrict_from_m0m_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1m, w1p, w2m, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(ii-2,jj,kk-1,ist_p00)) / (abs(fsten(ii-2,jj,kk-2,ist_p0p)) &
           &                                   +abs(fsten(ii-2,jj,kk-1,ist_p0p)) + eps)
      w1p = abs(fsten(ii-1,jj,kk-1,ist_p00)) / (abs(fsten(ii-1,jj,kk-2,ist_p0p)) &
           &                                   +abs(fsten(ii-1,jj,kk-1,ist_p0p)) + eps)
      w2m = abs(fsten(ii-1,jj,kk-2,ist_00p)) / (abs(fsten(ii-2,jj,kk-2,ist_p0p)) &
           &                                   +abs(fsten(ii-1,jj,kk-2,ist_p0p)) + eps)
      w2p = abs(fsten(ii-1,jj,kk-1,ist_00p)) / (abs(fsten(ii-2,jj,kk-1,ist_p0p)) &
           &                                   +abs(fsten(ii-1,jj,kk-1,ist_p0p)) + eps)
      wmm = abs(fsten(ii-2,jj,kk-2,ist_p0p)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(ii-1,jj,kk-2,ist_p0p)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(ii-2,jj,kk-1,ist_p0p)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(ii-1,jj,kk-1,ist_p0p)) * (1.d0 + w1p + w2p)
      r = wpp / (wmm+wpm+wmp+wpp+eps)
    end function restrict_from_m0m_to

    elemental function restrict_from_00m_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1, w2
      w1 = abs(fsten(ii,jj,kk-2,ist_00p))
      w2 = abs(fsten(ii,jj,kk-1,ist_00p))
      if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
         r = 0.5d0
      else
         r = w2 / (w1+w2)
      end if
    end function restrict_from_00m_to

    elemental function restrict_from_p0m_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1m, w1p, w2m, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(ii  ,jj,kk-1,ist_p00)) / (abs(fsten(ii  ,jj,kk-2,ist_p0p)) &
           &                                   +abs(fsten(ii  ,jj,kk-1,ist_p0p)) + eps)
      w1p = abs(fsten(ii+1,jj,kk-1,ist_p00)) / (abs(fsten(ii+1,jj,kk-2,ist_p0p)) &
           &                                   +abs(fsten(ii+1,jj,kk-1,ist_p0p)) + eps)
      w2m = abs(fsten(ii+1,jj,kk-2,ist_00p)) / (abs(fsten(ii  ,jj,kk-2,ist_p0p)) &
           &                                   +abs(fsten(ii+1,jj,kk-2,ist_p0p)) + eps)
      w2p = abs(fsten(ii+1,jj,kk-1,ist_00p)) / (abs(fsten(ii  ,jj,kk-1,ist_p0p)) &
           &                                   +abs(fsten(ii+1,jj,kk-1,ist_p0p)) + eps)
      wmm = abs(fsten(ii  ,jj,kk-2,ist_p0p)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(ii+1,jj,kk-2,ist_p0p)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(ii  ,jj,kk-1,ist_p0p)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(ii+1,jj,kk-1,ist_p0p)) * (1.d0 + w1p + w2p)
      r = wmp / (wmm+wpm+wmp+wpp+eps)
    end function restrict_from_p0m_to

    elemental function restrict_from_mpm_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r
      r = 1.d0
      r = r      + abs(fsten(ii-1,jj+1,kk-1,ist_p00)) / &
           &     ( abs(fsten(ii-1,jj  ,kk-2,ist_ppp)) &
           &     + abs(fsten(ii-1,jj+1,kk-2,ist_ppp)) &
           &     + abs(fsten(ii-1,jj  ,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj+1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj  ,kk-1,ist_0p0)) / &
           &     ( abs(fsten(ii-2,jj  ,kk-2,ist_ppp)) &
           &     + abs(fsten(ii-1,jj  ,kk-2,ist_ppp)) &
           &     + abs(fsten(ii-2,jj  ,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj  ,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj+1,kk-1,ist_00p)) / &
           &     ( abs(fsten(ii-2,jj  ,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj  ,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-2,jj+1,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj+1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj  ,kk-1,ist_pp0)) / &
           &     ( abs(fsten(ii-1,jj  ,kk-2,ist_ppp)) &
           &     + abs(fsten(ii-1,jj  ,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj+1,kk-1,ist_p0p)) / &
           &     ( abs(fsten(ii-1,jj  ,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj+1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj  ,kk-1,ist_0pp)) / &
           &     ( abs(fsten(ii-2,jj  ,kk-1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj  ,kk-1,ist_ppp)) + eps)
      r = r * abs(fsten(ii-1,jj  ,kk-1,ist_ppp)) * fsten(ii-1,jj+1,kk-1,ist_inv)
    end function restrict_from_mpm_to

    elemental function restrict_from_0pm_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1m, w1p, w2m, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(ii,jj  ,kk-1,ist_0p0)) / (abs(fsten(ii,jj  ,kk-2,ist_0pp)) &
           &                                   +abs(fsten(ii,jj  ,kk-1,ist_0pp)) + eps)
      w1p = abs(fsten(ii,jj+1,kk-1,ist_0p0)) / (abs(fsten(ii,jj+1,kk-2,ist_0pp)) &
           &                                   +abs(fsten(ii,jj+1,kk-1,ist_0pp)) + eps)
      w2m = abs(fsten(ii,jj+1,kk-2,ist_00p)) / (abs(fsten(ii,jj  ,kk-2,ist_0pp)) &
           &                                   +abs(fsten(ii,jj+1,kk-2,ist_0pp)) + eps)
      w2p = abs(fsten(ii,jj+1,kk-1,ist_00p)) / (abs(fsten(ii,jj  ,kk-1,ist_0pp)) &
           &                                   +abs(fsten(ii,jj+1,kk-1,ist_0pp)) + eps)
      wmm = abs(fsten(ii,jj  ,kk-2,ist_0pp)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(ii,jj+1,kk-2,ist_0pp)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(ii,jj  ,kk-1,ist_0pp)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(ii,jj+1,kk-1,ist_0pp)) * (1.d0 + w1p + w2p)
      r = wmp / (wmm+wpm+wmp+wpp+eps)
    end function restrict_from_0pm_to

    elemental function restrict_from_ppm_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r
      r = 1.d0
      r = r      + abs(fsten(ii  ,jj+1,kk-1,ist_p00)) / &
           &     ( abs(fsten(ii  ,jj  ,kk-2,ist_ppp)) &
           &     + abs(fsten(ii  ,jj+1,kk-2,ist_ppp)) &
           &     + abs(fsten(ii  ,jj  ,kk-1,ist_ppp)) &
           &     + abs(fsten(ii  ,jj+1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii+1,jj  ,kk-1,ist_0p0)) / &
           &     ( abs(fsten(ii  ,jj  ,kk-2,ist_ppp)) &
           &     + abs(fsten(ii+1,jj  ,kk-2,ist_ppp)) &
           &     + abs(fsten(ii  ,jj  ,kk-1,ist_ppp)) &
           &     + abs(fsten(ii+1,jj  ,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii+1,jj+1,kk-1,ist_00p)) / &
           &     ( abs(fsten(ii  ,jj  ,kk-1,ist_ppp)) &
           &     + abs(fsten(ii+1,jj  ,kk-1,ist_ppp)) &
           &     + abs(fsten(ii  ,jj+1,kk-1,ist_ppp)) &
           &     + abs(fsten(ii+1,jj+1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii  ,jj  ,kk-1,ist_pp0)) / &
           &     ( abs(fsten(ii  ,jj  ,kk-2,ist_ppp)) &
           &     + abs(fsten(ii  ,jj  ,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii  ,jj+1,kk-1,ist_p0p)) / &
           &     ( abs(fsten(ii  ,jj  ,kk-1,ist_ppp)) &
           &     + abs(fsten(ii  ,jj+1,kk-1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii+1,jj  ,kk-1,ist_0pp)) / &
           &     ( abs(fsten(ii  ,jj  ,kk-1,ist_ppp)) &
           &     + abs(fsten(ii+1,jj  ,kk-1,ist_ppp)) + eps)
      r = r * abs(fsten(ii  ,jj  ,kk-1,ist_ppp)) * fsten(ii+1,jj+1,kk-1,ist_inv)
    end function restrict_from_ppm_to

    elemental function restrict_from_mm0_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1m, w1p, w2m, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(ii-2,jj-1,kk,ist_p00)) / (abs(fsten(ii-2,jj-2,kk,ist_pp0)) &
           &                                   +abs(fsten(ii-2,jj-1,kk,ist_pp0)) + eps)
      w1p = abs(fsten(ii-1,jj-1,kk,ist_p00)) / (abs(fsten(ii-1,jj-2,kk,ist_pp0)) &
           &                                   +abs(fsten(ii-1,jj-1,kk,ist_pp0)) + eps)
      w2m = abs(fsten(ii-1,jj-2,kk,ist_0p0)) / (abs(fsten(ii-2,jj-2,kk,ist_pp0)) &
           &                                   +abs(fsten(ii-1,jj-2,kk,ist_pp0)) + eps)
      w2p = abs(fsten(ii-1,jj-1,kk,ist_0p0)) / (abs(fsten(ii-2,jj-1,kk,ist_pp0)) &
           &                                   +abs(fsten(ii-1,jj-1,kk,ist_pp0)) + eps)
      wmm = abs(fsten(ii-2,jj-2,kk,ist_pp0)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(ii-1,jj-2,kk,ist_pp0)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(ii-2,jj-1,kk,ist_pp0)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(ii-1,jj-1,kk,ist_pp0)) * (1.d0 + w1p + w2p)
      r = wpp / (wmm+wpm+wmp+wpp+eps)
    end function restrict_from_mm0_to

    elemental function restrict_from_0m0_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1, w2
      w1 = abs(fsten(ii,jj-2,kk,ist_0p0))
      w2 = abs(fsten(ii,jj-1,kk,ist_0p0))
      if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
         r = 0.5d0
      else
         r = w2 / (w1+w2)
      end if
    end function restrict_from_0m0_to

    elemental function restrict_from_pm0_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1m, w1p, w2m, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(ii  ,jj-1,kk,ist_p00)) / (abs(fsten(ii  ,jj-2,kk,ist_pp0)) &
           &                                   +abs(fsten(ii  ,jj-1,kk,ist_pp0)) + eps)
      w1p = abs(fsten(ii+1,jj-1,kk,ist_p00)) / (abs(fsten(ii+1,jj-2,kk,ist_pp0)) &
           &                                   +abs(fsten(ii+1,jj-1,kk,ist_pp0)) + eps)
      w2m = abs(fsten(ii+1,jj-2,kk,ist_0p0)) / (abs(fsten(ii  ,jj-2,kk,ist_pp0)) &
           &                                   +abs(fsten(ii+1,jj-2,kk,ist_pp0)) + eps)
      w2p = abs(fsten(ii+1,jj-1,kk,ist_0p0)) / (abs(fsten(ii  ,jj-1,kk,ist_pp0)) &
           &                                   +abs(fsten(ii+1,jj-1,kk,ist_pp0)) + eps)
      wmm = abs(fsten(ii  ,jj-2,kk,ist_pp0)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(ii+1,jj-2,kk,ist_pp0)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(ii  ,jj-1,kk,ist_pp0)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(ii+1,jj-1,kk,ist_pp0)) * (1.d0 + w1p + w2p)
      r = wmp / (wmm+wpm+wmp+wpp+eps)
    end function restrict_from_pm0_to

    elemental function restrict_from_m00_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1, w2
      w1 = abs(fsten(ii-2,jj,kk,ist_p00))
      w2 = abs(fsten(ii-1,jj,kk,ist_p00))
      if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
         r = 0.5d0
      else
         r = w2 / (w1+w2)
      end if
    end function restrict_from_m00_to

    elemental function restrict_from_000_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r
      r = 1.d0
    end function restrict_from_000_to

    elemental function restrict_from_p00_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1, w2
      w1 = abs(fsten(ii  ,jj,kk,ist_p00))
      w2 = abs(fsten(ii+1,jj,kk,ist_p00))
      if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
         r = 0.5d0
      else
         r = w1 / (w1+w2)
      end if
    end function restrict_from_p00_to

    elemental function restrict_from_mp0_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1m, w1p, w2m, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(ii-2,jj+1,kk,ist_p00)) / (abs(fsten(ii-2,jj  ,kk,ist_pp0)) &
           &                                   +abs(fsten(ii-2,jj+1,kk,ist_pp0)) + eps)
      w1p = abs(fsten(ii-1,jj+1,kk,ist_p00)) / (abs(fsten(ii-1,jj  ,kk,ist_pp0)) &
           &                                   +abs(fsten(ii-1,jj+1,kk,ist_pp0)) + eps)
      w2m = abs(fsten(ii-1,jj  ,kk,ist_0p0)) / (abs(fsten(ii-2,jj  ,kk,ist_pp0)) &
           &                                   +abs(fsten(ii-1,jj  ,kk,ist_pp0)) + eps)
      w2p = abs(fsten(ii-1,jj+1,kk,ist_0p0)) / (abs(fsten(ii-2,jj+1,kk,ist_pp0)) &
           &                                   +abs(fsten(ii-1,jj+1,kk,ist_pp0)) + eps)
      wmm = abs(fsten(ii-2,jj  ,kk,ist_pp0)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(ii-1,jj  ,kk,ist_pp0)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(ii-2,jj+1,kk,ist_pp0)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(ii-1,jj+1,kk,ist_pp0)) * (1.d0 + w1p + w2p)
      r = wpm / (wmm+wpm+wmp+wpp+eps)
    end function restrict_from_mp0_to

    elemental function restrict_from_0p0_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1, w2
      w1 = abs(fsten(ii,jj  ,kk,ist_0p0))
      w2 = abs(fsten(ii,jj+1,kk,ist_0p0))
      if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
         r = 0.5d0
      else
         r = w1 / (w1+w2)
      end if
    end function restrict_from_0p0_to

    elemental function restrict_from_pp0_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1m, w1p, w2m, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(ii  ,jj+1,kk,ist_p00)) / (abs(fsten(ii  ,jj  ,kk,ist_pp0)) &
           &                                   +abs(fsten(ii  ,jj+1,kk,ist_pp0)) + eps)
      w1p = abs(fsten(ii+1,jj+1,kk,ist_p00)) / (abs(fsten(ii+1,jj  ,kk,ist_pp0)) &
           &                                   +abs(fsten(ii+1,jj+1,kk,ist_pp0)) + eps)
      w2m = abs(fsten(ii+1,jj  ,kk,ist_0p0)) / (abs(fsten(ii  ,jj  ,kk,ist_pp0)) &
           &                                   +abs(fsten(ii+1,jj  ,kk,ist_pp0)) + eps)
      w2p = abs(fsten(ii+1,jj+1,kk,ist_0p0)) / (abs(fsten(ii  ,jj+1,kk,ist_pp0)) &
           &                                   +abs(fsten(ii+1,jj+1,kk,ist_pp0)) + eps)
      wmm = abs(fsten(ii  ,jj  ,kk,ist_pp0)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(ii+1,jj  ,kk,ist_pp0)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(ii  ,jj+1,kk,ist_pp0)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(ii+1,jj+1,kk,ist_pp0)) * (1.d0 + w1p + w2p)
      r = wmm / (wmm+wpm+wmp+wpp+eps)
    end function restrict_from_pp0_to

    elemental function restrict_from_mmp_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r
      r = 1.d0
      r = r      + abs(fsten(ii-1,jj-1,kk+1,ist_p00)) / &
           &     ( abs(fsten(ii-1,jj-2,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-2,kk+1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk+1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj-1,kk+1,ist_0p0)) / &
           &     ( abs(fsten(ii-2,jj-1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-2,jj-1,kk+1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk+1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj-1,kk  ,ist_00p)) / &
           &     ( abs(fsten(ii-2,jj-2,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-2,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-2,jj-1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk  ,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj-1,kk+1,ist_pp0)) / &
           &     ( abs(fsten(ii-1,jj-1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk+1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj-1,kk  ,ist_p0p)) / &
           &     ( abs(fsten(ii-1,jj-2,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk  ,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj-1,kk  ,ist_0pp)) / &
           &     ( abs(fsten(ii-2,jj-1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj-1,kk  ,ist_ppp)) + eps)
      r = r * abs(fsten(ii-1,jj-1,kk  ,ist_ppp)) * fsten(ii-1,jj-1,kk+1,ist_inv)
    end function restrict_from_mmp_to

    elemental function restrict_from_0mp_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1m, w1p, w2m, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(ii,jj-2,kk+1,ist_0p0)) / (abs(fsten(ii,jj-2,kk  ,ist_0pp)) &
           &                                   +abs(fsten(ii,jj-2,kk+1,ist_0pp)) + eps)
      w1p = abs(fsten(ii,jj-1,kk+1,ist_0p0)) / (abs(fsten(ii,jj-1,kk  ,ist_0pp)) &
           &                                   +abs(fsten(ii,jj-1,kk+1,ist_0pp)) + eps)
      w2m = abs(fsten(ii,jj-1,kk  ,ist_00p)) / (abs(fsten(ii,jj-2,kk  ,ist_0pp)) &
           &                                   +abs(fsten(ii,jj-1,kk  ,ist_0pp)) + eps)
      w2p = abs(fsten(ii,jj-1,kk+1,ist_00p)) / (abs(fsten(ii,jj-2,kk+1,ist_0pp)) &
           &                                   +abs(fsten(ii,jj-1,kk+1,ist_0pp)) + eps)
      wmm = abs(fsten(ii,jj-2,kk  ,ist_0pp)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(ii,jj-1,kk  ,ist_0pp)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(ii,jj-2,kk+1,ist_0pp)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(ii,jj-1,kk+1,ist_0pp)) * (1.d0 + w1p + w2p)
      r = wpm / (wmm+wpm+wmp+wpp+eps)
    end function restrict_from_0mp_to

    elemental function restrict_from_pmp_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r
      r = 1.d0
      r = r      + abs(fsten(ii  ,jj-1,kk+1,ist_p00)) / &
           &     ( abs(fsten(ii  ,jj-2,kk  ,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-2,kk+1,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-1,kk+1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii+1,jj-1,kk+1,ist_0p0)) / &
           &     ( abs(fsten(ii  ,jj-1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii+1,jj-1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-1,kk+1,ist_ppp)) &
           &     + abs(fsten(ii+1,jj-1,kk+1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii+1,jj-1,kk  ,ist_00p)) / &
           &     ( abs(fsten(ii  ,jj-2,kk  ,ist_ppp)) &
           &     + abs(fsten(ii+1,jj-2,kk  ,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii+1,jj-1,kk  ,ist_ppp)) + eps)
      r = r      + abs(fsten(ii  ,jj-1,kk+1,ist_pp0)) / &
           &     ( abs(fsten(ii  ,jj-1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-1,kk+1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii  ,jj-1,kk  ,ist_p0p)) / &
           &     ( abs(fsten(ii  ,jj-2,kk  ,ist_ppp)) &
           &     + abs(fsten(ii  ,jj-1,kk  ,ist_ppp)) + eps)
      r = r      + abs(fsten(ii+1,jj-1,kk  ,ist_0pp)) / &
           &     ( abs(fsten(ii  ,jj-1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii+1,jj-1,kk  ,ist_ppp)) + eps)
      r = r * abs(fsten(ii  ,jj-1,kk  ,ist_ppp)) * fsten(ii+1,jj-1,kk+1,ist_inv)
    end function restrict_from_pmp_to

    elemental function restrict_from_m0p_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1m, w1p, w2m, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(ii-2,jj,kk+1,ist_p00)) / (abs(fsten(ii-2,jj,kk  ,ist_p0p)) &
           &                                   +abs(fsten(ii-2,jj,kk+1,ist_p0p)) + eps)
      w1p = abs(fsten(ii-1,jj,kk+1,ist_p00)) / (abs(fsten(ii-1,jj,kk  ,ist_p0p)) &
           &                                   +abs(fsten(ii-1,jj,kk+1,ist_p0p)) + eps)
      w2m = abs(fsten(ii-1,jj,kk  ,ist_00p)) / (abs(fsten(ii-2,jj,kk  ,ist_p0p)) &
           &                                   +abs(fsten(ii-1,jj,kk  ,ist_p0p)) + eps)
      w2p = abs(fsten(ii-1,jj,kk+1,ist_00p)) / (abs(fsten(ii-2,jj,kk+1,ist_p0p)) &
           &                                   +abs(fsten(ii-1,jj,kk+1,ist_p0p)) + eps)
      wmm = abs(fsten(ii-2,jj,kk  ,ist_p0p)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(ii-1,jj,kk  ,ist_p0p)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(ii-2,jj,kk+1,ist_p0p)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(ii-1,jj,kk+1,ist_p0p)) * (1.d0 + w1p + w2p)
      r = wpm / (wmm+wpm+wmp+wpp+eps)
    end function restrict_from_m0p_to

    elemental function restrict_from_00p_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1, w2
      w1 = abs(fsten(ii,jj,kk  ,ist_00p))
      w2 = abs(fsten(ii,jj,kk+1,ist_00p))
      if (w1 .eq. 0.d0 .and. w2 .eq. 0.d0) then
         r = 0.5d0
      else
         r = w1 / (w1+w2)
      end if
    end function restrict_from_00p_to

    elemental function restrict_from_p0p_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1m, w1p, w2m, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(ii  ,jj,kk+1,ist_p00)) / (abs(fsten(ii  ,jj,kk  ,ist_p0p)) &
           &                                   +abs(fsten(ii  ,jj,kk+1,ist_p0p)) + eps)
      w1p = abs(fsten(ii+1,jj,kk+1,ist_p00)) / (abs(fsten(ii+1,jj,kk  ,ist_p0p)) &
           &                                   +abs(fsten(ii+1,jj,kk+1,ist_p0p)) + eps)
      w2m = abs(fsten(ii+1,jj,kk  ,ist_00p)) / (abs(fsten(ii  ,jj,kk  ,ist_p0p)) &
           &                                   +abs(fsten(ii+1,jj,kk  ,ist_p0p)) + eps)
      w2p = abs(fsten(ii+1,jj,kk+1,ist_00p)) / (abs(fsten(ii  ,jj,kk+1,ist_p0p)) &
           &                                   +abs(fsten(ii+1,jj,kk+1,ist_p0p)) + eps)
      wmm = abs(fsten(ii  ,jj,kk  ,ist_p0p)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(ii+1,jj,kk  ,ist_p0p)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(ii  ,jj,kk+1,ist_p0p)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(ii+1,jj,kk+1,ist_p0p)) * (1.d0 + w1p + w2p)
      r = wmm / (wmm+wpm+wmp+wpp+eps)
    end function restrict_from_p0p_to

    elemental function restrict_from_mpp_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r
      r = 1.d0
      r = r      + abs(fsten(ii-1,jj+1,kk+1,ist_p00)) / &
           &     ( abs(fsten(ii-1,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj+1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj  ,kk+1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj+1,kk+1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj  ,kk+1,ist_0p0)) / &
           &     ( abs(fsten(ii-2,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-2,jj  ,kk+1,ist_ppp)) &
           &     + abs(fsten(ii-1,jj  ,kk+1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj+1,kk  ,ist_00p)) / &
           &     ( abs(fsten(ii-2,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-2,jj+1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj+1,kk  ,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj  ,kk+1,ist_pp0)) / &
           &     ( abs(fsten(ii-1,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj  ,kk+1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj+1,kk  ,ist_p0p)) / &
           &     ( abs(fsten(ii-1,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj+1,kk  ,ist_ppp)) + eps)
      r = r      + abs(fsten(ii-1,jj  ,kk  ,ist_0pp)) / &
           &     ( abs(fsten(ii-2,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii-1,jj  ,kk  ,ist_ppp)) + eps)
      r = r * abs(fsten(ii-1,jj  ,kk  ,ist_ppp)) * fsten(ii-1,jj+1,kk+1,ist_inv)
    end function restrict_from_mpp_to

    elemental function restrict_from_0pp_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r, w1m, w1p, w2m, w2p, wmm, wpm, wmp, wpp
      w1m = abs(fsten(ii,jj  ,kk+1,ist_0p0)) / (abs(fsten(ii,jj  ,kk  ,ist_0pp)) &
           &                                   +abs(fsten(ii,jj  ,kk+1,ist_0pp)) + eps)
      w1p = abs(fsten(ii,jj+1,kk+1,ist_0p0)) / (abs(fsten(ii,jj+1,kk  ,ist_0pp)) &
           &                                   +abs(fsten(ii,jj+1,kk+1,ist_0pp)) + eps)
      w2m = abs(fsten(ii,jj+1,kk  ,ist_00p)) / (abs(fsten(ii,jj  ,kk  ,ist_0pp)) &
           &                                   +abs(fsten(ii,jj+1,kk  ,ist_0pp)) + eps)
      w2p = abs(fsten(ii,jj+1,kk+1,ist_00p)) / (abs(fsten(ii,jj  ,kk+1,ist_0pp)) &
           &                                   +abs(fsten(ii,jj+1,kk+1,ist_0pp)) + eps)
      wmm = abs(fsten(ii,jj  ,kk  ,ist_0pp)) * (1.d0 + w1m + w2m)
      wpm = abs(fsten(ii,jj+1,kk  ,ist_0pp)) * (1.d0 + w1p + w2m)
      wmp = abs(fsten(ii,jj  ,kk+1,ist_0pp)) * (1.d0 + w1m + w2p)
      wpp = abs(fsten(ii,jj+1,kk+1,ist_0pp)) * (1.d0 + w1p + w2p)
      r = wmm / (wmm+wpm+wmp+wpp+eps)
    end function restrict_from_0pp_to

    elemental function restrict_from_ppp_to (ii,jj,kk) result(r)
      integer, intent(in) :: ii,jj,kk
      real(amrex_real) :: r
      r = 1.d0
      r = r      + abs(fsten(ii  ,jj+1,kk+1,ist_p00)) / &
           &     ( abs(fsten(ii  ,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii  ,jj+1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii  ,jj  ,kk+1,ist_ppp)) &
           &     + abs(fsten(ii  ,jj+1,kk+1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii+1,jj  ,kk+1,ist_0p0)) / &
           &     ( abs(fsten(ii  ,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii+1,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii  ,jj  ,kk+1,ist_ppp)) &
           &     + abs(fsten(ii+1,jj  ,kk+1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii+1,jj+1,kk  ,ist_00p)) / &
           &     ( abs(fsten(ii  ,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii+1,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii  ,jj+1,kk  ,ist_ppp)) &
           &     + abs(fsten(ii+1,jj+1,kk  ,ist_ppp)) + eps)
      r = r      + abs(fsten(ii  ,jj  ,kk+1,ist_pp0)) / &
           &     ( abs(fsten(ii  ,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii  ,jj  ,kk+1,ist_ppp)) + eps)
      r = r      + abs(fsten(ii  ,jj+1,kk  ,ist_p0p)) / &
           &     ( abs(fsten(ii  ,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii  ,jj+1,kk  ,ist_ppp)) + eps)
      r = r      + abs(fsten(ii+1,jj  ,kk  ,ist_0pp)) / &
           &     ( abs(fsten(ii  ,jj  ,kk  ,ist_ppp)) &
           &     + abs(fsten(ii+1,jj  ,kk  ,ist_ppp)) + eps)
      r = r * abs(fsten(ii  ,jj  ,kk  ,ist_ppp)) * fsten(ii+1,jj+1,kk+1,ist_inv)
    end function restrict_from_ppp_to

  end subroutine amrex_mlndlap_stencil_rap


#ifdef AMREX_USE_EB

  subroutine amrex_mlndlap_set_integral (lo, hi, intg, glo, ghi) &
       bind(c,name='amrex_mlndlap_set_integral')
    integer, dimension(3) :: lo, hi, glo, ghi
    real(amrex_real), intent(inout) :: intg(glo(1):hi(1),lo(2):hi(2),lo(3):hi(3),n_Sintg)
    integer :: i,j,k
    real(amrex_real), parameter :: offth = 1.d0/144.d0
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             intg(i,j,k,i_S_x    ) = 0.d0
             intg(i,j,k,i_S_y    ) = 0.d0
             intg(i,j,k,i_S_z    ) = 0.d0
             intg(i,j,k,i_S_x2   ) = twelfth
             intg(i,j,k,i_S_y2   ) = twelfth
             intg(i,j,k,i_S_z2   ) = twelfth
             intg(i,j,k,i_S_x_y  ) = 0.d0
             intg(i,j,k,i_S_x_z  ) = 0.d0
             intg(i,j,k,i_S_y_z  ) = 0.d0
             intg(i,j,k,i_S_x2_y ) = 0.d0
             intg(i,j,k,i_S_x2_z ) = 0.d0
             intg(i,j,k,i_S_x_y2 ) = 0.d0
             intg(i,j,k,i_S_y2_z ) = 0.d0
             intg(i,j,k,i_S_x_z2 ) = 0.d0
             intg(i,j,k,i_S_y_z2 ) = 0.d0
             intg(i,j,k,i_S_x2_y2) = offth
             intg(i,j,k,i_S_x2_z2) = offth
             intg(i,j,k,i_S_y2_z2) = offth
          end do
       end do
    end do
  end subroutine amrex_mlndlap_set_integral

  subroutine amrex_mlndlap_set_integral_eb (lo, hi, intg, glo, ghi, flag, flo, fhi, &
       vol, vlo, vhi, ax, axlo, axhi, ay, aylo, ayhi, az, azlo, azhi, bcen, blo, bhi) &
       bind(c,name='amrex_mlndlap_set_integral_eb')
    use amrex_ebcellflag_module, only : is_single_valued_cell, is_regular_cell, is_covered_cell
    integer, dimension(3) :: lo, hi, glo, ghi, flo, fhi, vlo, vhi, axlo, axhi, aylo, ayhi, &
         azlo, azhi, blo, bhi
    real(amrex_real), intent(inout) :: intg( glo(1): ghi(1), glo(2): ghi(2), glo(3): ghi(3),n_Sintg)
    real(amrex_real), intent(in   ) :: vol ( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3))
    real(amrex_real), intent(in   ) :: ax  (axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3))
    real(amrex_real), intent(in   ) :: ay  (aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3))
    real(amrex_real), intent(in   ) :: az  (azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))
    real(amrex_real), intent(in   ) :: bcen( blo(1): bhi(1), blo(2): bhi(2), blo(3): bhi(3),3)
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3))

    call amrex_mlndlap_set_integral(lo, hi, intg, glo, ghi)

  end subroutine amrex_mlndlap_set_integral_eb


  subroutine amrex_mlndlap_set_connection (lo, hi, conn, clo, chi, intg, glo, ghi, flag, flo, fhi, &
       vol, vlo, vhi) bind(c,name='amrex_mlndlap_set_connection')
    use amrex_ebcellflag_module, only : is_single_valued_cell, is_regular_cell, is_covered_cell
    integer, dimension(3) :: lo, hi, clo, chi, glo, ghi, flo, fhi, vlo, vhi
    real(amrex_real), intent(inout) :: conn( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3),n_conn)
    real(amrex_real), intent(inout) :: intg( glo(1): ghi(1), glo(2): ghi(2), glo(3): ghi(3),n_Sintg)
    real(amrex_real), intent(in   ) :: vol ( vlo(1): vhi(1), vlo(2): vhi(2), vlo(3): vhi(3))
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3))

    integer :: i,j,k
    real(amrex_real), parameter :: almostone = 1.d0 - 1.d2*epsilon(1._amrex_real)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (is_covered_cell(flag(i,j,k))) then
                conn(i,j,k,:) = 0.d0
             else if (is_regular_cell(flag(i,j,k)) .or. vol(i,j,k).ge.almostone) then
                conn(i,j,k,:) = 1.d0
             else

                ! Scaled by 9
                conn(i,j,k,i_c_xmym) = 0.5625d0*vol(i,j,k) &
                     + 2.25d0*(-intg(i,j,k,i_S_x ) - intg(i,j,k,i_S_y) &
                     &         +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2)) &
                     + 9.d0*(intg(i,j,k,i_S_x_y ) - intg(i,j,k,i_S_x2_y) &
                     &      -intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2))

                ! Scaled by 18
                conn(i,j,k,i_c_xmyb) = 1.125d0*vol(i,j,k) &
                     + 4.5d0*(-intg(i,j,k,i_S_x) + intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_y2)) &
                     + 18.d0*( intg(i,j,k,i_S_x_y2) - intg(i,j,k,i_S_x2_y2))

                ! Scaled by 9
                conn(i,j,k,i_c_xmyp) =  0.5625d0*vol(i,j,k) &
                     + 2.25d0*(-intg(i,j,k,i_S_x ) + intg(i,j,k,i_S_y) &
                     &         +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2)) &
                     + 9.d0*(-intg(i,j,k,i_S_x_y ) + intg(i,j,k,i_S_x2_y) &
                     &       -intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2))

                ! Scaled by 18
                conn(i,j,k,i_c_xbym) = 1.125d0*vol(i,j,k) &
                     + 4.5d0*(-intg(i,j,k,i_S_y) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2)) &
                     + 18.d0*(intg(i,j,k,i_S_x2_y) - intg(i,j,k,i_S_x2_y2))

                ! Scaled by 36
                conn(i,j,k,i_c_xbyb) = 2.25d0*vol(i,j,k) &
                     + 9.d0*(-intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_y2)) &
                     + 36.d0*intg(i,j,k,i_S_x2_y2)

                ! Scaled by 18
                conn(i,j,k,i_c_xbyp) =  1.125d0*vol(i,j,k) &
                     + 4.5d0*( intg(i,j,k,i_S_y) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2)) &
                     + 18.d0*(-intg(i,j,k,i_S_x2_y) - intg(i,j,k,i_S_x2_y2))

                ! Scaled by 9
                conn(i,j,k,i_c_xpym) = 0.5625d0*vol(i,j,k) &
                     + 2.25d0*( intg(i,j,k,i_S_x ) - intg(i,j,k,i_S_y) &
                     &         +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2)) &
                     + 9.d0*(-intg(i,j,k,i_S_x_y ) - intg(i,j,k,i_S_x2_y) &
                     &       +intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2))

                ! Scaled by 18
                conn(i,j,k,i_c_xpyb) = 1.125d0*vol(i,j,k) &
                     + 4.5d0*( intg(i,j,k,i_S_x) + intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_y2)) &
                     + 18.d0*(-intg(i,j,k,i_S_x_y2) - intg(i,j,k,i_S_x2_y2))

                ! Scaled by 9
                conn(i,j,k,i_c_xpyp) =  0.5625d0*vol(i,j,k) &
                     + 2.25d0*( intg(i,j,k,i_S_x ) + intg(i,j,k,i_S_y) &
                     &         +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_y2)) &
                     + 9.d0*( intg(i,j,k,i_S_x_y ) + intg(i,j,k,i_S_x2_y) &
                     &       +intg(i,j,k,i_S_x_y2) + intg(i,j,k,i_S_x2_y2))

                ! Scaled by 9
                conn(i,j,k,i_c_xmzm) = 0.5625d0*vol(i,j,k) &
                     + 2.25d0*(-intg(i,j,k,i_S_x) - intg(i,j,k,i_S_z) &
                     &         +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2)) &
                     + 9.d0*(intg(i,j,k,i_S_x_z) - intg(i,j,k,i_S_x2_z) &
                     &      -intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2))

                ! Scaled by 18
                conn(i,j,k,i_c_xmzb) = 1.125d0*vol(i,j,k) &
                     + 4.5d0*(-intg(i,j,k,i_S_x) + intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_z2)) &
                     + 18.d0*(intg(i,j,k,i_S_x_z2) - intg(i,j,k,i_S_x2_z2))

                ! Scaled by 9
                conn(i,j,k,i_c_xmzp) = 0.5625d0*vol(i,j,k) &
                     + 2.25d0*(-intg(i,j,k,i_S_x  ) + intg(i,j,k,i_S_z) &
                     &         +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2)) &
                     + 9.d0*(-intg(i,j,k,i_S_x_z  ) + intg(i,j,k,i_S_x2_z) &
                     &       -intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2))

                ! Scaled by 18
                conn(i,j,k,i_c_xbzm) = 1.125d0*vol(i,j,k) &
                     + 4.5d0*(-intg(i,j,k,i_S_z) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2)) &
                     + 18.d0*(intg(i,j,k,i_S_x2_z) - intg(i,j,k,i_S_x2_z2))

                ! Scaled by 18
                conn(i,j,k,i_c_xbzb) = 2.25d0*vol(i,j,k) &
                     + 9.d0*(-intg(i,j,k,i_S_x2) - intg(i,j,k,i_S_z2)) &
                     + 36.d0*intg(i,j,k,i_S_x2_z2)

                ! Scaled by 18
                conn(i,j,k,i_c_xbzp) = 1.125d0*vol(i,j,k) &
                     + 4.5d0*( intg(i,j,k,i_S_z) - intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2)) &
                     + 18.d0*(-intg(i,j,k,i_S_x2_z) - intg(i,j,k,i_S_x2_z2))

                ! Scaled by 9
                conn(i,j,k,i_c_xpzm) = 0.5625d0*vol(i,j,k) &
                     + 2.25d0*( intg(i,j,k,i_S_x ) - intg(i,j,k,i_S_z) &
                     &         +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2)) &
                     + 9.d0*(-intg(i,j,k,i_S_x_z ) - intg(i,j,k,i_S_x2_z) &
                     &       +intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2))

                ! Scaled by 18
                conn(i,j,k,i_c_xpzb) = 1.125d0*vol(i,j,k) &
                     + 4.5d0*( intg(i,j,k,i_S_x   ) + intg(i,j,k,i_S_x2   ) - intg(i,j,k,i_S_z2)) &
                     + 18.d0*(-intg(i,j,k,i_S_x_z2) - intg(i,j,k,i_S_x2_z2))

                ! Scaled by 9
                conn(i,j,k,i_c_xpzp) = 0.5625d0*vol(i,j,k) &
                     + 2.25d0*( intg(i,j,k,i_S_x ) + intg(i,j,k,i_S_z) &
                     &         +intg(i,j,k,i_S_x2) + intg(i,j,k,i_S_z2)) &
                     + 9.d0*( intg(i,j,k,i_S_x_z ) + intg(i,j,k,i_S_x2_z) &
                     &       +intg(i,j,k,i_S_x_z2) + intg(i,j,k,i_S_x2_z2))

                ! Scaled by 9
                conn(i,j,k,i_c_ymzm) = 0.5625d0*vol(i,j,k) &
                     + 2.25d0*(-intg(i,j,k,i_S_y) - intg(i,j,k,i_S_z) &
                     &         +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2)) &
                     + 9.d0*(intg(i,j,k,i_S_y_z) - intg(i,j,k,i_S_y2_z) &
                     &      -intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2))

                ! Scaled by 18
                conn(i,j,k,i_c_ymzb) = 1.125d0*vol(i,j,k) &
                     + 4.5d0*(-intg(i,j,k,i_S_y) + intg(i,j,k,i_S_y2) - intg(i,j,k,i_S_z2)) &
                     + 18.d0*(intg(i,j,k,i_S_y_z2) - intg(i,j,k,i_S_y2_z2))

                ! Scaled by 9
                conn(i,j,k,i_c_ymzp) = 0.5625d0*vol(i,j,k) &
                     + 2.25d0*(-intg(i,j,k,i_S_y ) + intg(i,j,k,i_S_z) &
                     &         +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2)) &
                     + 9.d0*(-intg(i,j,k,i_S_y_z ) + intg(i,j,k,i_S_y2_z) &
                     &       -intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2))

                ! Scaled by 18
                conn(i,j,k,i_c_ybzm) = 1.125d0*vol(i,j,k) &
                     + 4.5d0*(-intg(i,j,k,i_S_z) - intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2)) &
                     + 18.d0*(intg(i,j,k,i_S_y2_z) - intg(i,j,k,i_S_y2_z2))

                ! Scaled by 36
                conn(i,j,k,i_c_ybzb) = 2.25d0*vol(i,j,k) &
                     + 9.d0*(-intg(i,j,k,i_S_y2) - intg(i,j,k,i_S_z2)) &
                     + 36.d0*intg(i,j,k,i_S_y2_z2)

                ! Scaled by 18
                conn(i,j,k,i_c_ybzp) = 1.125d0*vol(i,j,k) &
                     + 4.5d0*( intg(i,j,k,i_S_z) - intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2)) &
                     + 18.d0*(-intg(i,j,k,i_S_y2_z) - intg(i,j,k,i_S_y2_z2))

                ! Scaled by 9
                conn(i,j,k,i_c_ypzm) = 0.5625d0*vol(i,j,k) &
                     + 2.25d0*( intg(i,j,k,i_S_y ) - intg(i,j,k,i_S_z) &
                     &         +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2)) &
                     + 9.d0*(-intg(i,j,k,i_S_y_z ) - intg(i,j,k,i_S_y2_z) &
                     &       +intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2))

                ! Scaled by 18
                conn(i,j,k,i_c_ypzb) = 1.125d0*vol(i,j,k) &
                     + 4.5d0*( intg(i,j,k,i_S_y   ) + intg(i,j,k,i_S_y2) - intg(i,j,k,i_S_z2)) &
                     + 18.d0*(-intg(i,j,k,i_S_y_z2) - intg(i,j,k,i_S_y2_z2))

                ! Scaled by 9
                conn(i,j,k,i_c_ypzp) = 0.5625d0*vol(i,j,k) &
                     + 2.25d0*( intg(i,j,k,i_S_y ) + intg(i,j,k,i_S_z) &
                     &         +intg(i,j,k,i_S_y2) + intg(i,j,k,i_S_z2)) &
                     + 9.d0*( intg(i,j,k,i_S_y_z ) + intg(i,j,k,i_S_y2_z) &
                     &       +intg(i,j,k,i_S_y_z2) + intg(i,j,k,i_S_y2_z2))

             end if
          end do
       end do
    end do

  end subroutine amrex_mlndlap_set_connection

  subroutine amrex_mlndlap_set_stencil_eb (lo, hi, sten, tlo, thi, sig, glo, ghi, &
       conn, clo, chi, dxinv) bind(c,name='amrex_mlndlap_set_stencil_eb')
    integer, dimension(3), intent(in) :: lo, hi, tlo, thi, glo, ghi, clo, chi
    real(amrex_real), intent(inout) :: sten(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),n_sten)
    real(amrex_real), intent(in   ) ::  sig(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3))
    real(amrex_real), intent(in   ) :: conn(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),n_conn)
    real(amrex_real), intent(in) :: dxinv(3)

    integer :: i,j,k
    real(amrex_real) :: facx, facy, facz
    
    facx = (1.d0/36.d0)*dxinv(1)*dxinv(1)
    facy = (1.d0/36.d0)*dxinv(2)*dxinv(2)
    facz = (1.d0/36.d0)*dxinv(3)*dxinv(3)

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! i+1,j,k
             sten(i,j,k,ist_p00) = ( &
                  sig(i,j  ,k  )*(4.d0*facx*conn(i,j  ,k  ,i_c_ymzm) - 2.d0*facy*conn(i,j  ,k  ,i_c_xbzm) - 2.d0*facz*conn(i,j  ,k  ,i_c_xbym) ) + &
                  sig(i,j-1,k  )*(4.d0*facx*conn(i,j-1,k  ,i_c_ypzm) - 2.d0*facy*conn(i,j-1,k  ,i_c_xbzm) - 2.d0*facz*conn(i,j-1,k  ,i_c_xbyp) ) + &
                  sig(i,j  ,k-1)*(4.d0*facx*conn(i,j  ,k-1,i_c_ymzp) - 2.d0*facy*conn(i,j  ,k-1,i_c_xbzp) - 2.d0*facz*conn(i,j  ,k-1,i_c_xbym) ) + &
                  sig(i,j-1,k-1)*(4.d0*facx*conn(i,j-1,k-1,i_c_ypzp) - 2.d0*facy*conn(i,j-1,k-1,i_c_xbzp) - 2.d0*facz*conn(i,j-1,k-1,i_c_xbyp) ) )

             ! i,j+1,k
             sten(i,j,k,ist_0p0) = ( &
                  sig(i  ,j,k  )*(-2.d0*facx*conn(i  ,j,k  ,i_c_ybzm) + 4.d0*facy*conn(i  ,j,k  ,i_c_xmzm) - 2.d0*facz*conn(i  ,j,k  ,i_c_xmyb) ) + &
                  sig(i-1,j,k  )*(-2.d0*facx*conn(i-1,j,k  ,i_c_ybzm) + 4.d0*facy*conn(i-1,j,k  ,i_c_xpzm) - 2.d0*facz*conn(i-1,j,k  ,i_c_xpyb) ) + &
                  sig(i  ,j,k-1)*(-2.d0*facx*conn(i  ,j,k-1,i_c_ybzp) + 4.d0*facy*conn(i  ,j,k-1,i_c_xmzp) - 2.d0*facz*conn(i  ,j,k-1,i_c_xmyb) ) + &
                  sig(i-1,j,k-1)*(-2.d0*facx*conn(i-1,j,k-1,i_c_ybzp) + 4.d0*facy*conn(i-1,j,k-1,i_c_xpzp) - 2.d0*facz*conn(i-1,j,k-1,i_c_xpyb) ) )

             ! i,j,k+1
             sten(i,j,k,ist_00p) = ( &
                  sig(i  ,j  ,k)*(-2.d0*facx*conn(i  ,j  ,k,i_c_ymzb) - 2.d0*facy*conn(i  ,j  ,k,i_c_xmzb) + 4.d0*facz*conn(i  ,j  ,k,i_c_xmym) ) + &
                  sig(i-1,j  ,k)*(-2.d0*facx*conn(i-1,j  ,k,i_c_ymzb) - 2.d0*facy*conn(i-1,j  ,k,i_c_xpzb) + 4.d0*facz*conn(i-1,j  ,k,i_c_xpym) ) + &
                  sig(i  ,j-1,k)*(-2.d0*facx*conn(i  ,j-1,k,i_c_ypzb) - 2.d0*facy*conn(i  ,j-1,k,i_c_xmzb) + 4.d0*facz*conn(i  ,j-1,k,i_c_xmyp) ) + &
                  sig(i-1,j-1,k)*(-2.d0*facx*conn(i-1,j-1,k,i_c_ypzb) - 2.d0*facy*conn(i-1,j-1,k,i_c_xpzb) + 4.d0*facz*conn(i-1,j-1,k,i_c_xpyp) ) )

             ! i+1,j+1,k
             sten(i,j,k,ist_pp0) = ( &
                  sig(i,j,k  )*(2.d0*facx*conn(i,j,k  ,i_c_ybzm) + 2.d0*facy*conn(i,j,k  ,i_c_xbzm) - facz*conn(i,j,k  ,i_c_xbyb) ) + &
                  sig(i,j,k-1)*(2.d0*facx*conn(i,j,k-1,i_c_ybzp) + 2.d0*facy*conn(i,j,k-1,i_c_xbzp) - facz*conn(i,j,k-1,i_c_xbyb) ) )
             
             ! i+1,j,k+1
             sten(i,j,k,ist_p0p) = ( &
                  sig(i,j,k  )*(2.d0*facx*conn(i,j,k  ,i_c_ymzb) - facy*conn(i,j,k  ,i_c_xbzb) + 2.d0*facz*conn(i,j,k  ,i_c_xbym) ) + &
                  sig(i,j-1,k)*(2.d0*facx*conn(i,j-1,k,i_c_ypzb) - facy*conn(i,j-1,k,i_c_xbzb) + 2.d0*facz*conn(i,j-1,k,i_c_xbyp) ) )

             ! i,j+1,k+1
             sten(i,j,k,ist_0pp) = ( &
                  sig(i  ,j,k)*(-facx*conn(i  ,j,k,i_c_ybzb) + 2.d0*facy*conn(i  ,j,k,i_c_xmzb) + 2.d0*facz*conn(i  ,j,k,i_c_xmyb) ) + &
                  sig(i-1,j,k)*(-facx*conn(i-1,j,k,i_c_ybzb) + 2.d0*facy*conn(i-1,j,k,i_c_xpzb) + 2.d0*facz*conn(i-1,j,k,i_c_xpyb) ) )

             ! i+1,j+1,k+1
             sten(i,j,k,ist_ppp) = sig(i,j,k) * (facx*conn(i,j,k,i_c_ybzb) + facy*conn(i,j,k,i_c_xbzb) + facz*conn(i,j,k,i_c_xbyb) ) 
          end do
       end do
    end do

  end subroutine amrex_mlndlap_set_stencil_eb

  subroutine amrex_mlndlap_divu_eb (lo, hi, rhs, rlo, rhi, vel, vlo, vhi, vfrac, flo, fhi, &
       intg, glo, ghi, msk, mlo, mhi, dxinv) &
       bind(c,name='amrex_mlndlap_divu_eb')
    integer, dimension(3), intent(in) :: lo, hi, rlo, rhi, vlo, vhi, flo, fhi, glo, ghi, mlo, mhi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(amrex_real), intent(in   ) :: vel(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),3)
    real(amrex_real), intent(in   ) :: vfrac(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(amrex_real), intent(in   ) :: intg(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),n_Sintg)
    integer, intent(in) :: msk(mlo(1):mhi(1),mlo(2):mhi(2),mlo(3):mhi(3))
    integer :: i,j,k
    real(amrex_real) :: facx, facy, facz

    facx = 0.25d0*dxinv(1)
    facy = 0.25d0*dxinv(2)
    facz = 0.25d0*dxinv(3)

    do    k = lo(3), hi(3)
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (msk(i,j,k) .ne. dirichlet) then

             rhs(i,j,k) = facx*( &
                  &              vel(i-1,j-1,k  ,1)*(    -vfrac(i-1,j-1,k  ) &
                  &                                  -2.d0*intg(i-1,j-1,k  ,i_S_y) &
                  &                                  +2.d0*intg(i-1,j-1,k  ,i_S_z) &
                  &                                  +4.d0*intg(i-1,j-1,k  ,i_S_y_z)) &
                  &             +vel(i  ,j-1,k  ,1)*(     vfrac(i  ,j-1,k  ) &
                  &                                  +2.d0*intg(i  ,j-1,k  ,i_S_y) &
                  &                                  -2.d0*intg(i  ,j-1,k  ,i_S_z) &
                  &                                  -4.d0*intg(i  ,j-1,k  ,i_S_y_z)) &
                  &             +vel(i-1,j  ,k  ,1)*(    -vfrac(i-1,j  ,k  ) &
                  &                                  +2.d0*intg(i-1,j  ,k  ,i_S_y) &
                  &                                  +2.d0*intg(i-1,j  ,k  ,i_S_z) &
                  &                                  -4.d0*intg(i-1,j  ,k  ,i_S_y_z)) &
                  &             +vel(i  ,j  ,k  ,1)*(     vfrac(i  ,j  ,k  ) &
                  &                                  -2.d0*intg(i  ,j  ,k  ,i_S_y) &
                  &                                  -2.d0*intg(i  ,j  ,k  ,i_S_z) &
                  &                                  +4.d0*intg(i  ,j  ,k  ,i_S_y_z)) &
                  &             +vel(i-1,j-1,k-1,1)*(    -vfrac(i-1,j-1,k-1) &
                  &                                  -2.d0*intg(i-1,j-1,k-1,i_S_y) &
                  &                                  -2.d0*intg(i-1,j-1,k-1,i_S_z) &
                  &                                  -4.d0*intg(i-1,j-1,k-1,i_S_y_z)) &
                  &             +vel(i  ,j-1,k-1,1)*(     vfrac(i  ,j-1,k-1) &
                  &                                  +2.d0*intg(i  ,j-1,k-1,i_S_y) &
                  &                                  +2.d0*intg(i  ,j-1,k-1,i_S_z) &
                  &                                  +4.d0*intg(i  ,j-1,k-1,i_S_y_z)) &
                  &             +vel(i-1,j  ,k-1,1)*(    -vfrac(i-1,j  ,k-1) &
                  &                                  +2.d0*intg(i-1,j  ,k-1,i_S_y) &
                  &                                  -2.d0*intg(i-1,j  ,k-1,i_S_z) &
                  &                                  +4.d0*intg(i-1,j  ,k-1,i_S_y_z)) &
                  &             +vel(i  ,j  ,k-1,1)*(     vfrac(i  ,j  ,k-1) &
                  &                                  -2.d0*intg(i  ,j  ,k-1,i_S_y) &
                  &                                  +2.d0*intg(i  ,j  ,k-1,i_S_z) &
                  &                                  -4.d0*intg(i  ,j  ,k-1,i_S_y_z)) ) &
                  &     + facy*( &
                  &              vel(i-1,j-1,k  ,2)*(    -vfrac(i-1,j-1,k  ) &
                  &                                  -2.d0*intg(i-1,j-1,k  ,i_S_x) &
                  &                                  +2.d0*intg(i-1,j-1,k  ,i_S_z) &
                  &                                  +4.d0*intg(i-1,j-1,k  ,i_S_x_z)) &
                  &             +vel(i  ,j-1,k  ,2)*(    -vfrac(i  ,j-1,k  ) &
                  &                                  +2.d0*intg(i  ,j-1,k  ,i_S_x) &
                  &                                  +2.d0*intg(i  ,j-1,k  ,i_S_z) &
                  &                                  -4.d0*intg(i  ,j-1,k  ,i_S_x_z)) &
                  &             +vel(i-1,j  ,k  ,2)*(     vfrac(i-1,j  ,k  ) &
                  &                                  +2.d0*intg(i-1,j  ,k  ,i_S_x) &
                  &                                  -2.d0*intg(i-1,j  ,k  ,i_S_z) &
                  &                                  -4.d0*intg(i-1,j  ,k  ,i_S_x_z)) &
                  &             +vel(i  ,j  ,k  ,2)*(     vfrac(i  ,j  ,k  ) &
                  &                                  -2.d0*intg(i  ,j  ,k  ,i_S_x) &
                  &                                  -2.d0*intg(i  ,j  ,k  ,i_S_z) &
                  &                                  +4.d0*intg(i  ,j  ,k  ,i_S_x_z)) &
                  &             +vel(i-1,j-1,k-1,2)*(    -vfrac(i-1,j-1,k-1) &
                  &                                  -2.d0*intg(i-1,j-1,k-1,i_S_x) &
                  &                                  -2.d0*intg(i-1,j-1,k-1,i_S_z) &
                  &                                  -4.d0*intg(i-1,j-1,k-1,i_S_x_z)) &
                  &             +vel(i  ,j-1,k-1,2)*(    -vfrac(i  ,j-1,k-1) &
                  &                                  +2.d0*intg(i  ,j-1,k-1,i_S_x) &
                  &                                  -2.d0*intg(i  ,j-1,k-1,i_S_z) &
                  &                                  +4.d0*intg(i  ,j-1,k-1,i_S_x_z)) &
                  &             +vel(i-1,j  ,k-1,2)*(     vfrac(i-1,j  ,k-1) &
                  &                                  +2.d0*intg(i-1,j  ,k-1,i_S_x) &
                  &                                  +2.d0*intg(i-1,j  ,k-1,i_S_z) &
                  &                                  +4.d0*intg(i-1,j  ,k-1,i_S_x_z)) &
                  &             +vel(i  ,j  ,k-1,2)*(     vfrac(i  ,j  ,k-1) &
                  &                                  -2.d0*intg(i  ,j  ,k-1,i_S_x) &
                  &                                  +2.d0*intg(i  ,j  ,k-1,i_S_z) &
                  &                                  -4.d0*intg(i  ,j  ,k-1,i_S_x_z)) ) &
                  &     + facz*( &
                  &              vel(i-1,j-1,k  ,3)*(     vfrac(i-1,j-1,k  ) &
                  &                                  +2.d0*intg(i-1,j-1,k  ,i_S_x) &
                  &                                  +2.d0*intg(i-1,j-1,k  ,i_S_y) &
                  &                                  +4.d0*intg(i-1,j-1,k  ,i_S_x_y)) &
                  &             +vel(i  ,j-1,k  ,3)*(     vfrac(i  ,j-1,k  ) &
                  &                                  -2.d0*intg(i  ,j-1,k  ,i_S_x) &
                  &                                  +2.d0*intg(i  ,j-1,k  ,i_S_y) &
                  &                                  -4.d0*intg(i  ,j-1,k  ,i_S_x_y)) &
                  &             +vel(i-1,j  ,k  ,3)*(     vfrac(i-1,j  ,k  ) &
                  &                                  +2.d0*intg(i-1,j  ,k  ,i_S_x) &
                  &                                  -2.d0*intg(i-1,j  ,k  ,i_S_y) &
                  &                                  -4.d0*intg(i-1,j  ,k  ,i_S_x_y)) &
                  &             +vel(i  ,j  ,k  ,3)*(     vfrac(i  ,j  ,k  ) &
                  &                                  -2.d0*intg(i  ,j  ,k  ,i_S_x) &
                  &                                  -2.d0*intg(i  ,j  ,k  ,i_S_y) &
                  &                                  +4.d0*intg(i  ,j  ,k  ,i_S_x_y)) &
                  &             +vel(i-1,j-1,k-1,3)*(    -vfrac(i-1,j-1,k-1) &
                  &                                  -2.d0*intg(i-1,j-1,k-1,i_S_x) &
                  &                                  -2.d0*intg(i-1,j-1,k-1,i_S_y) &
                  &                                  -4.d0*intg(i-1,j-1,k-1,i_S_x_y)) &
                  &             +vel(i  ,j-1,k-1,3)*(    -vfrac(i  ,j-1,k-1) &
                  &                                  +2.d0*intg(i  ,j-1,k-1,i_S_x) &
                  &                                  -2.d0*intg(i  ,j-1,k-1,i_S_y) &
                  &                                  +4.d0*intg(i  ,j-1,k-1,i_S_x_y)) &
                  &             +vel(i-1,j  ,k-1,3)*(    -vfrac(i-1,j  ,k-1) &
                  &                                  -2.d0*intg(i-1,j  ,k-1,i_S_x) &
                  &                                  +2.d0*intg(i-1,j  ,k-1,i_S_y) &
                  &                                  +4.d0*intg(i-1,j  ,k-1,i_S_x_y)) &
                  &             +vel(i  ,j  ,k-1,3)*(    -vfrac(i  ,j  ,k-1) &
                  &                                  +2.d0*intg(i  ,j  ,k-1,i_S_x) &
                  &                                  +2.d0*intg(i  ,j  ,k-1,i_S_y) &
                  &                                  -4.d0*intg(i  ,j  ,k-1,i_S_x_y)) )
          else
             rhs(i,j,k) = 0.d0
          end if
       end do
    end do
    end do

  end subroutine amrex_mlndlap_divu_eb

  subroutine amrex_mlndlap_mknewu_eb (lo, hi, u, ulo, uhi, p, plo, phi, sig, slo, shi, &
       vfrac, vlo, vhi, intg, glo, ghi, dxinv) bind(c,name='amrex_mlndlap_mknewu_eb')
    integer, dimension(3), intent(in) :: lo, hi, ulo, uhi, plo, phi, slo, shi, vlo, vhi, glo, ghi
    real(amrex_real), intent(in) :: dxinv(3)
    real(amrex_real), intent(inout) ::   u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),3)
    real(amrex_real), intent(in   ) ::   p(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3))
    real(amrex_real), intent(in   ) :: sig(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3))
    real(amrex_real), intent(in   )::vfrac(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
    real(amrex_real), intent(in   ) ::intg(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),n_Sintg)

    integer :: i, j, k
    real(amrex_real) :: dpdx, dpdy, dpdz, dpp_xy, dpp_xz, dpp_yz, dpp_xyz

    do    k = lo(3), hi(3)
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (vfrac(i,j,k) .eq. 0.d0) then
             u(i,j,k,1) = 0.d0
             u(i,j,k,2) = 0.d0
             u(i,j,k,3) = 0.d0
          else
             dpdx = 0.25d0*(-p(i,j,k  )+p(i+1,j,k  )-p(i,j+1,k  )+p(i+1,j+1,k  ) &
                            -p(i,j,k+1)+p(i+1,j,k+1)-p(i,j+1,k+1)+p(i+1,j+1,k+1))
             dpdy = 0.25d0*(-p(i,j,k  )-p(i+1,j,k  )+p(i,j+1,k  )+p(i+1,j+1,k  ) &
                            -p(i,j,k+1)-p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1))
             dpdz = 0.25d0*(-p(i,j,k  )-p(i+1,j,k  )-p(i,j+1,k  )-p(i+1,j+1,k  ) &
                            +p(i,j,k+1)+p(i+1,j,k+1)+p(i,j+1,k+1)+p(i+1,j+1,k+1))

             dpp_xy = (p(i+1,j+1,k+1) - p(i,j+1,k+1) - p(i+1,j,k+1) + p(i,j,k+1) &
                      +p(i+1,j+1,k  ) - p(i,j+1,k  ) - p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k)

             dpp_xz = (p(i+1,j+1,k+1) - p(i,j+1,k+1) + p(i+1,j,k+1) - p(i,j,k+1) &
                      -p(i+1,j+1,k  ) + p(i,j+1,k  ) - p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k)

             dpp_yz = (p(i+1,j+1,k+1) + p(i,j+1,k+1) - p(i+1,j,k+1) - p(i,j,k+1) &
                      -p(i+1,j+1,k  ) - p(i,j+1,k  ) + p(i+1,j,k  ) + p(i,j,k  ) ) / vfrac(i,j,k)

             dpp_xyz = (p(i+1,j+1,k+1) - p(i,j+1,k+1) - p(i+1,j,k+1) + p(i,j,k+1) &
                       -p(i+1,j+1,k  ) + p(i,j+1,k  ) + p(i+1,j,k  ) - p(i,j,k  ) ) / vfrac(i,j,k)

             u(i,j,k,1) = u(i,j,k,1) - sig(i,j,k)*dxinv(1)*(dpdx + .5d0*intg(i,j,k,i_S_y  )*dpp_xy + &
                                                                   .5d0*intg(i,j,k,i_S_z  )*dpp_xz + &
                                                                        intg(i,j,k,i_S_y_z)*dpp_xyz )
             u(i,j,k,2) = u(i,j,k,2) - sig(i,j,k)*dxinv(2)*(dpdy + .5d0*intg(i,j,k,i_S_x  )*dpp_xy + &
                                                                   .5d0*intg(i,j,k,i_S_z  )*dpp_yz + &
                                                                        intg(i,j,k,i_S_x_z)*dpp_xyz )
             u(i,j,k,3) = u(i,j,k,3) - sig(i,j,k)*dxinv(3)*(dpdz + .5d0*intg(i,j,k,i_S_x  )*dpp_xz + &
                                                                   .5d0*intg(i,j,k,i_S_y  )*dpp_yz + & 
                                                                        intg(i,j,k,i_S_x_y)*dpp_xyz )

          end if
       end do
    end do
    end do

  end subroutine amrex_mlndlap_mknewu_eb

#endif

end module amrex_mlnodelap_3d_module
