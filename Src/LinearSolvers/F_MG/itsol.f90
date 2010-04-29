module itsol_module

  use bl_types
  use multifab_module
  use stencil_module

  implicit none

  integer, private, parameter :: def_bicg_max_iter = 1000
  integer, private, parameter :: def_cg_max_iter = 1000
  real(dp_t), private, parameter :: ZERO = 0.0_dp_t
  real(dp_t), private, parameter :: ONE  = 1.0_dp_t

  private :: itsol_defect, itsol_precon
  private :: jacobi_precon_1d, jacobi_precon_2d, jacobi_precon_3d
  private :: nodal_precon_1d, nodal_precon_2d, nodal_precon_3d

contains

    subroutine jacobi_precon_1d(a, u, r, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)  :: a(:,0:)
      real(kind=dp_t), intent(inout) :: u(1-ng:)
      real(kind=dp_t), intent(in)  :: r(:)
      integer :: i, nx
      nx = size(a,dim=1)
      do i = 1, nx
         u(i) = r(i)/a(i,0)
      end do
    end subroutine jacobi_precon_1d

    subroutine jacobi_precon_2d(a, u, r, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)  :: a(:,:,0:)
      real(kind=dp_t), intent(inout) :: u(1-ng:,1-ng:)
      real(kind=dp_t), intent(in)  :: r(:,:)
      integer :: i, j, nx, ny
      ny = size(a,dim=2)
      nx = size(a,dim=1)
      do j = 1, ny
         do i = 1, nx
            u(i,j) = r(i,j)/a(i,j,0)
         end do
      end do
    end subroutine jacobi_precon_2d

    subroutine jacobi_precon_3d(a, u, r, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)  :: a(:,:,:,0:)
      real(kind=dp_t), intent(inout) :: u(1-ng:,1-ng:,1-ng:)
      real(kind=dp_t), intent(in)  :: r(:,:,:)
      integer i, j, k, nx, ny, nz
      nz = size(a,dim=3)
      ny = size(a,dim=2)
      nx = size(a,dim=1)
      !$OMP PARALLEL DO PRIVATE(j,i,k) IF(nz.ge.4)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               u(i,j,k) = r(i,j,k)/a(i,j,k,0)
            end do
         end do
      end do
      !$OMP END PARALLEL DO
    end subroutine jacobi_precon_3d

    subroutine nodal_precon_1d(a, u, r, mm, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)  :: a(:,0:)
      real(kind=dp_t), intent(inout) :: u(1-ng:)
      real(kind=dp_t), intent(in)  :: r(0:)
      integer, intent(in)  :: mm(:)
      integer :: i, nx
      nx = size(a,dim=1)
      do i = 1, nx
         if (.not. bc_dirichlet(mm(i),1,0)) &
            u(i) = r(i)/a(i,0)
      end do
    end subroutine nodal_precon_1d

    subroutine nodal_precon_2d(a, u, r, mm, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)  :: a(:,:,0:)
      real(kind=dp_t), intent(inout) :: u(1-ng:,1-ng:)
      real(kind=dp_t), intent(in)  :: r(0:,0:)
      integer, intent(in)  :: mm(:,:)
      integer :: i, j, nx, ny
      ny = size(a,dim=2)
      nx = size(a,dim=1)
      do j = 1, ny
         do i = 1, nx
            if (.not. bc_dirichlet(mm(i,j),1,0)) then
               u(i,j) = r(i,j)/a(i,j,0)
            end if
         end do
      end do
    end subroutine nodal_precon_2d

    subroutine nodal_precon_3d(a, u, r, mm, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)  :: a(:,:,:,0:)
      real(kind=dp_t), intent(inout) :: u(1-ng:,1-ng:,1-ng:)
      real(kind=dp_t), intent(in)  :: r(0:,0:,0:)
      integer, intent(in)  :: mm(:,:,:)
      integer :: i, j, k, nx, ny, nz
      nz = size(a,dim=3)
      ny = size(a,dim=2)
      nx = size(a,dim=1)
      !$OMP PARALLEL DO PRIVATE(j,i,k) IF(nz.ge.4)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (.not. bc_dirichlet(mm(i,j,k),1,0)) then
                  u(i,j,k) = r(i,j,k)/a(i,j,k,0)
               end if
            end do
         end do
      end do
      !$OMP END PARALLEL DO
    end subroutine nodal_precon_3d

    subroutine diag_init_cc_1d(a, ng_a, r, ng_r, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(lo(1)-ng_a:,0:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:   )

      integer         :: i, nc
      real(kind=dp_t) :: denom

      nc = size(a,dim=2)-1

      do i = lo(1),hi(1)
         denom = 1.d0 / a(i,0)
         r(i     ) = r(i     ) * denom
         a(i,1:nc) = a(i,1:nc) * denom
         a(i,0   ) = 1.d0
      end do

    end subroutine diag_init_cc_1d

    subroutine diag_init_cc_2d(a, ng_a, r, ng_r, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(lo(1)-ng_a:,lo(2)-ng_a:,0:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:,lo(2)-ng_r:   )

      integer         :: i, j, nc
      real(kind=dp_t) :: denom

      nc = size(a,dim=3)-1

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            denom = 1.d0 / a(i,j,0)
            r(i,j     ) = r(i,j     ) * denom
            a(i,j,1:nc) = a(i,j,1:nc) * denom
            a(i,j,0   ) = 1.d0
         end do
      end do

    end subroutine diag_init_cc_2d

    subroutine diag_init_cc_3d(a, ng_a, r, ng_r, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:,0:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:   )

      integer         :: i, j, k, nc
      real(kind=dp_t) :: denom

      nc = size(a,dim=4)-1

      !$OMP PARALLEL DO PRIVATE(j,i,k) IF(nz.ge.4)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               denom = 1.d0 / a(i,j,k,0)
               r(i,j,k     ) = r(i,j,k     ) * denom
               a(i,j,k,1:nc) = a(i,j,k,1:nc) * denom
               a(i,j,k,0   ) = 1.d0
            end do
         end do
      end do

      !$OMP END PARALLEL DO
    end subroutine diag_init_cc_3d

    subroutine diag_init_nd_1d(a, ng_a, r, ng_r, mm, ng_m, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r, ng_m
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(lo(1)-ng_a:,0:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:   )
      integer        , intent(inout)  :: mm(lo(1)-ng_m:   )

      integer         :: i, nc
      real(kind=dp_t) :: denom

      nc = size(a,dim=2)-1

      do i = lo(1),hi(1)+1
         if (.not. bc_dirichlet(mm(i),1,0)) then
            denom = 1.d0 / a(i,0)
            r(i     ) = r(i     ) * denom
            a(i,1:nc) = a(i,1:nc) * denom
            a(i,0   ) = 1.d0
         end if
      end do

    end subroutine diag_init_nd_1d

    subroutine diag_init_nd_2d(a, ng_a, r, ng_r, mm, ng_m, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r, ng_m
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(lo(1)-ng_a:,lo(2)-ng_a:,0:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:,lo(2)-ng_r:   )
      integer        , intent(inout)  :: mm(lo(1)-ng_m:,lo(2)-ng_m:   )

      integer         :: i, j, nc
      real(kind=dp_t) :: denom

      nc = size(a,dim=3)-1

      do j = lo(2),hi(2)+1
         do i = lo(1),hi(1)+1
            if (.not. bc_dirichlet(mm(i,j),1,0)) then
               denom = 1.d0 / a(i,j,0)
               r(i,j     ) = r(i,j     ) * denom
               a(i,j,1:nc) = a(i,j,1:nc) * denom
               a(i,j,0   ) = 1.d0
            end if
         end do
      end do

    end subroutine diag_init_nd_2d

    subroutine diag_init_nd_3d(a, ng_a, r, ng_r, mm, ng_m, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r, ng_m
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:,0:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:   )
      integer        , intent(inout)  :: mm(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:   )

      integer         :: i, j, k, nc, nz
      real(kind=dp_t) :: denom

      nc = size(a,dim=4)-1
      nz = size(a,dim=3)

      !$OMP PARALLEL DO PRIVATE(j,i,k) IF(nz.ge.4)
      do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)+1
               if (.not. bc_dirichlet(mm(i,j,k),1,0)) then
                  denom = 1.d0 / a(i,j,k,0)
                  r(i,j,k     ) = r(i,j,k     ) * denom
                  a(i,j,k,1:nc) = a(i,j,k,1:nc) * denom
                  a(i,j,k,0   ) = 1.d0
               end if
            end do
         end do
      end do

      !$OMP END PARALLEL DO
    end subroutine diag_init_nd_3d

  function itsol_converged(rr, uu, Anorm, bnorm, eps, abs_eps) result(r)
    use bl_prof_module
    type(multifab), intent(in) :: rr, uu
    real(dp_t), intent(in) :: Anorm, bnorm, eps
    real(dp_t), intent(in), optional :: abs_eps

    real(dp_t)             :: norm_rr, norm_uu
    logical :: r
    type(bl_prof_timer), save :: bpt

    call build(bpt, "its_converged")
    norm_rr = norm_inf(rr)
    norm_uu = norm_inf(uu)
    if (present(abs_eps)) then
!     r = (norm_rr <= eps*(Anorm*norm_uu + bnorm)) .or. &
!         (norm_rr <= epsilon(Anorm)*Anorm) .or. &
!         (norm_rr <= abs_eps)
      r = (norm_rr <= eps*(bnorm)) .or. &
          (norm_rr <= abs_eps)
    else
!     r = (norm_rr <= eps*(Anorm*norm_uu + bnorm)) .or. &
!         (norm_rr <= epsilon(Anorm)*Anorm)
      r = (norm_rr <= eps*(bnorm)) 
    endif
    call destroy(bpt)
  end function itsol_converged

  ! computes rr = aa * uu
  subroutine itsol_stencil_apply(aa, rr, uu, mm, uniform_dh)
    use bl_prof_module
    use stencil_nodal_module, only: stencil_apply_1d_nodal, stencil_apply_2d_nodal,  stencil_apply_3d_nodal
    type(multifab), intent(in) :: aa
    type(multifab), intent(inout) :: rr
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in) :: mm
    logical, intent(in),optional :: uniform_dh
    logical :: luniform_dh
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: ap(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)

    integer :: i, n, lo(rr%dim), hi(rr%dim)
    logical :: nodal_flag, lcross
    type(bl_prof_timer), save :: bpt

    call build(bpt, "its_stencil_apply")

    luniform_dh = .false. ; if ( present(uniform_dh) ) luniform_dh = uniform_dh

    lcross = ((ncomp(aa) == 5) .or. (ncomp(aa) == 7))

    call multifab_fill_boundary(uu, cross = lcross)

    nodal_flag = nodal_q(uu)

    do i = 1, rr%nboxes
       if ( multifab_remote(rr, i) ) cycle
       rp => dataptr(rr, i)
       up => dataptr(uu, i)
       ap => dataptr(aa, i)
       mp => dataptr(mm, i)
       lo = lwb(get_box(uu,i))
       hi = upb(get_box(uu,i))
       do n = 1, rr%nc
          select case(rr%dim)
          case (1)
             if ( .not. nodal_flag) then
                call stencil_apply_1d(ap(:,1,1,:), rp(:,1,1,n), rr%ng, up(:,1,1,n), uu%ng,  &
                                      mp(:,1,1,1), lo, hi)
             else
                call stencil_apply_1d_nodal(ap(:,1,1,:), rp(:,1,1,n), up(:,1,1,n),  &
                     mp(:,1,1,1), uu%ng)
             end if
          case (2)
             if ( .not. nodal_flag) then
                call stencil_apply_2d(ap(:,:,1,:), rp(:,:,1,n), rr%ng, up(:,:,1,n), uu%ng,  &
                     mp(:,:,1,1), lo, hi)
             else
                call stencil_apply_2d_nodal(ap(:,:,1,:), rp(:,:,1,n), up(:,:,1,n),  &
                     mp(:,:,1,1), uu%ng)
             end if
          case (3)
             if ( .not. nodal_flag) then
                call stencil_apply_3d(ap(:,:,:,:), rp(:,:,:,n), rr%ng, up(:,:,:,n), uu%ng,  &
                     mp(:,:,:,1))
             else
                call stencil_apply_3d_nodal(ap(:,:,:,:), rp(:,:,:,n), up(:,:,:,n),  &
                     mp(:,:,:,1), uu%ng, luniform_dh)
             end if
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine itsol_stencil_apply

  ! computes rr = aa * uu - rh
  subroutine itsol_defect(aa, rr, rh, uu, mm, uniform_dh)
    use bl_prof_module
    type(multifab), intent(inout) :: uu, rr
    type(multifab), intent(in) :: rh, aa
    type(imultifab), intent(in) :: mm
    logical, intent(in), optional :: uniform_dh
    type(bl_prof_timer), save :: bpt
    call build(bpt, "its_defect")
    call itsol_stencil_apply(aa, rr, uu, mm, uniform_dh)
    call saxpy(rr, rh, -1.0_dp_t, rr)
    call destroy(bpt)
  end subroutine itsol_defect

  subroutine itsol_BiCGStab_solve(aa, uu, rh, mm, eps, max_iter, verbose, &
       stat, singular_in, uniform_dh, nodal_mask)
    use bl_prof_module
    integer, intent(in) :: max_iter
    type(imultifab), intent(in) :: mm
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: rh
    type(multifab), intent(in) :: aa
    integer, intent(out), optional :: stat
    logical, intent(in), optional :: singular_in
    logical, intent(in), optional :: uniform_dh
    type(multifab), intent(in), optional :: nodal_mask
    real(kind=dp_t) :: eps
    type(layout) :: la
    integer, intent(in) :: verbose
    type(multifab) :: rr, rt, pp, ph, vv, tt, ss, sh
    real(kind=dp_t) :: rho_1, alpha, beta, omega, rho, Anorm, bnorm, rnorm, den
    real(dp_t) :: rho_hg, rho_orig, volume, tres0, small, norm_rr, norm_uu
    integer :: i
    integer :: cnt, ng_for_res
    logical :: nodal_solve
    logical :: singular
    type(bl_prof_timer), save :: bpt

    type(multifab) :: rh_local, aa_local

    call build(bpt, "its_BiCGStab_solve")

    if ( present(stat) ) stat = 0

    singular = .false.
    if ( present(singular_in) ) singular = singular_in

    ng_for_res = 0; if ( nodal_q(rh) ) ng_for_res = 1
    nodal_solve = .False.; if ( ng_for_res /= 0 ) nodal_solve = .TRUE.

    la = aa%la
    call multifab_build(rr, la, 1, ng_for_res, rh%nodal)
    call multifab_build(rt, la, 1, ng_for_res, rh%nodal)
    call multifab_build(pp, la, 1, ng_for_res, rh%nodal)
    call multifab_build(ph, la, 1, uu%ng     , rh%nodal)
    call multifab_build(vv, la, 1, ng_for_res, rh%nodal)
    call multifab_build(tt, la, 1, ng_for_res, rh%nodal)
    call multifab_build(sh, la, 1, uu%ng     , rh%nodal)
    call multifab_build(ss, la, 1, ng_for_res, rh%nodal)

    if ( nodal_solve ) then
       call setval(rr, ZERO, all=.true.)
       call setval(rt, ZERO, all=.true.)
       call setval(pp, ZERO, all=.true.)
       call setval(vv, ZERO, all=.true.)
       call setval(tt, ZERO, all=.true.)
       call setval(ss, ZERO, all=.true.)
    end if

    ! Use these for local preconditioning
    call multifab_build(rh_local, la, rh%nc, rh%ng, rh%nodal)
    call multifab_build(aa_local, la, aa%nc, aa%ng, aa%nodal)

    call copy(rh_local, 1, rh, 1, nc = rh%nc, ng = rh%ng)
    call copy(aa_local, 1, aa, 1, nc = aa%nc, ng = aa%ng)

    call diag_initialize(aa_local,rh_local,mm)

    call copy(ph, uu, ng = ph%ng)
    call copy(sh, uu, ng = sh%ng)

    cnt = 0
    ! compute rr = aa * uu - rh
    call itsol_defect(aa_local, rr, rh_local, uu, mm, uniform_dh); cnt = cnt + 1

    if (singular .and. nodal_solve) then
      call setval(ss,ONE)
      if (present(nodal_mask)) then
            rho = dot(rr, ss, nodal_mask)
         volume = dot(ss, ss, nodal_mask)
      else
            rho = dot(rr, ss)
         volume = dot(ss,ss)
      end if
!     print *,'SINGULAR ADJUSTMENT ',rho,' OVER ',volume 
      rho = rho / volume
      call saxpy(rr,-rho,ss)
      call setval(ss,ZERO,all=.true.)
    end if

    call copy(rt, rr)
    if (present(nodal_mask)) then
       rho = dot(rt, rr, nodal_mask)
    else 
       rho = dot(rt, rr)
    end if
    rho_orig = rho

    tres0 = norm_inf(rr)
    Anorm = stencil_norm(aa_local)
    small = epsilon(Anorm)
    bnorm = norm_inf(rh_local)

    if ( parallel_IOProcessor() .and. verbose > 0) then
       write(unit=*, fmt='("    BiCGStab: Initial error (error0) =        ",g15.8)') tres0
    end if

    i = 0
    if ( itsol_converged(rr, uu, Anorm, bnorm, eps) ) then
      if ( verbose > 0 ) then
         if ( tres0 < eps*bnorm ) then
            if ( parallel_IOProcessor() ) then
               write(unit=*, fmt='("    BiCGStab: Zero iterations: rnorm ",g15.8," < eps*bnorm ",g15.8)') &
                    tres0,eps*bnorm
            end if
        else
           norm_rr = norm_inf(rr)
           if ( norm_rr < epsilon(Anorm)*Anorm ) then
              if ( parallel_IOProcessor() ) then
                 write(unit=*, fmt='("    BiCGStab: Zero iterations: rnorm ",g15.8," < small*Anorm ",g15.8)') &
                      tres0,small*Anorm
              end if
           end if
        end if
     end if
     go to 100
    end if

    do i = 1, max_iter
       if (present(nodal_mask)) then
          rho = dot(rt, rr, nodal_mask)
       else 
          rho = dot(rt, rr)
       end if
       if ( i == 1 ) then
          call copy(pp, rr)
       else
          if ( rho_1 == ZERO ) then
             if ( present(stat) ) then
                call bl_warn("BiCGStab_SOLVE: failure 1")
                stat = 2
                goto 100
             end if
             call bl_error("BiCGStab: failure 1")
          end if
          if ( omega == ZERO ) then
             if ( present(stat) ) then
                call bl_warn("BiCGStab_SOLVE: failure 2")
                stat = 3
                goto 100
             end if
             call bl_error("BiCGStab: failure 2")
          end if
          beta = (rho/rho_1)*(alpha/omega)
          call saxpy(pp, -omega, vv)
          call saxpy(pp, rr, beta, pp)
       end if
       call itsol_precon(aa_local, ph, pp, mm, 0)
       call itsol_stencil_apply(aa_local, vv, ph, mm, uniform_dh); cnt = cnt + 1
       if (present(nodal_mask)) then
          den = dot(rt, vv, nodal_mask)
       else
          den = dot(rt, vv)
       end if 
       if ( den == ZERO ) then
          if ( present(stat) ) then
             call bl_warn("BICGSTAB_solve: breakdown in bicg, going with what I have")
             stat = 30
             goto 100
          endif
          call bl_error("BiCGStab: failure 3")
       end if
       alpha = rho/den
       call saxpy(uu, alpha, ph)
       call saxpy(ss, rr, -alpha, vv)
       rnorm = norm_inf(ss)
       if ( parallel_IOProcessor() .and. verbose > 1 ) then
          write(unit=*, fmt='("    BiCGStab: Half Iter        ",i4," rel. err. ",g15.8)') cnt/2, &
                             rnorm  /  (bnorm)
       end if
       if ( itsol_converged(ss, uu, Anorm, bnorm, eps) ) exit
       call itsol_precon(aa_local, sh, ss, mm,0)
       call itsol_stencil_apply(aa_local, tt, sh, mm, uniform_dh); cnt = cnt + 1

       if (present(nodal_mask)) then
          den = dot(tt, tt, nodal_mask)
       else
          den = dot(tt, tt)
       end if
       if ( den == ZERO ) then
          if ( present(stat) ) then
             call bl_warn("BICGSTAB_solve: breakdown in bicg, going with what I have")
             stat = 31
             goto 100
          endif
          call bl_error("BiCGStab: failure 3")
       end if
       if (present(nodal_mask)) then
          omega = dot(tt,ss,nodal_mask)/den
       else
          omega = dot(tt,ss)/den
       end if
       call saxpy(uu, omega, sh)
       call saxpy(rr, ss, -omega, tt)
       rnorm = norm_inf(rr)
       if ( parallel_IOProcessor() .and. verbose > 1) then
          write(unit=*, fmt='("    BiCGStab: Iteration        ",i4," rel. err. ",g15.8)') cnt/2, &
                             rnorm /  (bnorm)
       end if
       if ( .false. .and. nodal_solve ) then
          ! HACK, THIS IS USED TO MATCH THE HGPROJ STOPPING CRITERION
          call itsol_precon(aa_local, sh, rr, mm, 0)
          if (present(nodal_mask)) then
             rho_hg = dot(rr, sh, nodal_mask)
          else
             rho_hg = dot(rr, sh)
          end if
          if ( (abs(rho_hg) < rho_orig*eps) .or. &
              itsol_converged(rr, uu, Anorm, bnorm, eps) ) exit
       else
          if ( itsol_converged(rr, uu, Anorm, bnorm, eps) ) exit
       end if
       rho_1 = rho

    end do

    if ( verbose > 0 ) then
       if ( parallel_IOProcessor() ) then
          write(unit=*, fmt='("    BiCGStab: Final: Iteration  ", i3, " rel. err. ",g15.8)') cnt/2, &
               rnorm/ (bnorm)
       end if
       if ( rnorm < eps*bnorm ) then
          if ( parallel_IOProcessor() ) then
             write(unit=*, fmt='("    BiCGStab: Converged: rnorm ",g15.8," < eps*bnorm ",g15.8)') &
                  rnorm,eps*bnorm
          end if
       else
          norm_uu = norm_inf(uu)
          if ( rnorm < eps*Anorm*norm_uu ) then
             if ( parallel_IOProcessor() ) then
                write(unit=*, fmt='("    BiCGStab: Converged: rnorm ",g15.8," < eps*Anorm*sol_norm ",g15.8)') &
                     rnorm,eps*Anorm*norm_uu
             end if
          else if ( rnorm < epsilon(Anorm)*Anorm ) then
             if ( parallel_IOProcessor() ) then
                write(unit=*, fmt='("    BiCGStab: Converged: rnorm ",g15.8," < small*Anorm ",g15.8)') &
                     rnorm,small*Anorm
             end if
          end if
       end if
    end if

!    if (rnorm > bnorm) then
!       call setval(uu,ZERO,all=.true.)
!       if ( present(stat) ) stat = 1
!       if ( verbose > 0 .and.  parallel_IOProcessor() ) &
!          print *,'   BiCGStab: solution reset to zero'
!    end if

    if ( i > max_iter ) then
       if ( present(stat) ) then
          stat = 1
       else
          call bl_error("BiCGSolve: failed to converge");
       end if
    end if

    call destroy(rh_local)
    call destroy(aa_local)

100 continue

    call destroy(rr)
    call destroy(rt)
    call destroy(pp)
    call destroy(ph)
    call destroy(vv)
    call destroy(tt)
    call destroy(sh)
    call destroy(ss)

    call destroy(bpt)

  end subroutine itsol_BiCGStab_solve

  subroutine itsol_CG_Solve(aa, uu, rh, mm, eps, max_iter, verbose, &
                            stat, singular_in, uniform_dh, nodal_mask)
    use bl_prof_module
    integer, intent(in   )           :: max_iter, verbose
    integer, intent(  out), optional :: stat
    logical, intent(in   ), optional :: singular_in
    logical, intent(in   ), optional :: uniform_dh
    type(multifab), intent(in), optional :: nodal_mask

    type( multifab), intent(in) :: aa
    type( multifab), intent(inout) :: uu
    type( multifab), intent(in) :: rh
    type(imultifab), intent(in) :: mm

    real(dp_t), intent(in) :: eps
    type(multifab) :: rr, zz, pp, qq
    type(multifab) :: aa_local, rh_local
    real(kind = dp_t) :: rho_1, alpha, beta, Anorm, bnorm, rho, rnorm, den, tres0, small
    type(layout) :: la
    integer :: i, ng_for_res
    logical :: nodal_solve
    logical :: singular 
    integer :: cnt

    real(dp_t) :: rho_hg, rho_orig, volume, rho_hg_orig, norm_uu
    type(bl_prof_timer), save :: bpt

    call build(bpt, "its_CG_Solve")

    if ( present(stat) ) stat = 0

    singular = .false.
    if ( present(singular_in) ) singular = singular_in

    ng_for_res = 0; if ( nodal_q(rh) ) ng_for_res = 1
    nodal_solve = .FALSE.; if ( ng_for_res /= 0 ) nodal_solve = .TRUE.

    la = aa%la
    call multifab_build(rr, la, 1, ng_for_res, rh%nodal)
    call multifab_build(zz, la, 1, ng_for_res, rh%nodal)
    call multifab_build(pp, la, 1, uu%ng     , rh%nodal)
    call multifab_build(qq, la, 1, ng_for_res, rh%nodal)

    if ( nodal_solve ) then
       call setval(rr,ZERO,all=.true.)
       call setval(zz,ZERO,all=.true.)
       call setval(qq,ZERO,all=.true.)
    end if
    call setval(pp, ZERO, all=.true.)

    ! Use these for local preconditioning
    call multifab_build(rh_local, la, rh%nc, rh%ng, rh%nodal)
    call multifab_build(aa_local, la, aa%nc, aa%ng, aa%nodal)

    call copy(rh_local, 1, rh, 1, nc = rh%nc, ng = rh%ng)
    call copy(aa_local, 1, aa, 1, nc = aa%nc, ng = aa%ng)

    call diag_initialize(aa_local,rh_local,mm)

    cnt = 0
    ! compute rr = aa * uu - rh_local
    call itsol_defect(aa_local, rr, rh_local, uu, mm, uniform_dh)  
    cnt = cnt + 1

    if (singular .and. nodal_solve) then
      call setval(zz,ONE)
      if (present(nodal_mask)) then
            rho = dot(rr, zz, nodal_mask)
         volume = dot(zz,zz)
      else
            rho = dot(rr, zz)
         volume = dot(zz,zz)
      end if
!     print *,'SINGULAR ADJUSTMENT ',rho,' OVER ',volume 
      rho = rho / volume
      call saxpy(rr,-rho,zz)
      call setval(zz,ZERO,all=.true.)
    end if

    Anorm = stencil_norm(aa_local)
    small = epsilon(Anorm)
    bnorm = norm_inf(rh_local)
    tres0 = norm_inf(rr)

    if ( parallel_IOProcessor() .and. verbose > 0) then
       write(unit=*, fmt='("          CG: Initial error (error0) =        ",g15.8)') tres0
    end if

    i = 0
    if ( itsol_converged(rr, uu, Anorm, bnorm, eps) ) then
      if (parallel_IOProcessor() .and. verbose > 0) then
        if (tres0 < eps*bnorm) then
          write(unit=*, fmt='("          CG: Zero iterations: rnorm ",g15.8," < eps*bnorm ",g15.8)') &
                tres0,eps*bnorm
        else if (tres0 < epsilon(Anorm)*Anorm) then
          write(unit=*, fmt='("          CG: Zero iterations: rnorm ",g15.8," < small*Anorm ",g15.8)') &
                tres0,small*Anorm
        end if
      end if
      go to 100
    end if

    do i = 1, max_iter

       call itsol_precon(aa_local, zz, rr, mm, 0)
       if (present(nodal_mask)) then
          rho = dot(rr, zz, nodal_mask)
       else
          rho = dot(rr, zz)
       end if
       if ( i == 1 ) then
          call copy(pp, zz)
          rho_orig = rho

          call itsol_precon(aa_local, zz, rr, mm)
          if (present(nodal_mask)) then
             rho_hg_orig = dot(rr, zz, nodal_mask)
          else
             rho_hg_orig = dot(rr, zz)
          end if

       else
          if ( rho_1 == ZERO ) then
             if ( present(stat) ) then
                call bl_warn("CG_solve: failure 1")
                stat = 1
                goto 100
             end if
             call bl_error("CG_solve: failure 1")
          end if
          beta = rho/rho_1
          call saxpy(pp, zz, beta, pp)
       end if
       call itsol_stencil_apply(aa_local, qq, pp, mm, uniform_dh); cnt = cnt + 1
       if (present(nodal_mask)) then
          den = dot(pp, qq, nodal_mask)
       else
          den = dot(pp, qq)
       end if
       if ( den == ZERO ) then
          if ( present(stat) ) then
             call bl_warn("CG_solve: breakdown in solver, going with what I have")
             stat = 30
             goto 100
          end if
          call bl_error("CG_solve: failure 1")
       end if
       alpha = rho/den
       call saxpy(uu,   alpha, pp)
       call saxpy(rr, - alpha, qq)
       rnorm = norm_inf(rr)
       if ( parallel_IOProcessor() .and. verbose > 1) then
          write(unit=*, fmt='("          CG: Iteration        ",i4," rel. err. ",g15.8)') i, &
                             rnorm /  (bnorm)
       end if
       if ( .true. .and. nodal_solve ) then
          ! HACK, THIS IS USED TO MATCH THE HGPROJ STOPPING CRITERION
          call itsol_precon(aa_local, zz, rr, mm)
          if (present(nodal_mask)) then
             rho_hg = dot(rr, zz, nodal_mask)
          else
             rho_hg = dot(rr, zz)
          end if
          if ( (abs(rho_hg) < abs(rho_hg_orig)*eps) ) then
            exit
          end if
       else
          if ( itsol_converged(rr, uu, Anorm, bnorm, eps) ) exit
       end if
       rho_1 = rho

    end do

    if ( verbose > 0 ) then
       if ( parallel_IOProcessor() ) then
          write(unit=*, fmt='("          CG: Final: Iteration  ", i3, " rel. err. ",g15.8)') i, &
               rnorm/ (bnorm)
       end if
       if ( rnorm < eps*bnorm ) then
          if ( parallel_IOProcessor() ) then
             write(unit=*, fmt='("          CG: Converged: rnorm ",g15.8," < eps*bnorm ",g15.8)') &
                  rnorm,eps*bnorm
          end if
       else
          norm_uu = norm_inf(uu)
          if ( rnorm < eps*Anorm*norm_uu ) then
             if ( parallel_IOProcessor() ) then
                write(unit=*, fmt='("          CG: Converged: rnorm ",g15.8," < eps*Anorm*sol_norm ",g15.8)') &
                     rnorm,eps*Anorm*norm_uu
             end if
          else if ( rnorm < epsilon(Anorm)*Anorm ) then
             if ( parallel_IOProcessor() ) then
                write(unit=*, fmt='("          CG: Converged: rnorm ",g15.8," < small*Anorm ",g15.8)') &
                     rnorm,small*Anorm
             end if
          end if
       end if
    end if

!    if (rnorm > bnorm) call setval(uu,ZERO,all=.true.)

    if ( i > max_iter ) then
       if ( present(stat) ) then
          stat = 1
       else
          call bl_error("CG_solve: failed to converge");
       end if
    end if

100 continue

    call destroy(rr)
    call destroy(zz)
    call destroy(pp)
    call destroy(qq)

    call destroy(bpt)

    call destroy(aa_local)
    call destroy(rh_local)

  end subroutine itsol_CG_Solve

  subroutine itsol_precon(aa, uu, rh, mm, method)
    use bl_prof_module
    type(multifab), intent(in) :: aa
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: rh
    type(imultifab), intent(in) :: mm
    real(kind=dp_t), pointer, dimension(:,:,:,:) :: ap, up, rp
    integer, pointer, dimension(:,:,:,:) :: mp
    integer :: i, n
    integer, intent(in), optional :: method
    integer :: lm
    type(bl_prof_timer), save :: bpt

    call build(bpt, "its_precon")

    lm = 1; if ( present(method) ) lm = method

    select case (lm)
    case (0)
       call copy(uu, rh)
    case (1)
       do i = 1, rh%nboxes
          if ( multifab_remote(uu, i) ) cycle
          rp => dataptr(rh, i)
          up => dataptr(uu, i)
          ap => dataptr(aa, i)
          mp => dataptr(mm, i)
          do n = 1, uu%nc
             select case(uu%dim)
             case (1)
                if ( cell_centered_q(rh) ) then
                   call jacobi_precon_1d(ap(:,1,1,:), up(:,1,1,n), rp(:,1,1,n), uu%ng)
                else
                   call nodal_precon_1d(ap(:,1,1,:), up(:,1,1,n), rp(:,1,1,n), &
                                        mp(:,1,1,1),uu%ng)
                end if
             case (2)
                if ( cell_centered_q(rh) ) then
                   call jacobi_precon_2d(ap(:,:,1,:), up(:,:,1,n), rp(:,:,1,n), uu%ng)
                else
                   call nodal_precon_2d(ap(:,:,1,:), up(:,:,1,n), rp(:,:,1,n), &
                                        mp(:,:,1,1),uu%ng)
                end if
             case (3)
                if ( cell_centered_q(rh) ) then
                   call jacobi_precon_3d(ap(:,:,:,:), up(:,:,:,n), rp(:,:,:,n), uu%ng)
                else
                   call nodal_precon_3d(ap(:,:,:,:), up(:,:,:,n), rp(:,:,:,n), &
                                        mp(:,:,:,1),uu%ng)
                end if
             end select
          end do
       end do
    end select

    call destroy(bpt)

  end subroutine itsol_precon

  subroutine diag_initialize(aa, rh, mm)
    use bl_prof_module
    type( multifab), intent(in) :: aa
    type( multifab), intent(in) :: rh
    type(imultifab), intent(in) :: mm

    real(kind=dp_t), pointer, dimension(:,:,:,:) :: ap, rp
    integer        , pointer, dimension(:,:,:,:) :: mp
    integer                                      :: i
    integer                                      :: ng_a, ng_r, ng_m
    integer                                      :: lo(rh%dim),hi(rh%dim)
    type(bl_prof_timer), save                    :: bpt

    call build(bpt, "diag_initialize")

    ng_a = aa%ng
    ng_r = rh%ng
    ng_m = mm%ng

    do i = 1, rh%nboxes
       if ( multifab_remote(rh, i) ) cycle
       rp => dataptr(rh, i)
       ap => dataptr(aa, i)
       mp => dataptr(mm, i)
       lo = lwb(get_box(rh,i))
       hi = upb(get_box(rh,i))
       select case(rh%dim)
          case (1)
             if ( cell_centered_q(rh) ) then
                call diag_init_cc_1d(ap(:,1,1,:), ng_a, rp(:,1,1,1), ng_r, lo, hi)
             else
                call diag_init_nd_1d(ap(:,1,1,:), ng_a, rp(:,1,1,1), ng_r, mp(:,1,1,1), ng_m, lo, hi)
             end if
          case (2)
             if ( cell_centered_q(rh) ) then
                call diag_init_cc_2d(ap(:,:,1,:), ng_a, rp(:,:,1,1), ng_r, lo, hi)
             else
                call diag_init_nd_2d(ap(:,:,1,:), ng_a, rp(:,:,1,1), ng_r, mp(:,:,1,1), ng_m, lo, hi)
             end if
          case (3)
             if ( cell_centered_q(rh) ) then
                call diag_init_cc_3d(ap(:,:,:,:), ng_a, rp(:,:,:,1), ng_r, lo, hi)
             else
                call diag_init_nd_3d(ap(:,:,:,:), ng_a, rp(:,:,:,1), ng_r, mp(:,:,:,1), ng_m, lo, hi)
             end if
       end select
    end do

    call destroy(bpt)

  end subroutine diag_initialize

end module itsol_module

