module itsol_module

  use bl_types
  use multifab_module
  use stencil_module
  use stencil_nodal_module

  implicit none

  integer, private, parameter :: def_bicg_max_iter = 1000
  integer, private, parameter :: def_cg_max_iter = 1000
  real(dp_t), private, parameter :: ZERO = 0.0_dp_t
  real(dp_t), private, parameter :: ONE  = 1.0_dp_t

  private :: itsol_defect
  private :: itsol_precon
  private :: itsol_precon_st

contains

  !! ITSOL_BREAKDOWN: is supposed to detect 'bad' numbers
  !! that shouldn't be used in a CG/BICG
  !! Not Used
  function itsol_breakdown(alpha, how) result(r)
    logical :: r
    real(dp_t), intent(in) :: alpha
    integer, intent(out) :: how
    real(dp_t) :: beta

    r = .false.
    how = 0
    if ( alpha > ZERO ) then
       beta = ONE/alpha
       if ( beta <= ZERO ) then
          r = .true.
          how = 2
       end if
    else if ( alpha < ZERO ) then
       beta = ONE/alpha
       if ( beta >= ZERO ) then
          r = .true.
          how = 2
       end if
    else if ( alpha == ZERO ) then
       r = .true.
       how = 1
    else
       r = .true.
       how = 2
    end if

  end function itsol_breakdown

  function itsol_converged(rr, uu, Anorm, bnorm, eps) result(r)
    type(multifab), intent(in) :: rr, uu
    real(dp_t), intent(in) :: Anorm, bnorm, eps
    real(dp_t)             :: norm_rr, norm_uu
    logical :: r
    norm_rr = norm_inf(rr)
    norm_uu = norm_inf(uu)
    r = (norm_rr <= eps*(Anorm*norm_uu + bnorm)) .or. &
        (norm_rr <= epsilon(Anorm)*Anorm)
  end function itsol_converged

  subroutine itsol_stencil_apply(aa, rr, uu, mm)
    type(multifab), intent(in) :: aa
    type(multifab), intent(inout) :: rr
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in) :: mm
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: ap(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: i, n
    logical :: nodal_flag

    call multifab_fill_boundary(uu)

    nodal_flag = nodal_q(uu)

    do i = 1, rr%nboxes
       if ( multifab_remote(rr, i) ) cycle
       rp => dataptr(rr, i)
       up => dataptr(uu, i)
       ap => dataptr(aa, i)
       mp => dataptr(mm, i)
       do n = 1, rr%nc
          select case(rr%dim)
          case (1)
             if ( .not. nodal_flag) then
                call stencil_apply_1d(ap(:,1,1,:), rp(:,1,1,n), up(:,1,1,n),  &
                     mp(:,1,1,1), uu%ng)
             else
                call stencil_apply_1d_nodal(ap(:,1,1,:), rp(:,1,1,n), up(:,1,1,n),  &
                     mp(:,1,1,1), uu%ng)
             end if
          case (2)
             if ( .not. nodal_flag) then
                call stencil_apply_2d(ap(:,:,1,:), rp(:,:,1,n), up(:,:,1,n),  &
                     mp(:,:,1,1), uu%ng)
             else
                call stencil_apply_2d_nodal(ap(:,:,1,:), rp(:,:,1,n), up(:,:,1,n),  &
                     mp(:,:,1,1), uu%ng)
             end if
          case (3)
             if ( .not. nodal_flag) then
                call stencil_apply_3d(ap(:,:,:,:), rp(:,:,:,n), up(:,:,:,n),  &
                     mp(:,:,:,1), uu%ng)
             else
                call stencil_apply_3d_nodal(ap(:,:,:,:), rp(:,:,:,n), up(:,:,:,n),  &
                     mp(:,:,:,1), uu%ng)
             end if
          end select
       end do
    end do

  end subroutine itsol_stencil_apply

  subroutine itsol_defect(aa, rr, rh, uu, mm)
    type(multifab), intent(inout) :: uu, rr
    type(multifab), intent(in) :: rh, aa
    type(imultifab), intent(in) :: mm
    call itsol_stencil_apply(aa, rr, uu, mm)
    call saxpy(rr, rh, -1.0_dp_t, rr)
  end subroutine itsol_defect

  subroutine itsol_BiCGStab_solve(aa, uu, rh, mm, eps, max_iter, verbose, stat)
    integer, intent(in) :: max_iter
    type(imultifab), intent(in) :: mm
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: rh
    type(multifab), intent(in) :: aa
    integer, intent(out), optional :: stat
    real(kind=dp_t) :: eps
    type(layout) :: la
    integer, intent(in) :: verbose
    type(multifab) :: rr, rt, pp, ph, vv, tt, ss, sh
    real(kind=dp_t) :: rho_1, alpha, beta, omega, rho, Anorm, bnorm, nrm, den
    real(dp_t) :: rho_hg, rho_orig
    integer :: i
    integer :: cnt, ng_for_res
    logical :: nodal_solve

    if ( present(stat) ) stat = 0

    ng_for_res = 0; if ( nodal_q(rh) ) ng_for_res = 1
    nodal_solve = .False.; if ( ng_for_res /= 0 ) nodal_solve = .TRUE.

    la = aa%la
    call multifab_build(rr, la, 1, ng_for_res, rh%nodal)
    call multifab_build(rt, la, 1, ng_for_res, rh%nodal)
    call multifab_build(pp, la, 1, ng_for_res, rh%nodal)
    call multifab_build(ph, la, 1, 1,          rh%nodal)
    call multifab_build(vv, la, 1, ng_for_res, rh%nodal)
    call multifab_build(tt, la, 1, ng_for_res, rh%nodal)
    call multifab_build(sh, la, 1, 1,          rh%nodal)
    call multifab_build(ss, la, 1, ng_for_res, rh%nodal)

    if ( nodal_solve ) then
       call setval(rr, ZERO, all=.true.)
       call setval(rt, ZERO, all=.true.)
       call setval(pp, ZERO, all=.true.)
       call setval(vv, ZERO, all=.true.)
       call setval(tt, ZERO, all=.true.)
       call setval(ss, ZERO, all=.true.)
    end if

    call copy(ph, uu, all = .true.)
    call copy(sh, uu, all = .true.)

    cnt = 0
    call itsol_defect(aa, rr, rh, uu, mm); cnt = cnt + 1
    call copy(rt, rr)
    rho = dot(rt, rr)
    rho_orig = rho

    Anorm = stencil_norm(aa)
    bnorm = norm_inf(rh)
    if ( verbose > 1 ) then
       nrm = norm_inf(rr)
       if ( parallel_IOProcessor() ) then
          write(unit=*, fmt='(i3,": Anorm=",g15.8,", bnorm=",g15.8,", Rnorm=",g15.8)') cnt, Anorm, bnorm, nrm
       end if
    end if

    i = 0
    if ( itsol_converged(rr, uu, Anorm, bnorm, eps) ) goto 100

    do i = 1, max_iter
       rho = dot(rt, rr)
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
       call itsol_precon(aa, ph, pp, mm)
       call itsol_stencil_apply(aa, vv, ph, mm); cnt = cnt + 1
       den = dot(rt, vv)
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
       if ( verbose > 1 ) then
          nrm = norm_inf(ss)
          if ( parallel_IOProcessor() ) then
             write(unit=*, fmt='(i3,": Snorm=",g15.8)') cnt, nrm
          end if
       end if
       if ( itsol_converged(ss, uu, Anorm, bnorm, eps) ) exit
       call itsol_precon(aa, sh, ss, mm)
       call itsol_stencil_apply(aa, tt, sh, mm); cnt = cnt + 1

       den = dot(tt,tt)
       if ( den == ZERO ) then
          if ( present(stat) ) then
             call bl_warn("BICGSTAB_solve: breakdown in bicg, going with what I have")
             stat = 31
             goto 100
          endif
          call bl_error("BiCGStab: failure 3")
       end if
       omega = dot(tt,ss)/den
       call saxpy(uu, omega, sh)
       call saxpy(rr, ss, -omega, tt)
       if ( verbose > 1 ) then
          nrm = norm_inf(rr)
          if ( parallel_IOProcessor() ) then
             write(unit=*, fmt='(i3,": Rnorm=",g15.8)') cnt, nrm
          end if
       end if
       if ( .true. .and. nodal_solve ) then
          ! HACK, THIS IS USED TO MATCH THE HGPROJ STOPPING CRITERION
          call itsol_precon(aa, sh, rr, mm)
          rho_hg = dot(rr, sh)
          if ( (rho_hg < rho_orig*eps) .or. &
              itsol_converged(rr, uu, Anorm, bnorm, eps) ) exit
       else
          if ( itsol_converged(rr, uu, Anorm, bnorm, eps) ) exit
       end if
       rho_1 = rho
    end do

    if ( verbose > 0 ) then
       if ( parallel_IOProcessor() ) then
          write(unit=*, fmt='("BICGStab: iterations: ", i5)') cnt
       end if
    end if

100 continue
    call destroy(rr)
    call destroy(rt)
    call destroy(pp)
    call destroy(ph)
    call destroy(vv)
    call destroy(tt)
    call destroy(sh)
    call destroy(ss)

    if ( i > max_iter ) then
       if ( present(stat) ) then
          stat = 1
       else
          call bl_error("BiCGSolve: failed to converge");
       end if
    end if

  end subroutine itsol_BiCGStab_solve

  subroutine itsol_CG_Solve(aa, uu, rh, mm, eps, max_iter, verbose, stat)
    integer, intent(in) :: max_iter, verbose
    integer, intent(out), optional :: stat

    type( multifab), intent(in) :: aa
    type( multifab), intent(inout) :: uu
    type( multifab), intent(in) :: rh
    type(imultifab), intent(in) :: mm

    real(dp_t), intent(in) :: eps
    type(multifab) :: rr, zz, pp, qq
    real(kind = dp_t) :: rho_1, alpha, beta, Anorm, bnorm, rho, nrm, den
    type(layout) :: la
    integer :: i, ng_for_res
    logical :: nodal_solve
    integer :: cnt

    real(dp_t) :: rho_hg, rho_orig

    if ( present(stat) ) stat = 0

    ng_for_res = 0; if ( nodal_q(rh) ) ng_for_res = 1
    nodal_solve = .FALSE.; if ( ng_for_res /= 0 ) nodal_solve = .TRUE.

    la = aa%la
    call multifab_build(rr, la, 1, ng_for_res, rh%nodal)
    call multifab_build(zz, la, 1, ng_for_res, rh%nodal)
    call multifab_build(pp, la, 1, 1         , rh%nodal)
    call multifab_build(qq, la, 1, ng_for_res, rh%nodal)

    if ( nodal_solve ) then
       call setval(rr,ZERO,all=.true.)
       call setval(zz,ZERO,all=.true.)
       call setval(qq,ZERO,all=.true.)
    end if
    call setval(pp, ZERO, all=.true.)

    cnt = 0
    call itsol_defect(aa, rr, rh, uu, mm); cnt = cnt + 1

    Anorm = stencil_norm(aa)
    bnorm = norm_inf(rh)

    if ( verbose > 1 ) then
       nrm = norm_inf(rr)
       if ( parallel_IOProcessor() ) then
          write(unit=*, fmt='(i3,": Anorm=",g15.8,", bnorm=",g15.8,", Rnorm=",g15.8)') 0, Anorm, bnorm, nrm
       end if
    end if

    i = 0
    if ( itsol_converged(rr, uu, Anorm, bnorm, eps) ) goto 100

    do i = 1, max_iter
       call itsol_precon(aa, zz, rr, mm)
       rho = dot(rr, zz)
       if ( i == 1 ) then
          call copy(pp, zz)
          rho_orig = rho
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
       call itsol_stencil_apply(aa, qq, pp, mm); cnt = cnt + 1
       den = dot(pp, qq)
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
       if ( verbose > 1 ) then
          nrm = norm_inf(rr)
          if ( parallel_IOProcessor() ) then
             write(unit=*, fmt='(i3,": Rnorm=",g15.8)') i, nrm
          end if
       end if
       if ( .true. .and. nodal_solve ) then
          ! HACK, THIS IS USED TO MATCH THE HGPROJ STOPPING CRITERION
          call itsol_precon(aa, zz, rr, mm)
          rho_hg = dot(rr, zz)
          if ( (rho_hg < rho_orig*eps) .or. &
              itsol_converged(rr, uu, Anorm, bnorm, eps) ) exit
       else
          if ( itsol_converged(rr, uu, Anorm, bnorm, eps) ) exit
       end if
       rho_1 = rho
    end do

    if ( verbose > 0 ) then
       if ( parallel_IOProcessor() ) then
          write(unit=*, fmt='("CG: iterations: ", i5)') cnt
       end if
    end if

100 continue
    call destroy(rr)
    call destroy(zz)
    call destroy(pp)
    call destroy(qq)

    if ( i > max_iter ) then
       if ( present(stat) ) then
          stat = 1
       else
          call bl_error("CG_solve: failed to converge");
       end if
    end if

  end subroutine itsol_CG_Solve

  subroutine itsol_precon(aa, uu, rh, mm, method)
    type(multifab), intent(in) :: aa
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: rh
    type(imultifab), intent(in) :: mm
    real(kind=dp_t), pointer, dimension(:,:,:,:) :: ap, up, rp
    integer, pointer, dimension(:,:,:,:) :: mp
    integer :: i, n
    integer, intent(in), optional :: method
    integer :: lm
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
      !$OMP PARALLEL DO PRIVATE(j,i)
      do j = 1, ny
         do i = 1, nx
            u(i,j) = r(i,j)/a(i,j,0)
         end do
      end do
      !$OMP END PARALLEL DO
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
      !$OMP PARALLEL DO PRIVATE(j,i,k)
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
      !$OMP PARALLEL DO PRIVATE(j,i)
      do j = 1, ny
         do i = 1, nx
            if (.not. bc_dirichlet(mm(i,j),1,0)) then
               u(i,j) = r(i,j)/a(i,j,0)
            end if
         end do
      end do
      !$OMP END PARALLEL DO
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
      !$OMP PARALLEL DO PRIVATE(j,i,k)
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
  end subroutine itsol_precon

  subroutine itsol_BiCGStab_solve_st(st, uu, rh, eps, max_iter, verbose, stat)
    integer, intent(in) :: max_iter
    type(stencil), intent(in) :: st
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: rh
    integer, intent(out), optional :: stat
    real(kind=dp_t) :: eps
    type(layout) :: la
    integer, intent(in) :: verbose
    type(multifab) :: rr, rt, pp, ph, vv, tt, ss, sh
    real(kind=dp_t) :: rho_1, alpha, beta, omega, rho, Anorm, bnorm, nrm
    integer :: i
    integer :: cnt, ng_for_res
    logical :: nodal_solve

    if ( present(stat) ) stat = 0

    ng_for_res = 0; if ( nodal_q(rh) ) ng_for_res = 1
    nodal_solve = .False.; if ( ng_for_res /= 0 ) nodal_solve = .TRUE.

    la = st%ss%la
    call multifab_build(rr, la, 1, ng_for_res, rh%nodal)
    call multifab_build(rt, la, 1, ng_for_res, rh%nodal)
    call multifab_build(pp, la, 1, ng_for_res, rh%nodal)
    call multifab_build(ph, la, 1, 1,          rh%nodal)
    call multifab_build(vv, la, 1, ng_for_res, rh%nodal)
    call multifab_build(tt, la, 1, ng_for_res, rh%nodal)
    call multifab_build(sh, la, 1, 1,          rh%nodal)
    call multifab_build(ss, la, 1, ng_for_res, rh%nodal)

    if ( nodal_solve ) then
       call setval(rr, ZERO, all=.true.)
       call setval(rt, ZERO, all=.true.)
       call setval(pp, ZERO, all=.true.)
       call setval(vv, ZERO, all=.true.)
       call setval(tt, ZERO, all=.true.)
       call setval(ss, ZERO, all=.true.)
    end if

    call copy(ph, uu, all = .true.)
    call copy(sh, uu, all = .true.)

    cnt = 0
    call stencil_defect_st(st, rr, rh, uu); cnt = cnt + 1
    call copy(rt, rr)

    Anorm = stencil_norm_st(st)
    bnorm = norm_inf(rh)
    if ( verbose > 1 ) then
       nrm = norm_inf(rr)
       if ( parallel_IOProcessor() ) then
          write(unit=*, fmt='(i3,": Anorm=",g15.8,", bnorm=",g15.8,", Rnorm=",g15.8)') cnt, Anorm, bnorm, nrm
       end if
    end if

    i = 0
    if ( itsol_converged(rr, uu, Anorm, bnorm, eps) ) goto 100

    do i = 1, max_iter
       rho = dot(rt, rr)
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
       call itsol_precon_st(st, ph, pp)
       call stencil_apply_st(st, vv, ph); cnt = cnt + 1
       alpha = rho/dot(rt, vv)
       call saxpy(uu, alpha, ph)
       call saxpy(ss, rr, -alpha, vv)
       if ( verbose > 1 ) then
          nrm = norm_inf(ss)
          if ( parallel_IOProcessor() ) then
             write(unit=*, fmt='(i3,": Snorm=",g15.8)') cnt, nrm
          end if
       end if
       if ( itsol_converged(ss, uu, Anorm, bnorm, eps) ) exit
       call itsol_precon_st(st, sh, ss)
       call stencil_apply_st(st, tt, sh); cnt = cnt + 1

       omega = dot(tt,ss)/dot(tt,tt)
       call saxpy(uu, omega, sh)
       call saxpy(rr, ss, -omega, tt)
       if ( verbose > 1 ) then
          nrm = norm_inf(rr)
          if ( parallel_IOProcessor() ) then
             write(unit=*, fmt='(i3,": Rnorm=",g15.8)') cnt, nrm
          end if
       end if
       if ( itsol_converged(rr, uu, Anorm, bnorm, eps) ) exit
       rho_1 = rho
    end do

    if ( verbose > 0 ) then
       if ( parallel_IOProcessor() ) then
          write(unit=*, fmt='("BICGStab: iterations: ", i5)') cnt
       end if
    end if

100 continue
    call destroy(rr)
    call destroy(rt)
    call destroy(pp)
    call destroy(ph)
    call destroy(vv)
    call destroy(tt)
    call destroy(sh)
    call destroy(ss)

    if ( i > max_iter ) then
       if ( present(stat) ) then
          stat = 1
       else
          call bl_error("BiCGSolve: failed to converge");
       end if
    end if

  end subroutine itsol_BiCGStab_solve_st

  subroutine itsol_CG_Solve_st(st, uu, rh, eps, max_iter, verbose, stat)
    type(stencil), intent(in) :: st
    integer, intent(in) :: max_iter, verbose
    integer, intent(out), optional :: stat

    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: rh

    real(dp_t), intent(in) :: eps
    type(multifab) :: rr, zz, pp, qq
    real(kind = dp_t) :: rho_1, alpha, beta, Anorm, bnorm, rho, nrm
    type(layout) :: la
    integer :: i, ng_for_res
    logical :: nodal_solve
    integer :: cnt

    real(dp_t) :: rho_hg, rho_orig

    if ( present(stat) ) stat = 0

    ng_for_res = 0; if ( nodal_q(rh) ) ng_for_res = 1
    nodal_solve = .FALSE.; if ( ng_for_res /= 0 ) nodal_solve = .TRUE.

    la = st%ss%la
    call build(rr, la, 1, ng_for_res, rh%nodal)
    call build(zz, la, 1, ng_for_res, rh%nodal)
    call build(pp, la, 1, 1         , rh%nodal)
    call build(qq, la, 1, ng_for_res, rh%nodal)

    if ( nodal_solve ) then
       call setval(rr, ZERO, all=.true.)
       call setval(zz, ZERO, all=.true.)
       call setval(qq, ZERO, all=.true.)
    end if
    call setval(pp, ZERO, all=.true.)

    cnt = 0
    call stencil_defect_st(st, rr, rh, uu); cnt = cnt + 1

    Anorm = stencil_norm_st(st)
    bnorm = norm_inf(rh)

    if ( verbose > 1 ) then
       nrm = norm_inf(rr)
       if ( parallel_IOProcessor() ) then
          write(unit=*, fmt='(i3,": Anorm=",g15.8,", bnorm=",g15.8,", Rnorm=",g15.8)') 0, Anorm, bnorm, nrm
       end if
    end if

    i = 0
    if ( itsol_converged(rr, uu, Anorm, bnorm, eps) ) goto 100

    do i = 1, max_iter
       call itsol_precon_st(st, zz, rr)
       rho = dot(rr, zz)
       if ( i == 1 ) then
          call copy(pp, zz)
          rho_orig = rho
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
       call stencil_apply_st(st, qq, pp); cnt = cnt + 1
       alpha = rho/dot(pp, qq)
       call saxpy(uu,   alpha, pp)
       call saxpy(rr, - alpha, qq)
       if ( verbose > 1 ) then
          nrm = norm_inf(rr)
          if ( parallel_IOProcessor() ) then
             write(unit=*, fmt='(i3,": Rnorm=",g15.8)') i, nrm
          end if
       end if
       if ( .false. .and. nodal_solve ) then
          ! HACK, ONLY NEED THIS TO MATCH THE HGPROJ STOPPING CRITERION
          call itsol_precon_st(st, zz, rr)
          rho_hg = dot(rr, zz)
          if (rho_hg < rho_orig*eps) exit
       else
          if ( itsol_converged(rr, uu, Anorm, bnorm, eps) ) exit
       end if
       rho_1 = rho
    end do

    if ( verbose > 0 ) then
       if ( parallel_IOProcessor() ) then
          write(unit=*, fmt='("CG: iterations: ", i5)') cnt
       end if
    end if

100 continue
    call destroy(rr)
    call destroy(zz)
    call destroy(pp)
    call destroy(qq)

    if ( i > max_iter ) then
       if ( present(stat) ) then
          stat = 1
       else
          call bl_error("CG_solve: failed to converge");
       end if
    end if

  end subroutine itsol_CG_Solve_st

  subroutine itsol_precon_st(st, uu, rh, method)
    type(stencil), intent(in) :: st
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: rh
    real(kind=dp_t), pointer, dimension(:,:,:,:) :: ap, up, rp
    integer, pointer, dimension(:,:,:,:) :: mp
    integer :: i, n
    integer, intent(in), optional :: method
    integer :: lm
    lm = 1; if ( present(method) ) lm = method

    select case (lm)
    case (0)
       call copy(uu, rh)
    case (1)
       do i = 1, rh%nboxes
          if ( multifab_remote(uu, i) ) cycle
          rp => dataptr(rh, i)
          up => dataptr(uu, i)
          ap => dataptr(st%ss, i)
          mp => dataptr(st%mm, i)
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
      !$OMP PARALLEL DO PRIVATE(j,i)
      do j = 1, ny
         do i = 1, nx
            u(i,j) = r(i,j)/a(i,j,0)
         end do
      end do
      !$OMP END PARALLEL DO
    end subroutine jacobi_precon_2d
    subroutine jacobi_precon_3d(a, u, r, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)  :: a(:,:,:,0:)
      real(kind=dp_t), intent(inout) :: u(1-ng:,1-ng:,1-ng:)
      real(kind=dp_t), intent(in)  :: r(:,:,:)
      integer :: i, j, k, nx, ny, nz
      nz = size(a,dim=3)
      ny = size(a,dim=2)
      nx = size(a,dim=1)
      !$OMP PARALLEL DO PRIVATE(j,i,k)
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
      !$OMP PARALLEL DO PRIVATE(j,i)
      do j = 1, ny
         do i = 1, nx
            if (.not. bc_dirichlet(mm(i,j),1,0)) then
               u(i,j) = r(i,j)/a(i,j,0)
            end if
         end do
      end do
      !$OMP END PARALLEL DO
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
      !$OMP PARALLEL DO PRIVATE(j,i,k)
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
  end subroutine itsol_precon_st

end module itsol_module

