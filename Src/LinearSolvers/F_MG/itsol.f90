module itsol_module

  use bl_types
  use multifab_module
  use cc_stencil_module
  use cc_stencil_apply_module

  implicit none

  integer, private, parameter :: def_bicg_max_iter = 1000
  integer, private, parameter :: def_cg_max_iter   = 1000

  private :: itsol_defect, itsol_precon
  private :: jacobi_precon_1d, jacobi_precon_2d, jacobi_precon_3d
  private :: nodal_precon_1d, nodal_precon_2d, nodal_precon_3d

contains

    subroutine jacobi_precon_1d(a, u, r, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)    :: a(0:,:)
      real(kind=dp_t), intent(inout) :: u(1-ng:)
      real(kind=dp_t), intent(in)    :: r(:)
      integer :: i, nx
      nx = size(a,dim=2)
      do i = 1, nx
         u(i) = r(i)/a(0,i)
      end do
    end subroutine jacobi_precon_1d

    subroutine jacobi_precon_2d(a, u, r, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)    :: a(0:,:,:)
      real(kind=dp_t), intent(inout) :: u(1-ng:,1-ng:)
      real(kind=dp_t), intent(in)    :: r(:,:)
      integer :: i, j, nx, ny
      ny = size(a,dim=3)
      nx = size(a,dim=2)
      do j = 1, ny
         do i = 1, nx
            u(i,j) = r(i,j)/a(0,i,j)
         end do
      end do
    end subroutine jacobi_precon_2d

    subroutine jacobi_precon_3d(a, u, r, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)    :: a(0:,:,:,:)
      real(kind=dp_t), intent(inout) :: u(1-ng:,1-ng:,1-ng:)
      real(kind=dp_t), intent(in)    :: r(:,:,:)
      integer i, j, k, nx, ny, nz
      nz = size(a,dim=4)
      ny = size(a,dim=3)
      nx = size(a,dim=2)
      !$OMP PARALLEL DO PRIVATE(j,i,k) IF(nz.ge.7)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               u(i,j,k) = r(i,j,k)/a(0,i,j,k)
            end do
         end do
      end do
      !$OMP END PARALLEL DO
    end subroutine jacobi_precon_3d

    subroutine nodal_precon_1d(a, u, r, mm, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)  :: a(0:,:)
      real(kind=dp_t), intent(inout) :: u(1-ng:)
      real(kind=dp_t), intent(in)  :: r(0:)
      integer, intent(in)  :: mm(:)
      integer :: i, nx
      nx = size(a,dim=2)
      do i = 1, nx
         if (.not. bc_dirichlet(mm(i),1,0)) &
            u(i) = r(i)/a(0,i)
      end do
    end subroutine nodal_precon_1d

    subroutine nodal_precon_2d(a, u, r, mm, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)    :: a(0:,:,:)
      real(kind=dp_t), intent(inout) :: u(1-ng:,1-ng:)
      real(kind=dp_t), intent(in)    :: r(0:,0:)
      integer, intent(in)            :: mm(:,:)
      integer :: i, j, nx, ny
      ny = size(a,dim=3)
      nx = size(a,dim=2)
      do j = 1, ny
         do i = 1, nx
            if (.not. bc_dirichlet(mm(i,j),1,0)) then
               u(i,j) = r(i,j)/a(0,i,j)
            end if
         end do
      end do
    end subroutine nodal_precon_2d

    subroutine nodal_precon_3d(a, u, r, mm, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)    :: a(0:,:,:,:)
      real(kind=dp_t), intent(inout) :: u(1-ng:,1-ng:,1-ng:)
      real(kind=dp_t), intent(in)    :: r(0:,0:,0:)
      integer, intent(in)            :: mm(:,:,:)
      integer :: i, j, k, nx, ny, nz
      nz = size(a,dim=4)
      ny = size(a,dim=3)
      nx = size(a,dim=2)
      !$OMP PARALLEL DO PRIVATE(j,i,k) IF(nz.ge.7)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (.not. bc_dirichlet(mm(i,j,k),1,0)) then
                  u(i,j,k) = r(i,j,k)/a(0,i,j,k)
               end if
            end do
         end do
      end do
      !$OMP END PARALLEL DO
    end subroutine nodal_precon_3d

    subroutine diag_init_cc_1d(a, ng_a, r, ng_r, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(0:,lo(1)-ng_a:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:   )

      integer         :: i, nc
      real(kind=dp_t) :: denom

      nc = size(a,dim=1)-1

      ! Protect against divide by zero -- necessary for embedded boundary problems.
      do i = lo(1),hi(1)
         if (abs(a(0,i)) .gt. 0.d0) then
            denom = 1.d0 / a(0,i)
            r(i     ) = r(i     ) * denom
            a(1:nc,i) = a(1:nc,i) * denom
            a(0,i   ) = 1.d0
         end if
      end do

    end subroutine diag_init_cc_1d

    subroutine diag_init_cc_2d(a, ng_a, r, ng_r, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(0:,lo(1)-ng_a:,lo(2)-ng_a:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:,lo(2)-ng_r:   )

      integer         :: i, j, nc
      real(kind=dp_t) :: denom

      nc = size(a,dim=1)-1

      ! Protect against divide by zero -- necessary for embedded boundary problems.
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            if (abs(a(0,i,j)) .gt. 0.d0) then
               denom = 1.d0 / a(0,i,j)
               r(i,j     ) = r(i,j     ) * denom
               a(1:nc,i,j) = a(1:nc,i,j) * denom
               a(0,i,j   ) = 1.d0
            end if
         end do
      end do

    end subroutine diag_init_cc_2d

    subroutine diag_init_cc_3d(a, ng_a, r, ng_r, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(0:,lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:   )

      integer         :: i, j, k, nc
      real(kind=dp_t) :: denom

      nc = size(a,dim=1)-1

      !$OMP PARALLEL DO PRIVATE(j,i,k,denom) IF((hi(3)-lo(3)).ge.7)
      ! Protect against divide by zero -- necessary for embedded boundary problems.
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (abs(a(0,i,j,k)) .gt. 0.d0) then
                  denom = 1.d0 / a(0,i,j,k)
                  r(i,j,k     ) = r(i,j,k     ) * denom
                  a(1:nc,i,j,k) = a(1:nc,i,j,k) * denom
                  a(0,i,j,k   ) = 1.d0
               end if
            end do
         end do
      end do
      !$OMP END PARALLEL DO

    end subroutine diag_init_cc_3d

    subroutine diag_init_nd_1d(a, ng_a, r, ng_r, mm, ng_m, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r, ng_m
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(0:,lo(1)-ng_a:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:   )
      integer        , intent(inout)  :: mm(lo(1)-ng_m:   )

      integer         :: i, nc
      real(kind=dp_t) :: denom

      nc = size(a,dim=1)-1

      do i = lo(1),hi(1)+1
         if (.not. bc_dirichlet(mm(i),1,0)) then
            denom = 1.d0 / a(0,i)
            r(i     ) = r(i     ) * denom
            a(1:nc,i) = a(1:nc,i) * denom
            a(0,i   ) = 1.d0
         end if
      end do

    end subroutine diag_init_nd_1d

    subroutine diag_init_nd_2d(a, ng_a, r, ng_r, mm, ng_m, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r, ng_m
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(0:,lo(1)-ng_a:,lo(2)-ng_a:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:,lo(2)-ng_r:   )
      integer        , intent(inout)  :: mm(lo(1)-ng_m:,lo(2)-ng_m:   )

      integer         :: i, j, nc
      real(kind=dp_t) :: denom

      nc = size(a,dim=1)-1

      do j = lo(2),hi(2)+1
         do i = lo(1),hi(1)+1
            if (.not. bc_dirichlet(mm(i,j),1,0)) then
               denom = 1.d0 / a(0,i,j)
               r(i,j     ) = r(i,j     ) * denom
               a(1:nc,i,j) = a(1:nc,i,j) * denom
               a(0,i,j   ) = 1.d0
            end if
         end do
      end do

    end subroutine diag_init_nd_2d

    subroutine diag_init_nd_3d(a, ng_a, r, ng_r, mm, ng_m, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r, ng_m
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(0:,lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:   )
      integer        , intent(inout)  :: mm(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:   )

      integer         :: i, j, k, nc
      real(kind=dp_t) :: denom

      nc = size(a,dim=1)-1

      !$OMP PARALLEL DO PRIVATE(j,i,k,denom) IF((hi(3)-lo(3)).ge.6)
      do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)+1
               if (.not. bc_dirichlet(mm(i,j,k),1,0)) then
                  denom = 1.d0 / a(0,i,j,k)
                  r(i,j,k     ) = r(i,j,k     ) * denom
                  a(1:nc,i,j,k) = a(1:nc,i,j,k) * denom
                  a(0,i,j,k   ) = 1.d0
               end if
            end do
         end do
      end do
      !$OMP END PARALLEL DO

    end subroutine diag_init_nd_3d

  function itsol_converged(rr, uu, bnorm, eps, abs_eps, rrnorm) result(r)
    use bl_prof_module

    type(multifab), intent(in )           :: rr, uu
    real(dp_t),     intent(in )           :: bnorm, eps
    real(dp_t),     intent(in ), optional :: abs_eps
    real(dp_t),     intent(out), optional :: rrnorm

    real(dp_t) :: norm_rr, norm_uu, tnorms(2), rtnorms(2)
    logical    :: r

    type(bl_prof_timer), save :: bpt

    call build(bpt, "its_converged")
    !
    ! Elide a reduction by doing both reductions together.
    !
    tnorms(1) = norm_inf(rr, local=.true.)
    tnorms(2) = norm_inf(uu, local=.true.)

    call parallel_reduce(rtnorms, tnorms, MPI_MAX)

    norm_rr = rtnorms(1)
    norm_uu = rtnorms(2)

    if ( present(rrnorm) ) rrnorm = norm_rr

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

  ! Computes rr = aa * uu
  subroutine itsol_stencil_apply(aa, rr, uu, mm, stencil_type, lcross, uniform_dh)

    use bl_prof_module

    use nodal_stencil_module, only: stencil_apply_1d_nodal, &
                                    stencil_apply_2d_nodal, &
                                    stencil_apply_3d_nodal

    type(multifab), intent(in)    :: aa
    type(multifab), intent(inout) :: rr
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in)   :: mm
    integer, intent(in)           :: stencil_type
    logical, intent(in)           :: lcross
    logical, intent(in),optional  :: uniform_dh

    logical                       :: luniform_dh

    real(kind=dp_t), pointer :: rp(:,:,:,:), up(:,:,:,:), ap(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)

    integer :: i, n, lo(get_dim(rr)), hi(get_dim(rr)), dm
    logical :: nodal_flag
    type(bl_prof_timer), save :: bpt

    call build(bpt, "its_stencil_apply")

    luniform_dh = .false. ; if ( present(uniform_dh) ) luniform_dh = uniform_dh

    call multifab_fill_boundary(uu, cross = lcross)

    dm = get_dim(rr)

    nodal_flag = nodal_q(uu)

    do i = 1, nfabs(rr)
       rp => dataptr(rr, i)
       up => dataptr(uu, i)
       ap => dataptr(aa, i)
       mp => dataptr(mm, i)
       lo = lwb(get_box(uu,i))
       hi = upb(get_box(uu,i))
       do n = 1, ncomp(rr)
          select case(dm)
          case (1)
             if ( .not. nodal_flag) then
                call stencil_apply_1d(ap(:,:,1,1), rp(:,1,1,n), nghost(rr), up(:,1,1,n), nghost(uu),  &
                                      mp(:,1,1,1), lo, hi, stencil_type)
             else
                call stencil_apply_1d_nodal(ap(:,:,1,1), rp(:,1,1,n), up(:,1,1,n),  &
                     mp(:,1,1,1), nghost(uu), stencil_type)
             end if
          case (2)
             if ( .not. nodal_flag) then
                call stencil_apply_2d(ap(:,:,:,1), rp(:,:,1,n), nghost(rr), up(:,:,1,n), nghost(uu),  &
                     mp(:,:,1,1), lo, hi, stencil_type)
             else
                call stencil_apply_2d_nodal(ap(:,:,:,1), rp(:,:,1,n), up(:,:,1,n),  &
                     mp(:,:,1,1), nghost(uu), stencil_type)
             end if
          case (3)
             if ( .not. nodal_flag) then
                call stencil_apply_3d(ap(:,:,:,:), rp(:,:,:,n), nghost(rr), up(:,:,:,n), nghost(uu),  &
                                      mp(:,:,:,1), stencil_type)
             else
                call stencil_apply_3d_nodal(ap(:,:,:,:), rp(:,:,:,n), up(:,:,:,n),  &
                     mp(:,:,:,1), nghost(uu), stencil_type, luniform_dh)
             end if
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine itsol_stencil_apply

  ! computes rr = aa * uu - rh
  subroutine itsol_defect(ss, rr, rh, uu, mm, stencil_type, lcross, uniform_dh)
    use bl_prof_module
    type(multifab), intent(inout) :: uu, rr
    type(multifab), intent(in)    :: rh, ss
    type(imultifab), intent(in)   :: mm
    integer, intent(in)           :: stencil_type
    logical, intent(in)           :: lcross
    logical, intent(in), optional :: uniform_dh
    type(bl_prof_timer), save     :: bpt
    call build(bpt, "its_defect")
    call itsol_stencil_apply(ss, rr, uu, mm, stencil_type, lcross, uniform_dh)
    call saxpy(rr, rh, -1.0_dp_t, rr)
    call destroy(bpt)
  end subroutine itsol_defect

  subroutine itsol_BiCGStab_solve(aa, uu, rh, mm, eps, max_iter, verbose, stencil_type, lcross, &
       stat, singular_in, uniform_dh, nodal_mask)
    use bl_prof_module
    integer, intent(in) :: max_iter
    type(imultifab), intent(in) :: mm
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: rh
    type(multifab), intent(in) :: aa
    integer        , intent(in) :: stencil_type, verbose
    logical        , intent(in) :: lcross
    real(kind=dp_t), intent(in) :: eps

    integer, intent(out), optional :: stat
    logical, intent(in), optional :: singular_in
    logical, intent(in), optional :: uniform_dh
    type(multifab), intent(in), optional :: nodal_mask

    type(layout) :: la
    type(multifab) :: rr, rt, pp, ph, vv, tt, ss, rh_local, aa_local
    real(kind=dp_t) :: rho_1, alpha, beta, omega, rho, Anorm, bnorm, rnorm, den
    real(dp_t) :: rho_orig, volume, tres0, small, norm_rr, norm_uu
    real(dp_t) :: tnorms(3),rtnorms(3)
    integer :: i, cnt, ng_for_res
    logical :: nodal_solve, singular, nodal(get_dim(rh)), diag_inited
    real(dp_t), pointer :: pdst(:,:,:,:), psrc(:,:,:,:)
    type(bl_prof_timer), save :: bpt

    call build(bpt, "its_BiCGStab_solve")

    if ( present(stat) ) stat = 0

    singular    = .false.; if ( present(singular_in) ) singular    = singular_in
    ng_for_res  = 0;       if ( nodal_q(rh)          ) ng_for_res  = 1
    nodal_solve = .false.; if ( ng_for_res /= 0      ) nodal_solve = .true.

    nodal = nodal_flags(rh)

    la = get_layout(aa)

    call multifab_build(rr, la, 1, ng_for_res, nodal)
    call multifab_build(rt, la, 1, ng_for_res, nodal)
    call multifab_build(pp, la, 1, ng_for_res, nodal)
    call multifab_build(ph, la, 1, nghost(uu), nodal)
    call multifab_build(vv, la, 1, ng_for_res, nodal)
    call multifab_build(tt, la, 1, ng_for_res, nodal)
    call multifab_build(ss, la, 1, ng_for_res, nodal)
    !
    ! Use these for local preconditioning.
    !
    call multifab_build(rh_local, la, ncomp(rh), nghost(rh), nodal)
    call multifab_build(aa_local, la, ncomp(aa), nghost(aa), nodal_flags(aa), stencil = .true.)

    if ( nodal_solve ) then
       call setval(rr, ZERO, all=.true.)
       call setval(rt, ZERO, all=.true.)
       call setval(pp, ZERO, all=.true.)
       call setval(vv, ZERO, all=.true.)
       call setval(tt, ZERO, all=.true.)
       call setval(ss, ZERO, all=.true.)
    end if

    call copy(rh_local, 1, rh, 1, nc = ncomp(rh), ng = nghost(rh))
    !
    ! Copy aa -> aa_local; gotta do it by hand since it's a stencil multifab.
    !
    do i = 1, nfabs(aa)
       pdst => dataptr(aa_local, i)
       psrc => dataptr(aa      , i)
       call cpy_d(pdst, psrc)
    end do
    !
    ! Make sure to do singular adjustment *before* diagonalization.
    !
    if ( singular ) then
      call setval(ss,ONE)
      tnorms(1) = dot(rh_local, ss, nodal_mask, local = .true.)
      tnorms(2) = dot(      ss, ss, nodal_mask, local = .true.)
      call parallel_reduce(rtnorms(1:2), tnorms(1:2), MPI_SUM)
      rho    = rtnorms(1)
      volume = rtnorms(2)
      rho    = rho / volume
      if ( parallel_IOProcessor() .and. verbose > 0 ) then
         print *,'...singular adjustment to rhs: ',rho
      endif
      call saxpy(rh_local,-rho,ss)
      call setval(ss,ZERO,all=.true.)
    end if

    call diag_initialize(aa_local,rh_local,mm); diag_inited = .true.

    call copy(ph, uu, ng = nghost(ph))

    cnt = 0
    !
    ! Compute rr = aa * uu - rh.
    !
    call itsol_defect(aa_local, rr, rh_local, uu, mm, stencil_type, lcross, uniform_dh); cnt = cnt + 1

    call copy(rt, rr)
    rho      = dot(rt, rr, nodal_mask)
    rho_orig = rho
    !
    ! Elide some reductions by calculating local norms & then reducing all together.
    !
    tnorms(1) = norm_inf(rr,           local=.true.)
    tnorms(2) = norm_inf(rh_local,     local=.true.)
    tnorms(3) = stencil_norm(aa_local, local=.true.)

    call parallel_reduce(rtnorms, tnorms, MPI_MAX)

    tres0 = rtnorms(1)
    bnorm = rtnorms(2)
    Anorm = rtnorms(3)
    small = epsilon(Anorm)

    if ( parallel_IOProcessor() .and. verbose > 0 ) then
       if ( diag_inited ) then
          write(*,*) "   BiCGStab: A and rhs have been rescaled. So do the error."
       end if
       write(unit=*, fmt='("    BiCGStab: Initial error (error0) =        ",g15.8)') tres0
    end if 

    if ( itsol_converged(rr, uu, bnorm, eps, rrnorm=norm_rr) ) then
       if ( verbose > 0 ) then
          if ( tres0 < eps*bnorm ) then
             if ( parallel_IOProcessor() ) then
                write(unit=*, fmt='("    BiCGStab: Zero iterations: rnorm ",g15.8," < eps*bnorm ",g15.8)') &
                     tres0,eps*bnorm
             end if
          else
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

    rho_1 = ZERO

    do i = 1, max_iter
       rho = dot(rt, rr, nodal_mask)
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
       call itsol_stencil_apply(aa_local, vv, ph, mm, stencil_type, lcross, uniform_dh)
       cnt = cnt + 1
       den = dot(rt, vv, nodal_mask)
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
       if ( itsol_converged(ss, uu, bnorm, eps) ) exit
       call itsol_precon(aa_local, ph, ss, mm,0)
       call itsol_stencil_apply(aa_local, tt, ph, mm, stencil_type, lcross, uniform_dh) 
       cnt = cnt + 1
       !
       ! Elide a reduction here by calculating the two dot-products
       ! locally and then reducing them both in a single call.
       !
       tnorms(1) = dot(tt, tt, nodal_mask, local = .true.)
       tnorms(2) = dot(tt, ss, nodal_mask, local = .true.)

       call parallel_reduce(rtnorms(1:2), tnorms(1:2), MPI_SUM)

       den   = rtnorms(1)
       omega = rtnorms(2)

       if ( den == ZERO ) then
          if ( present(stat) ) then
             call bl_warn("BICGSTAB_solve: breakdown in bicg, going with what I have")
             stat = 31
             goto 100
          endif
          call bl_error("BiCGStab: failure 3")
       end if
       omega = omega/den
       call saxpy(uu, omega, ph)
       call saxpy(rr, ss, -omega, tt)
       rnorm = norm_inf(rr)
       if ( parallel_IOProcessor() .and. verbose > 1 ) then
          write(unit=*, fmt='("    BiCGStab: Iteration        ",i4," rel. err. ",g15.8)') cnt/2, &
               rnorm /  (bnorm)
       end if

       if ( itsol_converged(rr, uu, bnorm, eps) ) exit

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

    if ( rnorm > bnorm ) then
       call setval(uu,ZERO,all=.true.)
       if ( present(stat) ) stat = 1
       if ( verbose > 0 .and.  parallel_IOProcessor() ) &
            print *,'   BiCGStab: solution reset to zero'
    end if

    if ( i > max_iter ) then
       if ( present(stat) ) then
          stat = 1
       else
          call bl_error("BiCGSolve: failed to converge");
       end if
    end if

100 continue

    call destroy(rh_local)
    call destroy(aa_local)
    call destroy(rr)
    call destroy(rt)
    call destroy(pp)
    call destroy(ph)
    call destroy(vv)
    call destroy(tt)
    call destroy(ss)

    call destroy(bpt)

  end subroutine itsol_BiCGStab_solve

  subroutine itsol_CABiCGStab_solve(aa, uu, rh, mm, eps, max_iter, verbose, stencil_type, lcross, &
       stat, singular_in, uniform_dh, nodal_mask)
    use bl_prof_module
    integer,         intent(in   ) :: max_iter
    type(imultifab), intent(in   ) :: mm
    type(multifab),  intent(inout) :: uu
    type(multifab),  intent(in   ) :: rh
    type(multifab),  intent(in   ) :: aa
    integer        , intent(in   ) :: stencil_type, verbose
    logical        , intent(in   ) :: lcross
    real(kind=dp_t), intent(in   ) :: eps

    integer,        intent(out), optional :: stat
    logical,        intent(in ), optional :: singular_in
    logical,        intent(in ), optional :: uniform_dh
    type(multifab), intent(in ), optional :: nodal_mask

    type(layout)    :: la
    type(multifab)  :: rr, rt, pp, pr, ss, rh_local, aa_local,    ph, vv, tt
    real(kind=dp_t) :: rho_1, alpha, beta, omega, rho, Anorm, bnorm, rnorm, den
    real(dp_t)      :: rnorm0, small, norm_rr, norm_uu, delta, L2_norm_of_rt
    real(dp_t)      :: tnorms(3),rtnorms(3), L2_norm_of_resid, atime, gtime, time1, time2
    integer         :: i, m, cnt, ng_for_res
    logical         :: nodal_solve, singular, nodal(get_dim(rh))
    logical         :: BiCGStabFailed, BiCGStabConverged

    real(dp_t), pointer :: pdst(:,:,:,:), psrc(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    integer, parameter :: SSS_MAX = 4

    real(dp_t)  temp1(4*SSS_MAX+1)
    real(dp_t)  temp2(4*SSS_MAX+1)
    real(dp_t)  temp3(4*SSS_MAX+1)
    real(dp_t)     Tp(4*SSS_MAX+1, 4*SSS_MAX+1)
    real(dp_t)    Tpp(4*SSS_MAX+1, 4*SSS_MAX+1)
    real(dp_t)     aj(4*SSS_MAX+1)
    real(dp_t)     cj(4*SSS_MAX+1)
    real(dp_t)     ej(4*SSS_MAX+1)
    real(dp_t)   Tpaj(4*SSS_MAX+1)
    real(dp_t)   Tpcj(4*SSS_MAX+1)
    real(dp_t)  Tppaj(4*SSS_MAX+1)
    real(dp_t)     cG(4*SSS_MAX+1,4*SSS_MAX+1)     ! G in C++ code
    real(dp_t)     lG(4*SSS_MAX+1)                 ! g in C++ code
    real(dp_t)     Gg((4*SSS_MAX+1)*(4*SSS_MAX+2))

    call build(bpt, "its_CABiCGStab_solve")

    stop 'itsol_CABiCGStab_solve() not fully implemented'

    if ( present(stat) ) stat = 0

    singular    = .false.; if ( present(singular_in) ) singular    = singular_in
    ng_for_res  = 0;       if ( nodal_q(rh)          ) ng_for_res  = 1
    nodal_solve = .false.; if ( ng_for_res /= 0      ) nodal_solve = .true.

    la    = get_layout(aa)
    nodal = nodal_flags(rh)

    aj    = 0.0d0
    cj    = 0.0d0
    ej    = 0.0d0
    Tpaj  = 0.0d0
    Tpcj  = 0.0d0
    Tppaj = 0.0d0
    temp1 = 0.0d0
    temp2 = 0.0d0
    temp3 = 0.0d0

    call SetMonomialBasis(SSS_MAX)

    call multifab_build(rr, la, 1, ng_for_res, nodal)
    call multifab_build(rt, la, 1, ng_for_res, nodal)
    call multifab_build(pp, la, 1, ng_for_res, nodal)
    call multifab_build(ph, la, 1, nghost(uu), nodal)
    !
    ! Contains the matrix powers of pp[] and rr[].
    !
    ! First 2*SSS+1 components are powers of pp[].
    ! Next  2*SSS   components are powers of rr[].
    !
    call multifab_build(pr, la, 4*SSS_MAX+1, ng_for_res, nodal)
    !
    ! Use these for local preconditioning.
    !
    call multifab_build(rh_local, la, ncomp(rh), nghost(rh), nodal)
    call multifab_build(aa_local, la, ncomp(aa), nghost(aa), nodal_flags(aa), stencil = .true.)

    if ( nodal_solve ) then
       call setval(rr, ZERO, all = .true.)
       call setval(rt, ZERO, all = .true.)
       call setval(pp, ZERO, all = .true.)
       call setval(pr, ZERO, all = .true.)
    end if

    call copy(rh_local, 1, rh, 1, nc = ncomp(rh), ng = nghost(rh))
    !
    ! Copy aa -> aa_local; gotta do it by hand since it's a stencil multifab.
    !
    do i = 1, nfabs(aa)
       pdst => dataptr(aa_local, i)
       psrc => dataptr(aa      , i)
       call cpy_d(pdst, psrc)
    end do
    !
    ! Make sure to do singular adjustment *before* diagonalization.
    !
    if ( singular ) then
       call multifab_build(ss, la, 1, ng_for_res, nodal)
       call setval(ss,ONE)
       tnorms(1) = dot(rh_local, ss, nodal_mask, local = .true.)
       tnorms(2) = dot(      ss, ss, nodal_mask, local = .true.)
       call parallel_reduce(rtnorms(1:2), tnorms(1:2), MPI_SUM)
       rho = rtnorms(1) / rtnorms(2)
       if ( parallel_IOProcessor() .and. verbose > 0 ) then
          print *,'...singular adjustment to rhs: ', rho
       endif
       call saxpy(rh_local,-rho,ss)
       call destroy(ss)
    end if

    call diag_initialize(aa_local,rh_local,mm)

    call copy(ph, uu, ng = nghost(ph))

    cnt = 0
    !
    ! Compute rr = aa * uu - rh.
    !
    call itsol_defect(aa_local, rr, rh_local, uu, mm, stencil_type, lcross, uniform_dh); cnt = cnt + 1

    call copy(rt, rr); call copy(pp, rr)
    !
    ! Elide some reductions by calculating local norms & then reducing all together.
    !
    tnorms(1) = norm_inf(rr,           local = .true.)
    tnorms(2) = norm_inf(rh_local,     local = .true.)
    tnorms(3) = stencil_norm(aa_local, local = .true.)

    call parallel_reduce(rtnorms, tnorms, MPI_MAX)

    rnorm0 = rtnorms(1)
    bnorm  = rtnorms(2)
    Anorm  = rtnorms(3)

    small         = epsilon(Anorm)
    delta         = dot(rt, rr, nodal_mask)
    L2_norm_of_rt = sqrt(delta)

    if ( parallel_IOProcessor() .and. verbose > 0 ) then
       write(*,*) "   CABiCGStab: A and rhs have been rescaled. So has the error."
       write(unit=*, fmt='("    CABiCGStab: Initial error (error0) =        ",g15.8)') rnorm0
    end if 

    if ( itsol_converged(rr, uu, bnorm, eps, rrnorm=norm_rr) .or. delta .eq. 0.0d0 ) then
       if ( verbose > 0 ) then
          if ( rnorm0 < eps*bnorm ) then
             if ( parallel_IOProcessor() ) then
                write(unit=*, fmt='("    CABiCGStab: Zero iterations: rnorm ",g15.8," < eps*bnorm ",g15.8)') &
                     rnorm0,eps*bnorm
             end if
          else if ( norm_rr < epsilon(Anorm)*Anorm ) then
             if ( parallel_IOProcessor() ) then
                write(unit=*, fmt='("    CABiCGStab: Zero iterations: rnorm ",g15.8," < small*Anorm ",g15.8)') &
                     rnorm0,small*Anorm
             end if
          else if ( delta .eq. 0.0d0 ) then
             if ( parallel_IOProcessor() ) then
                write(unit=*, fmt='("    CABiCGStab: Zero iterations: delta == 0")')
             end if
          end if
       end if
       go to 100
    end if

    L2_norm_of_resid = 0

    BiCGStabFailed = .false. ; BiCGStabConverged = .false.

    atime = 0.0d0; gtime = 0.0d0

    do i = 1, max_iter

       call itsol_precon(aa_local, ph, pp, mm, 0)
       call itsol_stencil_apply(aa_local, vv, ph, mm, stencil_type, lcross, uniform_dh)

       cnt = cnt + 1


!       if ( itsol_converged(rr, uu, bnorm, eps) ) exit

    end do

    if ( verbose > 0 ) then
       if ( parallel_IOProcessor() ) then
          write(unit=*, fmt='("    CABiCGStab: Final: Iteration  ", i3, " rel. err. ",g15.8)') cnt/2, &
               rnorm/ (bnorm)
       end if
       if ( rnorm < eps*bnorm ) then
          if ( parallel_IOProcessor() ) then
             write(unit=*, fmt='("    CABiCGStab: Converged: rnorm ",g15.8," < eps*bnorm ",g15.8)') &
                  rnorm,eps*bnorm
          end if
       else
          norm_uu = norm_inf(uu)
          if ( rnorm < eps*Anorm*norm_uu ) then
             if ( parallel_IOProcessor() ) then
                write(unit=*, fmt='("    CABiCGStab: Converged: rnorm ",g15.8," < eps*Anorm*sol_norm ",g15.8)') &
                     rnorm,eps*Anorm*norm_uu
             end if
          else if ( rnorm < epsilon(Anorm)*Anorm ) then
             if ( parallel_IOProcessor() ) then
                write(unit=*, fmt='("    CABiCGStab: Converged: rnorm ",g15.8," < small*Anorm ",g15.8)') &
                     rnorm,small*Anorm
             end if
          end if
       end if
    end if

    if ( rnorm > bnorm ) then
       call setval(uu,ZERO,all=.true.)
       if ( present(stat) ) stat = 1
       if ( verbose > 0 .and.  parallel_IOProcessor() ) &
            print *,'   CABiCGStab: solution reset to zero'
    end if

    if ( i > max_iter ) then
       if ( present(stat) ) then
          stat = 1
       else
          call bl_error("CABiCGSolve: failed to converge");
       end if
    end if

100 continue

    call destroy(rh_local)
    call destroy(aa_local)
    call destroy(rr)
    call destroy(rt)
    call destroy(pp)
    call destroy(pr)
    call destroy(ph)

    call destroy(bpt)

  contains

    subroutine SetMonomialBasis(sss)

      integer, intent(in) :: sss

      Tp = 0.0d0

      do i = 0,2*sss-1
         Tp(i+1,i) = 1.0d0
      end do
      do i = 2*sss+1, 4*sss-1
         Tp(i+1,i) = 1.0d0
      end do

      Tpp = 0.0d0

      do i = 0,2*sss-2
         Tpp(i+2,i) = 1.0d0
      end do
      do i = 2*sss+1, 4*sss-2
         Tpp(i+2,i) = 1.0d0
      end do

    end subroutine SetMonomialBasis

  end subroutine itsol_CABiCGStab_solve

  subroutine itsol_CG_Solve(aa, uu, rh, mm, eps, max_iter, verbose, stencil_type, lcross, &
                            stat, singular_in, uniform_dh, nodal_mask)
    use bl_prof_module
    integer    , intent(in   ) :: max_iter, verbose, stencil_type
    logical    , intent(in   ) :: lcross
    real(dp_t) , intent(in   ) :: eps

    integer, intent(  out), optional :: stat
    logical, intent(in   ), optional :: singular_in
    logical, intent(in   ), optional :: uniform_dh
    type(multifab), intent(in), optional :: nodal_mask

    type( multifab), intent(in)    :: aa
    type( multifab), intent(inout) :: uu
    type( multifab), intent(in)    :: rh
    type(imultifab), intent(in)    :: mm

    type(multifab) :: rr, zz, pp, qq
    type(multifab) :: aa_local, rh_local
    real(kind = dp_t) :: rho_1, alpha, beta, Anorm, bnorm, rho, rnorm, den, tres0, small
    type(layout) :: la
    integer :: i, ng_for_res
    logical :: nodal_solve, nodal(get_dim(rh))
    logical :: singular 
    integer :: cnt
    real(dp_t), pointer :: pdst(:,:,:,:), psrc(:,:,:,:)
    real(dp_t) :: tnorms(3), rtnorms(3)

    real(dp_t) :: rho_orig, volume, rho_hg_orig, norm_uu
    type(bl_prof_timer), save :: bpt

    call build(bpt, "its_CG_Solve")

    if ( present(stat) ) stat = 0

    singular    = .false.; if ( present(singular_in) ) singular = singular_in
    ng_for_res  = 0;       if ( nodal_q(rh)          ) ng_for_res = 1
    nodal_solve = .false.; if ( ng_for_res /= 0      ) nodal_solve = .true.

    nodal = nodal_flags(rh)

    la = get_layout(aa)
    call multifab_build(rr, la, 1, ng_for_res, nodal)
    call multifab_build(zz, la, 1, ng_for_res, nodal)
    call multifab_build(pp, la, 1, nghost(uu), nodal)
    call multifab_build(qq, la, 1, ng_for_res, nodal)

    if ( nodal_solve ) then
       call setval(rr,ZERO,all=.true.)
       call setval(zz,ZERO,all=.true.)
       call setval(qq,ZERO,all=.true.)
    end if
    call setval(pp, ZERO, all=.true.)

    ! Use these for local preconditioning
    call multifab_build(rh_local, la, ncomp(rh), nghost(rh), nodal)

    call multifab_build(aa_local, la, ncomp(aa), nghost(aa), nodal_flags(aa), stencil = .true.)

    call copy(rh_local, 1, rh, 1, nc = ncomp(rh), ng = nghost(rh))

    ! Copy aa -> aa_local; gotta do it by hand since it's a stencil multifab.
    do i = 1, nfabs(aa)
       pdst => dataptr(aa_local, i)
       psrc => dataptr(aa      , i)
       call cpy_d(pdst, psrc)
    end do

    call diag_initialize(aa_local,rh_local,mm)

    cnt = 0
    ! compute rr = aa * uu - rh_local
    call itsol_defect(aa_local, rr, rh_local, uu, mm, stencil_type, lcross, uniform_dh)  
    cnt = cnt + 1

    if ( singular .and. nodal_solve ) then
      call setval(zz,ONE)
      rho = dot(rr, zz, nodal_mask)
      volume = dot(zz,zz)
      rho = rho / volume
      call saxpy(rr,-rho,zz)
      call setval(zz,ZERO,all=.true.)
    end if
    !
    ! Elide some reductions by calculating local norms & then reducing all together.
    !
    tnorms(1) = norm_inf(rr,           local=.true.)
    tnorms(2) = norm_inf(rh_local,     local=.true.)
    tnorms(3) = stencil_norm(aa_local, local=.true.)

    call parallel_reduce(rtnorms, tnorms, MPI_MAX)

    tres0 = rtnorms(1)
    bnorm = rtnorms(2)
    Anorm = rtnorms(3)

    small = epsilon(Anorm)

    if ( parallel_IOProcessor() .and. verbose > 0) then
       write(unit=*, fmt='("          CG: Initial error (error0) =        ",g15.8)') tres0
    end if

    i = 0
    if ( itsol_converged(rr, uu, bnorm, eps) ) then
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

    rho_1       = ZERO
    rho_hg_orig = ZERO

    do i = 1, max_iter

       call itsol_precon(aa_local, zz, rr, mm, 0)
       rho = dot(rr, zz, nodal_mask)
       if ( i == 1 ) then
          call copy(pp, zz)
          rho_orig = rho
          call itsol_precon(aa_local, zz, rr, mm)
          rho_hg_orig = dot(rr, zz, nodal_mask)
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
       call itsol_stencil_apply(aa_local, qq, pp, mm, stencil_type, lcross, uniform_dh) 
       cnt = cnt + 1
       den = dot(pp, qq, nodal_mask)
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
       if ( parallel_IOProcessor() .and. verbose > 1 ) then
          write(unit=*, fmt='("          CG: Iteration        ",i4," rel. err. ",g15.8)') i, &
                             rnorm /  (bnorm)
       end if

       if ( itsol_converged(rr, uu, bnorm, eps) ) exit

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
    integer :: i, n, dm
    integer, intent(in), optional :: method
    integer :: lm
    type(bl_prof_timer), save :: bpt

    call build(bpt, "its_precon")

    lm = 1; if ( present(method) ) lm = method

    dm = get_dim(uu)

    select case (lm)
    case (0)
       call copy(uu, rh)
    case (1)
       do i = 1, nfabs(rh)
          rp => dataptr(rh, i)
          up => dataptr(uu, i)
          ap => dataptr(aa, i)
          mp => dataptr(mm, i)
          do n = 1, ncomp(uu)
             select case(dm)
             case (1)
                if ( cell_centered_q(rh) ) then
                   call jacobi_precon_1d(ap(:,:,1,1), up(:,1,1,n), rp(:,1,1,n), nghost(uu))
                else
                   call nodal_precon_1d(ap(:,:,1,1), up(:,1,1,n), rp(:,1,1,n), &
                                        mp(:,1,1,1),nghost(uu))
                end if
             case (2)
                if ( cell_centered_q(rh) ) then
                   call jacobi_precon_2d(ap(:,:,:,1), up(:,:,1,n), rp(:,:,1,n), nghost(uu))
                else
                   call nodal_precon_2d(ap(:,:,:,1), up(:,:,1,n), rp(:,:,1,n), &
                                        mp(:,:,1,1),nghost(uu))
                end if
             case (3)
                if ( cell_centered_q(rh) ) then
                   call jacobi_precon_3d(ap(:,:,:,:), up(:,:,:,n), rp(:,:,:,n), nghost(uu))
                else
                   call nodal_precon_3d(ap(:,:,:,:), up(:,:,:,n), rp(:,:,:,n), &
                                        mp(:,:,:,1),nghost(uu))
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
    integer                                      :: i,dm
    integer                                      :: ng_a, ng_r, ng_m
    integer                                      :: lo(get_dim(rh)),hi(get_dim(rh))
    type(bl_prof_timer), save                    :: bpt

    call build(bpt, "diag_initialize")

    ng_a = nghost(aa)
    ng_r = nghost(rh)
    ng_m = nghost(mm)

    dm = get_dim(rh)

    do i = 1, nfabs(rh)
       rp => dataptr(rh, i)
       ap => dataptr(aa, i)
       mp => dataptr(mm, i)
       lo = lwb(get_box(rh,i))
       hi = upb(get_box(rh,i))
       select case(dm)
          case (1)
             if ( cell_centered_q(rh) ) then
                call diag_init_cc_1d(ap(:,:,1,1), ng_a, rp(:,1,1,1), ng_r, lo, hi)
             else
                call diag_init_nd_1d(ap(:,:,1,1), ng_a, rp(:,1,1,1), ng_r, mp(:,1,1,1), ng_m, lo, hi)
             end if
          case (2)
             if ( cell_centered_q(rh) ) then
                call diag_init_cc_2d(ap(:,:,:,1), ng_a, rp(:,:,1,1), ng_r, lo, hi)
             else
                call diag_init_nd_2d(ap(:,:,:,1), ng_a, rp(:,:,1,1), ng_r, mp(:,:,1,1), ng_m, lo, hi)
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

