module itsol_module

  use bl_constants_module
  use bl_types
  use bc_functions_module
  use multifab_module
  use cc_stencil_module
  use stencil_types_module
  use stencil_defect_module
  use stencil_util_module, only : stencil_multifab_copy, is_ibc_stencil

  implicit none

  integer, private, parameter :: def_bicg_max_iter = 1000
  integer, private, parameter :: def_cg_max_iter   = 1000

  private
  public :: itsol_bicgstab_solve, itsol_converged, itsol_cg_solve, itsol_cabicgstab_solve

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

    subroutine jacobi_precon_ibc_2d(a0, u, r, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)    :: a0
      real(kind=dp_t), intent(inout) :: u(1-ng:,1-ng:)
      real(kind=dp_t), intent(in)    :: r(:,:)
      integer :: i, j, nx, ny
      nx = size(r,dim=1)
      ny = size(r,dim=2)
      do j = 1, ny
         do i = 1, nx
            u(i,j) = r(i,j)*(one/a0)
         end do
      end do
    end subroutine jacobi_precon_ibc_2d

    subroutine jacobi_precon_3d(a, u, r, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)    :: a(0:,:,:,:)
      real(kind=dp_t), intent(inout) :: u(1-ng:,1-ng:,1-ng:)
      real(kind=dp_t), intent(in)    :: r(:,:,:)
      integer i, j, k, nx, ny, nz
      nz = size(a,dim=4)
      ny = size(a,dim=3)
      nx = size(a,dim=2)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               u(i,j,k) = r(i,j,k)/a(0,i,j,k)
            end do
         end do
      end do
    end subroutine jacobi_precon_3d

    subroutine jacobi_precon_ibc_3d(a0, u, r, ng)
      integer, intent(in) :: ng
      real(kind=dp_t), intent(in)    :: a0
      real(kind=dp_t), intent(inout) :: u(1-ng:,1-ng:,1-ng:)
      real(kind=dp_t), intent(in)    :: r(:,:,:)
      integer i, j, k, nx, ny, nz
      nx = size(r,dim=1)
      ny = size(r,dim=2)
      nz = size(r,dim=3)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               u(i,j,k) = r(i,j,k)*(one/a0)
            end do
         end do
      end do
    end subroutine jacobi_precon_ibc_3d

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
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (.not. bc_dirichlet(mm(i,j,k),1,0)) then
                  u(i,j,k) = r(i,j,k)/a(0,i,j,k)
               end if
            end do
         end do
      end do
    end subroutine nodal_precon_3d

    subroutine diag_init_cc_1d(a, ng_a, r, ng_r, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(0:,lo(1)-ng_a:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:   )

      integer         :: i, nc
      real(kind=dp_t) :: denom

      nc = size(a,dim=1)-1
      !
      ! Protect against divide by zero -- necessary for embedded boundary problems.
      !
      do i = lo(1),hi(1)
         if (abs(a(0,i)) .gt. zero) then
            denom = one / a(0,i)
            r(i     ) = r(i     ) * denom
            a(1:nc,i) = a(1:nc,i) * denom
            a(0,i   ) = one
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
      !
      ! Protect against divide by zero -- necessary for embedded boundary problems.
      !
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            if (abs(a(0,i,j)) .gt. zero) then
               denom = one / a(0,i,j)
               r(i,j     ) = r(i,j     ) * denom
               a(1:nc,i,j) = a(1:nc,i,j) * denom
               a(0,i,j   ) = one
            end if
         end do
      end do

    end subroutine diag_init_cc_2d

    subroutine diag_init_ibc_2d(a, r, ng_r, lo, hi)
      integer        , intent(in   )  :: ng_r
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(0:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:,lo(2)-ng_r:   )

      integer         :: i, j
      real(kind=dp_t) :: denom

      denom = one / a(0)
      a(0) = one
      a(1:) = a(1:) * denom

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            r(i,j) = r(i,j) * denom
         end do
      end do
    end subroutine diag_init_ibc_2d

    subroutine diag_init_cc_3d(a, ng_a, r, ng_r, lo, hi)
      integer        , intent(in   )  :: ng_a, ng_r
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(0:,lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:   )

      integer         :: i, j, k, nc
      real(kind=dp_t) :: denom

      nc = size(a,dim=1)-1
      !
      ! Protect against divide by zero -- necessary for embedded boundary problems.
      !
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (abs(a(0,i,j,k)) .gt. zero) then
                  denom = one / a(0,i,j,k)
                  r(i,j,k     ) = r(i,j,k     ) * denom
                  a(1:nc,i,j,k) = a(1:nc,i,j,k) * denom
                  a(0,i,j,k   ) = one
               end if
            end do
         end do
      end do

    end subroutine diag_init_cc_3d

    subroutine diag_init_ibc_3d(a, r, ng_r, lo, hi)
      integer        , intent(in   )  :: ng_r
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  ::  a(0:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:   )

      integer         :: i, j, k
      real(kind=dp_t) :: denom

      denom = one / a(0)
      a(0) = one
      a(1:) = a(1:) * denom

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               r(i,j,k) = r(i,j,k) * denom
            end do
         end do
      end do
    end subroutine diag_init_ibc_3d

    subroutine diag_init_nd_1d(sg, ng_sg, r, ng_r, mm, ng_m, lo, hi)
      integer        , intent(in   )  :: ng_sg, ng_r, ng_m
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  :: sg(lo(1)-ng_sg:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r :)
      integer        , intent(inout)  :: mm(lo(1)-ng_m :)

      integer         :: i
      real(kind=dp_t) :: ss0

      do i = lo(1),hi(1)+1
         if (.not. bc_dirichlet(mm(i),1,0)) then
              ss0 = -(sg(i)+sg(i-1))
             r(i) =  r(i) / ss0
         end if
      end do

    end subroutine diag_init_nd_1d

    subroutine diag_init_nd_2d(sg, ng_sg, r, ng_r, mm, ng_m, lo, hi, stencil_type)
      integer        , intent(in   )  :: ng_sg, ng_r, ng_m
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  :: sg(lo(1)-ng_sg:,lo(2)-ng_sg:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r :,lo(2)-ng_r :)
      integer        , intent(inout)  :: mm(lo(1)-ng_m :,lo(2)-ng_m :)
      integer        , intent(in   )  :: stencil_type

      integer         :: i, j
      real(kind=dp_t) :: ss0
      !
      ! NOTE NOTE : we only diagonalize the RHS here --
      !             we will diagnoalize the matrix in the apply routine itself
      !

      if (stencil_type .eq. ND_CROSS_STENCIL) then

         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)+1
               if (.not. bc_dirichlet(mm(i,j),1,0)) then
                  ss0 = -(sg(i-1,j-1)+sg(i,j-1)+sg(i-1,j)+sg(i,j))
                  r(i,j) =  r(i,j) / ss0
               end if
            end do
         end do

      else if (stencil_type .eq. ND_DENSE_STENCIL) then

         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)+1
               if (.not. bc_dirichlet(mm(i,j),1,0)) then
                  ss0 = -TWO3RD*(sg(i-1,j-1)+sg(i,j-1)+sg(i-1,j)+sg(i,j))
                  r(i,j) =  r(i,j) / ss0
               end if
            end do
         end do

      else if (stencil_type .eq. ND_VATER_STENCIL) then

         do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)+1
               if (.not. bc_dirichlet(mm(i,j),1,0)) then
                  ss0 = -THREE4TH*(sg(i-1,j-1)+sg(i,j-1)+sg(i-1,j)+sg(i,j))
                  r(i,j) =  r(i,j) / ss0
               end if
            end do
         end do

     else
        call bl_error("diag_init_nd_2d: dont know this stencil_type")
     end if

    end subroutine diag_init_nd_2d

    subroutine diag_init_nd_3d(sg, ng_sg, r, ng_r, mm, ng_m, lo, hi, stencil_type)
      integer        , intent(in   )  :: ng_sg, ng_r, ng_m
      integer        , intent(in   )  :: lo(:),hi(:)
      real(kind=dp_t), intent(inout)  :: sg(lo(1)-ng_sg:,lo(2)-ng_sg:,lo(3)-ng_sg:)
      real(kind=dp_t), intent(inout)  ::  r(lo(1)-ng_r :,lo(2)-ng_r :,lo(3)-ng_r :)
      integer        , intent(inout)  :: mm(lo(1)-ng_m :,lo(2)-ng_m :,lo(3)-ng_m :)
      integer        , intent(in   )  :: stencil_type

      integer         :: i, j, k
      real(kind=dp_t) :: ss0
      !
      ! NOTE NOTE : we only diagonalize the RHS here --
      !             we will diagnolize the matrix in the apply routine itself
      !

      if (stencil_type .eq. ND_CROSS_STENCIL) then

         do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)+1
         do i = lo(1),hi(1)+1
            if (.not. bc_dirichlet(mm(i,j,k),1,0)) then
               ss0 = -( sg(i-1,j-1,k-1) + sg(i-1,j  ,k-1) &
                       +sg(i  ,j-1,k-1) + sg(i  ,j  ,k-1) &
                       +sg(i-1,j-1,k  ) + sg(i-1,j  ,k  ) &
                       +sg(i  ,j-1,k  ) + sg(i  ,j  ,k  )) * THREE
               r(i,j,k) =  r(i,j,k) / ss0
            end if
         end do
         end do
         end do

      else if (stencil_type .eq. ND_DENSE_STENCIL) then

         do k = lo(3),hi(3)+1
         do j = lo(2),hi(2)+1
         do i = lo(1),hi(1)+1
            if (.not. bc_dirichlet(mm(i,j,k),1,0)) then
                   ! The extra factor of 4 accounts for the fact that fac = 1/(4*dx*dx) 
                   !      to be compatible with the cross stencil
               ss0 =  -( sg(i-1,j-1,k-1) + sg(i,j-1,k-1) &
                        +sg(i-1,j  ,k-1) + sg(i,j  ,k-1) &
                        +sg(i-1,j-1,k  ) + sg(i,j-1,k  ) &
                        +sg(i-1,j  ,k  ) + sg(i,j  ,k  ) ) * FOUR3RD
               r(i,j,k) =  r(i,j,k) / ss0
            end if
         end do
         end do
         end do

     else if (stencil_type .eq. ND_VATER_STENCIL) then
        call bl_error("diag_init_nd_3d: ND_VATER_STENCIL not implemented in 3-d")

     else
        call bl_error("diag_init_nd_3d: dont know this stencil_type")
     end if

    end subroutine diag_init_nd_3d

  function itsol_converged(rr, bnorm, eps, abs_eps, rrnorm, comm) result(r)
    use bl_prof_module

    type(multifab), intent(in )           :: rr
    real(dp_t),     intent(in )           :: bnorm, eps
    real(dp_t),     intent(in ), optional :: abs_eps
    real(dp_t),     intent(out), optional :: rrnorm
    integer,        intent(in ), optional :: comm

    real(dp_t) :: norm_rr
    logical    :: r

    if ( present(comm) ) then
       norm_rr = norm_inf(rr,comm=comm)
    else
       norm_rr = norm_inf(rr)
    end if

    if ( present(rrnorm) ) rrnorm = norm_rr

    if ( present(abs_eps) ) then
      r = (norm_rr <= eps*bnorm) .or. (norm_rr <= abs_eps)
    else
      r = (norm_rr <= eps*bnorm) 
    endif

  end function itsol_converged

  subroutine itsol_bicgstab_solve(aa, uu, rh, mm, eps, max_iter, verbose, stencil_type, lcross, &
                                  stat, singular_in, uniform_dh, nodal_mask, comm_in)

    use bl_prof_module

    integer,         intent(in   ) :: max_iter
    type(imultifab), intent(in   ) :: mm
    type(multifab),  intent(inout) :: uu
    type(multifab),  intent(in   ) :: rh
    type(multifab),  intent(in   ) :: aa
    integer,         intent(in   ) :: stencil_type, verbose
    logical,         intent(in   ) :: lcross
    real(kind=dp_t), intent(in   ) :: eps

    integer,         intent(out), optional :: stat
    logical,         intent(in ), optional :: singular_in
    logical,         intent(in ), optional :: uniform_dh
    type(multifab),  intent(in ), optional :: nodal_mask
    integer,         intent(in ), optional :: comm_in

    type(layout)    :: la
    type(multifab)  :: rr, rt, pp, ph, vv, tt, ss, rh_local, aa_local
    real(kind=dp_t) :: rho_1, alpha, beta, omega, rho, bnorm, rnorm, den
    real(dp_t)      :: tres0, tnorms(2),rtnorms(2)
    integer         :: i, cnt, ng_for_res, comm
    logical         :: nodal_solve, singular, nodal(get_dim(rh)), ioproc

    type(bl_prof_timer), save :: bpt

    if ( present(stat) ) stat = 0

    singular    = .false. ; if ( present(singular_in) ) singular    = singular_in
    ng_for_res  = 0       ; if ( nodal_q(rh)          ) ng_for_res  = 1
    nodal_solve = .false. ; if ( ng_for_res /= 0      ) nodal_solve = .true.

    comm = parallel_communicator() ; if ( present(comm_in) ) comm = comm_in

    call build(bpt, "its_BiCGStab_solve")

    if ( comm == parallel_null_communicator() ) then
       call destroy(bpt)
       return
    end if

    la     = get_layout(aa)
    nodal  = nodal_flags(rh)
    ioproc = parallel_IOProcessor(comm = comm)

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
    ! Copy aa -> aa_local; gotta to call the special copy
    !
    call stencil_multifab_copy(aa_local, aa)
    !
    ! Make sure to do singular adjustment *before* diagonalization.
    !
    if ( singular ) then
      call setval(ss,ONE)
      tnorms(1) = dot(rh_local, ss, nodal_mask, local = .true.)
      tnorms(2) = dot(      ss, ss, nodal_mask, local = .true.)
      call parallel_reduce(rtnorms, tnorms, MPI_SUM, comm = comm)
      rho = rtnorms(1) / rtnorms(2)
      if ( ioproc .and. verbose > 0 ) then
         print *,'   ...singular adjustment to rhs: ', rho
      endif
      call sub_sub(rh_local, rho)
      call setval(ss,ZERO,all=.true.)
    end if

    call diag_initialize(aa_local,rh_local,mm,stencil_type)

    call copy(ph, uu, ng = nghost(ph))

    cnt = 0
    !
    ! Compute rr = aa_local * uu - rh_local
    !
    if ( cell_centered_q(rr) ) then
        call compute_defect(aa_local, rr, rh_local, uu, mm, stencil_type, lcross, &
                            uniform_dh, bottom_solver=.true.) 
    else
        call compute_defect(aa_local, rr, rh_local, uu, mm, stencil_type, lcross, &
                            uniform_dh, bottom_solver=.true.,diagonalize=.true.) 
    end if
    cnt = cnt + 1

    call copy(rt, rr)

    rho = dot(rt, rr, nodal_mask, comm = comm)
    !
    ! Elide some reductions by calculating local norms & then reducing all together.
    !
    tnorms(1) = norm_inf(rr,       local = .true.) 
    tnorms(2) = norm_inf(rh_local, local = .true.)

    call parallel_reduce(rtnorms, tnorms, MPI_MAX, comm = comm)

    tres0 = rtnorms(1)
    bnorm = rtnorms(2)

    if ( ioproc .and. verbose > 0 ) then
       write(*,*) "   BiCGStab: A and rhs have been rescaled. So has the error."
       write(unit=*, fmt='("    BiCGStab: Initial error (error0) =        ",g15.8)') tres0
    end if 

    if ( itsol_converged(rr, bnorm, eps, comm = comm) ) then
       if ( ioproc .and. verbose > 0 ) then
          if ( tres0 < eps*bnorm ) then
             write(unit=*, fmt='("    BiCGStab: Zero iterations: rnorm ",g15.8," < eps*bnorm ",g15.8)') tres0,eps*bnorm
          end if
       end if
       go to 100
    end if

    rho_1 = ZERO

    do i = 1, max_iter
       rho = dot(rt, rr, nodal_mask, comm = comm)
       if ( i == 1 ) then
          call copy(pp, rr)
       else
          if ( rho_1 == ZERO ) then
             if ( present(stat) ) then
                call bl_warn("BiCGStab_SOLVE: failure 1"); stat = 2; goto 100
             end if
             call bl_error("BiCGStab: failure 1")
          end if
          if ( omega == ZERO ) then
             if ( present(stat) ) then
                call bl_warn("BiCGStab_SOLVE: failure 2"); stat = 3; goto 100
             end if
             call bl_error("BiCGStab: failure 2")
          end if
          beta = (rho/rho_1)*(alpha/omega)
          call saxpy(pp, -omega, vv)
          call rescale(pp, beta)
          call plus_plus(pp, rr)
       end if
       call copy(ph,pp)
      
       if ( cell_centered_q(ph) ) then
           call stencil_apply(aa_local, vv, ph, mm, stencil_type, lcross, uniform_dh, bottom_solver=.true.)
       else
           call stencil_apply(aa_local, vv, ph, mm, stencil_type, lcross, uniform_dh, bottom_solver=.true., &
                              diagonalize=.true.)
       end if

       cnt = cnt + 1
       den = dot(rt, vv, nodal_mask, comm = comm)
       if ( den == ZERO ) then
          if ( present(stat) ) then
             call bl_warn("BICGSTAB_solve: breakdown in bicg, going with what I have"); stat = 30; goto 100
          endif
          call bl_error("BiCGStab: failure 3")
       end if
       alpha = rho/den
       call saxpy(uu, alpha, ph)
       call saxpy(ss, rr, -alpha, vv)
       if ( verbose > 1 ) then
          rnorm = norm_inf(ss, comm = comm)
          if ( ioproc ) then
             write(unit=*, fmt='("    BiCGStab: Half Iter        ",i4," rel. err. ",g15.8)') cnt/2, rnorm/bnorm
          end if
       end if
       if ( itsol_converged(ss, bnorm, eps, rrnorm = rnorm, comm = comm) ) exit
       call copy(ph,ss)
       if ( cell_centered_q(ph) ) then
           call stencil_apply(aa_local, tt, ph, mm, stencil_type, lcross, uniform_dh, bottom_solver=.true.)
       else
           call stencil_apply(aa_local, tt, ph, mm, stencil_type, lcross, uniform_dh, bottom_solver=.true., &
                              diagonalize=.true.)
       end if
       cnt = cnt + 1
       !
       ! Elide a reduction here by calculating the two dot-products
       ! locally and then reducing them both in a single call.
       !
       tnorms(1) = dot(tt, tt, nodal_mask, local = .true.)
       tnorms(2) = dot(tt, ss, nodal_mask, local = .true.)

       call parallel_reduce(rtnorms, tnorms, MPI_SUM, comm = comm)

       den   = rtnorms(1)
       omega = rtnorms(2)

       if ( den == ZERO ) then
          if ( present(stat) ) then
             call bl_warn("BICGSTAB_solve: breakdown in bicg, going with what I have"); stat = 31; goto 100
          endif
          call bl_error("BiCGStab: failure 3")
       end if
       omega = omega/den
       call saxpy(uu, omega, ph)
       call saxpy(rr, ss, -omega, tt)
       if ( verbose > 1 ) then
          rnorm = norm_inf(rr, comm = comm)
          if ( ioproc ) then
             write(unit=*, fmt='("    BiCGStab: Iteration        ",i4," rel. err. ",g15.8)') cnt/2, rnorm/bnorm
          end if
       end if
       if ( itsol_converged(rr, bnorm, eps, rrnorm = rnorm, comm = comm) ) exit
       rho_1 = rho
    end do

    if ( ioproc .and. verbose > 0 ) then
       write(unit=*, fmt='("    BiCGStab: Final: Iteration  ", i3, " rel. err. ",g15.8)') cnt/2, rnorm/bnorm
       if ( rnorm < eps*bnorm ) then
          write(unit=*, fmt='("    BiCGStab: Converged: rnorm ",g15.8," < eps*bnorm ",g15.8)') rnorm,eps*bnorm
       end if
    end if

    if ( rnorm > bnorm ) then
       call setval(uu,ZERO,all=.true.)
       if ( present(stat) ) stat = 1
       if ( verbose > 0 .and. ioproc ) print *,'   BiCGStab: solution reset to zero'
    end if

    if ( i > max_iter ) then
       if ( present(stat) ) then
          stat = 1
       else
          call bl_error("BiCGSolve: failed to converge");
       end if
    end if

100 continue

    call multifab_destroy(rh_local)
    call multifab_destroy(aa_local)
    call multifab_destroy(rr)
    call multifab_destroy(rt)
    call multifab_destroy(pp)
    call multifab_destroy(ph)
    call multifab_destroy(vv)
    call multifab_destroy(tt)
    call multifab_destroy(ss)

    call destroy(bpt)

  end subroutine itsol_bicgstab_solve
  !
  ! This is a slightly simplified version of the BLAS 2 routine.
  !
  subroutine dgemv(alpha,a,x,beta,y,m,n)

    integer,    intent(in   ) :: m,n
    real(dp_t), intent(in   ) :: a(m,n),x(n),alpha,beta
    real(dp_t), intent(inout) :: y(m)
    !
    !  dgemv  performs 
    !
    !     y := alpha*A*x + beta*y
    !
    !  where alpha and beta are scalars, x and y are vectors and A is an  m x n matrix.
    !
    !  The vector and matrix arguments are not referenced when N = 0, or M = 0
    !
    !  -- Written on 22-October-1986.
    !     Jack Dongarra, Argonne National Lab.
    !     Jeremy Du Croz, Nag Central Office.
    !     Sven Hammarling, Nag Central Office.
    !     Richard Hanson, Sandia National Labs.
    !
    integer    :: i,j
    real(dp_t) :: temp
    !
    ! Quick return if possible.
    !
    if ( (m.eq.0) .or. (n.eq.0) .or. ((alpha.eq.zero).and.(beta.eq.one)) ) return
    !
    ! Start the operations. In this version the elements of A are
    ! accessed sequentially with one pass through A.
    !
    ! First form y := beta*y.
    !
    if ( beta.ne.one ) then
       if ( beta.eq.zero ) then
          do i = 1,m
             y(i) = zero
          end do
       else
          do i = 1,m
             y(i) = beta*y(i)
          end do
       end if
    end if
    if ( alpha.eq.zero ) return
    !
    ! Now form y := alpha*a*x + y.
    !
    do j = 1,n
       if ( x(j).ne.zero ) then
          temp = alpha*x(j)
          do i = 1,m
             y(i) = y(i) + temp*a(i,j)
          end do
       end if
    end do

  end subroutine dgemv

  subroutine itsol_CABiCGStab_solve(aa, uu, rh, mm, eps, max_iter, verbose, stencil_type, lcross, &
       stat, singular_in, uniform_dh, nodal_mask, comm_in)
    use bl_prof_module
    integer,         intent(in   ) :: max_iter
    type(imultifab), intent(in   ) :: mm
    type(multifab),  intent(inout) :: uu
    type(multifab),  intent(in   ) :: rh
    type(multifab),  intent(in   ) :: aa
    integer        , intent(in   ) :: stencil_type, verbose
    logical        , intent(in   ) :: lcross
    real(kind=dp_t), intent(in   ) :: eps

    integer,         intent(out), optional :: stat
    logical,         intent(in ), optional :: singular_in
    logical,         intent(in ), optional :: uniform_dh
    type(multifab),  intent(in ), optional :: nodal_mask
    integer,         intent(in ), optional :: comm_in

    type(layout)    :: la
    type(multifab)  :: rr, rt, pp, pr, ss, rh_local, aa_local, ph, tt
    real(kind=dp_t) :: alpha, beta, omega, rho, bnorm
    real(dp_t)      :: rnorm0, delta, delta_next, L2_norm_of_rt
    real(dp_t)      :: nrms(2),rnrms(2), L2_norm_of_resid, L2_norm_of_r
    integer         :: i, m, niters, ng_for_res, nit, comm
    logical         :: nodal_solve, singular, nodal(get_dim(rh))
    logical         :: BiCGStabFailed, BiCGStabConverged, ioproc
    real(dp_t)      :: g_dot_Tpaj, omega_numerator, omega_denominator, L2_norm_of_s

    type(bl_prof_timer), save :: bpt

    integer, parameter :: SSS = 4

    real(dp_t)  temp1(4*SSS+1)
    real(dp_t)  temp2(4*SSS+1)
    real(dp_t)  temp3(4*SSS+1)
    real(dp_t)     Tp(4*SSS+1, 4*SSS+1)
    real(dp_t)    Tpp(4*SSS+1, 4*SSS+1)
    real(dp_t)     aj(4*SSS+1)
    real(dp_t)     cj(4*SSS+1)
    real(dp_t)     ej(4*SSS+1)
    real(dp_t)   Tpaj(4*SSS+1)
    real(dp_t)   Tpcj(4*SSS+1)
    real(dp_t)  Tppaj(4*SSS+1)
    real(dp_t)      G(4*SSS+1, 4*SSS+1)
    real(dp_t)     gg(4*SSS+1)

    if ( present(stat) ) stat = 0

    singular    = .false. ; if ( present(singular_in) ) singular    = singular_in
    ng_for_res  = 0       ; if ( nodal_q(rh)          ) ng_for_res  = 1
    nodal_solve = .false. ; if ( ng_for_res /= 0      ) nodal_solve = .true.

    comm = parallel_communicator() ; if ( present(comm_in) ) comm = comm_in

    call build(bpt, "its_CABiCGStab_solve")

    if ( comm == parallel_null_communicator() ) then
       call destroy(bpt)
       return
    end if

    la     = get_layout(aa)
    nodal  = nodal_flags(rh)
    ioproc = parallel_IOProcessor(comm = comm)

    aj    = zero
    cj    = zero
    ej    = zero
    Tpaj  = zero
    Tpcj  = zero
    Tppaj = zero
    temp1 = zero
    temp2 = zero
    temp3 = zero

    call SetMonomialBasis()

    call multifab_build(rr, la, 1, ng_for_res, nodal)
    call multifab_build(rt, la, 1, ng_for_res, nodal)
    call multifab_build(pp, la, 1, ng_for_res, nodal)
    call multifab_build(tt, la, 2, ng_for_res, nodal)
    call multifab_build(ph, la, 2, nghost(uu), nodal)
    !
    ! Contains the matrix powers of pp[] and rr[].
    !
    ! First 2*SSS+1 components are powers of pp[].
    ! Next  2*SSS   components are powers of rr[].
    !
    call multifab_build(PR, la, 4*SSS+1, 0, nodal)
    !
    ! Use these for local preconditioning.
    !
    call multifab_build(rh_local, la, ncomp(rh), nghost(rh), nodal)
    call multifab_build(aa_local, la, ncomp(aa), nghost(aa), nodal_flags(aa), stencil = .true.)

    call copy(rh_local, 1, rh, 1, nc = ncomp(rh), ng = nghost(rh))

    call copy(ph, 1, uu, 1, 1, ng = nghost(ph))
    call copy(ph, 2, uu, 1, 1, ng = nghost(ph))
    !
    ! Copy aa -> aa_local; gotta to call the special copy
    !
    call stencil_multifab_copy(aa_local, aa)
    !
    ! Make sure to do singular adjustment *before* diagonalization.
    !
    if ( singular ) then
       call multifab_build(ss, la, 1, ng_for_res, nodal)
       call setval(ss,ONE)
       nrms(1) = dot(rh_local, ss, nodal_mask, local = .true.)
       nrms(2) = dot(      ss, ss, nodal_mask, local = .true.)
       call parallel_reduce(rnrms, nrms, MPI_SUM, comm = comm)
       rho = rnrms(1) / rnrms(2)
       if ( ioproc .and. verbose > 0 ) then
          print *,'   ...singular adjustment to rhs: ', rho
       endif
       call sub_sub(rh_local,rho)
       call multifab_destroy(ss)
    end if

    call diag_initialize(aa_local,rh_local,mm,stencil_type)

    !if (contains_nan(ph)) then; print*, '*** Got NaNs @ 1'; stop; endif

    !
    ! Compute rr = aa_local * uu - rh_local
    !
    if ( cell_centered_q(rt) ) then
        call compute_defect(aa_local, rt, rh_local, uu, mm, stencil_type, lcross, &
                            uniform_dh, bottom_solver=.true.)
    else
        call compute_defect(aa_local, rt, rh_local, uu, mm, stencil_type, lcross, &
                            uniform_dh, bottom_solver=.true.,diagonalize=.true.)
    end if

    call copy(rr,rt); call copy(pp,rt)
    !
    ! Elide some reductions by calculating local norms & then reducing all together.
    !
    nrms(1) = norm_inf(rt,       local = .true.)
    nrms(2) = norm_inf(rh_local, local = .true.)

    call parallel_reduce(rnrms, nrms, MPI_MAX, comm = comm)

    rnorm0 = rnrms(1)
    bnorm  = rnrms(2)

    delta         = dot(rt,rt,nodal_mask, comm = comm)
    L2_norm_of_rt = sqrt(delta)

    if ( ioproc .and. verbose > 0 ) then
       write(*,*) "   CABiCGStab: A and rhs have been rescaled. So has the error."
       write(unit=*, fmt='("    CABiCGStab: Initial error (error0) =        ",g15.8)') rnorm0
    end if 

    if ( itsol_converged(rr, bnorm, eps, comm = comm) .or. (delta.eq.zero) ) then
       if ( ioproc .and. verbose > 0 ) then
          if ( rnorm0 < eps*bnorm ) then
             write(unit=*, fmt='("    CABiCGStab: Zero iterations: rnorm ",g15.8," < eps*bnorm ",g15.8)') rnorm0,eps*bnorm
          else if ( delta .eq. zero ) then
             write(unit=*, fmt='("    CABiCGStab: Zero iterations: delta == 0")')
          end if
       end if
       go to 100
    end if

    L2_norm_of_resid = 0

    BiCGStabFailed = .false. ; BiCGStabConverged = .false.

    niters = 0; m = 1

    do while (m <= max_iter .and. (.not. BiCGStabFailed) .and. (.not. BiCGStabConverged))

       !if (contains_nan(ph)) then; print*, '*** Got NaNs @ 2'; stop; endif

       !
       ! Compute the matrix powers on pp[] & rr[] (monomial basis).
       ! The 2*SSS+1 powers of pp[] followed by the 2*SSS powers of rr[].
       !
       call copy(PR,1,pp,1,1,0)
       call copy(ph,1,pp,1,1,0)

       call copy(PR,2*SSS+2,rr,1,1,0)
       call copy(ph,2,      rr,1,1,0)

       do i = 2, 2*SSS
          !
          ! apply the stencil to pp & rr at the same time to cut down on comm time.
          !
          if ( cell_centered_q(ph) ) then
              call stencil_apply(aa_local, tt, ph, mm, stencil_type, lcross, uniform_dh, bottom_solver=.true.)
          else
              call stencil_apply(aa_local, tt, ph, mm, stencil_type, lcross, uniform_dh, bottom_solver=.true., &
                                 diagonalize=.true.)
          end if

          !if (contains_nan(ph)) then; print*, '*** Got NaNs @ 3'; stop; endif

          call copy(PR,i,tt,1,1,0)
          call copy(ph,1,tt,1,1,0)

          call copy(PR,2*SSS+1+i,tt,2,1,0)
          call copy(ph,2,        tt,2,1,0)
       end do
       !
       ! And the final power of pp[].
       !
       if ( cell_centered_q(ph) ) then
          call stencil_apply(aa_local, tt, ph, mm, stencil_type, lcross, uniform_dh, bottom_solver=.true.)
       else
          call stencil_apply(aa_local, tt, ph, mm, stencil_type, lcross, uniform_dh, bottom_solver=.true., &
                             diagonalize=.true.)
       end if

       call copy(PR,2*SSS+1,tt,1,1,0)

       !if (contains_nan(ph)) then; print*, '*** Got NaNs @ 4'; stop; endif

       call BuildGramMatrix()

       aj = 0; aj(1)       = 1
       cj = 0; cj(2*SSS+2) = 1
       ej = 0

       do nit = 1, SSS
          call dgemv(one,  Tp, aj, zero,  Tpaj, 4*SSS+1, 4*SSS+1)
          call dgemv(one,  Tp, cj, zero,  Tpcj, 4*SSS+1, 4*SSS+1)
          call dgemv(one, Tpp, aj, zero, Tppaj, 4*SSS+1, 4*SSS+1)

          g_dot_Tpaj = dot_product(gg,Tpaj)

          if ( g_dot_Tpaj == zero ) then
             if ( ioproc .and. verbose > 0 ) &
                  print*, "CGSolver_CABiCGStab: g_dot_Tpaj == 0, nit = ", nit
             BiCGStabFailed = .true. ; exit
          end if

          alpha = delta / g_dot_Tpaj

          if ( is_an_inf(alpha) ) then
             if ( verbose > 1 .and. ioproc ) &
                  print*, "CGSolver_CABiCGStab: alpha == inf, nit = ", nit
             BiCGStabFailed = .true. ; exit
          end if

          temp1 = Tpcj - alpha * Tppaj
          call dgemv(one, G, temp1, zero, temp2, 4*SSS+1, 4*SSS+1)
          temp3 = cj - alpha * Tpaj

          omega_numerator   = dot_product(temp3, temp2)
          omega_denominator = dot_product(temp1, temp2)
          !
          ! NOTE: omega_numerator/omega_denominator can be 0/x or 0/0, but should never be x/0.
          !
          ! If omega_numerator==0, and ||s||==0, then convergence, x=x+alpha*aj.
          ! If omega_numerator==0, and ||s||!=0, then stabilization breakdown.
          !
          ! Partial update of ej must happen before the check on omega to ensure forward progress !!!
          !
          ej = ej + alpha * aj
          !
          ! ej has been updated so consider that we've done an iteration since
          ! even if we break out of the loop we'll be able to update "uu".
          !
          niters = niters + 1
          !
          ! Calculate the norm of Saad's vector 's' to check intra s-step convergence.
          !
          temp1 = cj - alpha * Tpaj

          call dgemv(one, G, temp1, zero, temp2, 4*SSS+1, 4*SSS+1)

          L2_norm_of_s = dot_product(temp1,temp2)

          L2_norm_of_resid = zero; if ( L2_norm_of_s > 0 ) L2_norm_of_resid = sqrt(L2_norm_of_s)

          if ( L2_norm_of_resid < eps*L2_norm_of_rt ) then
             if ( verbose > 1 .and. (L2_norm_of_resid .eq. zero) .and. ioproc ) &
                  print*, "CGSolver_CABiCGStab: L2 norm of s: ", L2_norm_of_s
             BiCGStabConverged = .true. ; exit
          end if

          if ( omega_denominator .eq. zero ) then
             if ( verbose > 1 .and. ioproc ) &
                print*, "CGSolver_CABiCGStab: omega_denominator == 0, nit = ", nit
             BiCGStabFailed = .true. ; exit
          end if

          omega = omega_numerator / omega_denominator

          if ( verbose > 1 .and. ioproc ) then
             if ( omega .eq. zero  ) print*, "CGSolver_CABiCGStab: omega == 0, nit = ", nit
             if ( is_an_inf(omega) ) print*, "CGSolver_CABiCGStab: omega == inf, nit = ", nit
          end if

          if ( omega .eq. zero ) then
             BiCGStabFailed = .true. ; exit
          end if
          if ( is_an_inf(omega) ) then
             BiCGStabFailed = .true. ; exit
          end if
          !
          ! Complete the update of ej & cj now that omega is known to be ok.
          !
          ej = ej +  omega          * cj
          ej = ej - (omega * alpha) * Tpaj
          cj = cj -  omega          * Tpcj
          cj = cj -          alpha  * Tpaj
          cj = cj + (omega * alpha) * Tppaj
          !
          ! Do an early check of the residual to determine convergence.
          !
          call dgemv(one, G, cj, zero, temp1, 4*SSS+1, 4*SSS+1)
          !
          ! sqrt( (cj,Gcj) ) == L2 norm of the intermediate residual in exact arithmetic.
          ! However, finite precision can lead to the norm^2 being < 0 (Jim Demmel).
          ! If cj_dot_Gcj < 0 we flush to zero and consider ourselves converged.
          !
          L2_norm_of_r = dot_product(cj,temp1)

          L2_norm_of_resid = zero; if ( L2_norm_of_r > 0 ) L2_norm_of_resid = sqrt(L2_norm_of_r)

          if ( L2_norm_of_resid < eps*L2_norm_of_rt ) then
             if ( verbose > 1 .and. (L2_norm_of_resid .eq. zero) .and. ioproc ) &
                  print*, "CGSolver_CABiCGStab: L2_norm_of_r: ", L2_norm_of_r
             BiCGStabConverged = .true. ; exit
          end if

          delta_next = dot_product(gg,cj)

          if ( verbose > 1 .and. ioproc ) then
             if ( delta_next .eq. zero  ) print*, "CGSolver_CABiCGStab: delta == 0, nit = ", nit
             if ( is_an_inf(delta_next) ) print*, "CGSolver_CABiCGStab: delta == inf, nit = ", nit
          end if
          if ( delta_next .eq. zero ) then
             BiCGStabFailed = .true. ; exit
          end if
          if ( is_an_inf(delta_next) ) then
             BiCGStabFailed = .true. ; exit
          end if

          beta = (delta_next/delta)*(alpha/omega)

          if ( verbose > 1 .and. ioproc ) then
             if ( beta .eq. zero  ) print*, "CGSolver_CABiCGStab: beta == 0, nit = ", nit
             if ( is_an_inf(beta) ) print*, "CGSolver_CABiCGStab: beta == inf, nit = ", nit
          end if
          if ( beta .eq. zero ) then
             BiCGStabFailed = .true. ; exit
          end if
          if ( is_an_inf(beta) ) then
             BiCGStabFailed = .true. ; exit
          end if

          aj = cj +          beta  * aj
          aj = aj - (omega * beta) * Tpaj

          delta = delta_next
       end do
       !
       ! Update iterates.
       !
       do i = 1,4*SSS+1
          call saxpy(uu,1,ej(i),PR,i,1)
       end do

       call copy(pp,1,PR,1,1)
       call mult_mult(pp,aj(1))

       do i = 2,4*SSS+1
          call saxpy(pp,1,aj(i),PR,i,1)
       end do

       call copy(rr,1,PR,1,1)
       call mult_mult(rr,cj(1))

       do i = 2,4*SSS+1
          call saxpy(rr,1,cj(i),PR,i,1)
       end do

       if ( (.not. BiCGStabFailed) .and. (.not. BiCGStabConverged) ) m = m + SSS
    end do

    if ( ioproc .and. verbose > 0 ) then
       write(unit=*, fmt='("    CABiCGStab: Final: Iteration  ", i3, " rel. err. ",g15.8)') niters, L2_norm_of_resid
       if ( BiCGStabConverged ) then
          write(unit=*, fmt='("    CABiCGStab: Converged: rnorm ",g15.8," < eps*bnorm ",g15.8)') &
               L2_norm_of_resid,eps*L2_norm_of_rt
       end if
    end if

    if ( L2_norm_of_resid > L2_norm_of_rt ) then
       call setval(uu,ZERO,all=.true.)
       if ( present(stat) ) stat = 1
       if ( ioproc .and. verbose > 0 ) then
          print *,'   CABiCGStab: solution reset to zero'
       end if
    end if

    if ( m > max_iter .or. BiCGStabFailed ) then
       if ( present(stat) ) then
          stat = 1
       else
          call bl_error("CABiCGSolve: failed to converge");
       end if
    end if

100 continue

    call multifab_destroy(rh_local)
    call multifab_destroy(aa_local)
    call multifab_destroy(rr)
    call multifab_destroy(rt)
    call multifab_destroy(pp)
    call multifab_destroy(tt)
    call multifab_destroy(ph)
    call multifab_destroy(pr)

    call destroy(bpt)

  contains

    subroutine SetMonomialBasis ()

      Tp = zero

      do i = 1,2*SSS
         Tp(i+1,i) = one
      end do
      do i = 2*SSS+2, 4*SSS
         Tp(i+1,i) = one
      end do

      Tpp = zero

      do i = 1,2*SSS-1
         Tpp(i+2,i) = one
      end do
      do i = 2*SSS+2, 4*SSS-1
         Tpp(i+2,i) = one
      end do

    end subroutine SetMonomialBasis

    subroutine BuildGramMatrix ()

      integer, parameter :: Nrows = 4*SSS+1, Ncols = 4*SSS+2

      integer    :: mm, nn, cnt
      real(dp_t) :: Gram(Nrows,Ncols)
      real(dp_t) :: tmp((Ncols*(Ncols+1))/2-1), rtmp((Ncols*(Ncols+1))/2-1)

      !$OMP PARALLEL DO PRIVATE(mm,nn) SCHEDULE(static,1)
      do mm = 1, Nrows
         do nn = mm, Nrows
            Gram(mm,nn) = dot(PR, mm, PR, nn, nodal_mask = nodal_mask, local = .true.)
         end do
         Gram(mm,Ncols) = dot(PR, mm, rt,  1, nodal_mask = nodal_mask, local = .true.)
      end do
      !$OMP END PARALLEL DO
      !
      ! Put upper triangle into "tmp".
      !
      cnt = 1
      do mm = 1, Nrows
         do nn = mm, Nrows
            tmp(cnt) = Gram(mm,nn)
            cnt = cnt + 1
         end do
         tmp(cnt) = Gram(mm,Ncols)
         cnt = cnt + 1
      end do
      !
      ! Reduce upper triangle into "rtmp"
      !
      call parallel_reduce(rtmp, tmp, MPI_SUM, comm = comm)
      !
      ! Fill in upper triangle with "rtmp".
      !
      cnt = 1
      do mm = 1, Nrows
         do nn = mm, Nrows
            Gram(mm,nn) = rtmp(cnt)
            cnt = cnt + 1
         end do
         Gram(mm,Ncols) = rtmp(cnt)
         cnt = cnt + 1
      end do
      !
      ! Then fill in strict lower triangle using symmetry.
      !
      do mm = 1, Nrows
         do nn = 1, mm-1
            Gram(mm,nn) = Gram(nn,mm)
         end do
      end do
      !
      ! Form G[][] and g[] from Gram[][].
      !
      G(1:Nrows,1:Nrows) = Gram(1:Nrows,1:Nrows)
      !
      ! Last column goes to g[].
      !
      gg = Gram(:,Ncols)

    end subroutine BuildGramMatrix

  end subroutine itsol_CABiCGStab_solve

  subroutine itsol_cg_solve(aa, uu, rh, mm, eps, max_iter, verbose, stencil_type, lcross, &
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
    real(kind = dp_t) :: rho_1, alpha, beta, bnorm, rho, rnorm, den, tres0
    type(layout) :: la
    integer :: i, ng_for_res
    logical :: nodal_solve, nodal(get_dim(rh))
    logical :: singular 
    integer :: cnt
    real(dp_t) :: nrms(2), rnrms(2)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "its_CG_Solve")

    if ( present(stat) ) stat = 0

    singular    = .false. ; if ( present(singular_in) ) singular = singular_in
    ng_for_res  = 0       ; if ( nodal_q(rh)          ) ng_for_res = 1
    nodal_solve = .false. ; if ( ng_for_res /= 0      ) nodal_solve = .true.

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
    !
    ! Copy aa -> aa_local; gotta to call the special copy
    !
    call stencil_multifab_copy(aa_local, aa)

    cnt = 0
    !
    ! Compute rr = aa_local * uu - rh_local
    !
    if ( cell_centered_q(rr) ) then
        call compute_defect(aa_local, rr, rh_local, uu, mm, stencil_type, lcross, &
                            uniform_dh, bottom_solver=.true.)
    else 
        call compute_defect(aa_local, rr, rh_local, uu, mm, stencil_type, lcross, &
                            uniform_dh, bottom_solver=.true.,diagonalize=.true.)
    end if
    cnt = cnt + 1

    if ( singular .and. nodal_solve ) then
      call setval(zz,ONE)
      rho = dot(rr, zz, nodal_mask) / dot(zz,zz)
      call sub_sub(rr,rho)
      call setval(zz,ZERO,all=.true.)
    end if
    !
    ! Elide some reductions by calculating local norms & then reducing all together.
    !
    nrms(1) = norm_inf(rr,       local=.true.)
    nrms(2) = norm_inf(rh_local, local=.true.)

    call parallel_reduce(rnrms, nrms, MPI_MAX)

    tres0 = rnrms(1)
    bnorm = rnrms(2)

    if ( parallel_IOProcessor() .and. verbose > 0) then
       write(unit=*, fmt='("          CG: Initial error (error0) =        ",g15.8)') tres0
    end if

    i = 0
    if ( itsol_converged(rr, bnorm, eps) ) then
       if (parallel_IOProcessor() .and. verbose > 0) then
          if (tres0 < eps*bnorm) then
             write(unit=*, fmt='("          CG: Zero iterations: rnorm ",g15.8," < eps*bnorm ",g15.8)') tres0,eps*bnorm
          end if
       end if
       go to 100
    end if

    rho_1 = ZERO

    do i = 1, max_iter
       call copy(zz,rr)
       call itsol_precon(aa_local, zz, rr, mm)
       rho = dot(rr, zz, nodal_mask)
       if ( i == 1 ) then
          call copy(pp, zz)
       else
          if ( rho_1 == ZERO ) then
             if ( present(stat) ) then
                call bl_warn("CG_solve: failure 1"); stat = 1; goto 100
             end if
             call bl_error("CG_solve: failure 1")
          end if
          beta = rho/rho_1
          call rescale(pp, beta)
          call plus_plus(pp, zz)
       end if
       if ( cell_centered_q(pp) ) then
           call stencil_apply(aa_local, qq, pp, mm, stencil_type, lcross, uniform_dh, bottom_solver=.true.)
       else
           call stencil_apply(aa_local, qq, pp, mm, stencil_type, lcross, uniform_dh, bottom_solver=.true., &
                              diagonalize=.true.)
       end if
       cnt = cnt + 1
       den = dot(pp, qq, nodal_mask)
       if ( den == ZERO ) then
          if ( present(stat) ) then
             call bl_warn("CG_solve: breakdown in solver, going with what I have"); stat = 30; goto 100
          end if
          call bl_error("CG_solve: failure 1")
       end if
       alpha = rho/den
       call saxpy(uu,  alpha, pp)
       call saxpy(rr, -alpha, qq)
       if ( verbose > 1 ) then
          rnorm = norm_inf(rr)
          if ( parallel_IOProcessor() ) then
             write(unit=*, fmt='("          CG: Iteration        ",i4," rel. err. ",g15.8)') i,rnorm/bnorm
          end if
       end if
       if ( itsol_converged(rr, bnorm, eps, rrnorm = rnorm) ) exit
       rho_1 = rho
    end do

    if ( parallel_IOProcessor() .and. verbose > 0 ) then
       write(unit=*, fmt='("          CG: Final: Iteration  ", i3, " rel. err. ",g15.8)') i, rnorm/bnorm
       if ( rnorm < eps*bnorm ) then
          write(unit=*, fmt='("          CG: Converged: rnorm ",g15.8," < eps*bnorm ",g15.8)') rnorm,eps*bnorm
       end if
    end if

    if ( i > max_iter ) then
       if ( present(stat) ) then
          stat = 1
       else
          call bl_error("CG_solve: failed to converge");
       end if
    end if

100 continue

    call multifab_destroy(rr)
    call multifab_destroy(zz)
    call multifab_destroy(pp)
    call multifab_destroy(qq)

    call destroy(bpt)

    call multifab_destroy(aa_local)
    call multifab_destroy(rh_local)

  end subroutine itsol_cg_solve

  subroutine itsol_precon(aa, uu, rh, mm, method)
    use bl_prof_module
    type(multifab), intent(in) :: aa
    type(multifab), intent(inout) :: uu
    type(multifab), intent(in) :: rh
    type(imultifab), intent(in) :: mm
    real(kind=dp_t), pointer, dimension(:,:,:,:) :: ap, up, rp
    integer, pointer, dimension(:,:,:,:) :: mp
    integer :: i, n, dm, ngu
    integer, intent(in), optional :: method
    integer :: lm
    type(bl_prof_timer), save :: bpt

    call build(bpt, "its_precon")

    lm = 1; if ( present(method) ) lm = method

    dm = get_dim(uu)
    ngu = nghost(uu)

    select case (lm)
    case (0)
       call copy(uu, rh)
    case (1)
       !$OMP PARALLEL DO PRIVATE(i,n,rp,up,ap,mp)
       do i = 1, nfabs(rh)
          rp => dataptr(rh, i)
          up => dataptr(uu, i)
          ap => dataptr(aa, i)
          mp => dataptr(mm, i)

          do n = 1, ncomp(uu)
             select case(dm)
             case (1)
                if ( cell_centered_q(rh) ) then
                   call jacobi_precon_1d(ap(:,:,1,1), up(:,1,1,n), rp(:,1,1,n), ngu)
                else
                   call nodal_precon_1d(ap(:,:,1,1), up(:,1,1,n), rp(:,1,1,n), &
                                        mp(:,1,1,1),ngu)
                end if
             case (2)
                if (is_ibc_stencil(aa,i)) then
                   call jacobi_precon_ibc_2d(ap(1,1,1,1), up(:,:,1,n), rp(:,:,1,n), ngu)
                else if ( cell_centered_q(rh) ) then
                   call jacobi_precon_2d(ap(:,:,:,1), up(:,:,1,n), rp(:,:,1,n), ngu)
                else
                   call nodal_precon_2d(ap(:,:,:,1), up(:,:,1,n), rp(:,:,1,n), &
                                        mp(:,:,1,1),ngu)
                end if
             case (3)
                if (is_ibc_stencil(aa,i)) then
                   call jacobi_precon_ibc_3d(ap(1,1,1,1), up(:,:,:,n), rp(:,:,:,n), ngu)
                else if ( cell_centered_q(rh) ) then
                   call jacobi_precon_3d(ap(:,:,:,:), up(:,:,:,n), rp(:,:,:,n), ngu)
                else
                   call nodal_precon_3d(ap(:,:,:,:), up(:,:,:,n), rp(:,:,:,n), &
                                        mp(:,:,:,1),ngu)
                end if
             end select
          end do
       end do
       !$OMP END PARALLEL DO
    end select

    call destroy(bpt)

  end subroutine itsol_precon

  subroutine diag_initialize(aa, rh, mm, stencil_type)
    use bl_prof_module
    type( multifab), intent(in) :: aa
    type( multifab), intent(in) :: rh
    type(imultifab), intent(in) :: mm
    integer        , intent(in) :: stencil_type

    real(kind=dp_t), pointer, dimension(:,:,:,:) :: ap, rp
    integer        , pointer, dimension(:,:,:,:) :: mp
    integer                                      :: i,dm
    integer                                      :: ng_a, ng_r, ng_m
    integer                                      :: lo(get_dim(rh)),hi(get_dim(rh))

    ! do NOT add bl_prof_timer in this subroutine becuase not every MPI rank calls this

    ng_a = nghost(aa)
    ng_r = nghost(rh)
    ng_m = nghost(mm)

    dm = get_dim(rh)

    !$OMP PARALLEL DO PRIVATE(i,rp,ap,mp,lo,hi)
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
                call diag_init_nd_1d(ap(1,:,1,1), ng_a, rp(:,1,1,1), ng_r, mp(:,1,1,1), ng_m, &
                                     lo, hi)
             end if
          case (2)
             if (is_ibc_stencil(aa,i)) then
                call diag_init_ibc_2d(ap(:,1,1,1), rp(:,:,1,1), ng_r, lo, hi)
             else if ( cell_centered_q(rh) ) then
                call diag_init_cc_2d(ap(:,:,:,1), ng_a, rp(:,:,1,1), ng_r, lo, hi)
             else
                call diag_init_nd_2d(ap(1,:,:,1), ng_a, rp(:,:,1,1), ng_r, mp(:,:,1,1), ng_m, &
                                     lo, hi, stencil_type)
             end if
          case (3)
             if (is_ibc_stencil(aa,i)) then
                call diag_init_ibc_3d(ap(:,1,1,1), rp(:,:,:,1), ng_r, lo, hi)
             else if ( cell_centered_q(rh) ) then
                call diag_init_cc_3d(ap(:,:,:,:), ng_a, rp(:,:,:,1), ng_r, lo, hi)
             else
                call diag_init_nd_3d(ap(1,:,:,:), ng_a, rp(:,:,:,1), ng_r, mp(:,:,:,1), ng_m, &
                                     lo, hi, stencil_type)
             end if
       end select
    end do
    !$OMP END PARALLEL DO

  end subroutine diag_initialize

end module itsol_module

