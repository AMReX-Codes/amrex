module ml_solve_module

   use bl_types
   use mg_module
   use bndry_reg_module
   use ml_layout_module
   use multifab_module

   implicit none

   private

   public :: ml_cc_solve, ml_nd_solve

contains

   subroutine ml_cc_solve(mla,mgt,rh,full_soln,fine_flx,ref_ratio,do_diagnostics,rel_eps_in,abs_eps_in)

      use ml_cc_module , only : ml_cc

      type(ml_layout), intent(in   ) :: mla
      type(mg_tower ), intent(inout) :: mgt(:)
      type(multifab ), intent(inout) :: rh(:)
      type(multifab ), intent(inout) :: full_soln(:)
      type(bndry_reg), intent(inout) :: fine_flx(2:)
      integer        , intent(in   ) :: ref_ratio(:,:)
      integer        , intent(in   ) :: do_diagnostics
      real(dp_t)     , intent(in   ), optional :: rel_eps_in
      real(dp_t)     , intent(in   ), optional :: abs_eps_in

      type(boxarray)  :: bac
      type(lmultifab) :: fine_mask(mla%nlevel)
      integer         :: i, dm, n, nlevs, mglev
      real(dp_t)      :: rel_eps
      real(dp_t)      :: abs_eps

      dm    = mla%dim
      nlevs = mla%nlevel

      rel_eps =  1.d-12; if ( present(rel_eps_in) ) rel_eps = rel_eps_in
      abs_eps = -1.d0  ; if ( present(abs_eps_in) ) abs_eps = abs_eps_in

      do n = nlevs, 1, -1
        call lmultifab_build(fine_mask(n), mla%la(n), 1, 0)
        call setval(fine_mask(n), val = .true., all = .true.)
      end do
      do n = nlevs-1, 1, -1
        call copy(bac, get_boxarray(mla%la(n+1)))
        call boxarray_coarsen(bac, ref_ratio(n,:))
        call setval(fine_mask(n), .false., bac)
        call destroy(bac)
      end do

! ****************************************************************************

      call ml_cc(mla,mgt,rh,full_soln,fine_mask,ref_ratio,do_diagnostics,rel_eps,abs_eps, &
                 need_grad_phi_in=.true.)

! ****************************************************************************

!   Put boundary conditions of soln in fine_flx to get correct grad(phi) at
!     crse-fine boundaries (after soln correctly interpolated in ml_cc)
    do n = 2,nlevs
       mglev = mgt(n)%nlevels
       do i = 1, dm
          call ml_fill_fine_fluxes(mgt(n)%ss(mglev), fine_flx(n)%bmf(i,0), &
                                   full_soln(n), mgt(n)%mm(mglev), -1, i, &
                                   mgt(n)%lcross)
          call ml_fill_fine_fluxes(mgt(n)%ss(mglev), fine_flx(n)%bmf(i,1), &
                                   full_soln(n), mgt(n)%mm(mglev),  1, i, &
                                   mgt(n)%lcross)
       end do
    end do

    do n = 1,nlevs
      call lmultifab_destroy(fine_mask(n))
    end do

   end subroutine ml_cc_solve

   subroutine ml_fill_fine_fluxes(ss, flux, uu, mm, face, dim, lcross)

     use bl_prof_module
     use cc_stencil_apply_module

     type(multifab) , intent(inout) :: flux
     type(multifab) , intent(in   ) :: ss
     type(multifab) , intent(inout) :: uu
     type(imultifab), intent(in   ) :: mm
     integer        , intent(in   ) :: face, dim
     logical        , intent(in   ) :: lcross

     integer :: i, n
     real(kind=dp_t), pointer :: fp(:,:,:,:)
     real(kind=dp_t), pointer :: up(:,:,:,:)
     real(kind=dp_t), pointer :: sp(:,:,:,:)
     integer        , pointer :: mp(:,:,:,:)
     integer :: ng
     type(bl_prof_timer), save :: bpt

     call build(bpt, "ml_fill_fine_fluxes")

     ng = nghost(uu)

     if ( ncomp(uu) /= ncomp(flux) ) then
        call bl_error("ML_FILL_FINE_FLUXES: uu%nc /= flux%nc")
     end if

     call multifab_fill_boundary(uu, cross = lcross)

     !$OMP PARALLEL DO PRIVATE(i,n,fp,up,sp,mp)
     do i = 1, nfabs(flux)
        fp => dataptr(flux, i)
        up => dataptr(uu, i)
        sp => dataptr(ss, i)
        mp => dataptr(mm, i)
        do n = 1, ncomp(uu)
           select case(get_dim(ss))
           case (1)
              call stencil_fine_flux_1d(sp(:,:,1,1), fp(:,1,1,n), up(:,1,1,n), &
                   mp(:,1,1,1), ng, face, dim)
           case (2)
              call stencil_fine_flux_2d(sp(:,:,:,1), fp(:,:,1,n), up(:,:,1,n), &
                   mp(:,:,1,1), ng, face, dim)
           case (3)
              call stencil_fine_flux_3d(sp(:,:,:,:), fp(:,:,:,n), up(:,:,:,n), &
                   mp(:,:,:,1), ng, face, dim)
           end select
        end do
     end do
     !$OMP END PARALLEL DO

     call destroy(bpt)

   end subroutine ml_fill_fine_fluxes

!
! ******************************************************************************************
!

   subroutine ml_nd_solve(mla,mgt,rh,full_soln,one_sided_ss,ref_ratio,do_diagnostics,&
                          rel_eps_in,abs_eps_in)

       use ml_nd_module, only : ml_nd

       type(ml_layout), intent(in   )           :: mla
       type(mg_tower) , intent(inout)           :: mgt(:)
       type(multifab) , intent(inout)           :: rh(:)
       type(multifab) , intent(inout)           :: full_soln(:)
       type(multifab) , intent(in   )           :: one_sided_ss(2:)
       integer        , intent(in   )           :: ref_ratio(:,:)
       integer        , intent(in   )           :: do_diagnostics 
       real(dp_t)     , intent(in   ), optional :: rel_eps_in
       real(dp_t)     , intent(in   ), optional :: abs_eps_in

       type(lmultifab) :: fine_mask(mla%nlevel)
       integer         :: nlevs, n, dm
       logical         :: nodal(get_dim(rh(mla%nlevel)))
       real(dp_t)      :: rel_eps,abs_eps

       rel_eps =  1.d-12; if ( present(rel_eps_in) ) rel_eps = rel_eps_in
       abs_eps = -1.d0  ; if ( present(abs_eps_in) ) abs_eps = abs_eps_in

       nlevs = mla%nlevel
       dm    = get_dim(rh(nlevs))
       nodal = .true.

!      We are only considering the dense stencils here (3 in 1d, 9 in 2d, 27 in 3d)

       do n = nlevs, 1, -1
          call lmultifab_build(fine_mask(n), mla%la(n), 1, 0, nodal)
          if ( n < nlevs ) then
             call create_nodal_mask(fine_mask(n), &
                                    mgt(n  )%mm(mgt(n  )%nlevels), &
                                    mgt(n+1)%mm(mgt(n+1)%nlevels), &
                                    ref_ratio(n,:))
          else
             call setval(fine_mask(n), val = .true., all = .true.)
          endif
       end do

       call ml_nd(mla,mgt,rh,full_soln,fine_mask,one_sided_ss,ref_ratio,do_diagnostics,&
                  rel_eps,abs_eps)
     
       do n = 1,nlevs
          call lmultifab_destroy(fine_mask(n))
       end do

   end subroutine ml_nd_solve

   subroutine create_nodal_mask(mask,mm_crse,mm_fine,ir)

     type(lmultifab), intent(inout) :: mask
     type(imultifab), intent(in   ) :: mm_crse
     type(imultifab), intent(in   ) :: mm_fine
     integer        , intent(in   ) :: ir(:)

     type(box)                        :: cbox, fbox, isect
     type(boxarray)                   :: ba
     type(layout)                     :: la
     logical,               pointer   :: mkp(:,:,:,:)
     integer                          :: loc(get_dim(mask)), lof(get_dim(mask)), i, j, k, dims(4), proc, dm, ii, jj
     integer,               pointer   :: cmp(:,:,:,:), fmp(:,:,:,:)
     integer                          :: lo(get_dim(mask)), hi(get_dim(mask))
     integer,               parameter :: tag = 1071
     type(box_intersector), pointer   :: bi(:)
     type(bl_prof_timer),   save      :: bpt

     if ( .not. nodal_q(mask) ) call bl_error('create_nodal_mask(): mask NOT nodal')

     call build(bpt, "create_nodal_mask")

     call setval(mask, .true.)
     !
     ! Need to build temporary layout with nodal boxarray for the intersection tests below.
     !
     call copy(ba, get_boxarray(get_layout(mask)))
     call boxarray_nodalize(ba, nodal_flags(mask))
     call build(la, ba, boxarray_bbox(ba), mapping = LA_LOCAL)  ! LA_LOCAL ==> bypass processor distribution calculation.
     call destroy(ba)
     !
     !   Note :          mm_fine is  in fine space
     !   Note : mask and mm_crse are in crse space
     !
     dims = 1

     dm = get_dim(mask)

     do i = 1, nboxes(mm_fine%la)

        fbox =  box_nodalize(get_box(mm_fine%la,i),mm_fine%nodal)
        bi   => layout_get_box_intersector(la, coarsen(fbox,ir))

        do k = 1, size(bi)
           j = bi(k)%i

           if ( remote(mask%la,j) .and. remote(mm_fine%la,i) ) cycle

           cbox  = box_nodalize(get_box(mask%la,j),mask%nodal)
           loc   = lwb(cbox)
           isect = bi(k)%bx
           lo    = lwb(isect)
           hi    = upb(isect)

           if ( local(mask%la,j) .and. local(mm_fine%la,i) ) then
              ii  =  local_index(mm_fine,i)
              jj  =  local_index(mask,j)
              lof =  lwb(fbox)
              mkp => dataptr(mask,jj)
              cmp => dataptr(mm_crse,jj)
              fmp => dataptr(mm_fine,ii)

              select case (dm)
              case (1)
                 call create_nodal_mask_1d(mkp(:,1,1,1),cmp(:,1,1,1),loc,fmp(:,1,1,1),lof,lo,hi,ir)
              case (2)
                 call create_nodal_mask_2d(mkp(:,:,1,1),cmp(:,:,1,1),loc,fmp(:,:,1,1),lof,lo,hi,ir)
              case (3)
                 call create_nodal_mask_3d(mkp(:,:,:,1),cmp(:,:,:,1),loc,fmp(:,:,:,1),lof,lo,hi,ir)
              end select
           else if ( local(mm_fine%la,i) ) then
              !
              ! Must send mm_fine.
              !
              ii    =  local_index(mm_fine,i)
              isect =  intersection(refine(isect,ir),fbox)
              fmp   => dataptr(mm_fine, ii, isect)
              proc  =  get_proc(get_layout(mask), j)
              call parallel_send(fmp, proc, tag)
           else if ( local(mask%la,j) ) then
              !
              ! Must receive mm_fine.
              !
              jj    =  local_index(mask,j)
              isect =  intersection(refine(isect,ir),fbox)
              lof   =  lwb(isect)
              mkp   => dataptr(mask,jj)
              cmp   => dataptr(mm_crse,jj)
              proc  =  get_proc(get_layout(mm_fine),i)

              dims(1:dm) = extent(isect)
              allocate(fmp(dims(1),dims(2),dims(3),1))
              call parallel_recv(fmp, proc, tag)

              select case (dm)
              case (1)
                 call create_nodal_mask_1d(mkp(:,1,1,1),cmp(:,1,1,1),loc,fmp(:,1,1,1),lof,lo,hi,ir)
              case (2)
                 call create_nodal_mask_2d(mkp(:,:,1,1),cmp(:,:,1,1),loc,fmp(:,:,1,1),lof,lo,hi,ir)
              case (3)
                 call create_nodal_mask_3d(mkp(:,:,:,1),cmp(:,:,:,1),loc,fmp(:,:,:,1),lof,lo,hi,ir)
              end select

              deallocate(fmp)
           end if
        end do
        deallocate(bi)
     end do

     call destroy(la)
     call destroy(bpt)

   end subroutine create_nodal_mask

   subroutine create_nodal_mask_1d(mask,mm_crse,loc,mm_fine,lof,lo,hi,ir)

     integer, intent(in   ) :: loc(:),lof(:)
     logical, intent(inout) ::    mask(loc(1):)
     integer, intent(in   ) :: mm_crse(loc(1):)
     integer, intent(in   ) :: mm_fine(lof(1):)
     integer, intent(in   ) :: lo(:),hi(:)
     integer, intent(in   ) :: ir(:)

     integer :: i,fi

     do i = lo(1),hi(1)
        fi = i*ir(1)
        if (.not.  bc_dirichlet(mm_fine(fi),1,0) .or. bc_dirichlet(mm_crse(i),1,0)) &
             mask(i) = .false.
     end do

   end subroutine create_nodal_mask_1d

   subroutine create_nodal_mask_2d(mask,mm_crse,loc,mm_fine,lof,lo,hi,ir)

     integer, intent(in   ) :: loc(:),lof(:)
     logical, intent(inout) ::    mask(loc(1):,loc(2):)
     integer, intent(in   ) :: mm_crse(loc(1):,loc(2):)
     integer, intent(in   ) :: mm_fine(lof(1):,lof(2):)
     integer, intent(in   ) :: lo(:),hi(:)
     integer, intent(in   ) :: ir(:)

     integer :: i,j,fi,fj

     do j = lo(2),hi(2)
        fj = j*ir(2)
        do i = lo(1),hi(1)
           fi = i*ir(1)
           if (.not.  bc_dirichlet(mm_fine(fi,fj),1,0) .or. bc_dirichlet(mm_crse(i,j),1,0)) &
                mask(i,j) = .false.
        end do
     end do

   end subroutine create_nodal_mask_2d

   subroutine create_nodal_mask_3d(mask,mm_crse,loc,mm_fine,lof,lo,hi,ir)

     integer, intent(in   ) :: loc(:),lof(:)
     logical, intent(inout) ::    mask(loc(1):,loc(2):,loc(3):)
     integer, intent(in   ) :: mm_crse(loc(1):,loc(2):,loc(3):)
     integer, intent(in   ) :: mm_fine(lof(1):,lof(2):,lof(3):)
     integer, intent(in   ) :: lo(:),hi(:)
     integer, intent(in   ) :: ir(:)

     integer :: i,j,k,fi,fj,fk

     !$OMP PARALLEL DO PRIVATE(i,j,k,fi,fj,fk)
     do k = lo(3),hi(3)
        fk = k*ir(3)
        do j = lo(2),hi(2)
           fj = j*ir(2)
           do i = lo(1),hi(1)
              fi = i*ir(1)
              if (.not. bc_dirichlet(mm_fine(fi,fj,fk),1,0) .or. bc_dirichlet(mm_crse(i,j,k),1,0) ) &
                   mask(i,j,k) = .false.
           end do
        end do
     end do
     !$OMP END PARALLEL DO

   end subroutine create_nodal_mask_3d

end module ml_solve_module
