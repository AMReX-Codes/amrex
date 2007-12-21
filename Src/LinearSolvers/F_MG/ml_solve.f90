module ml_solve_module

   use bl_types
   use mg_module
   use ml_layout_module
   use bndry_reg_module
   use multifab_module

   implicit none

   private

   public :: ml_cc_solve, ml_nd_solve

contains

   subroutine ml_cc_solve(mla,mgt,rh,full_soln,fine_flx,ref_ratio,do_diagnostics)

      use ml_util_module
      use ml_cc_module

      type(ml_layout), intent(inout) :: mla
      type(mg_tower ), intent(inout) :: mgt(:)
      type(multifab ), intent(inout) :: rh(:)
      type(multifab ), intent(inout) :: full_soln(:)
      type(bndry_reg), intent(inout) :: fine_flx(2:)
      integer        , intent(in   ) :: ref_ratio(:,:)
      integer        , intent(in   ) :: do_diagnostics

      type(boxarray)           :: bac
      type(lmultifab), pointer :: fine_mask(:) => Null()
      integer                  :: i, dm, n, nlevs, mglev
      real(dp_t)               :: eps

      dm    = mla%dim
      nlevs = mla%nlevel

      eps = 1.d-12

      allocate(fine_mask(nlevs))

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

      call ml_cc(mla,mgt,rh,full_soln,fine_mask,ref_ratio,do_diagnostics,eps,.true.)

! ****************************************************************************

!   Put boundary conditions of soln in fine_flx to get correct grad(phi) at
!     crse-fine boundaries (after soln correctly interpolated in ml_cc)
    do n = 2,nlevs
       mglev = mgt(n)%nlevels
       do i = 1, dm
          call ml_fill_fine_fluxes(mgt(n)%ss(mglev), fine_flx(n)%bmf(i,0), &
                                   full_soln(n), mgt(n)%mm(mglev), -1, i)
          call ml_fill_fine_fluxes(mgt(n)%ss(mglev), fine_flx(n)%bmf(i,1), &
                                   full_soln(n), mgt(n)%mm(mglev),  1, i)
       end do
    end do

    do n = 1,nlevs
      call lmultifab_destroy(fine_mask(n))
    end do
    deallocate(fine_mask)

   end subroutine ml_cc_solve

   subroutine ml_nd_solve(mla,mgt,rh,full_soln,one_sided_ss,ref_ratio,do_diagnostics,eps_in)

       use ml_nd_module

       type(ml_layout), intent(inout)           :: mla
       type(mg_tower) , intent(inout)           :: mgt(:)
       type(multifab) , intent(inout)           :: rh(:)
       type(multifab) , intent(inout)           :: full_soln(:)
       type(multifab) , intent(in   )           :: one_sided_ss(2:)
       integer        , intent(in   )           :: ref_ratio(:,:)
       integer        , intent(in   )           :: do_diagnostics 
       real(dp_t)     , intent(in   ), optional :: eps_in

       type(lmultifab), pointer :: fine_mask(:) => Null()
       integer                  :: nlevs, n, dm, lo(rh(mla%nlevel)%dim), hi(rh(mla%nlevel)%dim)
       logical                  :: nodal(rh(mla%nlevel)%dim)
       real(dp_t)               :: eps

       eps = 1.d-12; if ( present(eps_in) ) eps = eps_in

       nlevs = mla%nlevel
       dm    = rh(nlevs)%dim
       nodal = .true.

       allocate(fine_mask(nlevs))

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

       call ml_nd(mla,mgt,rh,full_soln,fine_mask,one_sided_ss,ref_ratio,do_diagnostics,eps)
     
       do n = 1,nlevs
          call lmultifab_destroy(fine_mask(n))
       end do

       deallocate(fine_mask)

   contains
     !
     ! TODO - cache the communication stuff in the mask's layout?
     !
     subroutine create_nodal_mask(mask,mm_crse,mm_fine,ir)

       type(lmultifab), intent(inout) :: mask
       type(imultifab), intent(in   ) :: mm_crse
       type(imultifab), intent(in   ) :: mm_fine
       integer        , intent(in   ) :: ir(:)

       type(box)          :: cbox, fbox, isect
       logical, pointer   :: mkp(:,:,:,:)
       integer            :: loc(mask%dim), lof(mask%dim), i, j, k, dims(4), proc
       integer, pointer   :: cmp(:,:,:,:), fmp(:,:,:,:)
       integer, parameter :: tag = 1071

       type(box_intersector), pointer :: bi(:)

       type(bl_prof_timer), save :: bpt

       call build(bpt, "create_nodal_mask")

       call setval(mask,.true.)

       !       Note :          mm_fine is  in fine space
       !       Note : mask and mm_crse are in crse space

       dims = 1

       do i = 1,mm_fine%nboxes

          fbox =  get_ibox(mm_fine,i)
          bi   => layout_get_box_intersector(mask%la, coarsen(fbox,ir))

          do k = 1, size(bi)
             j = bi(k)%i

             if ( remote(mask,j) .and. remote(mm_fine,i) ) cycle

             cbox  = get_ibox(mask,j)
             loc   = lwb(cbox)
             isect = bi(k)%bx
             lo    = lwb(isect)
             hi    = upb(isect)

             if ( local(mask,j) .and. local(mm_fine,i) ) then
                lof =  lwb(fbox)
                mkp => dataptr(mask,j)
                cmp => dataptr(mm_crse,j)
                fmp => dataptr(mm_fine,i)

                select case (dm)
                case (2)
                   call create_nodal_mask_2d(mkp(:,:,1,1),cmp(:,:,1,1),loc,fmp(:,:,1,1),lof,lo,hi,ir)
                case (3)
                   call create_nodal_mask_3d(mkp(:,:,:,1),cmp(:,:,:,1),loc,fmp(:,:,:,1),lof,lo,hi,ir)
                end select
             else if ( local(mm_fine,i) ) then
                !
                ! Must send mm_fine.
                !
                isect =  intersection(refine(isect,ir),get_ibox(mm_fine,i))
                fmp   => dataptr(mm_fine, i, isect)
                proc  =  get_proc(mask%la, j)
                call parallel_send(fmp, proc, tag)
             else if ( local(mask,j) ) then
                !
                ! Must receive mm_fine.
                !
                isect =  intersection(refine(isect,ir),get_ibox(mm_fine,i))
                lof   =  lwb(isect)
                mkp   => dataptr(mask,j)
                cmp   => dataptr(mm_crse,j)
                proc  =  get_proc(mm_fine%la,i)

                dims(1:dm) = extent(isect)
                allocate(fmp(dims(1),dims(2),dims(3),dims(4)))
                call parallel_recv(fmp, proc, tag)

                select case (dm)
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

       call destroy(bpt)

     end subroutine create_nodal_mask

     subroutine create_nodal_mask_2d(mask,mm_crse,loc,mm_fine,lof,lo,hi,ir)

       integer, intent(in   ) :: loc(:),lof(:)
       logical, intent(inout) ::    mask(loc(1):,loc(2):)
       integer, intent(inout) :: mm_crse(loc(1):,loc(2):)
       integer, intent(inout) :: mm_fine(lof(1):,lof(2):)
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
       integer, intent(inout) :: mm_crse(loc(1):,loc(2):,loc(3):)
       integer, intent(inout) :: mm_fine(lof(1):,lof(2):,lof(3):)
       integer, intent(in   ) :: lo(:),hi(:)
       integer, intent(in   ) :: ir(:)

       integer :: i,j,k,fi,fj,fk

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

     end subroutine create_nodal_mask_3d

   end subroutine ml_nd_solve

end module ml_solve_module
