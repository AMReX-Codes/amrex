module nodal_mask_module

  use nodal_stencil_module
  use ml_layout_module
  use bl_prof_module

  implicit none

  private

  public :: create_nodal_mask

  contains

  subroutine create_nodal_mask(n,mask,mm_crse,mm_fine,mla)

  integer        , intent(in   ) :: n
  type(lmultifab), intent(inout) :: mask
  type(imultifab), intent(inout) :: mm_crse
  type(imultifab), intent(inout) :: mm_fine
  type(ml_layout), intent(in   ) :: mla

  integer   :: dm
  integer   :: ir(get_dim(mask)),hi_fine(get_dim(mask)),hi_crse(get_dim(mask))
  type(box) :: cbox

  logical, pointer :: mkp(:,:,:,:)
  integer, pointer :: cmp(:,:,:,:),fmp(:,:,:,:)

  type(layout)     :: lacfine, mm_fine_la
  type(lmultifab)  :: cfmask

  integer :: lo(get_dim(mask)),hi(get_dim(mask)),lof(get_dim(mask)),j

  type(bl_prof_timer), save :: bpt

  call build(bpt, "create_nodal_mask")

  dm = get_dim(mask)

! Note :          mm_fine is  in fine space
! Note : mask and mm_crse are in crse space
  
  hi_fine = upb(layout_get_pd(mla%la(n+1))) + 1
  hi_crse = upb(layout_get_pd(mla%la(n  ))) + 1

  ir = hi_fine / hi_crse

  mm_fine_la = get_layout(mm_fine)

  call layout_build_coarse(lacfine, mm_fine_la, ir)

  call lmultifab_build(cfmask, lacfine, ncomp(mask), nghost(mask), nodal_flags(mask))

  call setval(  mask,.true.)
  call setval(cfmask,.true.)

  !$OMP PARALLEL DO PRIVATE(j,cbox,lo,hi,lof,mkp,fmp)
  do j = 1,nfabs(cfmask)
     cbox = get_ibox(cfmask,j)
     lo   = lwb(cbox)
     hi   = upb(cbox)
     lof  = lwb(get_ibox(mm_fine,j))
     mkp  => dataptr(cfmask,j)
     fmp  => dataptr(mm_fine,j)
     select case (dm)
     case (1)
        call set_crsefine_nodal_mask_1d(mkp(:,1,1,1),fmp(:,1,1,1),lof,lo,hi,ir)
     case (2)
        call set_crsefine_nodal_mask_2d(mkp(:,:,1,1),fmp(:,:,1,1),lof,lo,hi,ir)
     case (3)
        call set_crsefine_nodal_mask_3d(mkp(:,:,:,1),fmp(:,:,:,1),lof,lo,hi,ir)
     end select
  end do
  !$OMP END PARALLEL DO

  call copy(mask, 1, cfmask, 1, ncomp(mask))

  call lmultifab_destroy(cfmask)

  !$OMP PARALLEL DO PRIVATE(j,cbox,lo,hi,mkp,cmp)
  do j = 1,nfabs(mask)
     cbox =  get_ibox(mask,j)
     lo   =  lwb(cbox)
     hi   =  upb(cbox)
     mkp  => dataptr(mask,j)
     cmp  => dataptr(mm_crse,j)
     select case (dm)
     case (1)
        call set_crse_nodal_mask_1d(mkp(:,1,1,1),cmp(:,1,1,1),lo,hi)
     case (2)
        call set_crse_nodal_mask_2d(mkp(:,:,1,1),cmp(:,:,1,1),lo,hi)
     case (3)
        call set_crse_nodal_mask_3d(mkp(:,:,:,1),cmp(:,:,:,1),lo,hi)
     end select
  end do
  !$OMP END PARALLEL DO

  call destroy(bpt)

  end subroutine create_nodal_mask

  subroutine set_crsefine_nodal_mask_1d(mask,mm_fine,lof,lo,hi,ir)

       integer, intent(in   ) :: lo(:),hi(:),lof(:),ir(:)
       logical, intent(inout) :: mask(lo(1):)
       integer, intent(inout) :: mm_fine(lof(1):)

       integer :: i

       do i = lo(1),hi(1)
          if (.not.  bc_dirichlet(mm_fine(i*ir(1)),1,0) ) &
             mask(i) = .false.
       end do

  end subroutine set_crsefine_nodal_mask_1d

  subroutine set_crsefine_nodal_mask_2d(mask,mm_fine,lof,lo,hi,ir)

       integer, intent(in   ) :: lo(:),hi(:),lof(:),ir(:)
       logical, intent(inout) :: mask(lo(1):,lo(2):)
       integer, intent(inout) :: mm_fine(lof(1):,lof(2):)

       integer :: i,j

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (.not.  bc_dirichlet(mm_fine(i*ir(1),j*ir(2)),1,0)) &
                mask(i,j) = .false.
          end do
       end do

  end subroutine set_crsefine_nodal_mask_2d

  subroutine set_crsefine_nodal_mask_3d(mask,mm_fine,lof,lo,hi,ir)

       integer, intent(in   ) :: lo(:),hi(:),lof(:),ir(:)
       logical, intent(inout) :: mask(lo(1):,lo(2):,lo(3):)
       integer, intent(inout) :: mm_fine(lof(1):,lof(2):,lof(3):)

       integer :: i,j,k

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (.not. bc_dirichlet(mm_fine(i*ir(1),j*ir(2),k*ir(3)),1,0)) &
                   mask(i,j,k) = .false.
             end do
          end do
       end do

  end subroutine set_crsefine_nodal_mask_3d

  subroutine set_crse_nodal_mask_1d(mask,mm_crse,lo,hi)

       integer, intent(in   ) :: lo(:),hi(:)
       logical, intent(inout) ::    mask(lo(1):)
       integer, intent(inout) :: mm_crse(lo(1):)

       integer :: i

       do i = lo(1),hi(1)
          if ( bc_dirichlet(mm_crse(i),1,0 ) ) &
             mask(i) = .false.
       end do

  end subroutine set_crse_nodal_mask_1d

  subroutine set_crse_nodal_mask_2d(mask,mm_crse,lo,hi)

       integer, intent(in   ) :: lo(:),hi(:)
       logical, intent(inout) ::    mask(lo(1):,lo(2):)
       integer, intent(inout) :: mm_crse(lo(1):,lo(2):)

       integer :: i,j

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if ( bc_dirichlet(mm_crse(i,j),1,0 ) ) &
                mask(i,j) = .false.
          end do
       end do

  end subroutine set_crse_nodal_mask_2d

  subroutine set_crse_nodal_mask_3d(mask,mm_crse,lo,hi)

       integer, intent(in   ) :: lo(:),hi(:)
       logical, intent(inout) ::    mask(lo(1):,lo(2):,lo(3):)
       integer, intent(inout) :: mm_crse(lo(1):,lo(2):,lo(3):)

       integer :: i,j,k
       logical :: jface,kface

       do k = lo(3),hi(3)
          kface = .false. ; if ( (k.eq.lo(3)) .or. (k.eq.hi(3)) ) kface = .true.
          do j = lo(2),hi(2)
             jface = .false. ; if ( (j.eq.lo(2)) .or. (j.eq.hi(2)) ) jface = .true.
             do i = lo(1),hi(1)
                if ( jface .or. kface .or. (i.eq.lo(1)) .or. (i.eq.hi(1)) ) then
                   if ( bc_dirichlet(mm_crse(i,j,k),1,0) ) mask(i,j,k) = .false.
                end if
             end do
          end do
       end do

     end subroutine set_crse_nodal_mask_3d

end module nodal_mask_module
