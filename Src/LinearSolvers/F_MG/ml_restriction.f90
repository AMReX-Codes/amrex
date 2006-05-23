module ml_restriction_module

  use bl_types
  use multifab_module
  use mg_restriction_module

  implicit none

  real(dp_t), private, parameter :: ZERO = 0.0_dp_t
  real(dp_t), private, parameter :: ONE  = 1.0_dp_t

contains

  subroutine ml_cc_restriction_c(crse, cc, fine, cf, ir, nc)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer, intent(in)           :: cc, cf, ir(:)
    integer, intent(in), optional :: nc

    integer             :: i, n, lnc, lo(fine%dim), hi(fine%dim), lof(fine%dim)
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    type(layout)        :: lacfine
    type(multifab)      :: cfine

    lnc = 1; if ( present(nc) ) lnc = nc

    call layout_build_coarse(lacfine, fine%la, ir)

    call build(cfine, lacfine, nc = lnc, ng = 0)

    do i = 1, fine%nboxes
       if ( remote(fine, i) ) cycle
       lof = lwb(get_ibox(fine, i)) - fine%ng 
       lo  = lwb(get_ibox(cfine,i))
       hi  = upb(get_ibox(cfine,i))
       do n = 1, lnc
          fp => dataptr(fine,  i, n+cf-1, 1)
          cp => dataptr(cfine, i, n,      1)
          select case (cfine%dim)
          case (1)
             call cc_restriction_1d(cp(:,1,1,1), lo, fp(:,1,1,1), lof, lo, hi, ir)
          case (2)
             call cc_restriction_2d(cp(:,:,1,1), lo, fp(:,:,1,1), lof, lo, hi, ir)
          case (3)
             call cc_restriction_3d(cp(:,:,:,1), lo, fp(:,:,:,1), lof, lo, hi, ir)
          end select
       end do
    end do

    call copy(crse, cc, cfine, 1, lnc)

    call multifab_destroy(cfine)

  end subroutine ml_cc_restriction_c

  subroutine ml_cc_restriction(crse, fine, ir)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer,        intent(in)    :: ir(:)
    if ( crse%nc .ne. fine%nc ) then
       call bl_error('ml_cc_restriction: crse & fine must have same # of components')
    end if
    call ml_cc_restriction_c(crse, 1, fine, 1, ir, crse%nc)
  end subroutine ml_cc_restriction

  subroutine ml_edge_restriction(crse, fine, ir, n)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer,        intent(in)    :: ir(:)
    integer,        intent(in)    :: n

    integer             :: i
    integer             :: lo(fine%dim), hi(fine%dim), loc(fine%dim), lof(fine%dim)
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    type(layout)        :: lacfine
    type(multifab)      :: cfine

    call layout_build_coarse(lacfine, fine%la, ir)

    call multifab_build(cfine, lacfine, nc = crse%nc, ng = 0, nodal = crse%nodal)

    do i = 1, fine%nboxes
       if ( remote(fine,i) ) cycle
       lo  = lwb(get_ibox(cfine,i))
       hi  = upb(get_ibox(cfine,i))
       loc = lwb(get_pbox(cfine,i))
       lof = lwb(get_pbox(fine, i))
       fp  => dataptr(fine,  i)
       cp  => dataptr(cfine, i)
       select case (crse%dim)
       case (1)
          call edge_restriction_1d(cp(:,1,1,1), loc, fp(:,1,1,1), lof, lo, hi, ir)
       case (2)
          call edge_restriction_2d(cp(:,:,1,1), loc, fp(:,:,1,1), lof, lo, hi, ir, n)
       case (3)
          call edge_restriction_3d(cp(:,:,:,1), loc, fp(:,:,:,1), lof, lo, hi, ir, n)
       end select
    end do

    call copy(crse, cfine)

    call multifab_destroy(cfine)

  end subroutine ml_edge_restriction

  subroutine ml_restrict_copy_sum(out, in)
     use bl_types
     real(dp_t), intent(inout) :: out(:,:,:,:)
     real(dp_t), intent(in   ) ::  in(:,:,:,:)
     out = out + in 
  end subroutine ml_restrict_copy_sum

  subroutine ml_nodal_restriction(crse, fine, mm_fine, mm_crse, face_type, ir, inject, zero_only)
    type(multifab),  intent(inout) :: fine
    type(multifab),  intent(inout) :: crse
    type(imultifab), intent(in   ) :: mm_fine
    type(imultifab), intent(in   ) :: mm_crse
    integer,         intent(in)    :: ir(:)
    integer,         intent(in)    :: face_type(:,:,:)
    logical,         intent(in), optional :: inject
    logical,         intent(in), optional :: zero_only

    integer             :: i, j, n, id
    integer             :: lo (fine%dim), hi (fine%dim), loc(fine%dim), lof(fine%dim)
    integer             :: lom_fine(fine%dim), lom_crse(fine%dim)
    logical             :: linject, lzero_only
    type(box)           :: fbox, cbox, isect
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    integer,    pointer :: mp_fine(:,:,:,:), mp_crse(:,:,:,:)
    integer             :: rmode
    type(layout)        :: lacfine
    type(multifab)      :: cfine
    type(imultifab)     :: mm_cfine

    if ( crse%nc .ne. fine%nc ) then
       call bl_error('ml_nodal_restriction: crse & fine must have same # of components')
    end if

    linject    = .false. ; if ( present(inject   ) ) linject    = inject
    lzero_only = .false. ; if ( present(zero_only) ) lzero_only = zero_only

    call layout_build_coarse(lacfine, fine%la, ir)
    call  multifab_build(cfine,    lacfine, nc =    crse%nc, ng = 0, nodal =    crse%nodal)
    call imultifab_build(mm_cfine, lacfine, nc = mm_crse%nc, ng = 0, nodal = mm_crse%nodal)
    call copy(   cfine,    crse)
    call copy(mm_cfine, mm_crse)

    rmode = 0

    if ( .not. linject ) then
       do i = 1, fine%nboxes
          if ( remote(fine, i) ) cycle
          lo       = lwb(get_ibox(cfine,   i))
          hi       = upb(get_ibox(cfine,   i))
          loc      = lwb(get_pbox(cfine,   i))
          lom_fine = lwb(get_pbox(mm_fine, i))
          do n = 1, fine%nc
             cp      => dataptr(cfine,   i, n, 1)
             mp_fine => dataptr(mm_fine, i, n, 1)
             select case (fine%dim)
             case (1)
                call nodal_zero_1d(cp(:,1,1,1), loc, mp_fine(:,1,1,1), lom_fine, lo, hi, ir)
             case (2)
                call nodal_zero_2d(cp(:,:,1,1), loc, mp_fine(:,:,1,1), lom_fine, lo, hi, ir)
             case (3)
                call nodal_zero_3d(cp(:,:,:,1), loc, mp_fine(:,:,:,1), lom_fine, lo, hi, ir)
             end select
          end do
       end do
       call copy(crse, cfine)
       call setval(cfine, 0.0_dp_t)
    end if

    if ( .not. lzero_only ) then
    do i = 1, fine%nboxes
       if ( remote(fine, i) ) cycle
       lo       = lwb(get_ibox(cfine,   i))
       hi       = upb(get_ibox(cfine,   i))
       lof      = lwb(get_pbox(fine,    i))
       loc      = lwb(get_pbox(cfine,   i))
       lom_crse = lwb(get_pbox(mm_cfine,i))
       lom_fine = lwb(get_pbox(mm_fine, i))
       do n = 1, fine%nc
          cp      => dataptr(cfine,   i, n, 1)
          fp      => dataptr(fine,    i, n, 1)
          mp_crse => dataptr(mm_cfine,i, n, 1)
          mp_fine => dataptr(mm_fine, i, n, 1)
          select case (fine%dim)
          case (1)
             call nodal_restriction_1d(cp(:,1,1,1), loc, fp(:,1,1,1), lof, &
                  mp_fine(:,1,1,1), lom_fine, &
                  mp_crse(:,1,1,1), lom_crse, lo, hi, ir, linject, rmode)
          case (2)
             call nodal_restriction_2d(cp(:,:,1,1), loc, fp(:,:,1,1), lof, &
                  mp_fine(:,:,1,1), lom_fine, &
                  mp_crse(:,:,1,1), lom_crse, lo, hi, ir, linject, rmode)
          case (3)
             call nodal_restriction_3d(cp(:,:,:,1), loc, fp(:,:,:,1), lof, &
                  mp_fine(:,:,:,1), lom_fine, &
                  mp_crse(:,:,:,1), lom_crse, lo, hi, ir, linject, rmode)
          end select
       end do
    end do

      if ( linject ) then
         call multifab_copy(crse, cfine)
      else
         call multifab_copy(crse, cfine, filter = ml_restrict_copy_sum)
         do i = 1,crse%dim
           if (crse%la%lap%pmask(i)) call periodic_add_copy(crse,cfine,i)
         end do
      end if

    end if

    call destroy(mm_cfine)
    call destroy(cfine)

  end subroutine ml_nodal_restriction

  subroutine periodic_add_copy(dst,src,dir)

      type(multifab), intent(inout) :: dst
      type(multifab), intent(in   ) :: src
      integer       , intent(in   ) :: dir

      type(box)             :: domain,bxi,bxj,bx_lo,bx_hi
      real(dp_t), pointer   :: ap(:,:,:,:), bp(:,:,:,:)
      integer               :: i,j,proc,lod(MAX_SPACEDIM),hid(MAX_SPACEDIM)
      logical               :: nodal(dst%dim)
      integer,    parameter :: tag = 1111

      real(kind=dp_t), dimension(:,:,:,:), allocatable :: pt

      if ( dst%nc .ne. src%nc ) then
         call bl_error('periodic_add_copy: src & dst must have same # of components')
      end if

      nodal = .true.

      domain = box_nodalize(dst%la%lap%pd,nodal)

      do j = 1, dst%nboxes
         !
         ! Add values at hi end of domain to lo end.
         !
         bxj = get_ibox(dst,j)

         if (bxj%lo(dir) == domain%lo(dir)) then
            call box_set_upb_d(bxj,dir,domain%lo(dir))
            do i = 1, src%nboxes
               if ( remote(dst,j) .and. remote(src,i) ) cycle
               bxi = get_ibox(src,i)
               if (bxi%hi(dir) == domain%hi(dir)) then
                  call box_set_lwb_d(bxi,dir,domain%lo(dir))
                  call box_set_upb_d(bxi,dir,domain%lo(dir))
                  bx_lo = box_intersection(bxi,bxj)
                  if (.not. box_empty(bx_lo)) then
                     bx_hi = bx_lo
                     call box_set_lwb_d(bx_hi,dir,domain%hi(dir))
                     call box_set_upb_d(bx_hi,dir,domain%hi(dir))

                     if ( local(dst,j) .and. local(src,i) ) then
                        ap => dataptr(dst,j,bx_lo)
                        bp => dataptr(src,i,bx_hi)
                        ap =  ap + bp
                     else if ( local(src,i) ) then
                        !
                        ! We own src.  Got to send it to processor owning dst.
                        !
                        bp   => dataptr(src,i,bx_hi)
                        proc =  get_proc(dst%la,j)
                        call parallel_send(bp, proc, tag)
                     else
                        !
                        ! We own dst.  Got to get src from processor owning it.
                        !
                        lod = 1; hid = 1
                        lod(1:src%dim) = lwb(bx_hi); hid(1:src%dim) = upb(bx_hi)
                        proc = get_proc(src%la,i)
                        allocate(pt(lod(1):hid(1),lod(2):hid(2),lod(3):hid(3),1:dst%nc))
                        call parallel_recv(pt, proc, tag)
                        ap => dataptr(dst,j,bx_lo)
                        ap =  ap + pt
                        deallocate(pt)
                     end if

                  end if
               end if
            end do
         end if
         !
         ! Add values at lo end of domain to hi end.
         !
         bxj = get_ibox(dst,j)

         if (bxj%hi(dir) == domain%hi(dir)) then
            call box_set_lwb_d(bxj,dir,domain%hi(dir))
            do i = 1, src%nboxes
               if ( remote(dst,j) .and. remote(src,i) ) cycle
               bxi = get_ibox(src,i)
               if (bxi%lo(dir) == domain%lo(dir)) then
                  call box_set_lwb_d(bxi,dir,domain%hi(dir))
                  call box_set_upb_d(bxi,dir,domain%hi(dir))
                  bx_hi = box_intersection(bxi,bxj)
                  if (.not. box_empty(bx_hi)) then
                     bx_lo = bx_hi
                     call box_set_lwb_d(bx_lo,dir,domain%lo(dir))
                     call box_set_upb_d(bx_lo,dir,domain%lo(dir))

                     if ( local(dst,j) .and. local(src,i) ) then
                        ap => dataptr(dst,j,bx_hi)
                        bp => dataptr(src,i,bx_lo)
                        ap =  ap + bp
                     else if ( local(src,i) ) then
                        !
                        ! We own src.  Got to send it to processor owning dst.
                        !
                        bp   => dataptr(src,i,bx_lo)
                        proc =  get_proc(dst%la,j)
                        call parallel_send(bp, proc, tag)
                     else
                        !
                        ! We own dst.  Got to get src from processor owning it.
                        !
                        lod = 1; hid = 1
                        lod(1:src%dim) = lwb(bx_lo); hid(1:src%dim) = upb(bx_lo)
                        proc = get_proc(src%la,i)
                        allocate(pt(lod(1):hid(1),lod(2):hid(2),lod(3):hid(3),1:dst%nc))
                        call parallel_recv(pt, proc, tag)
                        ap => dataptr(dst,j,bx_hi)
                        ap =  ap + pt
                        deallocate(pt)
                     end if

                  end if
               end if
            end do
         end if

      end do
  
 end subroutine periodic_add_copy

 subroutine ml_restriction(crse, fine, mm_fine, mm_crse, face_type, ir, inject, zero_only)
  type(multifab),  intent(inout) :: fine
  type(multifab),  intent(inout) :: crse
  type(imultifab), intent(in   ) :: mm_fine
  type(imultifab), intent(in   ) :: mm_crse
  integer,         intent(in)    :: ir(:)
  integer,         intent(in)    :: face_type(:,:,:)
  logical,         intent(in), optional :: inject
  logical,         intent(in), optional :: zero_only
  if ( crse%nc .ne. fine%nc ) then
     call bl_error('ml_restriction: crse & fine must have same # of components')
  end if
  if ( nodal_q(fine) ) then
     call ml_nodal_restriction(crse, fine, mm_fine, mm_crse, face_type, ir, inject, zero_only)
  else
     call ml_cc_restriction(crse, fine, ir)
  end if
end subroutine ml_restriction

end module ml_restriction_module
