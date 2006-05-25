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
         call periodic_add_copy(crse,cfine)
      end if

    end if

    call destroy(mm_cfine)
    call destroy(cfine)

  end subroutine ml_nodal_restriction

  subroutine periodic_add_copy(dst,src)

    type(multifab), intent(inout) :: dst
    type(multifab), intent(in   ) :: src

    type(box)             :: domain,bxi,bxj,bx_to,bx_from
    type(box)             :: domain_edge_from, domain_edge_to
    real(dp_t), pointer   :: ap(:,:,:,:), bp(:,:,:,:)
    integer               :: i,j,dir,idir,jdir,kdir,proc,lo(MAX_SPACEDIM),hi(MAX_SPACEDIM)
    logical               :: nodal(dst%dim)
    integer               :: shift_vector(3)
    integer,    parameter :: tag = 1111

    real(kind=dp_t), dimension(:,:,:,:), allocatable :: pt

    if ( dst%nc .ne. src%nc ) then
       call bl_error('periodic_add_copy: src & dst must have same # of components')
    end if

    if ( all(dst%la%lap%pmask .eqv. .false.) ) return

    nodal  = .true.

    domain = box_nodalize(dst%la%lap%pd,nodal)

    do kdir = -1,1

       if ( dst%dim < 3  .and. kdir /= 0                                         ) cycle
       if ( dst%dim == 3 .and. (.not. dst%la%lap%pmask(dst%dim)) .and. kdir /= 0 ) cycle

       if ( dst%dim == 3 ) shift_vector(3) = kdir * (box_extent_d(domain,dst%dim) - 1)

       do jdir = -1,1

          if ( .not. dst%la%lap%pmask(2) .and. jdir /= 0 ) cycle

          do idir = -1,1

             if ( .not. dst%la%lap%pmask(1) .and. idir /= 0                      ) cycle
             if ( dst%dim == 2 .and. (idir == 0 .and. jdir == 0)                 ) cycle
             if ( dst%dim == 3 .and. (idir == 0 .and. jdir == 0 .and. kdir == 0) ) cycle

             shift_vector(1) = idir * (box_extent_d(domain,1) - 1)
             shift_vector(2) = jdir * (box_extent_d(domain,2) - 1)

             domain_edge_from = intersection(domain,shift(domain, shift_vector))
             domain_edge_to   = intersection(domain,shift(domain,-shift_vector))

             do j = 1, dst%nboxes
                !
                ! Add values from domain_edge_from side to domain_edge_to side
                !
                bxj = intersection(get_ibox(dst,j),domain_edge_to)

                if ( .not. empty(bxj) ) then
                   do i = 1, src%nboxes
                      if ( remote(dst,j) .and. remote(src,i) ) cycle
                      bxi = intersection(get_ibox(src,i),domain_edge_from)
                      if ( .not. empty(bxi) ) then
                         bxi     = shift(bxi,-shift_vector)
                         bx_from = box_intersection(bxi,bxj)
                         if ( .not. empty(bx_from) ) then
                            bx_to   = bx_from
                            bx_from = shift(bx_from,shift_vector)

                            ! print *,'ADDING FROM BOX ',bx_from%lo(1),bx_from%lo(2),bx_from%hi(1),bx_from%hi(2)
                            ! print *,'         TO BOX ',bx_to%lo(1)  ,bx_to%lo(2)  ,bx_to%hi(1)  ,bx_to%hi(2)

                            if ( local(dst,j) .and. local(src,i) ) then
                               ap => dataptr(dst,j,bx_to)
                               bp => dataptr(src,i,bx_from)
                               ap =  ap + bp
                            else if ( local(src,i) ) then
                               !
                               ! We own src.  Got to send it to processor owning dst.
                               !
                               bp   => dataptr(src,i,bx_from)
                               proc =  get_proc(dst%la,j)
                               call parallel_send(bp, proc, tag)
                            else
                               !
                               ! We own dst.  Got to get src from processor owning it.
                               !
                               lo = 1; hi = 1
                               lo(1:src%dim) = lwb(bx_from); hi(1:src%dim) = upb(bx_from)
                               proc = get_proc(src%la,i)
                               allocate(pt(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:dst%nc))
                               call parallel_recv(pt, proc, tag)
                               ap => dataptr(dst,j,bx_to)
                               ap =  ap + pt
                               deallocate(pt)
                            end if

                         end if
                      end if
                   end do
                end if
             end do
          end do
       end do
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
