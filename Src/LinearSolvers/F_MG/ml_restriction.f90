module ml_restriction_module

  use bl_types
  use bl_prof_module
  use multifab_module
  use mg_restriction_module

  implicit none

  real(dp_t), private, parameter :: ZERO = 0.0_dp_t
  real(dp_t), private, parameter :: ONE  = 1.0_dp_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
       lof = lwb(get_pbox(fine, i))
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

    call destroy(cfine)

    call multifab_fill_boundary(crse)

  end subroutine ml_cc_restriction_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ml_cc_restriction(crse, fine, ir)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer,        intent(in)    :: ir(:)
    if ( crse%nc .ne. fine%nc ) then
       call bl_error('ml_cc_restriction: crse & fine must have same # of components')
    end if
    call ml_cc_restriction_c(crse, 1, fine, 1, ir, crse%nc)
  end subroutine ml_cc_restriction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ml_edge_restriction_c(crse, cc, fine, cf, ir, face, nc)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer,        intent(in)    :: cc, cf, ir(:)
    integer,        intent(in)    :: face
    integer, intent(in), optional :: nc

    integer             :: i, n, lnc
    integer             :: lo(fine%dim), hi(fine%dim), loc(fine%dim), lof(fine%dim)
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    type(layout)        :: lacfine
    type(multifab)      :: cfine

    lnc = 1; if ( present(nc) ) lnc = nc

    call layout_build_coarse(lacfine, fine%la, ir)

    call multifab_build(cfine, lacfine, nc = crse%nc, ng = 0, nodal = crse%nodal)

    do i = 1, fine%nboxes
       if ( remote(fine,i) ) cycle
       lo  = lwb(get_ibox(cfine,i))
       hi  = upb(get_ibox(cfine,i))
       loc = lwb(get_pbox(cfine,i))
       lof = lwb(get_pbox(fine, i))
       do n = 1, lnc
          fp  => dataptr(fine,  i, n+cf-1, 1)
          cp  => dataptr(cfine, i, n,      1)
          select case (crse%dim)
          case (1)
             call edge_restriction_1d(cp(:,1,1,1), loc, fp(:,1,1,1), lof, lo, hi, ir)
          case (2)
             call edge_restriction_2d(cp(:,:,1,1), loc, fp(:,:,1,1), lof, lo, hi, ir, face)
          case (3)
             call edge_restriction_3d(cp(:,:,:,1), loc, fp(:,:,:,1), lof, lo, hi, ir, face)
          end select
       enddo
    end do

    call copy(crse, cc, cfine, 1, lnc)

    call destroy(cfine)

  end subroutine ml_edge_restriction_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ml_edge_restriction(crse, fine, ir, face)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer,        intent(in)    :: ir(:)
    integer,        intent(in)    :: face

    if ( crse%nc .ne. fine%nc ) then
       call bl_error('ml_edge_restriction: crse & fine must have same # of components')
    end if
    call ml_edge_restriction_c(crse, 1, fine, 1, ir, face, crse%nc)

  end subroutine ml_edge_restriction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ml_restrict_copy_sum(out, in)
     use bl_types
     real(dp_t), intent(inout) :: out(:,:,:,:)
     real(dp_t), intent(in   ) ::  in(:,:,:,:)
     integer                   :: i, j, k, n
     !
     ! out = out + in 
     !
     do n = 1, size(out,4)
        do k = 1, size(out,3)
           do j = 1, size(out,2)
              do i = 1, size(out,1)
                 out(i,j,k,n) = out(i,j,k,n) + in(i,j,k,n)
              end do
           end do
        end do
     end do
  end subroutine ml_restrict_copy_sum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ml_nodal_restriction(crse, fine, mm_fine, mm_crse, face_type, ir, inject, zero_only)
    type(multifab),  intent(inout)        :: crse
    type(multifab),  intent(inout)        :: fine
    type(imultifab), intent(in   )        :: mm_fine
    type(imultifab), intent(in   )        :: mm_crse
    integer,         intent(in)           :: ir(:)
    integer,         intent(in)           :: face_type(:,:,:)
    logical,         intent(in), optional :: inject
    logical,         intent(in), optional :: zero_only

    integer             :: i, n, rmode
    integer             :: lo (fine%dim), hi (fine%dim), loc(fine%dim), lof(fine%dim)
    integer             :: lom_fine(fine%dim), lom_crse(fine%dim)
    logical             :: linject, lzero_only
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    integer,    pointer :: mp_fine(:,:,:,:), mp_crse(:,:,:,:)
    type(layout)        :: lacfine
    type(multifab)      :: cfine
    type(imultifab)     :: mm_cfine

    if ( crse%nc .ne. fine%nc ) then
       call bl_error('ml_nodal_restriction: crse & fine must have same # of components')
    end if

    linject    = .false. ; if ( present(inject   ) ) linject    = inject
    lzero_only = .false. ; if ( present(zero_only) ) lzero_only = zero_only

    call layout_build_coarse(lacfine, fine%la, ir)
    call  multifab_build(cfine, lacfine, nc = crse%nc, ng = 0, nodal = crse%nodal)
    call copy(cfine, crse)

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
       rmode = 0
       call imultifab_build(mm_cfine, lacfine, nc = mm_crse%nc, ng = 0, nodal = mm_crse%nodal)
       call copy(mm_cfine, mm_crse)
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

       call destroy(mm_cfine)

       if ( linject ) then
          call multifab_copy(crse, cfine)
       else
          call multifab_copy(crse, cfine, filter = ml_restrict_copy_sum)
          call periodic_add_copy(crse,cfine,synced=.false.)
       end if
    end if

    call destroy(cfine)

  end subroutine ml_nodal_restriction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !
  ! TODO - cache the communication pattern here in the dst layout?
  !
  subroutine periodic_add_copy(dst,src,synced)

    type(multifab), intent(inout) :: dst
    type(multifab), intent(inout) :: src    ! Not changed except by layout_get_box_intersector()
    logical,        intent(in   ) :: synced
    !
    ! if ( synced )
    !
    !   This version assumes that src IS synced up on each edge to start with  - for example,
    !   if a node in grid A on the lo-x side has value a, and the equivalent node in grid B 
    !   but on the hi-x side has value b, and the equivalent node in grid C on the hi-x side
    !   has value b also, then the final value of the nodes in A,B and C will be (a+b)
    !
    ! else
    !
    !   This version assumes that src is NOT synced up on each edge to start with  - for example,
    !             if a node in grid A on the lo-x side has value a and 
    !   the equivalent node in grid B on the hi-x side has value b and
    !   the equivalent node in grid C on the hi-x side has value c ,
    !   then the final value of the node in each grid A, B and C will be (a+b+c)
    !
    type(multifab)        :: temp_dst
    type(box)             :: domain,bxi,bxj,bx_to,bx_from
    type(box)             :: domain_edge_src, domain_edge_dst
    real(dp_t), pointer   :: ap(:,:,:,:)
    real(dp_t), pointer   :: bp(:,:,:,:)
    integer               :: i,j,ii,jj,idir,jdir,kdir,proc,lo(MAX_SPACEDIM),hi(MAX_SPACEDIM),dm
    logical               :: nodal(dst%dim)
    integer               :: shift_vector(3)
    integer,  parameter   :: tag = 1111

    type(box_intersector), pointer :: bisrc(:), bidst(:)

    real(kind=dp_t), dimension(:,:,:,:), allocatable :: pt

    type(bl_prof_timer), save :: bpt

    if ( dst%nc .ne. src%nc ) then
       call bl_error('periodic_add_copy: src & dst must have same # of components')
    end if

    if ( all(dst%la%lap%pmask .eqv. .false.) ) return

    call build(bpt, "periodic_add_copy")

    dm     = dst%dim
    lo     = 1
    hi     = 1
    nodal  = .true.
    domain = box_nodalize(dst%la%lap%pd,nodal)

    if ( synced ) call multifab_build(temp_dst,dst%la,dst%nc,dst%ng,nodal)

    do kdir = -1,1

       if ( dm < 3  .and. kdir /= 0                                    ) cycle
       if ( dm == 3 .and. (.not. dst%la%lap%pmask(dm)) .and. kdir /= 0 ) cycle

       if ( dm == 3 ) shift_vector(3) = kdir * (extent(domain,dm) - 1)

       do jdir = -1,1

          if ( .not. dst%la%lap%pmask(2) .and. jdir /= 0 ) cycle

          shift_vector(2) = jdir * (extent(domain,2) - 1)

          do idir = -1,1

             if ( .not. dst%la%lap%pmask(1) .and. idir /= 0                 ) cycle
             if ( dm == 2 .and. (idir == 0 .and. jdir == 0)                 ) cycle
             if ( dm == 3 .and. (idir == 0 .and. jdir == 0 .and. kdir == 0) ) cycle

             shift_vector(1) = idir * (extent(domain,1) - 1)

             domain_edge_src = intersection(domain,shift(domain, shift_vector))
             domain_edge_dst = intersection(domain,shift(domain,-shift_vector))

             if ( synced ) call setval(temp_dst,ZERO)
             !
             ! Add values from domain_edge_src side to domain_edge_dst side
             !
             bidst => layout_get_box_intersector(dst%la, domain_edge_dst)

             do jj = 1, size(bidst)
                j     =  bidst(jj)%i
                bxj   =  bidst(jj)%bx
                bisrc => layout_get_box_intersector(src%la, domain_edge_src)
                do ii = 1, size(bisrc)
                   i = bisrc(ii)%i
                   if ( remote(dst,j) .and. remote(src,i) ) cycle
                   bxi     = shift(bisrc(ii)%bx,-shift_vector)
                   bx_from = intersection(bxi,bxj)
                   if ( empty(bx_from) ) cycle
                   bx_to   = bx_from
                   bx_from = shift(bx_from,shift_vector)
                   if ( local(dst,j) .and. local(src,i) ) then
                      if ( synced ) then
                         ap => dataptr(temp_dst,j,bx_to)
                         bp => dataptr(src,i,bx_from)
                         ap =  bp
                      else
                         ap => dataptr(dst,j,bx_to)
                         bp => dataptr(src,i,bx_from)
                         ap =  ap + bp
                      end if
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
                      lo(1:src%dim) = lwb(bx_from); hi(1:src%dim) = upb(bx_from)
                      proc = get_proc(src%la,i)
                      allocate(pt(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:dst%nc))
                      call parallel_recv(pt, proc, tag)
                      if ( synced ) then
                         ap => dataptr(temp_dst,j,bx_to)
                         ap =  pt
                      else
                         ap => dataptr(dst,j,bx_to)
                         ap =  ap + pt
                      end if
                      deallocate(pt)
                   end if
                end do
                deallocate(bisrc)
             end do
             deallocate(bidst)
             if ( synced ) call saxpy(dst,ONE,temp_dst)
          end do
       end do
    end do
 
    if ( synced ) call destroy(temp_dst)

    call destroy(bpt)

  end subroutine periodic_add_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ml_restriction(crse, fine, mm_fine, mm_crse, face_type, ir, inject, zero_only)
    type(multifab),  intent(inout) :: fine
    type(multifab),  intent(inout) :: crse
    type(imultifab), intent(in   ) :: mm_fine
    type(imultifab), intent(in   ) :: mm_crse
    integer,         intent(in)    :: ir(:)
    integer,         intent(in)    :: face_type(:,:,:)
    logical,         intent(in), optional :: inject
    logical,         intent(in), optional :: zero_only
    type(bl_prof_timer), save :: bpt
    if ( crse%nc .ne. fine%nc ) then
       call bl_error('ml_restriction: crse & fine must have same # of components')
    end if
    call build(bpt, "ml_restriction")
    if ( nodal_q(fine) ) then
       call ml_nodal_restriction(crse, fine, mm_fine, mm_crse, face_type, ir, inject, zero_only)
    else
       call ml_cc_restriction(crse, fine, ir)
    end if
    call destroy(bpt)
 end subroutine ml_restriction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ml_restriction_module
