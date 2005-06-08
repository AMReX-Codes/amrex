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

    call multifab_build(cfine, lacfine, nc = lnc, ng = 0)

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

  subroutine ml_nodal_restriction(crse, fine, mm_fine, mm_crse, face_type, ir, inject)
    type(multifab),  intent(inout) :: fine
    type(multifab),  intent(inout) :: crse
    type(imultifab), intent(in   ) :: mm_fine
    type(imultifab), intent(in   ) :: mm_crse
    integer,         intent(in)    :: ir(:)
    integer,         intent(in)    :: face_type(:,:,:)
    logical,         intent(in), optional :: inject

    integer             :: i, j, n, id
    integer             :: lo (fine%dim), hi (fine%dim), loc(fine%dim), lof(fine%dim)
    integer             :: lom_fine(fine%dim), lom_crse(fine%dim)
    logical             :: linject
    type(box)           :: fbox, cbox, isect
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    integer,    pointer :: mp_fine(:,:,:,:), mp_crse(:,:,:,:)
    integer             :: rmode
    type(layout)        :: lacfine
    type(multifab)      :: cfine
    type(imultifab)     :: mm_cfine

    if ( crse%nc .ne. fine%nc ) then
       call bl_error('ml_restriction: crse & fine must have same # of components')
    end if

    linject = .false. ; if ( present(inject) ) linject = inject

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
    end if

    call destroy(mm_cfine)
    call destroy(cfine)

  end subroutine ml_nodal_restriction

 subroutine ml_restriction(crse, fine, mm_fine, mm_crse, face_type, ir, inject)
  type(multifab),  intent(inout) :: fine
  type(multifab),  intent(inout) :: crse
  type(imultifab), intent(in   ) :: mm_fine
  type(imultifab), intent(in   ) :: mm_crse
  integer,         intent(in)    :: ir(:)
  integer,         intent(in)    :: face_type(:,:,:)
  logical,         intent(in), optional :: inject
  if ( crse%nc .ne. fine%nc ) then
     call bl_error('ml_restriction: crse & fine must have same # of components')
  end if
  if ( nodal_q(fine) ) then
     call ml_nodal_restriction(crse, fine, mm_fine, mm_crse, face_type, ir, inject)
  else
     call ml_cc_restriction(crse, fine, ir)
  end if
end subroutine ml_restriction

end module ml_restriction_module
