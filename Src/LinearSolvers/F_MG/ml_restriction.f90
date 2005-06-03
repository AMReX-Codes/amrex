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
  !
  ! Is this used?  If so it needs to be parallelized!!!
  !
  subroutine ml_edge_restriction(crse, fine, ir, n)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer,        intent(in)    :: ir(:)

    integer             :: i, j, n
    integer             :: lo(fine%dim), hi(fine%dim), loc(fine%dim), lof(fine%dim)
    type(box)           :: fbox, cbox
    real(dp_t), pointer :: fp(:,:,:,:)
    real(dp_t), pointer :: cp(:,:,:,:)

!    print *,' *** norm_l1(crse) before ', norm_l1(crse)
!    print *,' *** norm_l2(crse) before ', norm_l2(crse)
!    print *,' *** norm_inf(crse) before ', norm_inf(crse)

    do j = 1, crse%nboxes
       cbox = get_ibox(crse,j)
       loc = lwb(cbox) - crse%ng
       do i = 1, fine%nboxes
          fbox = get_ibox(fine,i)
          lof(:) = lwb(fbox) - fine%ng 
          fbox = box_coarsen_v(fbox,ir)
          if (box_intersects(fbox,cbox)) then
             lo(:) = lwb(box_intersection(cbox,fbox))
             hi(:) = upb(box_intersection(cbox,fbox))
             fp => dataptr(fine, i)
             cp => dataptr(crse, j)
             select case (crse%dim)
             case (1)
                call edge_restriction_1d(cp(:,1,1,1), loc, fp(:,1,1,1), lof, lo, hi, ir)
             case (2)
                call edge_restriction_2d(cp(:,:,1,1), loc, fp(:,:,1,1), lof, lo, hi, ir, n)
             case (3)
                call edge_restriction_3d(cp(:,:,:,1), loc, fp(:,:,:,1), lof, lo, hi, ir, n)
             end select
          end if
       end do
    end do

!    print *,' *** norm_l1(crse) after ', norm_l1(crse)
!    print *,' *** norm_l2(crse) after ', norm_l2(crse)
!    print *,' *** norm_inf(crse) after ', norm_inf(crse)

  end subroutine ml_edge_restriction
  !
  ! TODO -- this needs to be parallelized!!!
  !
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

    if ( crse%nc .ne. fine%nc ) then
       call bl_error('ml_restriction: crse & fine must have same # of components')
    end if

    linject = .false. ; if ( present(inject) ) linject = inject

    call multifab_fill_boundary(fine)

    rmode = 0

  do j = 1, crse%nboxes

    cbox = get_ibox(crse,j)
    loc = lwb(cbox) - crse%ng
    lom_crse = lwb(get_box(mm_crse,j)) - mm_crse%ng

!   Set to zero here on the interior of each fine grid so don't have to 
!      within nodal_restriction
    do i = 1, fine%nboxes
       fbox = get_ibox(fine,i)
       fbox = box_coarsen_v(fbox,ir)
       do id = 1,fbox%dim
          if (face_type(i,id,1) .ne. BC_NEU) fbox = grow(fbox,-1,id,-1)
          if (face_type(i,id,2) .ne. BC_NEU) fbox = grow(fbox,-1,id,+1)
       end do
       if (box_intersects(fbox,cbox)) then
          call setval(crse%fbs(j), ZERO, box_intersection(fbox,cbox))
       end if
    end do

    do i = 1, fine%nboxes

      fbox = get_ibox(fine,i)
      lof(:) = lwb(fbox) - fine%ng 
      fbox = box_coarsen_v(fbox,ir)

      if (box_intersects(fbox,cbox)) then
        lo(:) = lwb(box_intersection(cbox,fbox))
        hi(:) = upb(box_intersection(cbox,fbox))

        fp      => dataptr(fine   ,i)
        mp_fine => dataptr(mm_fine,i)
        lom_fine(:) = lwb(get_box(mm_fine,i)) - mm_fine%ng

        do n = 1, fine%nc
           cp      => dataptr(crse   ,j, n, 1)
           mp_crse => dataptr(mm_crse,j, n, 1)
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
      end if
    end do
  end do

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
