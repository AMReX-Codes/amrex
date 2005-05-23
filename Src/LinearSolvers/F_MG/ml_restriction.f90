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
    integer, intent(in)           :: cc, cf
    integer, intent(in), optional :: nc
    integer, intent(in)           :: ir(:)

    integer             :: i, j, n, lnc, dm, lo(fine%dim), hi(fine%dim)
    integer             :: loc(fine%dim), lof(fine%dim)
    type(box)           :: fbox, cbox
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)

    lnc = 1; if ( present(nc) ) lnc = nc

    dm = crse%dim
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
             do n = 0, lnc - 1
                fp => dataptr(fine, i, cf + n, 1)
                cp => dataptr(crse, j, cc + n, 1)
                select case (dm)
                case (1)
                   call cc_restriction_1d(cp(:,1,1,1), loc, fp(:,1,1,1), lof, lo, hi, ir)
                case (2)
                   call cc_restriction_2d(cp(:,:,1,1), loc, fp(:,:,1,1), lof, lo, hi, ir)
                case (3)
                   call cc_restriction_3d(cp(:,:,:,1), loc, fp(:,:,:,1), lof, lo, hi, ir)
                end select
             end do
          end if
       end do
    end do
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

  subroutine ml_edge_restriction(crse, fine, ir)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer,        intent(in)    :: ir(:)

    integer             :: i, j, n, dm
    integer             :: lo(fine%dim), hi(fine%dim), loc(fine%dim), lof(fine%dim)
    type(box)           :: fbox, cbox
    real(dp_t), pointer :: fp(:,:,:,:)
    real(dp_t), pointer :: cp(:,:,:,:)

    dm = crse%dim
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
             do n = 1, dm
                hi(n) = hi(n)+1
                select case (dm)
                case (1)
                   call edge_restriction_1d(cp(:,1,1,n), loc, fp(:,1,1,n), lof, lo, hi, ir)
                case (2)
                   call edge_restriction_2d(cp(:,:,1,n), loc, fp(:,:,1,n), lof, lo, hi, ir, n)
                case (3)
                   call edge_restriction_3d(cp(:,:,:,n), loc, fp(:,:,:,n), lof, lo, hi, ir, n)
                end select
                hi(n) = hi(n)-1
             end do
          end if
       end do
    end do
  end subroutine ml_edge_restriction

 subroutine ml_restriction(crse, fine, mm_fine, mm_crse, face_type, ir, inject)
  type(multifab),  intent(inout) :: fine
  type(multifab),  intent(inout) :: crse
  type(imultifab), intent(in   ) :: mm_fine
  type(imultifab), intent(in   ) :: mm_crse
  integer,         intent(in)    :: ir(:)
  integer,         intent(in)    :: face_type(:,:,:)
  logical,         intent(in), optional :: inject

  integer             :: i, j, n, id, dm
  integer             :: lo (fine%dim), hi (fine%dim), loc(fine%dim), lof(fine%dim)
  integer             :: lom_fine(fine%dim), lom_crse(fine%dim)
  logical             :: local_inject, nodal_flag
  type(box)           :: fbox, cbox
  real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
  integer,    pointer :: mp_fine(:,:,:,:), mp_crse(:,:,:,:)

  integer :: mg_restriction_mode

  nodal_flag = nodal_q(fine)
  dm = crse%dim
  local_inject = .false. ; if (present(inject)) local_inject = inject

  if ( nodal_flag ) call multifab_fill_boundary(fine)

  mg_restriction_mode = 0
  
  do j = 1, crse%nboxes

    cbox = get_ibox(crse,j)
    loc = lwb(cbox) - crse%ng
    lom_crse = lwb(get_box(mm_crse,j)) - mm_crse%ng

!   Set to zero here on the interior of each fine grid so don't have to 
!      within nodal_restriction
    if ( nodal_flag ) then
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
    end if

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

        cp      => dataptr(crse   ,j)
        mp_crse => dataptr(mm_crse,j)

        do n = 1, 1
          select case (dm)
          case (1)
             if ( .not. nodal_flag ) then
               call cc_restriction_1d(cp(:,1,1,n), loc, fp(:,1,1,n), lof, lo, hi, ir)
             else
               call nodal_restriction_1d(cp(:,1,1,n), loc, fp(:,1,1,n), lof, &
                    mp_fine(:,1,1,1), lom_fine, &
                    mp_crse(:,1,1,1), lom_crse, lo, hi, ir, local_inject, &
                    mg_restriction_mode)
             end if
          case (2)
             if ( .not. nodal_flag ) then
               call cc_restriction_2d(cp(:,:,1,n), loc, fp(:,:,1,n), lof, lo, hi, ir)
             else
               call nodal_restriction_2d(cp(:,:,1,n), loc, fp(:,:,1,n), lof, &
                    mp_fine(:,:,1,1), lom_fine, &
                    mp_crse(:,:,1,1), lom_crse, lo, hi, ir, local_inject, &
                    mg_restriction_mode)
             end if
          case (3)
             if ( .not. nodal_flag ) then
               call cc_restriction_3d(cp(:,:,:,n), loc, fp(:,:,:,n), lof, lo, hi, ir)
             else
               call nodal_restriction_3d(cp(:,:,:,n), loc, fp(:,:,:,n), lof, &
                    mp_fine(:,:,:,1), lom_fine, &
                    mp_crse(:,:,:,1), lom_crse, lo, hi, ir, local_inject, &
                    mg_restriction_mode)
             end if
          end select
        end do
      end if
    end do
  end do

 end subroutine ml_restriction

end module ml_restriction_module
