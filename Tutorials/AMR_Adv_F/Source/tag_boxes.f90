module tag_boxes_module

  use multifab_module
  use bl_error_module

  implicit none 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MUST set this to .true. if tagging uses ghost cells (e.g., tagging on gradient). !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical, save :: tagging_needs_ghost_cells = .false.

contains

  subroutine tag_boxes(mf,tagboxes,dx,lev,aux_tag_mf)

    type( multifab)         , intent(in   ) :: mf
    type(lmultifab)         , intent(inout) :: tagboxes
    real(dp_t)              , intent(in   ) :: dx
    integer                 , intent(in   ) :: lev
    type(multifab), optional, intent(in   ) :: aux_tag_mf
    ! aux_tag_mf allows user to pass in additional multifabs for tagging logic
    
    ! local variables
    real(kind = dp_t), pointer :: mfp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i
    integer           :: lo(get_dim(mf)), hi(get_dim(mf))
    integer           :: tlo(4), mflo(4)
    type(mfiter)      :: mfi
    type(box)         :: bx

    if (present(aux_tag_mf)) then
       call bl_error("tag_boxes.f90: aux_tag_mf passed to tag_boxes without implementation")
    end if

    !$omp parallel private(mfp,tp,i,lo,hi,mfi,bx,tlo,mflo)
    call mfiter_build(mfi,tagboxes,.true.)
    do while(next_tile(mfi,i))
       bx = get_tilebox(mfi)
       lo =  lwb(bx)
       hi =  upb(bx)

       mfp => dataptr(mf, i)
       tp  => dataptr(tagboxes, i)

       mflo = lbound(mfp)
       tlo = lbound(tp)

       select case (get_dim(mf))
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),tlo(1:2),mfp(:,:,1,1),mflo(1:2),lo,hi,dx,lev)
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),tlo(1:3),mfp(:,:,:,1),mflo(1:3),lo,hi,dx,lev)
       end select
    end do
    !$omp end parallel

  end subroutine tag_boxes

  subroutine tag_boxes_2d(tagbox,tlo,mf,mflo,lo,hi,dx,lev)

    integer          , intent(in   ) :: lo(2),hi(2),tlo(2), mflo(2)
    logical          , intent(inout) :: tagbox( tlo(1):, tlo(2):)
    real(kind = dp_t), intent(in   ) ::     mf(mflo(1):,mflo(2):)
    real(dp_t)       , intent(in   ) :: dx
    integer          , intent(in   ) :: lev

    ! local variables
    integer :: i,j

    ! initially say that we do not want to tag any cells for refinement
    tagbox(lo(1):hi(1),lo(2):hi(2)) = .false.

    select case(lev)
    case (1)
       ! tag all boxes where the first component of mf >= 1.01
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (mf(i,j) .gt. 1.01_dp_t) then
                tagbox(i,j) = .true.
             end if
          end do
       enddo
    case (2)
       ! for level 2 tag all boxes where the first component of mf >= 1.1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (mf(i,j) .gt. 1.1_dp_t) then
                tagbox(i,j) = .true.
             end if
          end do
       end do
    case default
       ! for level 3 or greater tag all boxes where the first component of mf >= 1.5
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (mf(i,j) .gt. 1.5_dp_t) then
                tagbox(i,j) = .true.
             end if
          end do
       end do
    end select

  end subroutine tag_boxes_2d

  subroutine tag_boxes_3d(tagbox,tlo,mf,mflo,lo,hi,dx,lev)

    integer          , intent(in   ) :: lo(3),hi(3),tlo(3),mflo(3)
    logical          , intent(inout) :: tagbox( tlo(1):, tlo(2):, tlo(3):)
    real(kind = dp_t), intent(in   ) ::     mf(mflo(1):,mflo(2):,mflo(3):)
    real(dp_t)       , intent(in   ) :: dx
    integer          , intent(in   ) :: lev

    ! local variables
    integer :: i,j,k

    ! initially say that we do not want to tag any cells for refinement
    tagbox(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = .false.

    select case(lev)
    case (1)
       ! tag all boxes where the first component of mf >= 1.01
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j,k) .gt. 1.01_dp_t) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          enddo
       end do
    case (2)
       ! for level 2 tag all boxes where the first component of mf >= 1.1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j,k) .gt. 1.1_dp_t) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          end do
       end do
    case default
       ! for level 3 or greater tag all boxes where the first component of mf >= 1.5
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j,k) .gt. 1.5_dp_t) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          end do
       end do
    end select

  end subroutine tag_boxes_3d

end module tag_boxes_module
