module ml_prolongation_module

  use bl_types
  use multifab_module

  implicit none

  real (kind = dp_t), private, parameter :: ZERO  = 0.0_dp_t
  real (kind = dp_t), private, parameter :: ONE   = 1.0_dp_t
  real (kind = dp_t), private, parameter :: TWO   = 2.0_dp_t
  real (kind = dp_t), private, parameter :: THREE = 3.0_dp_t
  real (kind = dp_t), private, parameter :: FOUR  = 4.0_dp_t
  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t
  real (kind = dp_t), private, parameter :: FOURTH= 0.25_dp_t
  real (kind = dp_t), private, parameter :: EIGHTH= 0.125_dp_t
  real (kind = dp_t), private, parameter :: THREE_EIGHTHS= 0.375_dp_t

contains

 subroutine ml_prolongation(fine, crse, fine_domain, ir)
  type(multifab), intent(inout) :: fine
  type(multifab), intent(in   ) :: crse
  type(box) :: fine_domain, fbox, cbox, fbox_grown, cbox_refined
  integer :: lo (fine%dim), hi (fine%dim)
  integer :: loc(fine%dim)
  integer :: lof(fine%dim)
  integer, intent(in) :: ir(:)
  integer :: i, j, n

  integer :: dm
  real(kind=dp_t), pointer :: fp(:,:,:,:)
  real(kind=dp_t), pointer :: cp(:,:,:,:)

  dm = crse%dim

  do j = 1, crse%nboxes

    cbox = get_ibox(crse,j)
    loc = lwb(cbox) - crse%ng

    if ( nodal_q(fine) ) then
      cbox_refined = make_box(ir*lwb(cbox),ir*upb(cbox))
!     do idim=1,dm
!       cbox%lo(i) = ir(i)*cbox%lo(i)
!       cbox%hi(i) = ir(i)*cbox%hi(i)
!     end do
    else
      cbox_refined = box_refine_v(cbox,ir)
    end if

    do i = 1, fine%nboxes

      fbox   = get_ibox(fine,i)
      lof = lwb(fbox) - fine%ng 

!     Here we grow the fbox so we can fill its ghost cells with values
!          which will be used as fixed boundary values in the fine level relaxation
!     Note that the interpolation must be done in fine index space because
!       otherwise we might try to write values outside of the fine grid array.

      if ( nodal_q(fine) ) then
        fbox_grown = fbox
      else
        fbox_grown = box_intersection(grow(fbox,1),fine_domain)
      end if

      if (box_intersects(fbox_grown,cbox_refined)) then
        lo = lwb(box_intersection(cbox_refined,fbox_grown))
        hi = upb(box_intersection(cbox_refined,fbox_grown))
        fp => dataptr(fine, i)
        cp => dataptr(crse, j)
        do n = 1, 1
          select case (dm)
          case (1)
             if ( nodal_q(fine) ) then
                call ml_prolongation_1d_nodal(fp(:,1,1,n), lof, &
                                              cp(:,1,1,n), loc, &
                                              lo, hi, ir)
             else
                call ml_prolongation_1d_cc(fp(:,1,1,n), lof, &
                                           cp(:,1,1,n), loc, &
                                           lo, hi, ir)
             end if
          case (2)
             if ( nodal_q(fine) ) then
                call ml_prolongation_2d_nodal(fp(:,:,1,n), lof, &
                                              cp(:,:,1,n), loc, &
                                              lo, hi, ir)
             else
                call ml_prolongation_2d_cc(fp(:,:,1,n), lof, &
                                           cp(:,:,1,n), loc, &
                                           lo, hi, ir)
             end if
          case (3)
             if ( nodal_q(fine) ) then
                call ml_prolongation_3d_nodal(fp(:,:,:,n), lof, &
                                              cp(:,:,:,n), loc, &
                                              lo, hi, ir)
             else
                call ml_prolongation_3d_cc(fp(:,:,:,n), lof, &
                                           cp(:,:,:,n), loc, &
                                           lo, hi, ir)
             end if
          end select
        end do
      end if
    end do
  end do
 end subroutine ml_prolongation

  subroutine ml_prolongation_1d_cc(ff, lof, cc, loc, lo, hi, ir)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: ff(lof(1):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):)
    integer, intent(in) :: ir(:)
    integer i, ic

    do i = lo(1),hi(1)
       ic = i / ir(1) 
       ff(i) = ff(i) + cc(ic)
    end do

  end subroutine ml_prolongation_1d_cc

  subroutine ml_prolongation_2d_cc(ff, lof, cc, loc, lo, hi, ir)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: ff(lof(1):,lof(2):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):,loc(2):)
    integer, intent(in) :: ir(:)
    integer i, j, ic, jc

    do j = lo(2),hi(2)
       jc = j / ir(2) 
       do i = lo(1),hi(1)
          ic = i / ir(1) 
          ff(i,j) = ff(i,j) + cc(ic,jc)
       end do
    end do

  end subroutine ml_prolongation_2d_cc

  subroutine ml_prolongation_3d_cc(ff, lof, cc, loc, lo, hi, ir)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: ff(lof(1):,lof(2):,lof(3):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):,loc(2):,loc(3):)
    integer, intent(in) :: ir(:)
    integer i, j, k, ic, jc, kc

    do k = lo(3),hi(3)
       kc = k / ir(3)
       do j = lo(2),hi(2)
          jc = j / ir(2)
          do i = lo(1),hi(1)
             ic = i / ir(1)
             ff(i,j,k) = ff(i,j,k) + cc(ic,jc,kc)
          end do
       end do
    end do

  end subroutine ml_prolongation_3d_cc

  subroutine ml_prolongation_1d_nodal(ff, lof, cc, loc, lo, hi, ir)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: ff(lof(1):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):)
    integer, intent(in) :: ir(:)
    integer i, ic, l
    real (kind = dp_t) :: fac_left, fac_rght

    do i = lo(1),hi(1),ir(1)
       ic = i / ir(1) 
!      ff(i) = ff(i) + cc(ic)
       ff(i) =         cc(ic)
    end do

    do l = 1, ir(1)-1
       fac_rght = real(l,kind=dp_t) / real(ir(1),kind=dp_t)
       fac_left = ONE - fac_rght
       do i = lo(1), hi(1)-1, ir(1)
          ic = i / ir(1) 
!         ff(i+l) = ff(i+l) + fac_left*cc(ic) + fac_rght*cc(ic+1)
          ff(i+l) =           fac_left*cc(ic) + fac_rght*cc(ic+1)
       end do
    end do

  end subroutine ml_prolongation_1d_nodal

  subroutine ml_prolongation_2d_nodal(ff, lof, cc, loc, lo, hi, ir)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: ff(lof(1):,lof(2):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):,loc(2):)
    integer, intent(in) :: ir(:)
    integer i, j, ic, jc, l, m
    real (kind = dp_t) :: fac_left, fac_rght
    real (kind = dp_t) :: temp(lof(1):lof(1)+size(ff,dim=1)-1,&
                               lof(2):lof(2)+size(ff,dim=2)-1)

    do j = lo(2),hi(2),ir(2)
       jc = j / ir(2) 
       do i = lo(1),hi(1),ir(1)
          ic = i / ir(1) 
          temp(i,j) = cc(ic,jc)
       end do
    end do

!   Interpolate at fine nodes between coarse nodes in the i-direction only
    do j = lo(2),hi(2),ir(2)
       do l = 1, ir(1)-1
          fac_rght = real(l,kind=dp_t) / real(ir(1),kind=dp_t)
          fac_left = ONE - fac_rght
          do i = lo(1),hi(1)-1,ir(1)
             temp(i+l,j) = fac_left*temp(i,j) + fac_rght*temp(i+ir(1),j)
          end do
       end do
    end do

!   Interpolate in the j-direction using previously interpolated "temp"
    do m = 1, ir(2)-1
      fac_rght = real(m,kind=dp_t) / real(ir(2),kind=dp_t)
      fac_left = ONE - fac_rght
      do j = lo(2),hi(2)-1,ir(2)
      do i = lo(1),hi(1)
          temp(i,j+m) = fac_left*temp(i,j      ) &
                       +fac_rght*temp(i,j+ir(2))
      end do
      end do
    end do

    do j = lo(2),hi(2)
    do i = lo(1),hi(1)
!     ff(i,j) = ff(i,j) + temp(i,j)
      ff(i,j) =           temp(i,j)
    end do
    end do

  end subroutine ml_prolongation_2d_nodal

  subroutine ml_prolongation_3d_nodal(ff, lof, cc, loc, lo, hi, ir)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: ff(lof(1):,lof(2):,lof(3):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):,loc(2):,loc(3):)
    integer, intent(in) :: ir(:)
    integer i, j, k, ic, jc, kc, l, m, n
    real (kind = dp_t) :: fac_left, fac_rght
    real (kind = dp_t) :: temp(lof(1):lof(1)+size(ff,dim=1)-1,&
                               lof(2):lof(2)+size(ff,dim=2)-1,&
                               lof(3):lof(3)+size(ff,dim=3)-1)

!   Interpolate at coarse node locations only
    do k = lo(3),hi(3),ir(3)
       kc = k / ir(3) 
       do j = lo(2),hi(2),ir(2)
          jc = j / ir(2) 
          do i = lo(1),hi(1),ir(1)
             ic = i / ir(1) 
             temp(i,j,k) = cc(ic,jc,kc)
          end do
       end do
    end do

!   Interpolate at fine nodes between coarse nodes in the i-direction only
    do k = lo(3),hi(3),ir(3)
    do j = lo(2),hi(2),ir(2)
       do l = 1, ir(1)-1
          fac_rght = real(l,kind=dp_t) / real(ir(1),kind=dp_t)
          fac_left = ONE - fac_rght
          do i = lo(1),hi(1)-1,ir(1)
             temp(i+l,j,k) = fac_left*temp(i,j,k) + fac_rght*temp(i+ir(1),j,k)
          end do
       end do
    end do
    end do

!   Interpolate in the j-direction.
    do k = lo(3),hi(3),ir(3)
       do m = 1, ir(2)-1
          fac_rght = real(m,kind=dp_t) / real(ir(2),kind=dp_t)
          fac_left = ONE - fac_rght
          do j = lo(2),hi(2)-1,ir(2)
          do i = lo(1),hi(1)
              temp(i,j+m,k) = fac_left*temp(i,j      ,k) &
                             +fac_rght*temp(i,j+ir(2),k)
          end do
          end do
       end do
    end do

!   Interpolate in the k-direction.
    do n = 1, ir(3)-1
       fac_rght = real(n,kind=dp_t) / real(ir(3),kind=dp_t)
       fac_left = ONE - fac_rght
       do k = lo(3),hi(3)-1,ir(3)
       do j = lo(2),hi(2)
       do i = lo(1),hi(1)
           temp(i,j,k+n) = fac_left*temp(i,j,k) &
                          +fac_rght*temp(i,j,k+ir(3))
       end do
       end do
       end do
    end do

    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)
!     ff(i,j,k) = ff(i,j,k) + temp(i,j,k)
      ff(i,j,k) =             temp(i,j,k)
    end do
    end do
    end do

  end subroutine ml_prolongation_3d_nodal

  subroutine ml_interp_bcs(fine, crse, fine_domain, ir, side)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(in   ) :: crse
    type(box), intent(in) :: fine_domain
    integer, intent(in) :: ir(:)
    type(box) :: fbox, cbox, fbox_grown, cbox_refined
    integer :: lo (fine%dim), hi (fine%dim)
    integer :: loc(fine%dim), hic(fine%dim)
    integer :: lof(fine%dim)
    integer :: dm, side
    integer :: i, n, nc
    integer :: dir,face

    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)

    dm = fine%dim
    nc = fine%nc
    if ( fine%nc /= crse%nc ) then
       call bl_error("ML_INTERP_BCS: fine%nc /= crse%nc")
    end if

    if (side == 1) then
       dir = 1
       face = 1
    else if (side == -1) then
       dir = 1
       face = -1
    else if (side == 2) then
       dir = 2
       face = 1
    else if (side == -2) then
       dir = 2
       face = -1
    else if (side == 3) then
       dir = 3
       face = 1
    else if (side == -3) then
       dir = 3
       face = -1
    end if

    do i = 1, fine%nboxes

       cbox = get_ibox(crse,i)
       loc = lwb(cbox) - crse%ng
       hic = upb(cbox) + crse%ng

       cbox_refined = box_refine_v(cbox,ir)

       fbox   = get_ibox(fine,i)
       lof = lwb(fbox) - fine%ng 

       if ( .not. nodal_q(crse) ) then
          fbox_grown = box_grow_n_d_f(fbox,1,dir,face)
          fbox_grown = box_intersection(fbox_grown,fine_domain)
       else
          fbox_grown = fbox
       end if

       if (box_intersects(fbox_grown,cbox_refined)) then
          lo = lwb(box_intersection(cbox_refined,fbox_grown))
          hi = upb(box_intersection(cbox_refined,fbox_grown))
          fp => dataptr(fine, i)
          cp => dataptr(crse, i)
          do n = 1, nc 
             select case (dm)
             case (1)
                if ( cell_centered_q(crse) ) then
                   call ml_interp_bcs_1d(&
                        fp(:,1,1,n), lof, &
                        cp(:,1,1,n), loc, hic, &
                        lo, hi, ir, side)
                else if ( nodal_q(crse) ) then
                   call bl_error("ML_INTERP_BCS: nodal 1d not provided")
                end if
             case (2)
                if ( cell_centered_q(crse) ) then
                   call ml_interp_bcs_2d(&
                        fp(:,:,1,n), lof, cp(:,:,1,n), loc, hic, &
                        lo, hi, ir, side)
                else if ( nodal_q(crse) ) then
                   call ml_interp_bcs_2d_nodal(&
                        fp(:,:,1,n), lof, cp(:,:,1,n), loc, hic, &
                        lo, hi, ir, side)
                end if
             case (3)
                if ( cell_centered_q(crse) ) then
                   call ml_interp_bcs_3d(&
                        fp(:,:,:,n), lof, cp(:,:,:,n), loc, hic, &
                        lo, hi, ir, side)
                else if ( nodal_q(crse) )  then
                   call ml_interp_bcs_3d_nodal( &
                        fp(:,:,:,n), lof, cp(:,:,:,n), loc, hic, &
                        lo, hi, ir, side)
                end if
             end select
          end do
       end if
    end do

    call multifab_fill_boundary(fine)

  end subroutine ml_interp_bcs

  subroutine ml_interp_bcs_1d(ff, lof, cc, loc, hic, lo, hi, ir, side)
    integer, intent(in) :: loc(:), hic(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: ff(lof(1):)
    real (kind = dp_t), intent(in) :: cc(loc(1):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)
    integer i, ic

    i  = lo(1)
    ic = i / ir(1)
    ff(i) = cc(ic)

  end subroutine ml_interp_bcs_1d

  subroutine ml_interp_bcs_2d(ff, lof, cc, loc, hic, lo, hi, ir, side)
    integer, intent(in) :: loc(:), hic(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: ff(lof(1):,lof(2):)
    real (kind = dp_t), intent(inout) :: cc(loc(1):,loc(2):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)
    integer i, j, ic, jc, m
    real (kind = dp_t) :: xloc(0:3)
    real (kind = dp_t) :: first_der, second_der

    if (ir(1) == 2) then
       xloc(0) = -FOURTH
       xloc(1) =  FOURTH
    else
       xloc(0) = -THREE_EIGHTHS
       xloc(1) = -EIGHTH
       xloc(2) =  EIGHTH
       xloc(3) =  THREE_EIGHTHS
    end if

!   This is the same interpolation as done in IAMR.
    if (side == 1 .or. side == -1) then
      i  = lo(1)
      ic = i / ir(1)
      do jc = lo(2)/ir(2), hi(2)/ir(2)
        first_der  = ZERO
        second_der = ZERO
        if ( jc > loc(2) .and. jc < hic(2) ) then
          first_der  = HALF * (cc(ic,jc+1)-cc(ic,jc-1))
          second_der = (cc(ic,jc+1)+cc(ic,jc-1)-TWO*cc(ic,jc))
        else if ( jc > loc(2) ) then
          if (jc > loc(2)+1) then
            first_der  = HALF*( THREE*cc(ic,jc)-FOUR*cc(ic,jc-1)+cc(ic,jc-2))
            second_der = cc(ic,jc-2) - TWO*cc(ic,jc-1) + cc(ic,jc)
          else
            first_der  =         (cc(ic,jc  )-cc(ic,jc-1))
            second_der = ZERO
          end if
        else if ( jc < hic(2) ) then
          if ( jc < hic(2)-1 ) then
            first_der  = HALF*(-THREE*cc(ic,jc)+FOUR*cc(ic,jc+1)-cc(ic,jc+2))
            second_der = cc(ic,jc+2) - TWO*cc(ic,jc+1) + cc(ic,jc)
          else
            first_der  =        (cc(ic,jc+1)-cc(ic,jc  ))
            second_der = ZERO
          end if
        end if

        j = ir(2)*jc
        do m = 0, ir(2)-1
          ff(i,j+m) = cc(ic,jc) + xloc(m)*first_der + HALF*xloc(m)**2*second_der
        end do
      end do

    else if ( side == 2 .or. side == -2 ) then

      j  = lo(2)
      jc = j / ir(2)
      do ic = lo(1)/ir(1), hi(1)/ir(1)
        first_der  = ZERO
        second_der = ZERO
        if ( ic > loc(1) .and. ic < hic(1) ) then
          first_der  = HALF * (cc(ic+1,jc)-cc(ic-1,jc))
          second_der = (cc(ic+1,jc)+cc(ic-1,jc)-TWO*cc(ic,jc))
        else if ( ic > loc(1) ) then
          if ( ic > loc(1)+1 ) then
            first_der  = HALF*( THREE*cc(ic,jc)-FOUR*cc(ic-1,jc)+cc(ic-2,jc))
            second_der = cc(ic-2,jc) - TWO*cc(ic-1,jc) + cc(ic,jc)
          else
            first_der  = (cc(ic,jc)-cc(ic-1,jc))
            second_der = ZERO
          end if
        else if ( ic < hic(1) ) then
          if (ic < hic(1)-1) then
            first_der  = HALF*(-THREE*cc(ic,jc)+FOUR*cc(ic+1,jc)-cc(ic+2,jc))
            second_der = cc(ic+2,jc) - TWO*cc(ic+1,jc) + cc(ic,jc)
          else
            first_der  = (cc(ic+1,jc)-cc(ic,jc  ))
            second_der = ZERO
          end if
        end if

        i = ir(1)*ic
        do m = 0,ir(1)-1
          ff(i+m,j) = cc(ic,jc) + xloc(m)*first_der + HALF*xloc(m)**2*second_der
        end do
      end do

    end if

  end subroutine ml_interp_bcs_2d

  subroutine ml_interp_bcs_2d_nodal(ff, lof, cc, loc, hic, lo, hi, ir, side)
    integer, intent(in) :: loc(:), hic(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: ff(lof(1):,lof(2):)
    real (kind = dp_t), intent(inout) :: cc(loc(1):,loc(2):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)
    integer i, j, ic, jc, m
    real (kind = dp_t) :: xloc(0:3)
    real (kind = dp_t) :: first_der

    if (ir(1) == 2) then
       xloc(0) = ZERO
       xloc(1) = HALF
    else
       xloc(0) = ZERO
       xloc(1) = FOURTH
       xloc(2) = HALF
       xloc(3) = THREE*FOURTH
    end if

!   This is the same interpolation as done in IAMR.
    if (side == 1 .or. side == -1) then
       i  = lo(1)
       ic = i / ir(1)
       do jc = lo(2)/ir(2), hi(2)/ir(2)-1
          j = ir(2)*jc
          first_der = cc(ic,jc+1)-cc(ic,jc)
          do m = 0,ir(2)-1
             ff(i,j+m) = cc(ic,jc) + xloc(m)*first_der
          end do
          ff(i,hi(2)) = cc(ic,hi(2)/ir(2))
       end do

    else if (side == 2 .or. side == -2) then

       j  = lo(2)
       jc = j / ir(2)
       do ic = lo(1)/ir(1), hi(1)/ir(1)-1
          i = ir(1)*ic
          first_der = cc(ic+1,jc)-cc(ic,jc)
          do m = 0,ir(1)-1
             ff(i+m,j) = cc(ic,jc) + xloc(m)*first_der
          end do
          ff(hi(1),j) = cc(hi(1)/ir(1),jc)
       end do

    end if

  end subroutine ml_interp_bcs_2d_nodal

  subroutine ml_interp_bcs_3d(ff, lof, cc, loc, hic, lo, hi, ir, side)
    integer, intent(in) :: loc(:), hic(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: ff(lof(1):,lof(2):,lof(3):)
    real (kind = dp_t), intent(in) :: cc(loc(1):,loc(2):,loc(3):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)
    integer i, j, k, ic, jc, kc, m, n
    real (kind = dp_t) :: xloc(0:3)
    real (kind = dp_t) :: xder,x2der,yder,y2der,zder,z2der
    real (kind = dp_t) :: xyder,yzder,xzder

    if (ir(1) == 2) then
       xloc(0) = -FOURTH
       xloc(1) =  FOURTH
    else
       xloc(0) = -THREE_EIGHTHS
       xloc(1) = -EIGHTH
       xloc(2) =  EIGHTH
       xloc(3) =  THREE_EIGHTHS
    end if

!   This is the same interpolation as done in IAMR.
    if (side == 1 .or. side == -1) then
      i  = lo(1)
      ic = i / ir(1)

      do kc = lo(3)/ir(3), hi(3)/ir(3)
      do jc = lo(2)/ir(2), hi(2)/ir(2)
        if (jc > loc(2) .and. jc < hic(2)) then
          yder  = HALF  * (cc(ic,jc+1,kc) - cc(ic,jc-1,kc))
          y2der = HALF  * (cc(ic,jc+1,kc) + cc(ic,jc-1,kc) - TWO*cc(ic,jc,kc))
        else if (jc < hic(2)) then
          yder  = cc(ic,jc+1,kc) - cc(ic,jc,kc)
          y2der = ZERO
        else if (jc > loc(2)) then
          yder  = cc(ic,jc,kc) - cc(ic,jc-1,kc)
          y2der = ZERO
        end if

        if (kc > loc(3) .and. kc < hic(3)) then
          zder  = HALF  * (cc(ic,jc,kc+1) - cc(ic,jc,kc-1))
          z2der = HALF  * (cc(ic,jc,kc+1) + cc(ic,jc,kc-1) - TWO*cc(ic,jc,kc))
        else if (kc < hic(3)) then
          zder  = cc(ic,jc,kc+1) - cc(ic,jc,kc)
          z2der = ZERO
        else if (kc > loc(3)) then
          zder  = cc(ic,jc,kc) - cc(ic,jc,kc-1)
          z2der = ZERO
        end if

        if (jc > loc(2) .and. jc < hic(2) .and. kc > loc(3) .and. kc < hic(3)) then
          yzder = FOURTH * (cc(ic,jc+1,kc+1) + cc(ic,jc-1,kc-1) &
                           -cc(ic,jc-1,kc+1) - cc(ic,jc+1,kc-1) )
        else
          yzder = ZERO
        end if

        j = ir(2)*jc
        k = ir(3)*kc
        do n = 0,ir(3)-1
        do m = 0,ir(2)-1
          ff(i,j+m,k+n) = cc(ic,jc,kc) + xloc(m)*yder + xloc(n)*zder &
                                       + xloc(m)**2*y2der + xloc(n)**2*z2der &
                                       + xloc(m)*xloc(n)*yzder
        end do
        end do
      end do
      end do

    else if (side == 2 .or. side == -2) then

      j  = lo(2)
      jc = j / ir(2)
      do kc = lo(3)/ir(3), hi(3)/ir(3)
      do ic = lo(1)/ir(1), hi(1)/ir(1)
        if (ic > loc(1) .and. ic < hic(1)) then
          xder  = HALF  * (cc(ic+1,jc,kc) - cc(ic-1,jc,kc))
          x2der = HALF  * (cc(ic+1,jc,kc) + cc(ic-1,jc,kc) - TWO*cc(ic,jc,kc))
        else if (ic < hic(1)) then
          xder  = cc(ic+1,jc,kc) - cc(ic,jc,kc)
          x2der = ZERO
        else if (ic > loc(1)) then
          xder  = cc(ic,jc,kc) - cc(ic-1,jc,kc)
          x2der = ZERO
        end if

        if (kc > loc(3) .and. kc < hic(3)) then
          zder  = HALF  * (cc(ic,jc,kc+1) - cc(ic,jc,kc-1))
          z2der = HALF  * (cc(ic,jc,kc+1) + cc(ic,jc,kc-1) - TWO*cc(ic,jc,kc))
        else if (kc < hic(3)) then
          zder  = cc(ic,jc,kc+1) - cc(ic,jc,kc)
          z2der = ZERO
        else if (kc > loc(3)) then
          zder  = cc(ic,jc,kc) - cc(ic,jc,kc-1)
          z2der = ZERO
        end if

        if (ic > loc(1) .and. ic < hic(1) .and. kc > loc(3) .and. kc < hic(3)) then
          xzder = FOURTH * (cc(ic+1,jc,kc+1) + cc(ic-1,jc,kc-1) &
                           -cc(ic-1,jc,kc+1) - cc(ic+1,jc,kc-1) )
        else
          xzder = ZERO
        end if

        i = ir(1)*ic
        k = ir(3)*kc
        do n = 0,ir(3)-1
        do m = 0,ir(1)-1
          ff(i+m,j,k+n) = cc(ic,jc,kc) + xloc(m)*xder + xloc(n)*zder &
                                       + xloc(m)**2*x2der + xloc(n)**2*z2der &
                                       + xloc(m)*xloc(n)*xzder
        end do
        end do
      end do
      end do

    else if (side == 3 .or. side == -3) then

      k  = lo(3)
      kc = k / ir(3)
      do jc = lo(2)/ir(2), hi(2)/ir(2)
      do ic = lo(1)/ir(1), hi(1)/ir(1)
        if (jc > loc(2) .and. jc < hic(2)) then
          yder  = HALF  * (cc(ic,jc+1,kc) - cc(ic,jc-1,kc))
          y2der = HALF  * (cc(ic,jc+1,kc) + cc(ic,jc-1,kc) - TWO*cc(ic,jc,kc))
        else if (jc < hic(2)) then
          yder  = cc(ic,jc+1,kc) - cc(ic,jc,kc)
          y2der = ZERO
        else if (jc > loc(2)) then
          yder  = cc(ic,jc,kc) - cc(ic,jc-1,kc)
          y2der = ZERO
        end if

        if (ic > loc(1) .and. ic < hic(1)) then
          xder  = HALF  * (cc(ic+1,jc,kc) - cc(ic-1,jc,kc))
          x2der = HALF  * (cc(ic+1,jc,kc) + cc(ic-1,jc,kc) - TWO*cc(ic,jc,kc))
        else if (ic < hic(1)) then
          xder  = cc(ic+1,jc,kc) - cc(ic,jc,kc)
          x2der = ZERO
        else if (ic > loc(1)) then
          xder  = cc(ic,jc,kc) - cc(ic-1,jc,kc)
          x2der = ZERO
        end if

        if (jc > loc(2) .and. jc < hic(2) .and. ic > loc(1) .and. ic < hic(1)) then
          xyder = FOURTH * (cc(ic+1,jc+1,kc) + cc(ic-1,jc-1,kc) &
                           -cc(ic+1,jc-1,kc) - cc(ic-1,jc+1,kc) )
        else
          xyder = ZERO
        end if

        j = ir(2)*jc
        i = ir(1)*ic
        do n = 0,ir(1)-1
        do m = 0,ir(2)-1
          ff(i+n,j+m,k) = cc(ic,jc,kc) + xloc(m)*yder + xloc(n)*xder &
                                       + xloc(m)**2*y2der + xloc(n)**2*x2der &
                                       + xloc(m)*xloc(n)*xyder
        end do
        end do
      end do
      end do

    end if

  end subroutine ml_interp_bcs_3d

  subroutine ml_interp_bcs_3d_nodal(ff, lof, cc, loc, hic, lo, hi, ir, side)
    integer, intent(in) :: loc(:), hic(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: ff(lof(1):,lof(2):,lof(3):)
    real (kind = dp_t), intent(in) :: cc(loc(1):,loc(2):,loc(3):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)
    integer i, j, k, ic, jc, kc, m, n
    real (kind = dp_t) :: xloc(0:3)
    real (kind = dp_t) :: xder, yder, zder
    real (kind = dp_t) :: xyder,yzder,xzder

    if (ir(1) == 2) then
       xloc(0) = ZERO
       xloc(1) = HALF
    else
       xloc(0) = ZERO
       xloc(1) = FOURTH
       xloc(2) = HALF
       xloc(3) = THREE*FOURTH
    end if

!   This is the same interpolation as done in IAMR.
    if (side == 1 .or. side == -1) then
      i  = lo(1)
      ic = i / ir(1)

      do kc = lo(3)/ir(3), hi(3)/ir(3)-1
         do jc = lo(2)/ir(2), hi(2)/ir(2)-1

           yder  = cc(ic,jc+1,kc) - cc(ic,jc,kc)
           zder  = cc(ic,jc,kc+1) - cc(ic,jc,kc)
           yzder = cc(ic,jc,kc) + cc(ic,jc+1,kc+1) - cc(ic,jc+1,kc) - cc(ic,jc,kc+1)

           j = ir(2)*jc
           k = ir(3)*kc
           do n = 0,ir(3)-1
           do m = 0,ir(2)-1
             ff(i,j+m,k+n) = cc(ic,jc,kc) + xloc(m)*yder + xloc(n)*zder &
                                          + xloc(m)*xloc(n)*yzder
           end do
           end do
         end do

         jc = hi(2)/ir(2)
         j  = hi(2)
         zder  = cc(ic,jc,kc+1) - cc(ic,jc,kc)
         do m = 0,ir(3)-1
             ff(i,j,k+m) = cc(ic,jc,kc) + xloc(m)*zder
         end do

      end do

      kc = hi(3)/ir(3)
      k  = hi(3)
      do jc = lo(2)/ir(2), hi(2)/ir(2)-1
        yder  = cc(ic,jc+1,kc) - cc(ic,jc,kc)
        do m = 0,ir(2)-1
          ff(i,j+m,k) = cc(ic,jc,kc) + xloc(m)*yder
        end do
      end do

      ff(i,hi(2),hi(3)) = cc(ic,hi(2)/ir(2),hi(3)/ir(3))

    else if (side == 2 .or. side == -2) then

      j  = lo(2)
      jc = j / ir(2)

      do kc = lo(3)/ir(3), hi(3)/ir(3)-1
         do ic = lo(1)/ir(1), hi(1)/ir(1)-1

           xder  = cc(ic+1,jc,kc) - cc(ic,jc,kc)
           zder  = cc(ic,jc,kc+1) - cc(ic,jc,kc)
           xzder = cc(ic,jc,kc) + cc(ic+1,jc,kc+1) - cc(ic+1,jc,kc) - cc(ic,jc,kc+1)

           i = ir(1)*ic
           k = ir(3)*kc
           do n = 0,ir(3)-1
           do m = 0,ir(1)-1
             ff(i+m,j,k+n) = cc(ic,jc,kc) + xloc(m)*xder + xloc(n)*zder &
                                          + xloc(m)*xloc(n)*xzder
           end do
           end do
         end do

         ic = hi(1)/ir(1)
         i  = hi(1)
         zder  = cc(ic,jc,kc+1) - cc(ic,jc,kc)
         do m = 0,ir(3)-1
             ff(i,j,k+m) = cc(ic,jc,kc) + xloc(m)*zder
         end do
      end do

      kc = hi(3)/ir(3)
      k  = hi(3)
      do ic = lo(1)/ir(1), hi(1)/ir(1)-1
        xder  = cc(ic+1,jc,kc) - cc(ic,jc,kc)
        do m = 0,ir(1)-1
          ff(i+m,j,k) = cc(ic,jc,kc) + xloc(m)*xder
        end do
      end do

    else if (side == 3 .or. side == -3) then

      k  = lo(3)
      kc = k / ir(3)

      do jc = lo(2)/ir(2), hi(2)/ir(2)-1
         do ic = lo(1)/ir(1), hi(1)/ir(1)-1

           xder  = cc(ic+1,jc,kc) - cc(ic,jc,kc)
           yder  = cc(ic,jc+1,kc) - cc(ic,jc,kc)
           xyder = cc(ic,jc,kc) + cc(ic+1,jc+1,kc) - cc(ic+1,jc,kc) - cc(ic,jc+1,kc)

           i = ir(1)*ic
           j = ir(2)*jc
           do n = 0,ir(2)-1
           do m = 0,ir(1)-1
             ff(i+m,j+n,k) = cc(ic,jc,kc) + xloc(m)*xder + xloc(n)*yder &
                                          + xloc(m)*xloc(n)*xyder
           end do
           end do
         end do

         ic = hi(1)/ir(1)
         i  = hi(1)
         yder  = cc(ic,jc+1,kc) - cc(ic,jc,kc)
         do m = 0,ir(2)-1
             ff(i,j+m,k) = cc(ic,jc,kc) + xloc(m)*yder
         end do
      end do

      jc = hi(2)/ir(2)
      j  = hi(2)
      do ic = lo(1)/ir(1), hi(1)/ir(1)-1
        xder  = cc(ic+1,jc,kc) - cc(ic,jc,kc)
        do m = 0,ir(1)-1
          ff(i+m,j,k) = cc(ic,jc,kc) + xloc(m)*xder
        end do
      end do

    end if

  end subroutine ml_interp_bcs_3d_nodal

end module ml_prolongation_module
