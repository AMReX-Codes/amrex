module ml_interface_stencil_module

  use bl_types
  use multifab_module
  use stencil_module

  implicit none

  real (kind = dp_t), private, parameter :: ZERO  = 0.0_dp_t
  real (kind = dp_t), private, parameter :: ONE   = 1.0_dp_t
  real (kind = dp_t), private, parameter :: TWO   = 2.0_dp_t
  real (kind = dp_t), private, parameter :: THREE = 3.0_dp_t
  real (kind = dp_t), private, parameter :: FOUR  = 4.0_dp_t
  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t
  real (kind = dp_t), private, parameter :: THIRD = 0.3333333333333333_dp_t
  real (kind = dp_t), private, parameter :: FOURTH= 0.25_dp_t
  real (kind = dp_t), private, parameter :: EIGHTH= 0.125_dp_t

!   interface ml_stencil_interface
!      module procedure ml_interface_1d
!      module procedure ml_interface_2d
!      module procedure ml_interface_3d
!   end interface

contains

 subroutine ml_interface(res, flux, crse, ss, crse_domain, face, dim)
  type(multifab), intent(inout) :: res
  type(multifab), intent(in   ) :: flux
  type(multifab), intent(in   ) :: crse
  type(multifab), intent(in   ) :: ss
  type(box), intent(in) :: crse_domain
  type(box) :: rbox, fbox, cbox, sbox
  integer :: lo (res%dim), hi (res%dim)
  integer :: loc(res%dim)
  integer :: lof(res%dim)
  integer :: lor(res%dim)
  integer :: los(res%dim)
  integer :: lo_dom(res%dim), hi_dom(res%dim)
  integer :: dm
  integer :: face, dim
  integer :: dir
  integer :: i, j, n

  real(kind=dp_t), pointer :: rp(:,:,:,:)
  real(kind=dp_t), pointer :: fp(:,:,:,:)
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  real(kind=dp_t), pointer :: sp(:,:,:,:)

  dm = res%dim
  dir = dim

  lo_dom = lwb(crse_domain)
  hi_dom = upb(crse_domain)

  do j = 1, crse%nboxes

    cbox = get_ibox(crse,j)
    loc = lwb(cbox) - crse%ng

    rbox = get_ibox(res,j)
    lor = lwb(rbox) - res%ng

    sbox = get_ibox(ss,j)
    los = lwb(sbox) - ss%ng

    do i = 1, flux%nboxes

       fbox   = get_ibox(flux,i)
       lof = lwb(fbox)

       if (box_contains(crse_domain,fbox)) then
        if (.not. empty(box_intersection(cbox,fbox))) then

          lo = lwb(box_intersection(cbox,fbox))
          hi = upb(box_intersection(cbox,fbox))
          fp => dataptr(flux, i)
          cp => dataptr(crse, j)
          rp => dataptr(res, j)
          sp => dataptr(ss, j)
          do n = 1, 1
           select case (dm)
           case (1)
             call ml_interface_1d(rp(:,1,1,n), lor, &
                                  fp(:,1,1,n), lof, &
                                  cp(:,1,1,n), loc, &
                                  sp(:,1,1,:), los, &
                                  lo, hi, face, dim)
            case (2)
             call ml_interface_2d(rp(:,:,1,n), lor, &
                                  fp(:,:,1,n), lof, &
                                  cp(:,:,1,n), loc, &
                                  sp(:,:,1,:), los, &
                                  lo, hi, face, dim)
            case (3)
             call ml_interface_3d(rp(:,:,:,n), lor, &
                                  fp(:,:,:,n), lof, &
                                  cp(:,:,:,n), loc, &
                                  sp(:,:,:,:), los, &
                                  lo, hi, face, dim)
           end select
          end do

        end if
       end if

    end do
  end do
 end subroutine ml_interface 
subroutine ml_interface_c(res, cr, flux, cf, crse, ss, crse_domain, face, dim)
  type(multifab), intent(inout) :: res
  type(multifab), intent(in   ) :: flux
  type(multifab), intent(in   ) :: crse
  type(multifab), intent(in   ) :: ss
  integer, intent(in) :: cr, cf
  type(box), intent(in) :: crse_domain
  type(box) :: rbox, fbox, cbox, sbox
  integer :: lo (res%dim), hi (res%dim)
  integer :: loc(res%dim)
  integer :: lof(res%dim)
  integer :: lor(res%dim)
  integer :: los(res%dim)
  integer :: lo_dom(res%dim), hi_dom(res%dim)
  integer :: dm
  integer :: face, dim
  integer :: dir
  integer :: i, j, n

  real(kind=dp_t), pointer :: rp(:,:,:,:)
  real(kind=dp_t), pointer :: fp(:,:,:,:)
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  real(kind=dp_t), pointer :: sp(:,:,:,:)

  dm = res%dim
  dir = dim

  lo_dom = lwb(crse_domain)
  hi_dom = upb(crse_domain)

  do j = 1, crse%nboxes

    cbox = get_ibox(crse,j)
    loc = lwb(cbox) - crse%ng

    rbox = get_ibox(res,j)
    lor = lwb(rbox) - res%ng

    sbox = get_ibox(ss,j)
    los = lwb(sbox) - ss%ng

    do i = 1, flux%nboxes

       fbox   = get_ibox(flux,i)
       lof = lwb(fbox)

       if (box_contains(crse_domain,fbox)) then
        if (.not. empty(box_intersection(cbox,fbox))) then

          lo = lwb(box_intersection(cbox,fbox))
          hi = upb(box_intersection(cbox,fbox))
          fp => dataptr(flux, i, cf)
          cp => dataptr(crse, j, cr)
          rp => dataptr(res, j, cr)
          sp => dataptr(ss, j)
           select case (dm)
           case (1)
             call ml_interface_1d(rp(:,1,1,1), lor, &
                                  fp(:,1,1,1), lof, &
                                  cp(:,1,1,1), loc, &
                                  sp(:,1,1,:), los, &
                                  lo, hi, face, dim)
            case (2)
             call ml_interface_2d(rp(:,:,1,1), lor, &
                                  fp(:,:,1,1), lof, &
                                  cp(:,:,1,1), loc, &
                                  sp(:,:,1,:), los, &
                                  lo, hi, face, dim)
            case (3)
             call ml_interface_3d(rp(:,:,:,1), lor, &
                                  fp(:,:,:,1), lof, &
                                  cp(:,:,:,1), loc, &
                                  sp(:,:,:,:), los, &
                                  lo, hi, face, dim)
           end select

        end if
       end if

    end do
  end do
 end subroutine ml_interface_c

  subroutine ml_interface_1d(res, lor, fine_flux, lof, cc, loc, &
                             ss , los, lo, hi, face, dim)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: res(lor(1):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):)
    real (kind = dp_t), intent(in   ) :: ss(los(1):,0:)
    integer, intent(in) :: face, dim

    integer i
    real (kind = dp_t) :: crse_flux

    i = lo(1)

!   Lo side
    if (face == -1) then
      crse_flux = ss(i,1)*(cc(i)-cc(i+1))
      res(i) = res(i) + fine_flux(i) - crse_flux

!   Hi side
    else if (face == 1) then
      crse_flux = ss(i,2)*(cc(i)-cc(i-1))
      res(i) = res(i) + fine_flux(i) - crse_flux
    end if

  end subroutine ml_interface_1d

  subroutine ml_interface_2d(res, lor, fine_flux, lof, cc, loc, &
                             ss , los, lo, hi, face, dim)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):)
    real (kind = dp_t), intent(in   ) ::        cc(loc(1):,loc(2):)
    real (kind = dp_t), intent(in   ) ::        ss(los(1):,los(2):,0:)
    integer, intent(in) :: face, dim

    integer i, j
    real (kind = dp_t) :: crse_flux

    !   Hi i side
    if ( dim == 1 ) then
       if (face == 1) then
          i = lo(1)
          do j = lo(2),hi(2)
             crse_flux = ss(i,j,2)*(cc(i,j)-cc(i-1,j))
             res(i,j) = res(i,j) + fine_flux(i,j) - crse_flux
          end do
          !   Lo i side
       else if (face == -1) then
          i = lo(1)
          do j = lo(2),hi(2)
             crse_flux = ss(i,j,1)*(cc(i,j)-cc(i+1,j))
             res(i,j) = res(i,j) + fine_flux(i,j) - crse_flux
          end do
       end if
    else if ( dim == 2 ) then
       !   Hi j side
       if (face == 1) then
          j = lo(2)
          do i = lo(1),hi(1)
             crse_flux = ss(i,j,4)*(cc(i,j)-cc(i,j-1))
             res(i,j) = res(i,j) + fine_flux(i,j) - crse_flux
          end do
          !   Lo j side
       else if (face == -2) then
          j = lo(2)
          do i = lo(1),hi(1)
             crse_flux = ss(i,j,3)*(cc(i,j)-cc(i,j+1))
             res(i,j) = res(i,j) + fine_flux(i,j) - crse_flux
          end do
       end if
    end if
  end subroutine ml_interface_2d

  subroutine ml_interface_3d(res, lor, fine_flux, lof, cc, loc, &
                             ss , los, lo, hi, face, dim)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):,lor(3):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):,lof(3):)
    real (kind = dp_t), intent(in   ) ::        cc(loc(1):,loc(2):,loc(3):)
    real (kind = dp_t), intent(in   ) ::        ss(los(1):,los(2):,los(3):,0:)
    integer, intent(in) :: face, dim

    integer i, j, k
    real (kind = dp_t) :: crse_flux

!   Hi i side
    if ( dim == 1 ) then
       if (face == 1) then
      i = lo(1)
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
        crse_flux = ss(i,j,k,2)*(cc(i,j,k)-cc(i-1,j,k))
        res(i,j,k) = res(i,j,k) + fine_flux(i,j,k) - crse_flux
      end do
      end do
!   Lo i side
    else if (face == -1) then
      i = lo(1)
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
        crse_flux = ss(i,j,k,1)*(cc(i,j,k)-cc(i+1,j,k))
        res(i,j,k) = res(i,j,k) + fine_flux(i,j,k) - crse_flux
      end do
      end do
   end if
!   Hi j side
   else if ( dim ==  2 )  then
    if (face == 1) then
      j = lo(2)
      do k = lo(3),hi(3)
      do i = lo(1),hi(1)
        crse_flux = ss(i,j,k,4)*(cc(i,j,k)-cc(i,j-1,k))
        res(i,j,k) = res(i,j,k) + fine_flux(i,j,k) - crse_flux
      end do
      end do
!   Lo j side
    else if (face == -1) then
      j = lo(2)
      do k = lo(3),hi(3)
      do i = lo(1),hi(1)
        crse_flux = ss(i,j,k,3)*(cc(i,j,k)-cc(i,j+1,k))
        res(i,j,k) = res(i,j,k) + fine_flux(i,j,k) - crse_flux
      end do
      end do
   end if
      else if ( dim == 3 ) then
!   Hi k side
      if (face == 1) then
      k = lo(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        crse_flux = ss(i,j,k,6)*(cc(i,j,k)-cc(i,j,k-1))
        res(i,j,k) = res(i,j,k) + fine_flux(i,j,k) - crse_flux
      end do
      end do
!   Lo k side
    else if (face == -1) then
      k = lo(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        crse_flux = ss(i,j,k,5)*(cc(i,j,k)-cc(i,j,k+1))
        res(i,j,k) = res(i,j,k) + fine_flux(i,j,k) - crse_flux
      end do
      end do
    end if
    end if
  end subroutine ml_interface_3d

 subroutine ml_crse_contrib(res, flux, crse, ss, mm, crse_domain, ir, side)
  type(multifab), intent(inout) :: res
  type(multifab), intent(inout) :: flux
  type(multifab), intent(in   ) :: crse
  type(multifab), intent(in   ) :: ss
  type(imultifab),intent(in   ) :: mm
  integer, intent(in) :: ir(:)
  type(box), intent(in) :: crse_domain
  type(box) :: rbox, fbox, cbox, sbox, mbox
  integer :: lo (res%dim), hi (res%dim)
  integer :: loc(res%dim)
  integer :: lof(res%dim), hif(res%dim)
  integer :: lor(res%dim)
  integer :: los(res%dim)
  integer :: lom(res%dim)
  integer :: lo_dom(res%dim), hi_dom(res%dim)
  integer :: side
  integer :: dir
  integer :: i, j, n

  integer :: dm
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  real(kind=dp_t), pointer :: fp(:,:,:,:)
  real(kind=dp_t), pointer :: cp(:,:,:,:)
  real(kind=dp_t), pointer :: sp(:,:,:,:)
  integer,         pointer :: mp(:,:,:,:)

  dm  = res%dim
  dir = iabs(side)

  lo_dom = lwb(crse_domain)
  hi_dom = upb(crse_domain)+1

  do j = 1, crse%nboxes

    cbox = get_ibox(crse,j)
    loc = lwb(cbox) - crse%ng

    rbox = get_ibox(res,j)
    lor = lwb(rbox) - res%ng

    sbox = get_ibox(ss,j)
    los = lwb(sbox) - ss%ng

    do i = 1, flux%nboxes

      fbox   = get_ibox(flux,i)
      lof = lwb(fbox)
      hif = upb(fbox)

      mbox = get_ibox(mm,i)
      lom = lwb(mbox) - mm%ng

      if ((.not. (lof(dir) == lo_dom(dir) .or. lof(dir) == hi_dom(dir))) .and. & 
          box_intersects(cbox,fbox)) then
        lo(:) = lwb(box_intersection(cbox,fbox))
        hi(:) = upb(box_intersection(cbox,fbox))

        fp => dataptr(flux, i)
        mp => dataptr(mm  , i)

        cp => dataptr(crse, j)
        rp => dataptr(res , j)
        sp => dataptr(ss  , j)
        do n = 1, 1
           select case (dm)
           case (1)
              call ml_interface_1d_nodal(rp(:,1,1,n), lor, &
                                         fp(:,1,1,n), lof, hif, &
                                         cp(:,1,1,n), loc, &
                                         sp(:,1,1,:), los, &
                                         lo, hi, ir, side)
          case (2)
              call ml_interface_2d_nodal(rp(:,:,1,n), lor, &
                                         fp(:,:,1,n), lof, hif, &
                                         cp(:,:,1,n), loc, &
                                         sp(:,:,1,:), los, &
                                         mp(:,:,1,1), lom, &
                                         lo, hi, ir, side)
          case (3)
              call ml_interface_3d_nodal(rp(:,:,:,n), lor, &
                                         fp(:,:,:,n), lof, hif, &
                                         cp(:,:,:,n), loc, &
                                         sp(:,:,:,:), los, &
                                         mp(:,:,:,1), lom, &
                                         lo, hi, ir, side)
          end select
        end do
      end if
    end do
  end do
 end subroutine ml_crse_contrib

  subroutine ml_interface_1d_nodal(res, lor, fine_flux, lof, hif, cc, loc, &
                                   ss , los, lo, hi, ir, side)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lof(:), hif(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: res(lor(1):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):)
    real (kind = dp_t), intent(in   ) :: ss(los(1):,0:)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)

    integer i
    real (kind = dp_t) :: crse_flux

    i = lo(1)

!   Lo side
    if (side == -1) then
      crse_flux = ss(i,1)*(cc(i+1)-cc(i))
      res(i) = res(i) - fine_flux(i) + crse_flux

!   Hi side
    else if (side == 1) then
      crse_flux = ss(i,2)*(cc(i-1)-cc(i))
      res(i) = res(i) - fine_flux(i) + crse_flux
    end if

  end subroutine ml_interface_1d_nodal

  subroutine ml_interface_2d_nodal(res, lor, fine_flux, lof, hif, cc, loc, &
                             ss , los, mm, lom, lo, hi, ir, side)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lom(:)
    integer, intent(in) :: lof(:), hif(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):)
    real (kind = dp_t), intent(in   ) ::        cc(loc(1):,loc(2):)
    real (kind = dp_t), intent(in   ) ::        ss(los(1):,los(2):,0:)
    integer           , intent(in   ) ::        mm(lom(1):,lom(2):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)

    integer i, j
    real (kind = dp_t) :: crse_flux

    i = lo(1)
    j = lo(2)

!   NOTE: THESE STENCILS ONLY WORK FOR DX == DY.

!   NOTE: MM IS ON THE FINE GRID, NOT THE CRSE

!   Lo i side
    if (side == -1) then

      do j = lo(2),hi(2)

        if (j == lof(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),2,-1)) then
          crse_flux = (ss(i,j,8)*(cc(i+1,j+1) + HALF*cc(i+1,j) + &
                                 HALF * cc(i,j+1) - TWO*cc(i,j))) * HALF
        else if (j == hif(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),2,+1)) then
          crse_flux = (ss(i,j,3)*(cc(i+1,j-1) + HALF*cc(i+1,j) + &
                                 HALF * cc(i,j-1) - TWO*cc(i,j))) * HALF
        else
          crse_flux = ss(i,j,8)*(cc(i+1,j+1) + HALF*cc(i+1,j) + &
                                 HALF * cc(i,j+1) - TWO*cc(i,j)) &
                     +ss(i,j,3)*(cc(i+1,j-1) + HALF*cc(i+1,j) + &
                                 HALF * cc(i,j-1) - TWO*cc(i,j))
        end if

        res(i,j) = res(i,j) + crse_flux

        if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then
           res(i,j) = res(i,j) + fine_flux(i,j)
        end if

      end do

!   Hi i side
    else if (side ==  1) then
      do j = lo(2),hi(2)

        if (j == lof(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),2,-1)) then
          crse_flux = (ss(i,j,6)*(cc(i-1,j+1) + HALF*cc(i-1,j) + &
                                 HALF * cc(i,j+1) - TWO*cc(i,j))) * HALF
        else if (j == hif(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),2,+1)) then
          crse_flux = (ss(i,j,1)*(cc(i-1,j-1) + HALF*cc(i-1,j) + &
                                 HALF * cc(i,j-1) - TWO*cc(i,j))) * HALF
        else
          crse_flux = ss(i,j,6)*(cc(i-1,j+1) + HALF*cc(i-1,j) + &
                                 HALF * cc(i,j+1) - TWO*cc(i,j)) &
                     +ss(i,j,1)*(cc(i-1,j-1) + HALF*cc(i-1,j) + &
                                 HALF * cc(i,j-1) - TWO*cc(i,j))
        end if

        res(i,j) = res(i,j) + crse_flux

        if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) &
           res(i,j) = res(i,j) + fine_flux(i,j)

      end do

!   Lo j side
    else if (side == -2) then

      do i = lo(1),hi(1)

        if (i == lof(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),1,-1)) then
          crse_flux = (ss(i,j,8)*(cc(i+1,j+1) + HALF*cc(i+1,j) + &
                                 HALF * cc(i,j+1) - TWO*cc(i,j))) * HALF
        else if (i == hif(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),1,+1)) then
          crse_flux = (ss(i,j,6)*(cc(i-1,j+1) + HALF*cc(i-1,j) + &
                                 HALF * cc(i,j+1) - TWO*cc(i,j))) * HALF
        else
          crse_flux = ss(i,j,8)*(cc(i+1,j+1) + HALF*cc(i+1,j) + &
                                 HALF * cc(i,j+1) - TWO*cc(i,j)) &
                     +ss(i,j,6)*(cc(i-1,j+1) + HALF*cc(i-1,j) + &
                                 HALF * cc(i,j+1) - TWO*cc(i,j))
        end if

        res(i,j) = res(i,j) + crse_flux

        if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) &
           res(i,j) = res(i,j) + fine_flux(i,j)

      end do

!   Hi j side
    else if (side ==  2) then

      do i = lo(1),hi(1)

        if (i == lof(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),1,-1)) then
          crse_flux = (ss(i,j,3)*(cc(i+1,j-1) + HALF*cc(i+1,j) + &
                                 HALF * cc(i,j-1) - TWO*cc(i,j)) ) * HALF
        else if (i == hif(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),1,+1)) then
          crse_flux = (ss(i,j,1)*(cc(i-1,j-1) + HALF*cc(i-1,j) + &
                                 HALF * cc(i,j-1) - TWO*cc(i,j)) ) * HALF
        else
          crse_flux = ss(i,j,3)*(cc(i+1,j-1) + HALF*cc(i+1,j) + &
                                 HALF * cc(i,j-1) - TWO*cc(i,j)) &
                     +ss(i,j,1)*(cc(i-1,j-1) + HALF*cc(i-1,j) + &
                                 HALF * cc(i,j-1) - TWO*cc(i,j))
        end if

        res(i,j) = res(i,j) + crse_flux

        if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) &
           res(i,j) = res(i,j) + fine_flux(i,j)

      end do
    end if

  end subroutine ml_interface_2d_nodal

  subroutine ml_interface_3d_nodal(res, lor, fine_flux, lof, hif, cc, loc, &
                             ss, los, mm, lom, lo, hi, ir, side)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lof(:), hif(:)
    integer, intent(in) :: lom(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):,lor(3):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):,lof(3):)
    real (kind = dp_t), intent(in   ) ::        cc(loc(1):,loc(2):,loc(3):)
    real (kind = dp_t), intent(in   ) ::        ss(los(1):,los(2):,los(3):,0:)
    integer           , intent(in   ) ::        mm(lom(1):,lom(2):,lom(3):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)

    integer i, j, k
    integer ioff,joff,koff
    integer sig_mm,sig_mp,sig_pm,sig_pp
    real (kind = dp_t) :: crse_flux
    real (kind = dp_t) :: cell_mm,cell_mp,cell_pm,cell_pp

    i = lo(1)
    j = lo(2)
    k = lo(3)

    crse_flux = ZERO

!   NOTE: THESE STENCILS ONLY WORK FOR DX == DY.

!   NOTE: MM IS ON THE FINE GRID, NOT THE CRSE

!   Lo/Hi i side
    if (side == -1 .or. side == 1) then
      i = lo(1)
      if (side == -1) then
         ioff = i+1
         sig_mm =  3
         sig_pm =  8
         sig_mp = 15
         sig_pp = 20
      else
         ioff = i-1
         sig_mm =  1
         sig_pm =  6
         sig_mp = 13
         sig_pp = 18
      end if

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)

        cell_mm = ss(i,j,k,sig_mm)*(cc(ioff,j-1,k-1) + cc(ioff,j-1,k  ) &
                                   +cc(ioff,j  ,k-1) + cc(i  ,j-1,k-1) &
                             - FOUR*cc(i  ,j  ,k) )
        cell_pm = ss(i,j,k,sig_pm)*(cc(ioff,j+1,k-1) + cc(ioff,j+1,k  ) &
                                   +cc(ioff,j  ,k-1) + cc(i  ,j+1,k-1) &
                             - FOUR*cc(i  ,j  ,k) )
        cell_mp = ss(i,j,k,sig_mp)*(cc(ioff,j-1,k+1) + cc(ioff,j-1,k  ) &
                                   +cc(ioff,j  ,k+1) + cc(i  ,j-1,k+1) &
                             - FOUR*cc(i  ,j  ,k) )
        cell_pp = ss(i,j,k,sig_pp)*(cc(ioff,j+1,k+1) + cc(ioff,j+1,k  ) &
                                   +cc(ioff,j  ,k+1) + cc(i  ,j+1,k+1) &
                             - FOUR*cc(i  ,j  ,k) )

        crse_flux = zero

        if (k == lof(3) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),3,-1)) then
           if (j == lof(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),2,-1)) then
             crse_flux = THIRD*cell_pp 
           else if (j == hif(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),2,+1)) then
             crse_flux = THIRD*cell_mp
           else
             crse_flux = HALF*(cell_pp + cell_mp)
            end if
        else if (k == hif(3) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),3,+1)) then
           if (j == lof(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),2,-1)) then
             crse_flux = THIRD*cell_pm 
           else if (j == hif(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),2,+1)) then
             crse_flux = THIRD*cell_mm 
           else
             crse_flux = HALF*(cell_pm  + cell_mm)
           end if
        else 
           if (j == lof(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),2,-1)) then
             crse_flux = HALF*(cell_pm  + cell_pp)
           else if (j == hif(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),2,+1)) then
             crse_flux = HALF*(cell_mm  + cell_mp)
           else
             crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
           end if
        end if

        res(i,j,k) = res(i,j,k) + crse_flux
  
        if (bc_dirichlet(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,0)) then
           res(i,j,k) = res(i,j,k) + fine_flux(i,j,k)
        end if
      end do
      end do
!   Lo/Hi j side
    else if (side == -2 .or. side == 2) then
      j = lo(2)
      if (side == -2) then
         joff = j+1
         sig_mm =  6
         sig_pm =  8
         sig_mp = 18
         sig_pp = 20
      else
         joff = j-1
         sig_mm =  1
         sig_pm =  3
         sig_mp = 13
         sig_pp = 15
      end if
      do k = lo(3),hi(3)
      do i = lo(1),hi(1)

        cell_mm = ss(i,j,k,sig_mm)*(cc(i-1,joff,k-1) + cc(i-1,joff,k  ) &
                                   +cc(i  ,joff,k-1) + cc(i-1,j   ,k-1) &
                             - FOUR*cc(i  ,j  ,k) )
        cell_pm = ss(i,j,k,sig_pm)*(cc(i+1,joff,k-1) + cc(i+1,joff,k  ) &
                                   +cc(i  ,joff,k-1) + cc(i+1,j   ,k-1) &
                             - FOUR*cc(i  ,j  ,k) )
        cell_mp = ss(i,j,k,sig_mp)*(cc(i-1,joff,k+1) + cc(i-1,joff,k  ) &
                                   +cc(i  ,joff,k+1) + cc(i-1,j   ,k+1) &
                             - FOUR*cc(i  ,j  ,k) )
        cell_pp = ss(i,j,k,sig_pp)*(cc(i+1,joff,k+1) + cc(i+1,joff,k  ) &
                                   +cc(i  ,joff,k+1) + cc(i+1,j   ,k+1) &
                             - FOUR*cc(i  ,j  ,k) )

        if (k == lof(3) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),3,-1)) then
           if (i == lof(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,-1)) then
             crse_flux = THIRD*cell_pp 
           else if (i == hif(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,+1)) then
             crse_flux = THIRD*cell_mp
           else
             crse_flux = HALF*(cell_pp + cell_mp)
            end if
        else if (k == hif(3) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),3,+1)) then
           if (i == lof(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,-1)) then
             crse_flux = THIRD*cell_pm 
           else if (i == hif(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,+1)) then
             crse_flux = THIRD*cell_mm 
           else
             crse_flux = HALF*(cell_pm  + cell_mm)
           end if
        else 
           if (i == lof(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,-1)) then
             crse_flux = HALF*(cell_pm  + cell_pp)
           else if (i == hif(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,+1)) then
             crse_flux = HALF*(cell_mm  + cell_mp)
           else
             crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
           end if
        end if

        res(i,j,k) = res(i,j,k) + crse_flux
  
        if (bc_dirichlet(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,0)) then
           res(i,j,k) = res(i,j,k) + fine_flux(i,j,k)
        end if

      end do
      end do
!   Lo/Hi k side
    else if (side == -3 .or. side == 3) then
      k = lo(3)
      if (side == -3) then
         koff = k+1
         sig_mm = 13
         sig_pm = 15
         sig_mp = 18
         sig_pp = 20
      else
         koff = k-1
         sig_mm =  1
         sig_pm =  3
         sig_mp =  6
         sig_pp =  8
      end if

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

        cell_mm = ss(i,j,k,sig_mm)*(cc(i-1,j-1,koff) + cc(i-1,j  ,koff) &
                                   +cc(i  ,j-1,koff) + cc(i-1,j-1,k   ) &
                             - FOUR*cc(i  ,j  ,k) )
        cell_pm = ss(i,j,k,sig_pm)*(cc(i+1,j-1,koff) + cc(i+1,j  ,koff) &
                                   +cc(i  ,j-1,koff) + cc(i+1,j-1,k   ) &
                             - FOUR*cc(i  ,j  ,k) )
        cell_mp = ss(i,j,k,sig_mp)*(cc(i-1,j+1,koff) + cc(i-1,j  ,koff) &
                                   +cc(i  ,j+1,koff) + cc(i-1,j+1,k   ) &
                             - FOUR*cc(i  ,j  ,k) )
        cell_pp = ss(i,j,k,sig_pp)*(cc(i+1,j+1,koff) + cc(i+1,j  ,koff) &
                                   +cc(i  ,j+1,koff) + cc(i+1,j+1,k   ) &
                             - FOUR*cc(i  ,j  ,k) )

        if (j == lof(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),2,-1)) then
           if (i == lof(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,-1)) then
             crse_flux = THIRD*cell_pp 
           else if (i == hif(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,+1)) then
             crse_flux = THIRD*cell_mp
           else
             crse_flux = HALF*(cell_pp + cell_mp)
            end if
        else if (j == hif(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),2,+1)) then
           if (i == lof(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,-1)) then
             crse_flux = THIRD*cell_pm 
           else if (i == hif(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,+1)) then
             crse_flux = THIRD*cell_mm 
           else
             crse_flux = HALF*(cell_pm  + cell_mm)
           end if
        else 
           if (i == lof(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,-1)) then
             crse_flux = HALF*(cell_pm  + cell_pp)
           else if (i == hif(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,+1)) then
             crse_flux = HALF*(cell_mm  + cell_mp)
           else
             crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
           end if
        end if

        res(i,j,k) = res(i,j,k) + crse_flux

        if (bc_dirichlet(mm(ir(1)*i,ir(2)*j,ir(3)*k),1,0)) then
           res(i,j,k) = res(i,j,k) + fine_flux(i,j,k)
        end if
      end do
      end do
    end if

  end subroutine ml_interface_3d_nodal

end module ml_interface_stencil_module
