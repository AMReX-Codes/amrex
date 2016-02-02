module cc_interface_stencil_module

  use bl_types
  use layout_module
  use multifab_module
  use bc_functions_module
  use bl_constants_module

  implicit none

  private

  public :: ml_interface, ml_interface_c

contains

  subroutine ml_interface(res, flux, crse, ss, crse_domain, facemap, indxmap, efactor)
    use bl_prof_module
    type(multifab), intent(inout) :: res
    type(multifab), intent(in   ) :: flux
    type(multifab), intent(in   ) :: crse
    type(multifab), intent(in   ) :: ss
    type(box),      intent(in   ) :: crse_domain
    integer,        intent(in   ) :: facemap(:), indxmap(:)
    real(kind=dp_t),intent(in   ) :: efactor
    type(bl_prof_timer), save :: bpt
    call build(bpt, "ml_interf")
    call ml_interface_c(res, 1, flux, 1, crse, ss, crse_domain, facemap, indxmap, efactor)
    call destroy(bpt)
  end subroutine ml_interface

  subroutine ml_interface_c(res, cr, flux, cf, crse, ss, crse_domain, facemap, indxmap, efactor)
    use bl_prof_module
    use stencil_util_module, only : is_ibc_stencil

    type(multifab), intent(inout) :: res
    type(multifab), intent(in   ) :: flux
    type(multifab), intent(in   ) :: crse
    type(multifab), intent(in   ) :: ss
    integer,        intent(in   ) :: cr, cf
    type(box),      intent(in   ) :: crse_domain
    integer,        intent(in   ) :: facemap(:), indxmap(:)
    real(kind=dp_t),intent(in   ) :: efactor

    type(box)      :: fbox
    integer        :: lor(get_dim(res)), los(get_dim(res)), i, j, dm, face, dim, idim
    integer        :: lo (get_dim(res)), hi (get_dim(res)), loc(get_dim(res))
    logical        :: pmask(get_dim(res))
    real(kind=dp_t), pointer :: rp(:,:,:,:), fp(:,:,:,:), cp(:,:,:,:), sp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_interf_c")

    call bl_assert(.not.(any(flux%nodal) .or. any(crse%nodal)), 'ml_interface_c(): cannot have nodal')

    pmask = get_pmask(get_layout(res))
    dm    = get_dim(res)

    !$omp parallel private(fbox,lor,los,i,j,face,dim,idim,lo,hi,loc,rp,fp,cp,sp)
    do idim = 1, dm
       !$omp do schedule(dynamic)
       do i = 1, nfabs(flux)

          dim = abs(facemap(i))

          if (dim .eq. idim) then
             j = indxmap(i)
             face = sign(1, facemap(i))
             
             fbox = get_ibox(flux,i)
             if (pmask(dim) .and. (.not. contains(crse_domain,fbox)) ) then
                if ( face .eq. -1 ) then
                   fbox = shift(fbox,  extent(crse_domain,dim), dim)
                else
                   fbox = shift(fbox, -extent(crse_domain,dim), dim)
                end if
             end if

             fp   => dataptr(flux, i, cf)
             rp   => dataptr(res , j, cr)
             cp   => dataptr(crse, j, cr)
             sp   => dataptr(ss  , j)
             
             lo = lwb(fbox)
             hi = upb(fbox)             
             loc = lwb(get_pbox(crse,j))
             lor = lwb(get_pbox(res,j))

             if (is_ibc_stencil(ss,j)) then
                select case (dm)
                case (2)
                   call ml_interface_ibc_2d(rp(:,:,1,1), lor, cp(:,:,1,1), loc, &
                        fp(:,:,1,1), lo, lo, hi, face, dim, efactor, sp(dim+1,1,1,1))
                case (3)
                   call ml_interface_ibc_3d(rp(:,:,:,1), lor, cp(:,:,:,1), loc, &
                        fp(:,:,:,1), lo, lo, hi, face, dim, efactor, sp(dim+1,1,1,1))
                end select
             else
                los  =  lwb(get_pbox(ss,j))
                select case (dm)
                case (1)
                   call ml_interface_1d(rp(:,1,1,1), lor, cp(:,1,1,1), loc, &
                        sp(:,:,1,1), los, fp(:,1,1,1), lo, lo, face, efactor)
                case (2)
                   call ml_interface_2d(rp(:,:,1,1), lor, cp(:,:,1,1), loc, &
                        sp(:,:,:,1), los, fp(:,:,1,1), lo, lo, hi, face, dim, efactor)
                case (3)
                   call ml_interface_3d(rp(:,:,:,1), lor, cp(:,:,:,1), loc, &
                        sp(:,:,:,:), los, fp(:,:,:,1), lo, lo, hi, face, dim, efactor)
                end select
             end if
          end if
       end do
       !$omp end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine ml_interface_c


  subroutine ml_interface_1d(res, lor, cc, loc, ss , los, fine_flux, lof, lo, face, efactor)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lof(:) 
    integer, intent(in) :: lo(:)
    real (kind = dp_t), intent(inout) :: res(lor(1):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):)
    real (kind = dp_t), intent(in   ) :: ss(0:,los(1):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):)
    integer, intent(in) :: face
    real(kind=dp_t), intent(in) :: efactor

    integer :: i
    real (kind = dp_t) :: crse_flux

    i = lo(1)

    !   Lo side
    if (face == -1) then
       crse_flux = ss(1,i)*(cc(i)-cc(i+1))
       res(i) = res(i) + efactor*(fine_flux(i)-crse_flux)

       !   Hi side
    else if (face == 1) then
       crse_flux = ss(2,i)*(cc(i)-cc(i-1))
       res(i) = res(i) + efactor*(fine_flux(i)-crse_flux)
    end if

  end subroutine ml_interface_1d

  subroutine ml_interface_2d(res, lor, cc, loc, ss , los, fine_flux, lof, &
       lo, hi, face, dim, efactor)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: res(lor(1):,lor(2):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):,loc(2):)
    real (kind = dp_t), intent(in   ) :: ss(0:,los(1):,los(2):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):)
    integer, intent(in) :: face, dim
    real(kind=dp_t), intent(in) :: efactor

    integer :: i, j, ns
    real (kind = dp_t) :: crse_flux
    real (kind = dp_t), parameter :: FIFTEEN16TH = 15.0_dp_t/16.0_dp_t
 
    ns = size(ss,1)
!    ns = size(ss,dim=1)  ! some version of Intel does not like this

    !   Hi i side
    if ( dim == 1 ) then
       if (face == 1) then
          i = lo(1)
          do j = lo(2),hi(2)
             if (ns.eq.7) then
                crse_flux = ss(2,i,j)*(cc(i,j)-cc(i-1,j))
             else if (ns.eq.9) then
                crse_flux = &
                   FIFTEEN16TH*ss(2,i,j)*(cc(i-1,j)-cc(i  ,j)) &
                 +             ss(1,i,j)*(cc(i-2,j)-cc(i+1,j))

             endif
             res(i,j) = res(i,j) + efactor*(fine_flux(i,j)-crse_flux)
          end do
          !   Lo i side
       else if (face == -1) then
          i = lo(1)
          do j = lo(2),hi(2)
             if (ns.eq.7) then
                crse_flux = ss(1,i,j)*(cc(i,j)-cc(i+1,j))
             else if (ns.eq.9) then
                crse_flux = &
                   FIFTEEN16TH*ss(3,i,j)*(cc(i+1,j)-cc(i  ,j)) &
                 +             ss(4,i,j)*(cc(i+2,j)-cc(i-1,j))

             endif
             res(i,j) = res(i,j) + efactor*(fine_flux(i,j)-crse_flux)
          end do
       end if
    else if ( dim == 2 ) then
       !   Hi j side
       if (face == 1) then
          j = lo(2)
          do i = lo(1),hi(1)
             if (ns.eq.7) then
                crse_flux = ss(4,i,j)*(cc(i,j)-cc(i,j-1))
             else if (ns.eq.9) then
                crse_flux = &
                   FIFTEEN16TH*ss(6,i,j)*(cc(i,j-1)-cc(i,j  )) &
                 +             ss(5,i,j)*(cc(i,j-2)-cc(i,j+1))

             endif
             res(i,j) = res(i,j) + efactor*(fine_flux(i,j)-crse_flux)
          end do
          !   Lo j side
       else if (face == -1) then
          j = lo(2)
          do i = lo(1),hi(1)
             if (ns.eq.7) then
                crse_flux = ss(3,i,j)*(cc(i,j)-cc(i,j+1))
             else if (ns.eq.9) then
                crse_flux = &
                   FIFTEEN16TH*ss(7,i,j)*(cc(i,j+1)-cc(i,j  )) &
                 +             ss(8,i,j)*(cc(i,j+2)-cc(i,j-1))
             endif
             res(i,j) = res(i,j) + efactor*(fine_flux(i,j)-crse_flux)
          end do
       end if
    end if
  end subroutine ml_interface_2d

  subroutine ml_interface_3d(res, lor, cc, loc, ss , los, fine_flux, lof, &
       lo, hi, face, dim, efactor)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):,lor(3):)
    real (kind = dp_t), intent(in   ) ::        cc(loc(1):,loc(2):,loc(3):)
    real (kind = dp_t), intent(in   ) ::        ss(0:,los(1):,los(2):,los(3):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):,lof(3):)
    integer, intent(in) :: face, dim
    real(kind=dp_t), intent(in) :: efactor

    integer :: i, j, k

    !   Hi i side
    if ( dim == 1 ) then
       if (face == 1) then
          i = lo(1)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                res(i,j,k) = res(i,j,k) + efactor*( fine_flux(i,j,k) &
                     - ss(2,i,j,k)*(cc(i,j,k)-cc(i-1,j,k)) )
             end do
          end do
          !   Lo i side
       else if (face == -1) then
          i = lo(1)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                res(i,j,k) = res(i,j,k) + efactor*( fine_flux(i,j,k) &
                     - ss(1,i,j,k)*(cc(i,j,k)-cc(i+1,j,k)) )
             end do
          end do
       end if
       !   Hi j side
    else if ( dim ==  2 )  then
       if (face == 1) then
          j = lo(2)
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)
                res(i,j,k) = res(i,j,k) + efactor*( fine_flux(i,j,k) &
                     - ss(4,i,j,k)*(cc(i,j,k)-cc(i,j-1,k)) )
             end do
          end do
          !   Lo j side
       else if (face == -1) then
          j = lo(2)
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)
                res(i,j,k) = res(i,j,k) + efactor*( fine_flux(i,j,k) &
                     - ss(3,i,j,k)*(cc(i,j,k)-cc(i,j+1,k)) )
             end do
          end do
       end if
    else if ( dim == 3 ) then
       !   Hi k side
       if (face == 1) then
          k = lo(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                res(i,j,k) = res(i,j,k) + efactor*( fine_flux(i,j,k) &
                     - ss(6,i,j,k)*(cc(i,j,k)-cc(i,j,k-1)) )
             end do
          end do
          !   Lo k side
       else if (face == -1) then
          k = lo(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                res(i,j,k) = res(i,j,k) + efactor*( fine_flux(i,j,k) &
                     - ss(5,i,j,k)*(cc(i,j,k)-cc(i,j,k+1)) )
             end do
          end do
       end if
    end if

  end subroutine ml_interface_3d

  subroutine ml_interface_ibc_2d(res, lor, cc, loc, fine_flux, lof, lo, hi, face, dim, efactor, ss)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: res(lor(1):,lor(2):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):,loc(2):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):)
    integer, intent(in) :: face, dim
    real(kind=dp_t), intent(in) :: efactor
    real (kind = dp_t), intent(in) :: ss

    integer :: i, j

    if (dim .eq. 1) then

       i = lo(1)
       do j = lo(2),hi(2)
          res(i,j) = res(i,j) + efactor*( fine_flux(i,j) &
               - ss*(cc(i,j)-cc(i-face,j)) )
       end do

    else 

       j = lo(2)
       do i = lo(1),hi(1)
          res(i,j) = res(i,j) + efactor*( fine_flux(i,j) &
               - ss*(cc(i,j)-cc(i,j-face)) )
       end do

    end if
  end subroutine ml_interface_ibc_2d

  subroutine ml_interface_ibc_3d(res, lor, cc, loc, fine_flux, lof, lo, hi, face, dim, efactor, ss)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: res(lor(1):,lor(2):,lor(3):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):,loc(2):,loc(3):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):,lof(3):)
    integer, intent(in) :: face, dim
    real(kind=dp_t), intent(in) :: efactor
    real (kind = dp_t), intent(in) :: ss

    integer :: i, j, k

    if (dim .eq. 1) then

       i = lo(1)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             res(i,j,k) = res(i,j,k) + efactor*( fine_flux(i,j,k) &
                  - ss*(cc(i,j,k)-cc(i-face,j,k)) )
          end do
       end do

    else if (dim .eq. 2) then

       j = lo(2)
       do k = lo(3),hi(3)
          do i = lo(1),hi(1)
             res(i,j,k) = res(i,j,k) + efactor*( fine_flux(i,j,k) &
                  - ss*(cc(i,j,k)-cc(i,j-face,k)) )
          end do
       end do

    else 

       k = lo(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             res(i,j,k) = res(i,j,k) + efactor*( fine_flux(i,j,k) &
                  - ss*(cc(i,j,k)-cc(i,j,k-face)) )
          end do
       end do

    end if
  end subroutine ml_interface_ibc_3d

end module cc_interface_stencil_module
