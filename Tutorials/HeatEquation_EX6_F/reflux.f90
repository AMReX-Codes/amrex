module reflux_module

  use layout_module
  use multifab_module
  use multifab_fill_ghost_module
  use bndry_reg_module
  use bl_constants_module

  implicit none

  private

  public :: reflux, bndry_reg_add_fine_flx

contains

  subroutine reflux(phi_new,flux,bndry_flx,crse_domain,dx)
 
    type(multifab) , intent(inout) :: phi_new
    type(multifab) , intent(inout) :: flux(:)
    type(bndry_reg), intent(inout) :: bndry_flx
    type(box)      , intent(in   ) :: crse_domain
    real(kind=dp_t), intent(in   ) :: dx

    integer         :: i
    real(kind=dp_t) :: fac

    fac = 1.d0 / dx
 
    do i = 1, get_dim(phi_new)
       call reflux_i(phi_new, bndry_flx%obmf(i,0), flux(i), &
                     crse_domain, bndry_flx%ofacemap(:,i), bndry_flx%oindxmap(:,i), fac) 
    end do  
 
  end subroutine reflux

  subroutine reflux_i(phi, fine_flux, crse_flux, crse_domain, facemap, indxmap, efactor)
    use bl_prof_module

    type(multifab), intent(inout) :: phi
    type(multifab), intent(in   ) :: fine_flux
    type(multifab), intent(in   ) :: crse_flux
    type(box),      intent(in   ) :: crse_domain
    integer,        intent(in   ) :: facemap(:), indxmap(:)
    real(kind=dp_t),intent(in   ) :: efactor

    type(box)      :: fbox, cbox
    integer        :: lor(get_dim(phi)), i, j, dm, face, dim
    integer        :: lo (get_dim(phi)), hi (get_dim(phi)), loc(get_dim(phi))
    logical        :: pmask(get_dim(phi))
    real(kind=dp_t), pointer :: pp(:,:,:,:), fp(:,:,:,:), cp(:,:,:,:)

    pmask = get_pmask(get_layout(phi))
    dm    = get_dim(phi)

    !$OMP PARALLEL DO PRIVATE(fbox,cbox,lor,i,j,face,dim,lo,hi,loc,pp,fp,cp)
    do i = 1, nfabs(fine_flux)
       j = indxmap(i)
       dim = abs(facemap(i))
       face = sign(1, facemap(i))

       fbox = get_ibox(fine_flux,i)
       if (pmask(dim) .and. (.not. contains(crse_domain,fbox)) ) then
          if ( face .eq. -1 ) then
             fbox = shift(fbox,  extent(crse_domain,dim), dim)
          else
             fbox = shift(fbox, -extent(crse_domain,dim), dim)
          end if
       end if

       lo = lwb(fbox)
       hi = upb(fbox)
       fp => dataptr(fine_flux,i)

       cbox =  get_ibox(crse_flux,j)
       loc  =  lwb(get_pbox(crse_flux,j))
       lor  =  lwb(get_pbox(phi,j))
       pp   => dataptr(phi, j)
       cp   => dataptr(crse_flux, j)

       select case (dm)
       case (2)
          call reflux_2d(pp(:,:,1,1), lor, cp(:,:,1,1), loc, &
                fp(:,:,1,1), lo, lo, hi, face, dim, efactor)
       case (3)
          call reflux_3d(pp(:,:,:,1), lor, cp(:,:,:,1), loc, &
               fp(:,:,:,1), lo, lo, hi, face, dim, efactor)
       end select
    end do
    !$OMP END PARALLEL DO

  end subroutine reflux_i

  subroutine reflux_2d(phi, lor, crse_flux, loc, fine_flux, lof, &
       lo, hi, face, dim, efactor)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: phi(lor(1):,lor(2):)
    real (kind = dp_t), intent(in   ) :: crse_flux(loc(1):,loc(2):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):)
    integer, intent(in) :: face, dim
    real(kind=dp_t), intent(in) :: efactor

    integer :: i, j

    ! Note:  fine_flux is only one cell wide and is indexed by the coarse cell which is being changed 
    ! Note:  crse_flux is edge-based so we have to use the lo-side value when face =  1, and 
    !                                                  the hi-side value when face = -1

    if ( dim == 1 ) then
       i = lo(1)
       if (face .eq. -1) then
          do j = lo(2),hi(2)
             phi(i,j) = phi(i,j) + efactor*(fine_flux(i,j)-crse_flux(i+1,j))
          end do
       else
          do j = lo(2),hi(2)
             phi(i,j) = phi(i,j) - efactor*(fine_flux(i,j)-crse_flux(i,j))
          end do
       end if
    else if ( dim == 2 ) then
       j = lo(2)
       if (face .eq. -1) then
          do i = lo(1),hi(1)
             phi(i,j) = phi(i,j) + efactor*(fine_flux(i,j)-crse_flux(i,j+1))
          end do
       else
          do i = lo(1),hi(1)
             phi(i,j) = phi(i,j) - efactor*(fine_flux(i,j)-crse_flux(i,j))
          end do
       end if
    end if

  end subroutine reflux_2d

  subroutine reflux_3d(phi, lor, crse_flux, loc, fine_flux, lof, &
       lo, hi, face, dim, efactor)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::       phi(lor(1):,lor(2):,lor(3):)
    real (kind = dp_t), intent(in   ) :: crse_flux(loc(1):,loc(2):,loc(3):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):,lof(3):)
    integer, intent(in) :: face, dim
    real(kind=dp_t), intent(in) :: efactor

    integer :: i, j, k

    if ( dim == 1 ) then
       i = lo(1)
       if (face .eq. -1) then
          do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             phi(i,j,k) = phi(i,j,k) + efactor*(fine_flux(i,j,k)-crse_flux(i+1,j,k))
          end do
          end do
       else
          do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             phi(i,j,k) = phi(i,j,k) - efactor*(fine_flux(i,j,k)-crse_flux(i,j,k))
          end do
          end do
       end if
    else if ( dim == 2 ) then
       j = lo(2)
       if (face .eq. -1) then
          do k = lo(3),hi(3)
          do i = lo(1),hi(1)
             phi(i,j,k) = phi(i,j,k) + efactor*(fine_flux(i,j,k)-crse_flux(i,j+1,k))
          end do
          end do
       else
          do k = lo(3),hi(3)
          do i = lo(1),hi(1)
             phi(i,j,k) = phi(i,j,k) - efactor*(fine_flux(i,j,k)-crse_flux(i,j,k))
          end do
          end do
       end if
    else if ( dim == 3 ) then
       k = lo(3)
       if (face .eq. -1) then
          do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             phi(i,j,k) = phi(i,j,k) + efactor*(fine_flux(i,j,k)-crse_flux(i,j,k+1))
          end do
          end do
       else
          do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             phi(i,j,k) = phi(i,j,k) - efactor*(fine_flux(i,j,k)-crse_flux(i,j,k))
          end do
          end do
       end if
    end if

  end subroutine reflux_3d

  subroutine bndry_reg_add_fine_flx(bndry_flx, flux, ref_ratio)

    type(bndry_reg), intent(inout) :: bndry_flx
    type(multifab) , intent(inout) :: flux(:)
    integer        , intent(in   ) :: ref_ratio(:)

    integer :: i

    do i = 1, bndry_flx%dim
       call bndry_reg_add_fine_flx_fab( &
            bndry_flx%bmf(i,0), &
            flux(i), ref_ratio(i), &
            bndry_flx%facemap(:,i), bndry_flx%indxmap(:,i))
    end do

  end subroutine bndry_reg_add_fine_flx

  subroutine bndry_reg_add_fine_flx_fab(fine_flx, flux, ratio, &
                                        facemap, indxmap)

    use bl_prof_module

    type(multifab), intent(inout) :: fine_flx
    type(multifab), intent(in   ) :: flux
    integer       , intent(in   ) :: ratio
    integer       , intent(in   ) :: facemap(:)
    integer       , intent(in   ) :: indxmap(:)

    integer              :: i, j, dm, dim, face
    real(dp_t), pointer  :: fp(:,:,:,:)
    real(dp_t), pointer  :: up(:,:,:,:)
    type(box)            :: ba, bac
    integer, allocatable :: lo(:), hi(:), loc(:)

    dm = get_dim(flux)

    allocate(lo(dm),hi(dm),loc(dm))

    !$OMP PARALLEL DO PRIVATE(i,j,dim,face,fp,up,n)
    do i = 1, nfabs(fine_flx)
       j    = indxmap(i)
       dim  = abs(facemap(i))
       face = sign(1, facemap(i))

       ba = get_box(flux,j)
       lo = lwb(ba)
       hi = upb(ba)

       bac = get_box(fine_flx,i)
       loc  = lwb(bac)

       fp => dataptr(fine_flx, i)
       up => dataptr(flux, j)

       select case(dm)
       case (2)
          call flux_add_2d(fp(:,:,1,1), loc, up(:,:,1,1), lo, hi, &
                            ratio, face, dim)
       case (3)
          call flux_add_3d(fp(:,:,:,1), loc, up(:,:,:,1), lo, hi, &
                            ratio, face, dim)
          end select
    end do

    deallocate(lo,hi,loc)

  end subroutine bndry_reg_add_fine_flx_fab

  subroutine flux_add_2d(brf, loc, flux, lo, hi, ratio, face, dim)

    integer           , intent(in   ) ::  lo(:), hi(:)
    integer           , intent(in   ) :: loc(:)
    real (kind = dp_t), intent(inout) ::  brf(loc(1):,loc(2):)
    real (kind = dp_t), intent(in   ) :: flux(lo(1):,lo(2):)
    integer           , intent(in   ) :: ratio, face, dim

    integer    :: i,j,ic,jc
    real(dp_t) :: fac

    fac = 1.d0 / dble(ratio)

    ! Note that the bndryreg object, brf, is only one cell wide in the
    !     direction normal to the face of the grid it is associated with
    ! 
    ! The flux arrays, flux, are defined on the faces of every cell in each grid.

    if ( dim == 1 ) then

       if (face == -1) then 
          i = lo(1)
       else 
          i = hi(1)+1
       end if

       do j = lo(2), hi(2)
          jc = j/ratio
          brf(loc(1),jc) = brf(loc(1),jc) +  flux(i,j) * fac
       end do

    else if ( dim == 2 ) then

       if (face == -1) then 
          j = lo(2)
       else
          j = hi(2)+1
       end if

       do i = lo(1), hi(1)
          ic = i/ratio
          brf(ic,loc(2)) = brf(ic,loc(2)) + flux(i,j) * fac
       end do

    end if

  end subroutine flux_add_2d

  subroutine flux_add_3d(brf, loc, flux, lo, hi, ratio, face, dim)

    integer           , intent(in   ) ::  lo(:), hi(:)
    integer           , intent(in   ) :: loc(:)
    real (kind = dp_t), intent(inout) ::  brf(loc(1):,loc(2):,loc(3):)
    real (kind = dp_t), intent(in   ) :: flux( lo(1):, lo(2):, lo(3):)
    integer           , intent(in   ) :: ratio, face, dim

    integer    :: i,j,k,ic,jc,kc
    real(dp_t) :: fac

    fac = 1.d0 / dble(ratio*ratio)

    ! Note that the bndryreg object, brf, is only one cell wide in the
    !     direction normal to the face of the grid it is associated with
    ! 
    ! The flux arrays, flux, are defined on the faces of every cell in each grid.

    if ( dim == 1 ) then

       if (face == -1) then 
          i = lo(1)
       else 
          i = hi(1)+1
       end if

       do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          jc = j/ratio
          kc = k/ratio
          brf(loc(1),jc,kc) = brf(loc(1),jc,kc) + flux(i,j,k) * fac
       end do
       end do

    else if ( dim == 2 ) then

       if (face == -1) then 
          j = lo(2)
       else
          j = hi(2)+1
       end if

       do k = lo(3), hi(3)
       do i = lo(1), hi(1)
          kc = k/ratio
          ic = i/ratio
          brf(ic,loc(2),kc) = brf(ic,loc(2),kc) + flux(i,j,k) * fac
       end do
       end do

    else if ( dim == 3 ) then

       if (face == -1) then 
          k = lo(3)
       else 
          k = hi(3)+1
       end if

       do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          jc = j/ratio
          ic = i/ratio
          brf(ic,jc,loc(3)) = brf(ic,jc,loc(3)) + flux(i,j,k) * fac
       end do
       end do

    end if

  end subroutine flux_add_3d

end module reflux_module
