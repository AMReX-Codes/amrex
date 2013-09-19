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

  subroutine ml_interface(res, flux, crse, ss, crse_domain, face, dim, efactor)
    use bl_prof_module
    type(multifab), intent(inout) :: res
    type(multifab), intent(in   ) :: flux
    type(multifab), intent(in   ) :: crse
    type(multifab), intent(in   ) :: ss
    type(box),      intent(in   ) :: crse_domain
    integer,        intent(in   ) :: face, dim
    real(kind=dp_t),intent(in   ) :: efactor
    type(bl_prof_timer), save :: bpt
    call build(bpt, "ml_interf")
    call ml_interface_c(res, 1, flux, 1, crse, ss, crse_domain, face, dim, efactor)
    call destroy(bpt)
  end subroutine ml_interface

  subroutine ml_interface_c(res, cr, flux, cf, crse, ss, crse_domain, face, dim, efactor)
    use bl_prof_module
    use vector_i_module

    type(multifab), intent(inout) :: res
    type(multifab), intent(in   ) :: flux
    type(multifab), intent(in   ) :: crse
    type(multifab), intent(in   ) :: ss
    integer,        intent(in   ) :: cr, cf
    type(box),      intent(in   ) :: crse_domain
    integer,        intent(in   ) :: face, dim
    real(kind=dp_t),intent(in   ) :: efactor

    type(box)      :: fbox, cbox, isect
    integer        :: lor(get_dim(res)), los(get_dim(res)), i, j, k, shft, dm
    integer        :: lo (get_dim(res)), hi (get_dim(res)), loc(get_dim(res)), proc
    logical        :: pmask(get_dim(res))
    type(multifab) :: tflux
    type(list_box) :: bl
    type(boxarray) :: ba
    type(layout)   :: la
    type(vector_i) :: procmap, indxmap, shftmap

    type(box_intersector), pointer :: bi(:)

    real(kind=dp_t), pointer :: rp(:,:,:,:), fp(:,:,:,:), cp(:,:,:,:), sp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_interf_c")

    pmask = get_pmask(get_layout(res))
    !
    ! Build layout used in intersection test for below loop.
    !
    call copy(ba, get_boxarray(flux))

    !$OMP PARALLEL DO PRIVATE(i,fbox)
    do i = 1, nboxes(flux%la)
       fbox = box_nodalize(get_box(flux%la,i),flux%nodal)
       if ( pmask(dim) .and. (.not. contains(crse_domain,fbox)) ) then
          if ( face .eq. -1 ) then
             fbox = shift(fbox,  extent(crse_domain,dim), dim)
          else
             fbox = shift(fbox, -extent(crse_domain,dim), dim)
          end if
          call set_box(ba, i, fbox)
       end if
    end do
    !$OMP END PARALLEL DO

    call build(la, ba, boxarray_bbox(ba), mapping = LA_LOCAL)  ! LA_LOCAL ==> bypass processor distribution calculation.

    call destroy(ba)

    dm = get_dim(res)

    !$OMP PARALLEL DO PRIVATE(j,k,cbox,loc,lor,los,sp,rp,cp,bi,lo,hi)
    do j = 1, nfabs(crse)
       cbox =  get_ibox(crse,j)
       loc  =  lwb(get_pbox(crse,j))
       lor  =  lwb(get_pbox(res,j))
       los  =  lwb(get_pbox(ss,j))
       sp   => dataptr(ss, j)
       rp   => dataptr(res, j, cr)
       cp   => dataptr(crse, j, cr)
       bi   => layout_get_box_intersector(la, cbox)

       do k = 1, size(bi)
          lo = lwb(bi(k)%bx)
          hi = upb(bi(k)%bx)

          select case (dm)
          case (1)
             call ml_interface_1d_crse(rp(:,1,1,1), lor, cp(:,1,1,1), loc, sp(:,:,1,1), los, lo, face, efactor)
          case (2)
             call ml_interface_2d_crse(rp(:,:,1,1), lor, cp(:,:,1,1), loc, sp(:,:,:,1), los, lo, hi, face, dim, efactor)
          case (3)
             call ml_interface_3d_crse(rp(:,:,:,1), lor, cp(:,:,:,1), loc, sp(:,:,:,:), los, lo, hi, face, dim, efactor)
          end select
       end do

      deallocate(bi)
    end do
    !$OMP END PARALLEL DO
    !
    ! Build a multifab based on the intersections of flux with crse in such a way that each 
    ! intersecting box is owned by the same CPU as that owning the appropriate box in crse.
    !
    call build(procmap)
    call build(indxmap)
    call build(shftmap)

    do j = 1, nboxes(crse%la)
       cbox = box_nodalize(get_box(crse%la,j),crse%nodal)
       proc =  get_proc(crse%la,j)
       bi   => layout_get_box_intersector(la, cbox)

       do k = 1, size(bi)
          shft  = 0
          isect = bi(k)%bx
          fbox  = box_nodalize(get_box(flux%la,bi(k)%i),flux%nodal)

          if ( pmask(dim) .and. (.not. contains(crse_domain,fbox)) ) then
             !
             ! We need to remember the original flux box & whether or not it needs to be shifted.
             !
             shft = 1
             if ( face .eq. -1 ) then
                isect = shift(isect, -extent(crse_domain,dim), dim)
             else
                isect = shift(isect,  extent(crse_domain,dim), dim)
             end if
          end if

          call push_back(bl, isect)
          call push_back(procmap, proc)
          call push_back(indxmap, j)
          call push_back(shftmap, shft)
       end do

       deallocate(bi)
    end do

    call destroy(la)

    if ( empty(procmap) ) then
       !
       ! Nothing else to do ...
       !
       call destroy(bpt)
       call destroy(shftmap)
       call destroy(indxmap)
       call destroy(procmap)
       return
    end if

    call build(ba, bl, sort = .false.)
    call destroy(bl)
    call build(la, ba, boxarray_bbox(ba), pmask = pmask, explicit_mapping = dataptr(procmap, 1, size(procmap)))
    call destroy(ba)
    call build(tflux, la, nc = ncomp(flux), ng = 0)
    call copy(tflux, 1, flux, cf)  ! parallel copy

    !$OMP PARALLEL DO PRIVATE(i,j,cbox,fbox,lor,rp,lo,hi,fp)
    do i = 1, nfabs(tflux)
       j    = at(indxmap,global_index(tflux,i))
       cbox = box_nodalize(get_box(crse%la,j),crse%nodal)
       lor  =  lwb(get_pbox(res,local_index(res,j)))
       rp   => dataptr(res, local_index(res,j), cr)
       fbox =  get_ibox(tflux,i)

       if ( at(shftmap,global_index(tflux,i)) .eq. 1 ) then
          if (face .eq. -1) then
             fbox = shift(fbox,  extent(crse_domain,dim), dim)
          else
             fbox = shift(fbox, -extent(crse_domain,dim), dim)
          end if
       end if

       call bl_assert(intersects(cbox,fbox), 'ml_interface_c(): how did this happen?')

       lo  =  lwb(fbox)
       hi  =  upb(fbox)
       fp  => dataptr(tflux, i, 1)

       select case (dm)
       case (1)
          call ml_interface_1d_fine(rp(:,1,1,1), lor, fp(:,1,1,1), lo, lo, efactor)
       case (2)
          call ml_interface_2d_fine(rp(:,:,1,1), lor, fp(:,:,1,1), lo, lo, hi, efactor)
       case (3)
          call ml_interface_3d_fine(rp(:,:,:,1), lor, fp(:,:,:,1), lo, lo, hi, efactor)
       end select
    end do
    !$OMP END PARALLEL DO

    call destroy(shftmap)
    call destroy(indxmap)
    call destroy(procmap)
    call destroy(tflux)
    call destroy(la)
    call destroy(bpt)

    contains

      subroutine ml_interface_1d_fine(res, lor, fine_flux, lof, lo, efactor)
        integer, intent(in) :: lor(:)
        integer, intent(in) :: lof(:) 
        integer, intent(in) :: lo(:)
        real (kind = dp_t), intent(inout) :: res(lor(1):)
        real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):)
        real (kind = dp_t), intent(in   ) :: efactor

        integer :: i 

        i = lo(1)
        res(i) = res(i) + efactor*fine_flux(i)
      end subroutine ml_interface_1d_fine

      subroutine ml_interface_2d_fine(res, lor, fine_flux, lof, lo, hi, efactor)
        integer, intent(in) :: lor(:)
        integer, intent(in) :: lof(:)
        integer, intent(in) :: lo(:), hi(:)
        real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):)
        real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):)
        real (kind = dp_t), intent(in   ) :: efactor

        integer :: i, j

        do j = lo(2),hi(2)
           do i = lo(1),hi(1)
              res(i,j) = res(i,j) + efactor*fine_flux(i,j)
           end do
        end do
      end subroutine ml_interface_2d_fine

      subroutine ml_interface_3d_fine(res, lor, fine_flux, lof, lo, hi, efactor)
        integer, intent(in) :: lor(:)
        integer, intent(in) :: lof(:)
        integer, intent(in) :: lo(:), hi(:)
        real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):,lor(3):)
        real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):,lof(3):)
        real( kind = dp_t), intent(in   ) :: efactor

        integer :: i, j, k

        do k = lo(3),hi(3)
           do j = lo(2),hi(2)
              do i = lo(1),hi(1)
                 res(i,j,k) = res(i,j,k) + efactor*fine_flux(i,j,k)
              end do
           end do
        end do
      end subroutine ml_interface_3d_fine

  end subroutine ml_interface_c

  subroutine ml_interface_1d_crse(res, lor, cc, loc, &
       ss , los, lo, face, efactor)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lo(:)
    real (kind = dp_t), intent(inout) :: res(lor(1):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):)
    real (kind = dp_t), intent(in   ) :: ss(0:,los(1):)
    integer, intent(in) :: face
    real(kind=dp_t), intent(in) :: efactor

    integer :: i
    real (kind = dp_t) :: crse_flux

    i = lo(1)

    !   Lo side
    if (face == -1) then
       crse_flux = ss(1,i)*(cc(i)-cc(i+1))
       res(i) = res(i) - efactor*crse_flux

       !   Hi side
    else if (face == 1) then
       crse_flux = ss(2,i)*(cc(i)-cc(i-1))
       res(i) = res(i) - efactor*crse_flux
    end if

  end subroutine ml_interface_1d_crse

  subroutine ml_interface_2d_crse(res, lor, cc, loc, &
                                  ss , los, lo, hi, face, dim, efactor)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: res(lor(1):,lor(2):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):,loc(2):)
    real (kind = dp_t), intent(in   ) :: ss(0:,los(1):,los(2):)
    integer, intent(in) :: face, dim
    real(kind=dp_t), intent(in) :: efactor

    integer :: i, j, ns
    real (kind = dp_t) :: crse_flux
 
    ns = size(ss,dim=1)

    !   Hi i side
    if ( dim == 1 ) then
       if (face == 1) then
          i = lo(1)
          do j = lo(2),hi(2)
             if (ns.eq.7) then
                crse_flux = ss(2,i,j)*(cc(i,j)-cc(i-1,j))
             else if (ns.eq.9) then
                crse_flux = &
                   (15.d0/16.d0)*ss(2,i,j)*(cc(i-1,j)-cc(i  ,j)) &
                 +               ss(1,i,j)*(cc(i-2,j)-cc(i+1,j))

             endif
             res(i,j) = res(i,j) - efactor*crse_flux
          end do
          !   Lo i side
       else if (face == -1) then
          i = lo(1)
          do j = lo(2),hi(2)
             if (ns.eq.7) then
                crse_flux = ss(1,i,j)*(cc(i,j)-cc(i+1,j))
             else if (ns.eq.9) then
                crse_flux = &
                   (15.d0/16.d0)*ss(3,i,j)*(cc(i+1,j)-cc(i  ,j)) &
                 +               ss(4,i,j)*(cc(i+2,j)-cc(i-1,j))

             endif
             res(i,j) = res(i,j) - efactor*crse_flux
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
                   (15.d0/16.d0)*ss(6,i,j)*(cc(i,j-1)-cc(i,j  )) &
                 +               ss(5,i,j)*(cc(i,j-2)-cc(i,j+1))

             endif
             res(i,j) = res(i,j) - efactor*crse_flux
          end do
          !   Lo j side
       else if (face == -1) then
          j = lo(2)
          do i = lo(1),hi(1)
             if (ns.eq.7) then
                crse_flux = ss(3,i,j)*(cc(i,j)-cc(i,j+1))
             else if (ns.eq.9) then
                crse_flux = &
                   (15.d0/16.d0)*ss(7,i,j)*(cc(i,j+1)-cc(i,j  )) &
                 +               ss(8,i,j)*(cc(i,j+2)-cc(i,j-1))
             endif
             res(i,j) = res(i,j) - efactor*crse_flux
          end do
       end if
    end if
  end subroutine ml_interface_2d_crse

  subroutine ml_interface_3d_crse(res, lor, cc, loc, &
       ss , los, lo, hi, face, dim, efactor)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):,lor(3):)
    real (kind = dp_t), intent(in   ) ::        cc(loc(1):,loc(2):,loc(3):)
    real (kind = dp_t), intent(in   ) ::        ss(0:,los(1):,los(2):,los(3):)
    integer, intent(in) :: face, dim
    real(kind=dp_t), intent(in) :: efactor

    integer :: i, j, k

    !   Hi i side
    if ( dim == 1 ) then
       if (face == 1) then
          i = lo(1)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                res(i,j,k) = res(i,j,k) - efactor*(ss(2,i,j,k)*(cc(i,j,k)-cc(i-1,j,k)))
             end do
          end do
          !   Lo i side
       else if (face == -1) then
          i = lo(1)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                res(i,j,k) = res(i,j,k) - efactor*(ss(1,i,j,k)*(cc(i,j,k)-cc(i+1,j,k)))
             end do
          end do
       end if
       !   Hi j side
    else if ( dim ==  2 )  then
       if (face == 1) then
          j = lo(2)
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)
                res(i,j,k) = res(i,j,k) - efactor*(ss(4,i,j,k)*(cc(i,j,k)-cc(i,j-1,k)))
             end do
          end do
          !   Lo j side
       else if (face == -1) then
          j = lo(2)
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)
                res(i,j,k) = res(i,j,k) - efactor*(ss(3,i,j,k)*(cc(i,j,k)-cc(i,j+1,k)))
             end do
          end do
       end if
    else if ( dim == 3 ) then
       !   Hi k side
       if (face == 1) then
          k = lo(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                res(i,j,k) = res(i,j,k) - efactor*(ss(6,i,j,k)*(cc(i,j,k)-cc(i,j,k-1)))
             end do
          end do
          !   Lo k side
       else if (face == -1) then
          k = lo(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                res(i,j,k) = res(i,j,k) - efactor*(ss(5,i,j,k)*(cc(i,j,k)-cc(i,j,k+1)))
             end do
          end do
       end if
    end if

  end subroutine ml_interface_3d_crse

end module cc_interface_stencil_module
