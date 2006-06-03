module ml_interface_stencil_module

  use bl_types
  use multifab_module
  use stencil_module
  use bl_prof_module
  use bl_constants_module

  implicit none

contains

  subroutine ml_interface(res, flux, crse, ss, crse_domain, face, dim, efactor)
    type(multifab), intent(inout) :: res
    type(multifab), intent(in   ) :: flux
    type(multifab), intent(in   ) :: crse
    type(multifab), intent(in   ) :: ss
    type(box), intent(in) :: crse_domain
    integer, intent(in) :: face, dim
    real(kind=dp_t), intent(in) :: efactor
    type(bl_prof_timer), save :: bpt
    call build(bpt, "ml_interf")
    call ml_interface_c(res, 1, flux, 1, crse, ss, crse_domain, face, dim, efactor)
    call destroy(bpt)
  end subroutine ml_interface

  subroutine ml_interface_c(res, cr, flux, cf, crse, ss, crse_domain, face, dim, efactor)
    type(multifab), intent(inout) :: res
    type(multifab), intent(in   ) :: flux
    type(multifab), intent(in   ) :: crse
    type(multifab), intent(in   ) :: ss
    integer, intent(in) :: cr, cf
    type(box), intent(in) :: crse_domain
    integer, intent(in) :: face, dim
    real(kind=dp_t), intent(in) :: efactor

    type(box) :: rbox, fbox, cbox, sbox, isect
    integer :: lo (res%dim), hi (res%dim)
    integer :: loc(res%dim)
    integer :: lof(res%dim)
    integer :: lor(res%dim)
    integer :: los(res%dim)
    integer :: lo_dom(res%dim), hi_dom(res%dim)
    integer :: dm
    integer :: dir
    integer :: i, j

    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_interf_c")

    dm = res%dim
    dir = dim

    lo_dom = lwb(crse_domain)
    hi_dom = upb(crse_domain)

    do j = 1, crse%nboxes
       if ( remote(crse, j) ) cycle

       cbox = get_ibox(crse,j)
       loc = lwb(cbox) - crse%ng

       rbox = get_ibox(res,j)
       lor = lwb(rbox) - res%ng

       sbox = get_ibox(ss,j)
       los = lwb(sbox) - ss%ng

       cp => dataptr(crse, j, cr)
       rp => dataptr(res, j, cr)
       sp => dataptr(ss, j)

       do i = 1, flux%nboxes

          if ( remote(flux, i) ) cycle

          fbox = get_ibox(flux,i)
          lof = lwb(fbox)

          if ( contains(crse_domain,fbox) ) then
             isect = box_intersection(cbox,fbox)
             if ( empty(isect) ) cycle
             lo = lwb(isect)
             hi = upb(isect)
             fp => dataptr(flux, i, cf)
             select case (dm)
             case (1)
                call ml_interface_1d(rp(:,1,1,1), lor, &
                     fp(:,1,1,1), lof, &
                     cp(:,1,1,1), loc, &
                     sp(:,1,1,:), los, &
                     lo, hi, face, dim, efactor)
             case (2)
                call ml_interface_2d(rp(:,:,1,1), lor, &
                     fp(:,:,1,1), lof, &
                     cp(:,:,1,1), loc, &
                     sp(:,:,1,:), los, &
                     lo, hi, face, dim, efactor)
             case (3)
                call ml_interface_3d(rp(:,:,:,1), lor, &
                     fp(:,:,:,1), lof, &
                     cp(:,:,:,1), loc, &
                     sp(:,:,:,:), los, &
                     lo, hi, face, dim, efactor)
             end select
          end if
       end do
    end do
    call destroy(bpt)
  end subroutine ml_interface_c

  subroutine ml_interface_1d(res, lor, fine_flux, lof, cc, loc, &
       ss , los, lo, hi, face, dim, efactor)
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
    real(kind=dp_t), intent(in) :: efactor

    integer :: i
    real (kind = dp_t) :: crse_flux

    i = lo(1)

    !   Lo side
    if (face == -1) then
       crse_flux = ss(i,1)*(cc(i)-cc(i+1))
       res(i) = res(i) + efactor*(fine_flux(i) - crse_flux)

       !   Hi side
    else if (face == 1) then
       crse_flux = ss(i,2)*(cc(i)-cc(i-1))
       res(i) = res(i) + efactor*(fine_flux(i) - crse_flux)
    end if

  end subroutine ml_interface_1d

  subroutine ml_interface_2d(res, lor, fine_flux, lof, cc, loc, &
       ss , los, lo, hi, face, dim, efactor)
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
    real(kind=dp_t), intent(in) :: efactor

    integer :: i, j
    real (kind = dp_t) :: crse_flux

    !   Hi i side
    if ( dim == 1 ) then
       if (face == 1) then
          i = lo(1)
          do j = lo(2),hi(2)
             crse_flux = ss(i,j,2)*(cc(i,j)-cc(i-1,j))
             res(i,j) = res(i,j) + efactor*(fine_flux(i,j) - crse_flux)
          end do
          !   Lo i side
       else if (face == -1) then
          i = lo(1)
          do j = lo(2),hi(2)
             crse_flux = ss(i,j,1)*(cc(i,j)-cc(i+1,j))
             res(i,j) = res(i,j) + efactor*(fine_flux(i,j) - crse_flux)
          end do
       end if
    else if ( dim == 2 ) then
       !   Hi j side
       if (face == 1) then
          j = lo(2)
          do i = lo(1),hi(1)
             crse_flux = ss(i,j,4)*(cc(i,j)-cc(i,j-1))
             res(i,j) = res(i,j) + efactor*(fine_flux(i,j) - crse_flux)
          end do
          !   Lo j side
       else if (face == -1) then
          j = lo(2)
          do i = lo(1),hi(1)
             crse_flux = ss(i,j,3)*(cc(i,j)-cc(i,j+1))
             res(i,j) = res(i,j) + efactor*(fine_flux(i,j) - crse_flux)
          end do
       end if
    end if
  end subroutine ml_interface_2d

  subroutine ml_interface_3d(res, lor, fine_flux, lof, cc, loc, &
       ss , los, lo, hi, face, dim, efactor)
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
    real(kind=dp_t), intent(in) :: efactor

    integer :: i, j, k
    real (kind = dp_t) :: crse_flux

    !   Hi i side
    if ( dim == 1 ) then
       if (face == 1) then
          i = lo(1)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                crse_flux = ss(i,j,k,2)*(cc(i,j,k)-cc(i-1,j,k))
                res(i,j,k) = res(i,j,k) + efactor*(fine_flux(i,j,k) - crse_flux)
             end do
          end do
          !   Lo i side
       else if (face == -1) then
          i = lo(1)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                crse_flux = ss(i,j,k,1)*(cc(i,j,k)-cc(i+1,j,k))
                res(i,j,k) = res(i,j,k) + efactor*(fine_flux(i,j,k) - crse_flux)
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
                res(i,j,k) = res(i,j,k) + efactor*(fine_flux(i,j,k) - crse_flux)
             end do
          end do
          !   Lo j side
       else if (face == -1) then
          j = lo(2)
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)
                crse_flux = ss(i,j,k,3)*(cc(i,j,k)-cc(i,j+1,k))
                res(i,j,k) = res(i,j,k) + efactor*(fine_flux(i,j,k) - crse_flux)
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
                res(i,j,k) = res(i,j,k) + efactor*(fine_flux(i,j,k) - crse_flux)
             end do
          end do
          !   Lo k side
       else if (face == -1) then
          k = lo(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                crse_flux = ss(i,j,k,5)*(cc(i,j,k)-cc(i,j,k+1))
                res(i,j,k) = res(i,j,k) + efactor*(fine_flux(i,j,k) - crse_flux)
             end do
          end do
       end if
    end if
  end subroutine ml_interface_3d

  subroutine ml_crse_contrib(res, flux, crse, ss, mm_crse, mm_fine, crse_domain, ir, side)

    type(multifab), intent(inout) :: res
    type(multifab), intent(in   ) :: flux
    type(multifab), intent(in   ) :: crse
    type(multifab), intent(in   ) :: ss
    type(imultifab),intent(in   ) :: mm_crse
    type(imultifab),intent(in   ) :: mm_fine
    type(box),      intent(in   ) :: crse_domain
    integer,        intent(in   ) :: ir(:)
    integer                       :: side

    type(box) :: fbox, cbox, mbox, isect
    integer   :: lo (res%dim), hi (res%dim), loc(res%dim)
    integer   :: lof(res%dim), hif(res%dim), lor(res%dim), los(res%dim)
    integer   :: lomf(res%dim), lomc(res%dim), lo_dom(res%dim), hi_dom(res%dim)
    integer   :: lod(MAX_SPACEDIM), hid(MAX_SPACEDIM), loflux(res%dim), hiflux(res%dim)
    integer   :: dir, i, j, n, proc

    real(kind=dp_t), pointer   :: rp(:,:,:,:),fp(:,:,:,:),cp(:,:,:,:),sp(:,:,:,:)
    integer,         pointer   :: mp(:,:,:,:),mcp(:,:,:,:)
    integer,         parameter :: tag = 1103

    real(kind=dp_t), dimension(:,:,:,:), allocatable :: flxpt
    integer,         dimension(:,:,:,:), allocatable :: mmfpt

    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_crse_contrib")

    dir = iabs(side)

    lo_dom = lwb(crse_domain)
    hi_dom = upb(crse_domain)+1

    do j = 1, crse%nboxes

       cbox  = get_ibox(crse,j)
       loc   = lwb(get_pbox(crse,j))
       lomc  = lwb(get_pbox(mm_crse,j))
       lor   = lwb(get_pbox(res,j))
       los   = lwb(get_pbox(ss,j))

       do i = 1, flux%nboxes

          if ( remote(crse,j) .and. remote(flux,i) ) cycle

          fbox  = get_ibox(flux,i)
          isect = intersection(cbox,fbox)

          if ( empty(isect) ) cycle

          loflux = lwb(fbox)
          hiflux = upb(fbox)

          if ( ss%la%lap%pmask(dir) .or. (loflux(dir) /= lo_dom(dir) .and. loflux(dir) /= hi_dom(dir)) ) then

             lo = lwb(isect)
             hi = upb(isect)

             if ( local(crse,j) .and. local(flux,i) ) then

                lof  =  loflux
                hif  =  hiflux
                lomf =  lwb(get_pbox(mm_fine,i))
                fp   => dataptr(flux   , i)
                mp   => dataptr(mm_fine, i)
                cp   => dataptr(crse   , j)
                rp   => dataptr(res    , j)
                sp   => dataptr(ss     , j)
                mcp  => dataptr(mm_crse, j)

                select case (res%dim)
                case (1)
                   call ml_interface_1d_nodal(rp(:,1,1,1), lor, &
                        fp(:,1,1,1), lof, hif, &
                        cp(:,1,1,1), loc, &
                        sp(:,1,1,:), los, lo, hi, ir, side, loflux, hiflux)
                case (2)
                   call ml_interface_2d_nodal(rp(:,:,1,1), lor, &
                        fp(:,:,1,1), lof , hif, &
                        cp(:,:,1,1), loc , &
                        sp(:,:,1,:), los , &
                        mp(:,:,1,1), lomf, &
                        mcp(:,:,1,1), lomc, lo, hi, ir, side, loflux, hiflux)
                case (3)
                   call ml_interface_3d_nodal(rp(:,:,:,1), lor, &
                        fp(:,:,:,1), lof , hif, &
                        cp(:,:,:,1), loc , &
                        sp(:,:,:,:), los , &
                        mp(:,:,:,1), lomf, &
                        mcp(:,:,:,1), lomc, lo, hi, ir, side, loflux, hiflux)
                end select

             else if ( local(flux,i) ) then
                !
                ! Got to send flux & mm_fine to processor owning crse.
                !
                ! Need flux on isect; need mm_fine on isect*ir.
                !
                fp => dataptr(flux, i, isect)
                proc = get_proc(crse%la, j)
                call parallel_send(fp, proc, tag)
                isect%lo(1:res%dim) = isect%lo(1:res%dim) * ir(1:res%dim)
                isect%hi(1:res%dim) = isect%hi(1:res%dim) * ir(1:res%dim)
                mp => dataptr(mm_fine, i, isect)
                call parallel_send(mp, proc, tag)
             else
                !
                ! We own crse.  Got to get flux & mm_fine.
                !
                ! flux will be defined on isect; mm_fine on isect*ir.
                !
                lod = 1;                     hid = 1
                lod(1:res%dim) = lwb(isect); hid(1:res%dim) = upb(isect)
                allocate(flxpt(lod(1):hid(1),lod(2):hid(2),lod(3):hid(3),1:flux%nc))
                mbox = isect
                mbox%lo(1:res%dim) = mbox%lo(1:res%dim) * ir(1:res%dim)
                mbox%hi(1:res%dim) = mbox%hi(1:res%dim) * ir(1:res%dim)
                lod(1:res%dim) = lwb(mbox); hid(1:res%dim) = upb(mbox)
                allocate(mmfpt(lod(1):hid(1),lod(2):hid(2),lod(3):hid(3),1:mm_fine%nc))
                proc = get_proc(flux%la, i)
                call parallel_recv(flxpt, proc, tag)
                call parallel_recv(mmfpt, proc, tag)

                lof  =  lwb(isect)
                hif  =  upb(isect)
                lomf =  lwb(mbox)
                cp   => dataptr(crse   , j)
                rp   => dataptr(res    , j)
                sp   => dataptr(ss     , j)
                mcp  => dataptr(mm_crse, j)

                select case (res%dim)
                case (1)
                   call ml_interface_1d_nodal(rp(:,1,1,1), lor, &
                        flxpt(:,1,1,1), lof, hif, &
                        cp(:,1,1,1), loc, &
                        sp(:,1,1,:), los, lo, hi, ir, side, loflux, hiflux)
                case (2)
                   call ml_interface_2d_nodal(rp(:,:,1,1), lor, &
                        flxpt(:,:,1,1), lof , hif, &
                        cp(:,:,1,1), loc , &
                        sp(:,:,1,:), los , &
                        mmfpt(:,:,1,1), lomf, &
                        mcp(:,:,1,1), lomc, lo, hi, ir, side, loflux, hiflux)
                case (3)
                   call ml_interface_3d_nodal(rp(:,:,:,1), lor, &
                        flxpt(:,:,:,1), lof , hif, &
                        cp(:,:,:,1), loc , &
                        sp(:,:,:,:), los , &
                        mmfpt(:,:,:,1), lomf, &
                        mcp(:,:,:,1), lomc, lo, hi, ir, side, loflux, hiflux)
                end select

                deallocate(flxpt)
                deallocate(mmfpt)
             end if
          end if
       end do
    end do

    call destroy(bpt)

  end subroutine ml_crse_contrib

  subroutine ml_interface_1d_nodal(res, lor, fine_flux, lof, hif, cc, loc, &
       ss , los, lo, hi, ir, side, loflux, hiflux)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lof(:), hif(:), loflux(:), hiflux(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: res(lor(1):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):)
    real (kind = dp_t), intent(in   ) :: ss(los(1):,0:)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)

    integer :: i
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
       ss , los, mm_fine, lomf, mm_crse, lomc, lo, hi, ir, side, loflux, hiflux)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lomf(:)
    integer, intent(in) :: lomc(:)
    integer, intent(in) :: lof(:), hif(:), loflux(:), hiflux(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):)
    real (kind = dp_t), intent(in   ) ::        cc(loc(1):,loc(2):)
    real (kind = dp_t), intent(in   ) ::        ss(los(1):,los(2):,0:)
    integer           , intent(in   ) ::   mm_fine(lomf(1):,lomf(2):)
    integer           , intent(in   ) ::   mm_crse(lomc(1):,lomc(2):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)

    integer :: i, j
    real (kind = dp_t) :: crse_flux

    i = lo(1)
    j = lo(2)

    !   NOTE: THESE STENCILS ONLY WORK FOR DX == DY.

    if (side == -1) then

       do j = lo(2),hi(2)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0) .and. &
               (.not. bc_dirichlet(mm_crse(i,j),1,0))) then
             if (j == loflux(2) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,-1)) then
                crse_flux = (ss(i,j,8)*(cc(i+1,j+1) + HALF*cc(i+1,j) + &
                     HALF * cc(i,j+1) - TWO*cc(i,j))) * HALF
             else if (j == hiflux(2) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,+1)) then
                crse_flux = (ss(i,j,3)*(cc(i+1,j-1) + HALF*cc(i+1,j) + &
                     HALF * cc(i,j-1) - TWO*cc(i,j))) * HALF
             else
                crse_flux = ss(i,j,8)*(cc(i+1,j+1) + HALF*cc(i+1,j) + &
                     HALF * cc(i,j+1) - TWO*cc(i,j)) &
                     +ss(i,j,3)*(cc(i+1,j-1) + HALF*cc(i+1,j) + &
                     HALF * cc(i,j-1) - TWO*cc(i,j))
             end if
             res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
          end if
       end do

    else if (side ==  1) then

       do j = lo(2),hi(2)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0) .and. &
               (.not. bc_dirichlet(mm_crse(i,j),1,0))) then
             if (j == loflux(2) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,-1)) then
                crse_flux = (ss(i,j,6)*(cc(i-1,j+1) + HALF*cc(i-1,j) + &
                     HALF * cc(i,j+1) - TWO*cc(i,j))) * HALF
             else if (j == hiflux(2) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,+1)) then
                crse_flux = (ss(i,j,1)*(cc(i-1,j-1) + HALF*cc(i-1,j) + &
                     HALF * cc(i,j-1) - TWO*cc(i,j))) * HALF
             else
                crse_flux = ss(i,j,6)*(cc(i-1,j+1) + HALF*cc(i-1,j) + &
                     HALF * cc(i,j+1) - TWO*cc(i,j)) &
                     +ss(i,j,1)*(cc(i-1,j-1) + HALF*cc(i-1,j) + &
                     HALF * cc(i,j-1) - TWO*cc(i,j))
             end if
             res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
          end if
       end do

    else if (side == -2) then

       do i = lo(1),hi(1)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0) .and. &
               (.not. bc_dirichlet(mm_crse(i,j),1,0))) then
             if (i == loflux(1) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,-1)) then
                crse_flux = (ss(i,j,8)*(cc(i+1,j+1) + HALF*cc(i+1,j) + &
                     HALF * cc(i,j+1) - TWO*cc(i,j))) * HALF
             else if (i == hiflux(1) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,+1)) then
                crse_flux = (ss(i,j,6)*(cc(i-1,j+1) + HALF*cc(i-1,j) + &
                     HALF * cc(i,j+1) - TWO*cc(i,j))) * HALF
             else
                crse_flux = ss(i,j,8)*(cc(i+1,j+1) + HALF*cc(i+1,j) + &
                     HALF * cc(i,j+1) - TWO*cc(i,j)) &
                     +ss(i,j,6)*(cc(i-1,j+1) + HALF*cc(i-1,j) + &
                     HALF * cc(i,j+1) - TWO*cc(i,j))
             end if
             res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
          end if
       end do

    else if (side ==  2) then

       do i = lo(1),hi(1)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0) .and. &
               (.not. bc_dirichlet(mm_crse(i,j),1,0))) then
             if (i == loflux(1) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,-1)) then
                crse_flux = (ss(i,j,3)*(cc(i+1,j-1) + HALF*cc(i+1,j) + &
                     HALF * cc(i,j-1) - TWO*cc(i,j)) ) * HALF
             else if (i == hiflux(1) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,+1)) then
                crse_flux = (ss(i,j,1)*(cc(i-1,j-1) + HALF*cc(i-1,j) + &
                     HALF * cc(i,j-1) - TWO*cc(i,j)) ) * HALF
             else
                crse_flux = ss(i,j,3)*(cc(i+1,j-1) + HALF*cc(i+1,j) + &
                     HALF * cc(i,j-1) - TWO*cc(i,j)) &
                     +ss(i,j,1)*(cc(i-1,j-1) + HALF*cc(i-1,j) + &
                     HALF * cc(i,j-1) - TWO*cc(i,j))
             end if
             res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
          end if
       end do

    end if

  end subroutine ml_interface_2d_nodal

  subroutine ml_interface_3d_nodal(res, lor, fine_flux, lof, hif, cc, loc, &
       ss, los, mm_fine, lomf, mm_crse, lomc, lo, hi, ir, side, loflux, hiflux)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lof(:), hif(:), loflux(:), hiflux(:)
    integer, intent(in) :: lomf(:)
    integer, intent(in) :: lomc(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):,lor(3):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):,lof(3):)
    real (kind = dp_t), intent(in   ) ::        cc(loc(1):,loc(2):,loc(3):)
    real (kind = dp_t), intent(in   ) ::        ss(los(1):,los(2):,los(3):,0:)
    integer           , intent(in   ) ::   mm_fine(lomf(1):,lomf(2):,lomf(3):)
    integer           , intent(in   ) ::   mm_crse(lomc(1):,lomc(2):,lomc(3):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)

    integer :: i, j, k
    integer :: ioff, joff, koff
    integer :: sig_mm, sig_mp, sig_pm, sig_pp
    logical :: lo_i_not,lo_j_not,lo_k_not,hi_i_not,hi_j_not,hi_k_not
    real (kind = dp_t) :: crse_flux
    real (kind = dp_t) :: cell_mm, cell_mp, cell_pm, cell_pp

    i = lo(1)
    j = lo(2)
    k = lo(3)

    !   NOTE: THESE STENCILS ONLY WORK FOR DX == DY.

    !   NOTE: MM IS ON THE FINE GRID, NOT THE CRSE

    !   Lo/Hi i side
    if (side == -1 .or. side == 1) then

       if (side == -1) then
          ioff   = i+1
          sig_mm =  3
          sig_pm =  8
          sig_mp = 15
          sig_pp = 20
       else
          ioff   = i-1
          sig_mm =  1
          sig_pm =  6
          sig_mp = 13
          sig_pp = 18
       end if

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)

           if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,0) .and. &
                (.not.bc_dirichlet(mm_crse(i,j,k),1,0))) then

             lo_j_not = ( (j == loflux(2)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),2,-1) )
             hi_j_not = ( (j == hiflux(2)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),2,+1) )
             lo_k_not = ( (k == loflux(3)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),3,-1) )
             hi_k_not = ( (k == hiflux(3)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),3,+1) )

             cell_mm = ss(i,j,k,sig_mm)*(cc(ioff,j-1,k-1) + cc(ioff,j-1,k  ) &
                  +cc(ioff,j  ,k-1) + cc(i  ,j-1,k-1) - FOUR*cc(i  ,j  ,k) )
             cell_pm = ss(i,j,k,sig_pm)*(cc(ioff,j+1,k-1) + cc(ioff,j+1,k  ) &
                  +cc(ioff,j  ,k-1) + cc(i  ,j+1,k-1) - FOUR*cc(i  ,j  ,k) )
             cell_mp = ss(i,j,k,sig_mp)*(cc(ioff,j-1,k+1) + cc(ioff,j-1,k  ) &
                  +cc(ioff,j  ,k+1) + cc(i  ,j-1,k+1) - FOUR*cc(i  ,j  ,k) )
             cell_pp = ss(i,j,k,sig_pp)*(cc(ioff,j+1,k+1) + cc(ioff,j+1,k  ) &
                  +cc(ioff,j  ,k+1) + cc(i  ,j+1,k+1) - FOUR*cc(i  ,j  ,k) )

             crse_flux = zero

             if (lo_k_not) then
                if (lo_j_not) then
                   crse_flux = THIRD*cell_pp 
                else if (hi_j_not) then
                   crse_flux = THIRD*cell_mp
                else
                   crse_flux = HALF*(cell_pp + cell_mp)
                end if
             else if (hi_k_not) then
                if (lo_j_not) then
                   crse_flux = THIRD*cell_pm 
                else if (hi_j_not) then
                   crse_flux = THIRD*cell_mm 
                else
                   crse_flux = HALF*(cell_pm  + cell_mm)
                end if
             else 
                if (lo_j_not) then
                   crse_flux = HALF*(cell_pm  + cell_pp)
                else if (hi_j_not) then
                   crse_flux = HALF*(cell_mm  + cell_mp)
                else
                   crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
                end if
             end if

             res(i,j,k) = res(i,j,k) + crse_flux + fine_flux(i,j,k)
           end if
          end do
       end do

       !   Lo/Hi j side
    else if (side == -2 .or. side == 2) then

       if (side == -2) then
          joff   = j+1
          sig_mm =  6
          sig_pm =  8
          sig_mp = 18
          sig_pp = 20
       else
          joff   = j-1
          sig_mm =  1
          sig_pm =  3
          sig_mp = 13
          sig_pp = 15
       end if
       do k = lo(3),hi(3)
          do i = lo(1),hi(1)

           if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,0) .and. &
                (.not.bc_dirichlet(mm_crse(i,j,k),1,0))) then

             lo_i_not = ( (i == loflux(1)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,-1) )
             hi_i_not = ( (i == hiflux(1)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,+1) )
             lo_k_not = ( (k == loflux(3)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),3,-1) )
             hi_k_not = ( (k == hiflux(3)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),3,+1) )

             cell_mm = ss(i,j,k,sig_mm)*(cc(i-1,joff,k-1) + cc(i-1,joff,k  ) &
                  +cc(i  ,joff,k-1) + cc(i-1,j   ,k-1) - FOUR*cc(i  ,j  ,k) )
             cell_pm = ss(i,j,k,sig_pm)*(cc(i+1,joff,k-1) + cc(i+1,joff,k  ) &
                  +cc(i  ,joff,k-1) + cc(i+1,j   ,k-1) - FOUR*cc(i  ,j  ,k) )
             cell_mp = ss(i,j,k,sig_mp)*(cc(i-1,joff,k+1) + cc(i-1,joff,k  ) &
                  +cc(i  ,joff,k+1) + cc(i-1,j   ,k+1) - FOUR*cc(i  ,j  ,k) )
             cell_pp = ss(i,j,k,sig_pp)*(cc(i+1,joff,k+1) + cc(i+1,joff,k  ) &
                  +cc(i  ,joff,k+1) + cc(i+1,j   ,k+1) - FOUR*cc(i  ,j  ,k) )

             if (lo_k_not) then
                if (lo_i_not) then
                   crse_flux = THIRD*cell_pp 
                else if (hi_i_not) then
                   crse_flux = THIRD*cell_mp
                else
                   crse_flux = HALF*(cell_pp + cell_mp)
                end if
             else if (hi_k_not) then
                if (lo_i_not) then
                   crse_flux = THIRD*cell_pm 
                else if (hi_i_not) then
                   crse_flux = THIRD*cell_mm 
                else
                   crse_flux = HALF*(cell_pm  + cell_mm)
                end if
             else 
                if (lo_i_not) then
                   crse_flux = HALF*(cell_pm  + cell_pp)
                else if (hi_i_not) then
                   crse_flux = HALF*(cell_mm  + cell_mp)
                else
                   crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
                end if
             end if

             res(i,j,k) = res(i,j,k) + crse_flux + fine_flux(i,j,k)
           end if
          end do
       end do
       !   Lo/Hi k side
    else if (side == -3 .or. side == 3) then
       k = lo(3)
       if (side == -3) then
          koff   = k+1
          sig_mm = 13
          sig_pm = 15
          sig_mp = 18
          sig_pp = 20
       else
          koff   = k-1
          sig_mm =  1
          sig_pm =  3
          sig_mp =  6
          sig_pp =  8
       end if

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

           if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,0) .and. &
                (.not.bc_dirichlet(mm_crse(i,j,k),1,0))) then

             lo_i_not = ( (i == loflux(1)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,-1) )
             hi_i_not = ( (i == hiflux(1)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,+1) )
             lo_j_not = ( (j == loflux(2)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),2,-1) )
             hi_j_not = ( (j == hiflux(2)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),2,+1) )

             cell_mm = ss(i,j,k,sig_mm)*(cc(i-1,j-1,koff) + cc(i-1,j  ,koff) &
                  +cc(i  ,j-1,koff) + cc(i-1,j-1,k   ) - FOUR*cc(i  ,j  ,k) )
             cell_pm = ss(i,j,k,sig_pm)*(cc(i+1,j-1,koff) + cc(i+1,j  ,koff) &
                  +cc(i  ,j-1,koff) + cc(i+1,j-1,k   ) - FOUR*cc(i  ,j  ,k) )
             cell_mp = ss(i,j,k,sig_mp)*(cc(i-1,j+1,koff) + cc(i-1,j  ,koff) &
                  +cc(i  ,j+1,koff) + cc(i-1,j+1,k   ) - FOUR*cc(i  ,j  ,k) )
             cell_pp = ss(i,j,k,sig_pp)*(cc(i+1,j+1,koff) + cc(i+1,j  ,koff) &
                  +cc(i  ,j+1,koff) + cc(i+1,j+1,k   ) - FOUR*cc(i  ,j  ,k) )

             if (lo_j_not) then
                if (lo_i_not) then
                   crse_flux = THIRD*cell_pp 
                else if (hi_i_not) then
                   crse_flux = THIRD*cell_mp
                else
                   crse_flux = HALF*(cell_pp + cell_mp)
                end if
             else if (hi_j_not) then
                if (lo_i_not) then
                   crse_flux = THIRD*cell_pm 
                else if (hi_i_not) then
                   crse_flux = THIRD*cell_mm 
                else
                   crse_flux = HALF*(cell_pm  + cell_mm)
                end if
             else 
                if (lo_i_not) then
                   crse_flux = HALF*(cell_pm  + cell_pp)
                else if (hi_i_not) then
                   crse_flux = HALF*(cell_mm  + cell_mp)
                else
                   crse_flux = cell_mm  + cell_mp + cell_pm + cell_pp
                end if
             end if

             res(i,j,k) = res(i,j,k) + crse_flux + fine_flux(i,j,k)
           end if
          end do
       end do
    end if

  end subroutine ml_interface_3d_nodal

end module ml_interface_stencil_module
