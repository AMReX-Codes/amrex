module ml_interface_stencil_module

  use bl_types
  use layout_module
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
    integer :: dm
    integer :: i, j
    logical :: pmask(res%dim)

    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_interf_c")

    pmask = layout_get_pmask(res%la)

    dm = res%dim

    ! This loop only operates on crse data on the same crse layout -- should be
    !   easy to parallelize
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

          if (pmask(dim)) then
             isect = box_intersection(crse_domain,fbox)
             if ( empty(isect) ) then
                if (face .eq. -1) then
                   fbox = box_shift_d(fbox,extent(crse_domain,dim),dim)
                else
                   fbox = box_shift_d(fbox,-extent(crse_domain,dim),dim)
                end if
             end if
          end if

          lof = lwb(fbox)

          isect = box_intersection(cbox,fbox)
          if ( empty(isect) ) cycle
          lo = lwb(isect)
          hi = upb(isect)
          fp => dataptr(flux, i, cf)
          select case (dm)
          case (1)
             call ml_interface_1d_crse( &
                  rp(:,1,1,1), lor, &
                  cp(:,1,1,1), loc, &
                  sp(:,1,1,:), los, &
                  lo, hi, face, dim, efactor)
          case (2)
             call ml_interface_2d_crse( &
                  rp(:,:,1,1), lor, &
                  cp(:,:,1,1), loc, &
                  sp(:,:,1,:), los, &
                  lo, hi, face, dim, efactor)
          case (3)
             call ml_interface_3d_crse( &
                  rp(:,:,:,1), lor, &
                  cp(:,:,:,1), loc, &
                  sp(:,:,:,:), los, &
                  lo, hi, face, dim, efactor)
          end select
       end do
    end do

    ! This loop only uses flux registers from the fine layout to modify crse data on the crse layout 
    !   -- harder to parallelize
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

          if (pmask(dim)) then
             isect = box_intersection(crse_domain,fbox)
             if ( empty(isect) ) then
                if (face .eq. -1) then
                   fbox = box_shift_d(fbox,extent(crse_domain,dim),dim)
                else
                   fbox = box_shift_d(fbox,-extent(crse_domain,dim),dim)
                end if
             end if
          end if

          lof = lwb(fbox)

          isect = box_intersection(cbox,fbox)
          if ( empty(isect) ) cycle
          lo = lwb(isect)
          hi = upb(isect)
          fp => dataptr(flux, i, cf)
          select case (dm)
          case (1)
             call ml_interface_1d_fine( &
                  rp(:,1,1,1), lor, &
                  fp(:,1,1,1), lof, &
                  lo, hi, efactor)
          case (2)
             call ml_interface_2d_fine( &
                  rp(:,:,1,1), lor, &
                  fp(:,:,1,1), lof, &
                  lo, hi, efactor)
          case (3)
             call ml_interface_3d_fine( &
                  rp(:,:,:,1), lor, &
                  fp(:,:,:,1), lof, &
                  lo, hi, efactor)
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine ml_interface_c

  subroutine ml_interface_1d_crse(res, lor, cc, loc, &
       ss , los, lo, hi, face, dim, efactor)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: res(lor(1):)
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
       res(i) = res(i) - efactor*crse_flux

       !   Hi side
    else if (face == 1) then
       crse_flux = ss(i,2)*(cc(i)-cc(i-1))
       res(i) = res(i) - efactor*crse_flux
    end if

  end subroutine ml_interface_1d_crse

  subroutine ml_interface_1d_fine(res, lor, fine_flux, lof, lo, hi, efactor)

    integer, intent(in) :: lor(:)
    integer, intent(in) :: lof(:) 
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: res(lor(1):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):)
    real (kind = dp_t), intent(in   ) :: efactor

    integer :: i 

    i = lo(1)
    res(i) = res(i) + efactor*fine_flux(i)

  end subroutine ml_interface_1d_fine

  subroutine ml_interface_2d_crse(res, lor, cc, loc, &
                                  ss , los, lo, hi, face, dim, efactor)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):)
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
             res(i,j) = res(i,j) - efactor*crse_flux
          end do
          !   Lo i side
       else if (face == -1) then
          i = lo(1)
          do j = lo(2),hi(2)
             crse_flux = ss(i,j,1)*(cc(i,j)-cc(i+1,j))
             res(i,j) = res(i,j) - efactor*crse_flux
          end do
       end if
    else if ( dim == 2 ) then
       !   Hi j side
       if (face == 1) then
          j = lo(2)
          do i = lo(1),hi(1)
             crse_flux = ss(i,j,4)*(cc(i,j)-cc(i,j-1))
             res(i,j) = res(i,j) - efactor*crse_flux
          end do
          !   Lo j side
       else if (face == -1) then
          j = lo(2)
          do i = lo(1),hi(1)
             crse_flux = ss(i,j,3)*(cc(i,j)-cc(i,j+1))
             res(i,j) = res(i,j) - efactor*crse_flux
          end do
       end if
    end if
  end subroutine ml_interface_2d_crse

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

  subroutine ml_interface_3d_crse(res, lor, cc, loc, &
       ss , los, lo, hi, face, dim, efactor)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):,lor(3):)
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
                res(i,j,k) = res(i,j,k) - efactor*crse_flux
             end do
          end do
          !   Lo i side
       else if (face == -1) then
          i = lo(1)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                crse_flux = ss(i,j,k,1)*(cc(i,j,k)-cc(i+1,j,k))
                res(i,j,k) = res(i,j,k) - efactor*crse_flux
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
                res(i,j,k) = res(i,j,k) - efactor*crse_flux
             end do
          end do
          !   Lo j side
       else if (face == -1) then
          j = lo(2)
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)
                crse_flux = ss(i,j,k,3)*(cc(i,j,k)-cc(i,j+1,k))
                res(i,j,k) = res(i,j,k) - efactor*crse_flux
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
                res(i,j,k) = res(i,j,k) - efactor*crse_flux
             end do
          end do
          !   Lo k side
       else if (face == -1) then
          k = lo(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                crse_flux = ss(i,j,k,5)*(cc(i,j,k)-cc(i,j,k+1))
                res(i,j,k) = res(i,j,k) - efactor*crse_flux
             end do
          end do
       end if
    end if

  end subroutine ml_interface_3d_crse

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

  subroutine ml_crse_contrib_fancy(res, flux, crse, ss, mm_crse, mm_fine, crse_domain, ir, side)

    type(multifab), intent(inout) :: res
    type(multifab), intent(in   ) :: flux
    type(multifab), intent(in   ) :: crse
    type(multifab), intent(in   ) :: ss
    type(imultifab),intent(in   ) :: mm_crse
    type(imultifab),intent(in   ) :: mm_fine
    type(box),      intent(in   ) :: crse_domain
    integer,        intent(in   ) :: ir(:)
    integer                       :: side

    type(box) :: fbox, mbox, isect
    integer   :: lo (res%dim), hi (res%dim), loc(res%dim)
    integer   :: lof(res%dim), hif(res%dim), lor(res%dim), los(res%dim)
    integer   :: lomf(res%dim), lomc(res%dim)
    integer   :: lod(MAX_SPACEDIM), hid(MAX_SPACEDIM), loflux(res%dim), hiflux(res%dim)
    integer   :: i, j, ii, np

    real(kind=dp_t), pointer   :: rp(:,:,:,:),fp(:,:,:,:),cp(:,:,:,:),sp(:,:,:,:)
    integer,         pointer   :: mp(:,:,:,:),mcp(:,:,:,:)
    integer,         parameter :: tag = 1103

    real(kind=dp_t), dimension(:,:,:,:), allocatable :: flxpt
    integer,         dimension(:,:,:,:), allocatable :: mmfpt

    integer                 :: fsh(MAX_SPACEDIM+1), msh(MAX_SPACEDIM+1)
    integer,    allocatable :: rcnt(:), rdsp(:), scnt(:), sdsp(:)
    real(dp_t), allocatable :: g_snd_d(:), g_rcv_d(:)
    integer,    allocatable :: g_snd_i(:), g_rcv_i(:)
    type(fluxassoc)         :: fa

    np = parallel_nprocs()
    fa = layout_fluxassoc(crse%la, flux%la, crse%nodal, flux%nodal, side, crse_domain, ir)
    !
    ! Do all the local work.
    !
    do ii = 1, fa%flux%l_con%ncpy
       i     =  fa%flux%l_con%cpy(ii)%ns
       j     =  fa%flux%l_con%cpy(ii)%nd
       isect =  fa%flux%l_con%cpy(ii)%sbx
       lo    =  lwb(isect)
       hi    =  upb(isect)
       fbox  =  get_ibox(flux,i)
       lof   =  lwb(fbox)
       hif   =  upb(fbox)
       loc   =  lwb(get_pbox(crse,   j))
       lomc  =  lwb(get_pbox(mm_crse,j))
       lor   =  lwb(get_pbox(res,    j))
       los   =  lwb(get_pbox(ss,     j))
       lomf  =  lwb(get_pbox(mm_fine,i))
       fp    => dataptr(flux,        i)
       mp    => dataptr(mm_fine,     i)
       cp    => dataptr(crse,        j)
       rp    => dataptr(res,         j)
       sp    => dataptr(ss,          j)
       mcp   => dataptr(mm_crse,     j)

       select case (res%dim)
       case (1)
          call ml_interface_1d_nodal(rp(:,1,1,1), lor, &
               fp(:,1,1,1), lof, hif, cp(:,1,1,1), loc, &
               sp(:,1,1,:), los, lo, hi, ir, side, lof, hif)
       case (2)
          call ml_interface_2d_nodal(rp(:,:,1,1), lor, &
               fp(:,:,1,1), lof , hif, cp(:,:,1,1), loc , sp(:,:,1,:), los , &
               mp(:,:,1,1), lomf, mcp(:,:,1,1), lomc, lo, hi, ir, side, lof, hif)
       case (3)
          call ml_interface_3d_nodal(rp(:,:,:,1), lor, &
               fp(:,:,:,1), lof , hif, cp(:,:,:,1), loc , sp(:,:,:,:), los , &
               mp(:,:,:,1), lomf, mcp(:,:,:,1), lomc, lo, hi, ir, side, lof, hif)
       end select
    end do
    !
    ! Now send/recv the flux data
    !
    allocate(g_snd_d(fa%flux%r_con%svol))
    allocate(g_rcv_d(fa%flux%r_con%rvol))

    do i = 1, fa%flux%r_con%nsnd
       fp => dataptr(flux, fa%flux%r_con%snd(i)%ns, fa%flux%r_con%snd(i)%sbx)
       g_snd_d(1 + fa%flux%r_con%snd(i)%pv:fa%flux%r_con%snd(i)%av) = reshape(fp, fa%flux%r_con%snd(i)%s1)
    end do

    allocate(rcnt(0:np-1), rdsp(0:np-1), scnt(0:np-1), sdsp(0:np-1))

    rcnt = 0; scnt = 0; rdsp = 0; sdsp = 0

    do i = 1, fa%flux%r_con%nsp
       ii = fa%flux%r_con%str(i)%pr
       scnt(ii) = fa%flux%r_con%str(i)%sz
       sdsp(ii) = fa%flux%r_con%str(i)%pv
    end do
    do i = 1, fa%flux%r_con%nrp
       ii = fa%flux%r_con%rtr(i)%pr
       rcnt(ii) = fa%flux%r_con%rtr(i)%sz
       rdsp(ii) = fa%flux%r_con%rtr(i)%pv
    end do
    call parallel_alltoall(g_rcv_d, rcnt, rdsp, g_snd_d, scnt, sdsp)
    !
    ! Now send/recv mask data.
    !
    allocate(g_snd_i(fa%mask%r_con%svol))
    allocate(g_rcv_i(fa%mask%r_con%rvol))

    do i = 1, fa%mask%r_con%nsnd
       mp => dataptr(mm_fine, fa%mask%r_con%snd(i)%ns, fa%mask%r_con%snd(i)%sbx)
       g_snd_i(1 + fa%mask%r_con%snd(i)%pv:fa%mask%r_con%snd(i)%av) = reshape(mp, fa%mask%r_con%snd(i)%s1)
    end do

    rcnt = 0; scnt = 0; rdsp = 0; sdsp = 0

    do i = 1, fa%mask%r_con%nsp
       ii = fa%mask%r_con%str(i)%pr
       scnt(ii) = fa%mask%r_con%str(i)%sz
       sdsp(ii) = fa%mask%r_con%str(i)%pv
    end do
    do i = 1, fa%mask%r_con%nrp
       ii = fa%mask%r_con%rtr(i)%pr
       rcnt(ii) = fa%mask%r_con%rtr(i)%sz
       rdsp(ii) = fa%mask%r_con%rtr(i)%pv
    end do
    call parallel_alltoall(g_rcv_i, rcnt, rdsp, g_snd_i, scnt, sdsp)
    !
    ! Got all the remote data.  Use it.
    !
    do i = 1, fa%flux%r_con%nrcv
       j      =  fa%flux%r_con%rcv(i)%nd
       fsh    =  fa%flux%r_con%rcv(i)%sh
       msh    =  fa%mask%r_con%rcv(i)%sh
       isect  =  fa%flux%r_con%rcv(i)%sbx
       loflux =  lwb(fa%fbxs(i))
       hiflux =  upb(fa%fbxs(i))
       lof    =  lwb(isect)
       hif    =  upb(isect)
       lomf   =  lwb(fa%mask%r_con%rcv(i)%sbx)
       loc    =  lwb(get_pbox(crse,   j))
       lomc   =  lwb(get_pbox(mm_crse,j))
       lor    =  lwb(get_pbox(res,    j))
       los    =  lwb(get_pbox(ss,     j))
       cp     => dataptr(crse,        j)
       rp     => dataptr(res,         j)
       sp     => dataptr(ss,          j)
       mcp    => dataptr(mm_crse,     j)

       lod = 1;                     hid = 1
       lod(1:res%dim) = lwb(isect); hid(1:res%dim) = upb(isect)
       allocate(flxpt(lod(1):hid(1),lod(2):hid(2),lod(3):hid(3),1))
       flxpt = reshape(g_rcv_d(1 + fa%flux%r_con%rcv(i)%pv:fa%flux%r_con%rcv(i)%av), fsh)

       mbox = fa%mask%r_con%rcv(i)%sbx
       lod(1:res%dim) = lwb(mbox); hid(1:res%dim) = upb(mbox)
       allocate(mmfpt(lod(1):hid(1),lod(2):hid(2),lod(3):hid(3),1))
       mmfpt = reshape(g_rcv_i(1 + fa%mask%r_con%rcv(i)%pv:fa%mask%r_con%rcv(i)%av), msh)

       select case (res%dim)
       case (1)
          call ml_interface_1d_nodal(rp(:,1,1,1), lor, &
               flxpt(:,1,1,1), lof, hif, cp(:,1,1,1), loc, &
               sp(:,1,1,:), los, lof, hif, ir, side, loflux, hiflux)
       case (2)
          call ml_interface_2d_nodal(rp(:,:,1,1), lor, &
               flxpt(:,:,1,1), lof , hif, cp(:,:,1,1), loc , sp(:,:,1,:), los , &
               mmfpt(:,:,1,1), lomf, mcp(:,:,1,1), lomc, lof, hif, ir, side, loflux, hiflux)
       case (3)
          call ml_interface_3d_nodal(rp(:,:,:,1), lor, &
               flxpt(:,:,:,1), lof , hif, cp(:,:,:,1), loc , sp(:,:,:,:), los , &
               mmfpt(:,:,:,1), lomf, mcp(:,:,:,1), lomc, lof, hif, ir, side, loflux, hiflux)
       end select

       deallocate(flxpt)
       deallocate(mmfpt)
    end do

  end subroutine ml_crse_contrib_fancy

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
    type(bl_prof_timer), save     :: bpt
    call build(bpt, "ml_crse_contrib")
    call ml_crse_contrib_fancy(res, flux, crse, ss, mm_crse, mm_fine, crse_domain, ir, side)
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

    integer :: i, j, ioff, joff, sig_i, sig_j
    real (kind = dp_t) :: crse_flux

    i = lo(1)
    j = lo(2)

    !   NOTE: THESE STENCILS ONLY WORK FOR DX == DY.

    if (size(ss,dim=3) .eq. 9) then
    ! Dense stencil

      if (side == -1) then

       do j = lo(2),hi(2)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0) .and. &
               (.not. bc_dirichlet(mm_crse(i,j),1,0))) then

             if (j == loflux(2) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,-1)) then
                crse_flux = HALF * ss(i,j,8) * &
                 (cc(i+1,j+1) + HALF*cc(i+1,j) + HALF * cc(i,j+1) - TWO*cc(i,j))

             else if (j == hiflux(2) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,+1)) then
                crse_flux = HALF * ss(i,j,3) * &
                 (cc(i+1,j-1) + HALF*cc(i+1,j) + HALF * cc(i,j-1) - TWO*cc(i,j))

             else
                crse_flux = ss(i,j,8) * &
                 (cc(i+1,j+1) + HALF*cc(i+1,j) + HALF * cc(i,j+1) - TWO*cc(i,j)) &
                           +ss(i,j,3) * &
                 (cc(i+1,j-1) + HALF*cc(i+1,j) + HALF * cc(i,j-1) - TWO*cc(i,j))
             end if

             if (ir(1) .eq. 2) then
               crse_flux = crse_flux * 4.0_dp_t
             else if (ir(1) .eq. 4) then
               crse_flux = crse_flux * 16.0_dp_t
             end if

             res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
          end if
       end do

      else if (side ==  1) then

       do j = lo(2),hi(2)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0) .and. &
               (.not. bc_dirichlet(mm_crse(i,j),1,0))) then

             if (j == loflux(2) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,-1)) then
                crse_flux = HALF * ss(i,j,6) * &
                 (cc(i-1,j+1) + HALF*cc(i-1,j) + HALF * cc(i,j+1) - TWO*cc(i,j))

             else if (j == hiflux(2) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,+1)) then
                crse_flux = HALF * ss(i,j,1) * &
                 (cc(i-1,j-1) + HALF*cc(i-1,j) + HALF * cc(i,j-1) - TWO*cc(i,j))

             else
                crse_flux = ss(i,j,6) * &
                 (cc(i-1,j+1) + HALF*cc(i-1,j) + HALF * cc(i,j+1) - TWO*cc(i,j)) &
                           +ss(i,j,1) * &
                 (cc(i-1,j-1) + HALF*cc(i-1,j) + HALF * cc(i,j-1) - TWO*cc(i,j))
             end if

             if (ir(1) .eq. 2) then
               crse_flux = crse_flux * 4.0_dp_t
             else if (ir(1) .eq. 4) then
               crse_flux = crse_flux * 16.0_dp_t
             end if

             res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
          end if
       end do

      else if (side == -2) then

       do i = lo(1),hi(1)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0) .and. &
               (.not. bc_dirichlet(mm_crse(i,j),1,0))) then

             if (i == loflux(1) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,-1)) then
                crse_flux = HALF * ss(i,j,8) * &
                 (cc(i+1,j+1) + HALF*cc(i+1,j) + HALF * cc(i,j+1) - TWO*cc(i,j))

             else if (i == hiflux(1) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,+1)) then
                crse_flux = HALF * ss(i,j,6) * &
                 (cc(i-1,j+1) + HALF*cc(i-1,j) + HALF * cc(i,j+1) - TWO*cc(i,j))

             else
                crse_flux = ss(i,j,8) * & 
                 (cc(i+1,j+1) + HALF*cc(i+1,j) + HALF * cc(i,j+1) - TWO*cc(i,j)) &
                           +ss(i,j,6) * &
                 (cc(i-1,j+1) + HALF*cc(i-1,j) + HALF * cc(i,j+1) - TWO*cc(i,j))
             end if

             if (ir(2) .eq. 2) then
               crse_flux = crse_flux * 4.0_dp_t
             else if (ir(2) .eq. 4) then
               crse_flux = crse_flux * 16.0_dp_t
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

             if (ir(2) .eq. 2) then
               crse_flux = crse_flux * 4.0_dp_t
             else if (ir(2) .eq. 4) then
               crse_flux = crse_flux * 16.0_dp_t
             end if

             res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
          end if
       end do

      end if

    else if (size(ss,dim=3) .eq. 5) then
    ! Cross stencil

      if (side == -1 .or. side == 1) then

       if (side == -1) then
          ioff   = i+1
          sig_i  = 1
       else if (side == 1) then
          ioff   = i-1
          sig_i  = 2
       end if

       do j = lo(2),hi(2)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0) .and. &
               (.not. bc_dirichlet(mm_crse(i,j),1,0))) then

             if (j == loflux(2) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,-1)) then
                crse_flux = FOURTH*ss(i,j,sig_i) * (cc(ioff,j  )-cc(i,j)) &
                           +FOURTH*ss(i,j,3    ) * (cc(i   ,j+1)-cc(i,j)) 

             else if (j == hiflux(2) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,+1)) then
                crse_flux = FOURTH*ss(i,j,sig_i) * (cc(ioff,j  )-cc(i,j)) &
                           +FOURTH*ss(i,j,4    ) * (cc(i   ,j-1)-cc(i,j)) 

             else
                crse_flux =      ss(i,j,sig_i) * (cc(ioff,j  )-cc(i,j)) &
                           +HALF*ss(i,j,3    ) * (cc(i   ,j+1)-cc(i,j)) &
                           +HALF*ss(i,j,4    ) * (cc(i   ,j-1)-cc(i,j)) 
             end if

             if (ir(1) .eq. 2) then
               crse_flux = crse_flux * 4.0_dp_t
             else if (ir(1) .eq. 4) then
               crse_flux = crse_flux * 16.0_dp_t
             end if

             res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
          end if
       end do

      else if (side == -2 .or. side == 2) then

       if (side == -2) then
          joff   = j+1
          sig_j  = 3
       else if (side == 2) then
          joff   = j-1
          sig_j  = 4
       end if

       do i = lo(1),hi(1)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0) .and. &
               (.not. bc_dirichlet(mm_crse(i,j),1,0))) then

             if (i == loflux(1) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,-1)) then
                crse_flux = FOURTH*ss(i,j,sig_j) * (cc(i  ,joff)-cc(i,j)) &
                           +FOURTH*ss(i,j,1    ) * (cc(i+1,j   )-cc(i,j)) 

             else if (i == hiflux(1) .and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,+1)) then
                crse_flux = FOURTH*ss(i,j,sig_j) * (cc(i  ,joff)-cc(i,j)) &
                           +FOURTH*ss(i,j,2    ) * (cc(i-1,j   )-cc(i,j)) 

             else
                crse_flux =      ss(i,j,sig_j) * (cc(i  ,joff)-cc(i,j)) &
                           +HALF*ss(i,j,1    ) * (cc(i+1,j   )-cc(i,j)) &
                           +HALF*ss(i,j,2    ) * (cc(i-1,j   )-cc(i,j)) 
             end if
             if (ir(1) .eq. 2) then
               crse_flux = crse_flux * 4.0_dp_t
             else if (ir(1) .eq. 4) then
               crse_flux = crse_flux * 16.0_dp_t
             end if
             res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
          end if
       end do

      end if
    
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
    integer :: sig_i,sig_j,sig_k
    logical :: lo_i_not,lo_j_not,lo_k_not,hi_i_not,hi_j_not,hi_k_not
    real (kind = dp_t) :: crse_flux
    real (kind = dp_t) :: cell_mm, cell_mp, cell_pm, cell_pp

    i = lo(1)
    j = lo(2)
    k = lo(3)

    !   NOTE: THESE STENCILS ONLY WORK FOR DX == DY.
    !   NOTE: MM IS ON THE FINE GRID, NOT THE CRSE

    if ( (size(ss,dim=4) .eq. 27) .or. (size(ss,dim=4) .eq. 21) ) then
    ! Dense stencil

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

             if (ir(1) .eq. 2) then
               crse_flux = crse_flux * 8.0_dp_t
             else if (ir(1) .eq. 4) then
               crse_flux = crse_flux * 64.0_dp_t
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

             if (ir(2) .eq. 2) then
               crse_flux = crse_flux * 8.0_dp_t
             else if (ir(2) .eq. 4) then
               crse_flux = crse_flux * 64.0_dp_t
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

             if (ir(3) .eq. 2) then
               crse_flux = crse_flux * 8.0_dp_t
             else if (ir(3) .eq. 4) then
               crse_flux = crse_flux * 64.0_dp_t
             end if

             res(i,j,k) = res(i,j,k) + crse_flux + fine_flux(i,j,k)
           end if
          end do
       end do
    end if

    else if (size(ss,dim=4) .eq. 7) then
    ! Cross stencil

    !   Lo/Hi i side
    if (side == -1 .or. side == 1) then

       if (side == -1) then
          ioff   = i+1
          sig_i  = 1
       else if (side == 1) then
          ioff   = i-1
          sig_i  = 2
       end if

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)

           if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,0) .and. &
                (.not.bc_dirichlet(mm_crse(i,j,k),1,0))) then

             lo_j_not = ( (j == loflux(2)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),2,-1) )
             hi_j_not = ( (j == hiflux(2)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),2,+1) )
             lo_k_not = ( (k == loflux(3)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),3,-1) )
             hi_k_not = ( (k == hiflux(3)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),3,+1) )

             cell_mm = FOURTH*ss(i,j,k,sig_i)*(cc(ioff,j,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    4)*(cc(i,j-1,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    6)*(cc(i,j,k-1)-cc(i,j,k)) 

             cell_pm = FOURTH*ss(i,j,k,sig_i)*(cc(ioff,j,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    3)*(cc(i,j+1,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    6)*(cc(i,j,k-1)-cc(i,j,k)) 

             cell_mp = FOURTH*ss(i,j,k,sig_i)*(cc(ioff,j,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    4)*(cc(i,j-1,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    5)*(cc(i,j,k+1)-cc(i,j,k)) 

             cell_pp = FOURTH*ss(i,j,k,sig_i)*(cc(ioff,j,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    3)*(cc(i,j+1,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    5)*(cc(i,j,k+1)-cc(i,j,k)) 

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

             if (ir(1) .eq. 2) then
               crse_flux = crse_flux * 8.0_dp_t
             else if (ir(1) .eq. 4) then
               crse_flux = crse_flux * 64.0_dp_t
             end if

             res(i,j,k) = res(i,j,k) + crse_flux + fine_flux(i,j,k)

           end if
          end do
       end do

       !   Lo/Hi j side
    else if (side == -2 .or. side == 2) then

       if (side == -2) then
          joff   = j+1
          sig_j  = 3
       else
          joff   = j-1
          sig_j  = 4
       end if
       do k = lo(3),hi(3)
          do i = lo(1),hi(1)

           if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,0) .and. &
                (.not.bc_dirichlet(mm_crse(i,j,k),1,0))) then

             lo_i_not = ( (i == loflux(1)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,-1) )
             hi_i_not = ( (i == hiflux(1)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,+1) )
             lo_k_not = ( (k == loflux(3)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),3,-1) )
             hi_k_not = ( (k == hiflux(3)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),3,+1) )

             cell_mm = FOURTH*ss(i,j,k,sig_j)*(cc(i,joff,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    2)*(cc(i-1,j,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    6)*(cc(i,j,k-1)-cc(i,j,k)) 

             cell_pm = FOURTH*ss(i,j,k,sig_j)*(cc(i,joff,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    1)*(cc(i+1,j,k)-cc(i,j,k))& 
                      +FOURTH*ss(i,j,k,    6)*(cc(i,j,k-1)-cc(i,j,k)) 

             cell_mp = FOURTH*ss(i,j,k,sig_j)*(cc(i,joff,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    2)*(cc(i-1,j,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    5)*(cc(i,j,k+1)-cc(i,j,k)) 

             cell_pp = FOURTH*ss(i,j,k,sig_j)*(cc(i,joff,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    1)*(cc(i+1,j,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    5)*(cc(i,j,k+1)-cc(i,j,k)) 

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

             if (ir(2) .eq. 2) then
               crse_flux = crse_flux * 8.0_dp_t
             else if (ir(2) .eq. 4) then
               crse_flux = crse_flux * 64.0_dp_t
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
          sig_k  = 5
       else
          koff   = k-1
          sig_k  = 6
       end if

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

           if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,0) .and. &
                (.not.bc_dirichlet(mm_crse(i,j,k),1,0))) then

             lo_i_not = ( (i == loflux(1)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,-1) )
             hi_i_not = ( (i == hiflux(1)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,+1) )
             lo_j_not = ( (j == loflux(2)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),2,-1) )
             hi_j_not = ( (j == hiflux(2)).and. .not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),2,+1) )

             cell_mm = FOURTH*ss(i,j,k,sig_k)*(cc(i,j,koff)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    2)*(cc(i-1,j,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    4)*(cc(i,j-1,k)-cc(i,j,k)) 

             cell_pm = FOURTH*ss(i,j,k,sig_k)*(cc(i,j,koff)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    1)*(cc(i+1,j,k)-cc(i,j,k))& 
                      +FOURTH*ss(i,j,k,    4)*(cc(i,j-1,k)-cc(i,j,k)) 

             cell_mp = FOURTH*ss(i,j,k,sig_k)*(cc(i,j,koff)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    2)*(cc(i-1,j,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    3)*(cc(i,j+1,k)-cc(i,j,k)) 

             cell_pp = FOURTH*ss(i,j,k,sig_k)*(cc(i,j,koff)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    1)*(cc(i+1,j,k)-cc(i,j,k)) &
                      +FOURTH*ss(i,j,k,    3)*(cc(i,j+1,k)-cc(i,j,k)) 

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

             if (ir(3) .eq. 2) then
               crse_flux = crse_flux * 8.0_dp_t
             else if (ir(3) .eq. 4) then
               crse_flux = crse_flux * 64.0_dp_t
             end if

             res(i,j,k) = res(i,j,k) + crse_flux + fine_flux(i,j,k)
           end if
          end do
       end do
    end if

    end if

  end subroutine ml_interface_3d_nodal

end module ml_interface_stencil_module
