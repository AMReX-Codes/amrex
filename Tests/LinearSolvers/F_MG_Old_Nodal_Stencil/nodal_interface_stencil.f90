module nodal_interface_stencil_module

  use bl_types
  use layout_module
  use multifab_module
  use bc_functions_module
  use bl_constants_module

  implicit none

  private

  public :: ml_crse_contrib, ml_fine_contrib

contains

  subroutine ml_crse_contrib_fancy(res, flux, crse, ss, mm_crse, mm_fine, crse_domain, ir, side)

    type(multifab), intent(inout) :: res
    type(multifab), intent(inout) :: flux
    type(multifab), intent(in   ) :: crse
    type(multifab), intent(in   ) :: ss
    type(imultifab),intent(in   ) :: mm_crse
    type(imultifab),intent(in   ) :: mm_fine
    type(box),      intent(in   ) :: crse_domain
    integer,        intent(in   ) :: ir(:)
    integer                       :: side

    type(box) :: fbox, mbox, isect
    integer   :: lo (get_dim(res)), hi (get_dim(res)), loc(get_dim(res))
    integer   :: lof(get_dim(res)), hif(get_dim(res)), lor(get_dim(res)), los(get_dim(res))
    integer   :: lomf(get_dim(res)), lomc(get_dim(res))
    integer   :: lod(MAX_SPACEDIM), hid(MAX_SPACEDIM), loflux(get_dim(res)), hiflux(get_dim(res))
    integer   :: i, j, ii, np, av, dm

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

    dm = get_dim(res)
    np = parallel_nprocs()
    fa = layout_fluxassoc(crse%la, flux%la, nodal_flags(crse), nodal_flags(flux), side, crse_domain, ir)
    !
    ! Do all the local work.
    !
    do ii = 1, fa%flux%l_con%ncpy
       i     =  fa%flux%l_con%cpy(ii)%lns
       j     =  fa%flux%l_con%cpy(ii)%lnd
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

       select case (dm)
       case (1)
          call ml_interface_1d_nodal(rp(:,1,1,1), lor, &
               fp(:,1,1,1), lof, cp(:,1,1,1), loc, &
               sp(:,:,1,1), los, lo, side)
       case (2)
          call ml_interface_2d_nodal(rp(:,:,1,1), lor, &
               fp(:,:,1,1), lof , cp(:,:,1,1), loc , sp(:,:,:,1), los , &
               mp(:,:,1,1), lomf, mcp(:,:,1,1), lomc, lo, hi, ir, side, lof, hif)
       case (3)
          call ml_interface_3d_nodal(rp(:,:,:,1), lor, &
               fp(:,:,:,1), lof , cp(:,:,:,1), loc , sp(:,:,:,:), los , &
               mp(:,:,:,1), lomf, mcp(:,:,:,1), lomc, lo, hi, ir, side, lof, hif)
       end select
    end do

    if (np == 1) return
    !
    ! Now send/recv the flux data
    !
    allocate(g_snd_d(fa%flux%r_con%svol))
    allocate(g_rcv_d(fa%flux%r_con%rvol))

    do i = 1, fa%flux%r_con%nsnd
       fp => dataptr(flux, fa%flux%r_con%snd(i)%lns, fa%flux%r_con%snd(i)%sbx)
       call reshape_d_4_1(g_snd_d, 1 + fa%flux%r_con%snd(i)%pv, fp)
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
       mp => dataptr(mm_fine, flux,fa%mask%r_con%snd(i)%lns, fa%mask%r_con%snd(i)%sbx)
       call reshape_i_4_1(g_snd_i, 1 + fa%mask%r_con%snd(i)%pv, mp)
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
       j      =  fa%flux%r_con%rcv(i)%lnd
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
       lod(1:dm) = lwb(isect); hid(1:dm) = upb(isect)
       allocate(flxpt(lod(1):hid(1),lod(2):hid(2),lod(3):hid(3),1))
       call reshape_d_1_4(flxpt, g_rcv_d, 1 + fa%flux%r_con%rcv(i)%pv, fsh)

       mbox = fa%mask%r_con%rcv(i)%sbx
       lod(1:dm) = lwb(mbox); hid(1:dm) = upb(mbox)
       allocate(mmfpt(lod(1):hid(1),lod(2):hid(2),lod(3):hid(3),1))
       av = fa%mask%r_con%rcv(i)%pv + volume(mbox)
       call reshape_i_1_4(mmfpt, g_rcv_i, 1 + fa%mask%r_con%rcv(i)%pv, msh)

       select case (dm)
       case (1)
          call ml_interface_1d_nodal(rp(:,1,1,1), lor, &
               flxpt(:,1,1,1), lof, cp(:,1,1,1), loc, &
               sp(:,:,1,1), los, lof, side)
       case (2)
          call ml_interface_2d_nodal(rp(:,:,1,1), lor, &
               flxpt(:,:,1,1), lof , cp(:,:,1,1), loc , sp(:,:,:,1), los , &
               mmfpt(:,:,1,1), lomf, mcp(:,:,1,1), lomc, lof, hif, ir, side, loflux, hiflux)
       case (3)
          call ml_interface_3d_nodal(rp(:,:,:,1), lor, &
               flxpt(:,:,:,1), lof , cp(:,:,:,1), loc , sp(:,:,:,:), los , &
               mmfpt(:,:,:,1), lomf, mcp(:,:,:,1), lomc, lof, hif, ir, side, loflux, hiflux)
       end select

       deallocate(flxpt)
       deallocate(mmfpt)
    end do

  end subroutine ml_crse_contrib_fancy

  subroutine ml_crse_contrib(res, flux, crse, ss, mm_crse, mm_fine, crse_domain, ir, side)
    use bl_prof_module
    type(multifab), intent(inout) :: res
    type(multifab), intent(inout) :: flux
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

  subroutine ml_interface_1d_nodal(res, lor, fine_flux, lof, cc, loc, ss , los, lo, side)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lof(:)
    integer, intent(in) :: lo(:)
    real (kind = dp_t), intent(inout) :: res(lor(1):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):)
    real (kind = dp_t), intent(in   ) :: cc(loc(1):)
    real (kind = dp_t), intent(in   ) :: ss(0:,los(1):)
    integer, intent(in) :: side

    integer :: i
    real (kind = dp_t) :: crse_flux

    i = lo(1)

    !   Lo side
    if (side == -1) then
       crse_flux = ss(1,i)*(cc(i+1)-cc(i))
       res(i) = res(i) - fine_flux(i) + crse_flux
       !   Hi side
    else if (side == 1) then
       crse_flux = ss(2,i)*(cc(i-1)-cc(i))
       res(i) = res(i) - fine_flux(i) + crse_flux
    end if

  end subroutine ml_interface_1d_nodal

  subroutine ml_interface_2d_nodal(res, lor, fine_flux, lof, cc, loc, &
       ss , los, mm_fine, lomf, mm_crse, lomc, lo, hi, ir, side, loflux, hiflux)
    use cc_stencil_module
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lomf(:)
    integer, intent(in) :: lomc(:)
    integer, intent(in) :: lof(:), loflux(:), hiflux(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):)
    real (kind = dp_t), intent(in   ) ::        cc(loc(1):,loc(2):)
    real (kind = dp_t), intent(in   ) ::        ss(0:,los(1):,los(2):)
    integer           , intent(in   ) ::   mm_fine(lomf(1):,lomf(2):)
    integer           , intent(in   ) ::   mm_crse(lomc(1):,lomc(2):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)

    integer :: i, j, ioff, joff, sig_i, sig_j
    real (kind = dp_t) :: crse_flux
    logical llo,lhi

    i = lo(1)
    j = lo(2)

    ioff = 0; joff = 0; sig_i = 0; sig_j = 0

    !   NOTE: THESE STENCILS ONLY WORK FOR DX == DY.

    if (size(ss,dim=1) .eq. 9) then
    ! Dense stencil

      if (side == -1) then

       do j = lo(2),hi(2)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0)) then
             if (.not. bc_dirichlet(mm_crse(i,j),1,0)) then

                llo = .false.
                lhi = .false.

                if (j == loflux(2)) then
                   if (.not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,-1)) llo = .true.
                end if

                if (j == hiflux(2)) then
                   if (.not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,+1)) lhi = .true.
                end if
                
                if (llo) then
                   crse_flux = HALF * ss(8,i,j) * &
                        (cc(i+1,j+1) + HALF*cc(i+1,j) + HALF * cc(i,j+1) - TWO*cc(i,j))
                else if (lhi) then
                   crse_flux = HALF * ss(3,i,j) * &
                        (cc(i+1,j-1) + HALF*cc(i+1,j) + HALF * cc(i,j-1) - TWO*cc(i,j))
                else
                   crse_flux = ss(8,i,j) * &
                        (cc(i+1,j+1) + HALF*cc(i+1,j) + HALF * cc(i,j+1) - TWO*cc(i,j)) &
                        +ss(3,i,j) * &
                        (cc(i+1,j-1) + HALF*cc(i+1,j) + HALF * cc(i,j-1) - TWO*cc(i,j))
                end if

                if (ir(1) .eq. 2) then
                   crse_flux = crse_flux * 4.0_dp_t
                else if (ir(1) .eq. 4) then
                   crse_flux = crse_flux * 16.0_dp_t
                end if

                res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
             end if
          end if
       end do

      else if (side ==  1) then

       do j = lo(2),hi(2)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0)) then
             if (.not. bc_dirichlet(mm_crse(i,j),1,0)) then

                llo = .false.
                lhi = .false.

                if (j == loflux(2)) then
                   if (.not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,-1)) llo = .true.
                end if

                if (j == hiflux(2)) then
                   if (.not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,+1)) lhi = .true.
                end if

                if (llo) then
                   crse_flux = HALF * ss(6,i,j) * &
                        (cc(i-1,j+1) + HALF*cc(i-1,j) + HALF * cc(i,j+1) - TWO*cc(i,j))
                else if (lhi) then
                   crse_flux = HALF * ss(1,i,j) * &
                        (cc(i-1,j-1) + HALF*cc(i-1,j) + HALF * cc(i,j-1) - TWO*cc(i,j))
                else
                   crse_flux = ss(6,i,j) * &
                        (cc(i-1,j+1) + HALF*cc(i-1,j) + HALF * cc(i,j+1) - TWO*cc(i,j)) &
                        +ss(1,i,j) * &
                        (cc(i-1,j-1) + HALF*cc(i-1,j) + HALF * cc(i,j-1) - TWO*cc(i,j))
                end if

                if (ir(1) .eq. 2) then
                   crse_flux = crse_flux * 4.0_dp_t
                else if (ir(1) .eq. 4) then
                   crse_flux = crse_flux * 16.0_dp_t
                end if

                res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
             end if
          end if
       end do

      else if (side == -2) then

       do i = lo(1),hi(1)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0)) then
             if (.not. bc_dirichlet(mm_crse(i,j),1,0)) then

                llo = .false.
                lhi = .false.

                if (i == loflux(1)) then
                   if (.not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,-1)) llo = .true.
                end if

                if (i == hiflux(1)) then
                   if (.not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,+1)) lhi = .true.
                end if

                if (llo) then
                   crse_flux = HALF * ss(8,i,j) * &
                        (cc(i+1,j+1) + HALF*cc(i+1,j) + HALF * cc(i,j+1) - TWO*cc(i,j))
                else if (lhi) then
                   crse_flux = HALF * ss(6,i,j) * &
                        (cc(i-1,j+1) + HALF*cc(i-1,j) + HALF * cc(i,j+1) - TWO*cc(i,j))
                else
                   crse_flux = ss(8,i,j) * & 
                        (cc(i+1,j+1) + HALF*cc(i+1,j) + HALF * cc(i,j+1) - TWO*cc(i,j)) &
                        +ss(6,i,j) * &
                        (cc(i-1,j+1) + HALF*cc(i-1,j) + HALF * cc(i,j+1) - TWO*cc(i,j))
                end if

                if (ir(2) .eq. 2) then
                   crse_flux = crse_flux * 4.0_dp_t
                else if (ir(2) .eq. 4) then
                   crse_flux = crse_flux * 16.0_dp_t
                end if

                res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
             end if
          end if
       end do

    else if (side ==  2) then

       do i = lo(1),hi(1)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0)) then
             if (.not. bc_dirichlet(mm_crse(i,j),1,0)) then

                llo = .false.
                lhi = .false.

                if (i == loflux(1)) then
                   if (.not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,-1)) llo = .true.
                end if

                if (i == hiflux(1)) then
                   if (.not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,+1)) lhi = .true.
                end if

                if (llo) then
                   crse_flux = (ss(3,i,j)*(cc(i+1,j-1) + HALF*cc(i+1,j) + &
                        HALF * cc(i,j-1) - TWO*cc(i,j)) ) * HALF
                else if (lhi) then
                   crse_flux = (ss(1,i,j)*(cc(i-1,j-1) + HALF*cc(i-1,j) + &
                        HALF * cc(i,j-1) - TWO*cc(i,j)) ) * HALF
                else
                   crse_flux = ss(3,i,j)*(cc(i+1,j-1) + HALF*cc(i+1,j) + &
                        HALF * cc(i,j-1) - TWO*cc(i,j)) &
                        +ss(1,i,j)*(cc(i-1,j-1) + HALF*cc(i-1,j) + &
                        HALF * cc(i,j-1) - TWO*cc(i,j))
                end if

                if (ir(2) .eq. 2) then
                   crse_flux = crse_flux * 4.0_dp_t
                else if (ir(2) .eq. 4) then
                   crse_flux = crse_flux * 16.0_dp_t
                end if

                res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
             end if
          end if
       end do

      end if

   else if (size(ss,dim=1) .eq. 5) then
      !
      ! Cross stencil
      !
      if (side == -1 .or. side == 1) then

       if (side == -1) then
          ioff   = i+1
          sig_i  = 1
       else if (side == 1) then
          ioff   = i-1
          sig_i  = 2
       end if

       do j = lo(2),hi(2)
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0)) then
             if (.not. bc_dirichlet(mm_crse(i,j),1,0)) then

                llo = .false.
                lhi = .false.

                if (j == loflux(2)) then
                   if (.not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,-1)) llo = .true.
                   end if

                if (j == hiflux(2)) then
                   if (.not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),2,+1)) lhi = .true.
                end if

                if (llo) then
                   crse_flux = FOURTH*ss(sig_i,i,j) * (cc(ioff,j  )-cc(i,j)) &
                        +FOURTH*ss(3,i,j) * (cc(i   ,j+1)-cc(i,j)) 
                else if (lhi) then
                   crse_flux = FOURTH*ss(sig_i,i,j) * (cc(ioff,j  )-cc(i,j)) &
                        +FOURTH*ss(4,i,j) * (cc(i   ,j-1)-cc(i,j)) 
                else
                   crse_flux =      ss(sig_i,i,j) * (cc(ioff,j  )-cc(i,j)) &
                        +HALF*ss(3,i,j) * (cc(i   ,j+1)-cc(i,j)) &
                        +HALF*ss(4,i,j) * (cc(i   ,j-1)-cc(i,j)) 
                end if

                if (ir(1) .eq. 2) then
                   crse_flux = crse_flux * 4.0_dp_t
                else if (ir(1) .eq. 4) then
                   crse_flux = crse_flux * 16.0_dp_t
                end if

                res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
             end if
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
          if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j),1,0)) then
             if (.not. bc_dirichlet(mm_crse(i,j),1,0)) then

                llo = .false.
                lhi = .false.

                if (i == loflux(1)) then
                   if (.not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,-1)) llo = .true.
                end if

                if (i == hiflux(1)) then
                   if (.not. bc_neumann(mm_fine(ir(1)*i,ir(2)*j),1,+1)) lhi = .true.
                end if

                if (llo) then
                   crse_flux = FOURTH*ss(sig_j,i,j) * (cc(i  ,joff)-cc(i,j)) &
                        +FOURTH*ss(1,i,j) * (cc(i+1,j   )-cc(i,j)) 
                else if (lhi) then
                   crse_flux = FOURTH*ss(sig_j,i,j) * (cc(i  ,joff)-cc(i,j)) &
                        +FOURTH*ss(2,i,j) * (cc(i-1,j   )-cc(i,j)) 
                else
                   crse_flux =      ss(sig_j,i,j) * (cc(i  ,joff)-cc(i,j)) &
                        +HALF*ss(1,i,j) * (cc(i+1,j   )-cc(i,j)) &
                        +HALF*ss(2,i,j) * (cc(i-1,j   )-cc(i,j)) 
                end if
                if (ir(1) .eq. 2) then
                   crse_flux = crse_flux * 4.0_dp_t
                else if (ir(1) .eq. 4) then
                   crse_flux = crse_flux * 16.0_dp_t
                end if
                res(i,j) = res(i,j) + crse_flux + fine_flux(i,j)
             end if
          end if
       end do

      end if
    
    end if

  end subroutine ml_interface_2d_nodal

  subroutine ml_interface_3d_nodal(res, lor, fine_flux, lof, cc, loc, &
       ss, los, mm_fine, lomf, mm_crse, lomc, lo, hi, ir, side, loflux, hiflux)

    use nodal_stencil_module

    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: los(:)
    integer, intent(in) :: lof(:), loflux(:), hiflux(:)
    integer, intent(in) :: lomf(:)
    integer, intent(in) :: lomc(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::       res(lor(1):,lor(2):,lor(3):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):,lof(3):)
    real (kind = dp_t), intent(in   ) ::        cc(loc(1):,loc(2):,loc(3):)
    real (kind = dp_t), intent(in   ) ::        ss(0:,los(1):,los(2):,los(3):)
    integer           , intent(in   ) ::   mm_fine(lomf(1):,lomf(2):,lomf(3):)
    integer           , intent(in   ) ::   mm_crse(lomc(1):,lomc(2):,lomc(3):)
    integer, intent(in) :: side
    integer, intent(in) :: ir(:)

    integer :: i, j, k, ifine, jfine, kfine
    integer :: ioff, joff, koff
    integer :: sig_mm, sig_mp, sig_pm, sig_pp
    integer :: sig_i,sig_j,sig_k
    logical :: lo_i_not,lo_j_not,lo_k_not,hi_i_not,hi_j_not,hi_k_not
    real (kind = dp_t) :: crse_flux
    real (kind = dp_t) :: cell_mm, cell_mp, cell_pm, cell_pp

    sig_i = 0; ioff = 0
    !
    !   NOTE: THESE STENCILS ONLY WORK FOR DX == DY.
    !   NOTE: MM IS ON THE FINE GRID, NOT THE CRSE
    !
    if ( (size(ss,dim=1) .eq. 27) .or. (size(ss,dim=1) .eq. 21) ) then
       !
       ! Dense stencil
       !
       !   Lo/Hi i side
       !
       if (side == -1 .or. side == 1) then

          i     = lo(1)
          ifine = ir(1)*i

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

          !$OMP PARALLEL DO PRIVATE(j,k,lo_j_not,hi_j_not,lo_k_not,hi_k_not) &
          !$OMP PRIVATE(cell_mm,cell_pm,cell_mp,cell_pp,crse_flux,jfine,kfine)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                kfine = ir(3)*k
                jfine = ir(2)*j

                if (bc_dirichlet(mm_fine(ir(1)*i,ir(2)*j,ir(3)*k),1,0) .and. (.not.bc_dirichlet(mm_crse(i,j,k),1,0))) then

                   lo_j_not = .false.; hi_j_not = .false.; lo_k_not = .false.; hi_k_not = .false.

                   if (j == loflux(2) .and. (.not. bc_neumann(mm_fine(ifine,jfine,kfine),2,-1))) lo_j_not = .true.
                   if (j == hiflux(2) .and. (.not. bc_neumann(mm_fine(ifine,jfine,kfine),2,+1))) hi_j_not = .true.
                   if (k == loflux(3) .and. (.not. bc_neumann(mm_fine(ifine,jfine,kfine),3,-1))) lo_k_not = .true.
                   if (k == hiflux(3) .and. (.not. bc_neumann(mm_fine(ifine,jfine,kfine),3,+1))) hi_k_not = .true.

                   cell_mm = ss(sig_mm,i,j,k)*(cc(ioff,j-1,k-1) + cc(ioff,j-1,k  ) &
                        +cc(ioff,j  ,k-1) + cc(i  ,j-1,k-1) - FOUR*cc(i  ,j  ,k) )
                   cell_pm = ss(sig_pm,i,j,k)*(cc(ioff,j+1,k-1) + cc(ioff,j+1,k  ) &
                        +cc(ioff,j  ,k-1) + cc(i  ,j+1,k-1) - FOUR*cc(i  ,j  ,k) )
                   cell_mp = ss(sig_mp,i,j,k)*(cc(ioff,j-1,k+1) + cc(ioff,j-1,k  ) &
                        +cc(ioff,j  ,k+1) + cc(i  ,j-1,k+1) - FOUR*cc(i  ,j  ,k) )
                   cell_pp = ss(sig_pp,i,j,k)*(cc(ioff,j+1,k+1) + cc(ioff,j+1,k  ) &
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
          !$OMP END PARALLEL DO
          !
          !   Lo/Hi j side
          !
       else if (side == -2 .or. side == 2) then

          j     = lo(2)
          jfine = ir(2)*j

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

          !$OMP PARALLEL DO PRIVATE(i,k,lo_i_not,hi_i_not,lo_k_not,hi_k_not) &
          !$OMP PRIVATE(cell_mm,cell_pm,cell_mp,cell_pp,crse_flux,ifine,kfine)
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                kfine = ir(3)*k
                ifine = ir(1)*i

                if (bc_dirichlet(mm_fine(ifine,jfine,kfine),1,0) .and. (.not.bc_dirichlet(mm_crse(i,j,k),1,0))) then

                   lo_i_not = .false.; hi_i_not = .false.; lo_k_not = .false.; hi_k_not = .false.

                   if (i == loflux(1) .and. (.not. bc_neumann(mm_fine(ifine,jfine,kfine),1,-1))) lo_i_not = .true.
                   if (i == hiflux(1) .and. (.not. bc_neumann(mm_fine(ifine,jfine,kfine),1,+1))) hi_i_not = .true.
                   if (k == loflux(3) .and. (.not. bc_neumann(mm_fine(ifine,jfine,kfine),3,-1))) lo_k_not = .true.
                   if (k == hiflux(3) .and. (.not. bc_neumann(mm_fine(ifine,jfine,kfine),3,+1))) hi_k_not = .true.

                   cell_mm = ss(sig_mm,i,j,k)*(cc(i-1,joff,k-1) + cc(i-1,joff,k  ) &
                        +cc(i  ,joff,k-1) + cc(i-1,j   ,k-1) - FOUR*cc(i  ,j  ,k) )
                   cell_pm = ss(sig_pm,i,j,k)*(cc(i+1,joff,k-1) + cc(i+1,joff,k  ) &
                        +cc(i  ,joff,k-1) + cc(i+1,j   ,k-1) - FOUR*cc(i  ,j  ,k) )
                   cell_mp = ss(sig_mp,i,j,k)*(cc(i-1,joff,k+1) + cc(i-1,joff,k  ) &
                        +cc(i  ,joff,k+1) + cc(i-1,j   ,k+1) - FOUR*cc(i  ,j  ,k) )
                   cell_pp = ss(sig_pp,i,j,k)*(cc(i+1,joff,k+1) + cc(i+1,joff,k  ) &
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
          !$OMP END PARALLEL DO
          !
          !   Lo/Hi k side
          !
       else if (side == -3 .or. side == 3) then

          k     = lo(3)
          kfine = ir(3)*k

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

          !$OMP PARALLEL DO PRIVATE(i,j,lo_i_not,hi_i_not,lo_j_not,hi_j_not) &
          !$OMP PRIVATE(cell_mm,cell_pm,cell_mp,cell_pp,crse_flux,ifine,jfine)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                jfine = ir(2)*j
                ifine = ir(1)*i

                if (bc_dirichlet(mm_fine(ifine,jfine,kfine),1,0) .and. (.not.bc_dirichlet(mm_crse(i,j,k),1,0))) then

                   lo_i_not = .false.; hi_i_not = .false.; lo_j_not = .false.; hi_j_not = .false.

                   if (i == loflux(1) .and. (.not. bc_neumann(mm_fine(ifine,jfine,kfine),1,-1))) lo_i_not = .true.
                   if (i == hiflux(1) .and. (.not. bc_neumann(mm_fine(ifine,jfine,kfine),1,+1))) hi_i_not = .true.
                   if (j == loflux(2) .and. (.not. bc_neumann(mm_fine(ifine,jfine,kfine),2,-1))) lo_j_not = .true.
                   if (j == hiflux(2) .and. (.not. bc_neumann(mm_fine(ifine,jfine,kfine),2,+1))) hi_j_not = .true.

                   cell_mm = ss(sig_mm,i,j,k)*(cc(i-1,j-1,koff) + cc(i-1,j  ,koff) &
                        +cc(i  ,j-1,koff) + cc(i-1,j-1,k   ) - FOUR*cc(i  ,j  ,k) )
                   cell_pm = ss(sig_pm,i,j,k)*(cc(i+1,j-1,koff) + cc(i+1,j  ,koff) &
                        +cc(i  ,j-1,koff) + cc(i+1,j-1,k   ) - FOUR*cc(i  ,j  ,k) )
                   cell_mp = ss(sig_mp,i,j,k)*(cc(i-1,j+1,koff) + cc(i-1,j  ,koff) &
                        +cc(i  ,j+1,koff) + cc(i-1,j+1,k   ) - FOUR*cc(i  ,j  ,k) )
                   cell_pp = ss(sig_pp,i,j,k)*(cc(i+1,j+1,koff) + cc(i+1,j  ,koff) &
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
          !$OMP END PARALLEL DO
       end if

    else if (size(ss,dim=1) .eq. 7) then
       !
       ! Cross stencil
       !
       !   Lo/Hi i side
       !
       if (side == -1 .or. side == 1) then

          i     = lo(1)
          ifine = ir(1)*i

          if (side == -1) then
             ioff   = i+1
             sig_i  = 1
          else if (side == 1) then
             ioff   = i-1
             sig_i  = 2
          end if

          !$OMP PARALLEL DO PRIVATE(j,k,lo_j_not,hi_j_not,lo_k_not,hi_k_not) &
          !$OMP PRIVATE(cell_mm,cell_pm,cell_mp,cell_pp,crse_flux,jfine,kfine)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)

                kfine = ir(3)*k
                jfine = ir(2)*j

                if (bc_dirichlet(mm_fine(ifine,jfine,kfine),1,0) .and. (.not. bc_dirichlet(mm_crse(i,j,k),1,0))) then

                   lo_j_not = .false.; hi_j_not = .false.; lo_k_not = .false.; hi_k_not = .false.

                   if (j == loflux(2)) then
                      if (.not. bc_neumann(mm_fine(ifine,jfine,kfine),2,-1)) lo_j_not = .true.
                   end if
                   if (j == hiflux(2)) then
                      if (.not. bc_neumann(mm_fine(ifine,jfine,kfine),2,+1)) hi_j_not = .true.
                   end if
                   if (k == loflux(3)) then
                      if (.not. bc_neumann(mm_fine(ifine,jfine,kfine),3,-1)) lo_k_not = .true.
                   end if
                   if (k == hiflux(3)) then
                      if (.not. bc_neumann(mm_fine(ifine,jfine,kfine),3,+1)) hi_k_not = .true.
                   end if

                   cell_mm = FOURTH*ss(sig_i,i,j,k)*(cc(ioff,j,k)-cc(i,j,k)) &
                        +FOURTH*ss(4,i,j,k)*(cc(i,j-1,k)-cc(i,j,k)) &
                        +FOURTH*ss(6,i,j,k)*(cc(i,j,k-1)-cc(i,j,k)) 

                   cell_pm = FOURTH*ss(sig_i,i,j,k)*(cc(ioff,j,k)-cc(i,j,k)) &
                        +FOURTH*ss(3,i,j,k)*(cc(i,j+1,k)-cc(i,j,k)) &
                        +FOURTH*ss(6,i,j,k)*(cc(i,j,k-1)-cc(i,j,k)) 

                   cell_mp = FOURTH*ss(sig_i,i,j,k)*(cc(ioff,j,k)-cc(i,j,k)) &
                        +FOURTH*ss(4,i,j,k)*(cc(i,j-1,k)-cc(i,j,k)) &
                        +FOURTH*ss(5,i,j,k)*(cc(i,j,k+1)-cc(i,j,k)) 

                   cell_pp = FOURTH*ss(sig_i,i,j,k)*(cc(ioff,j,k)-cc(i,j,k)) &
                        +FOURTH*ss(3,i,j,k)*(cc(i,j+1,k)-cc(i,j,k)) &
                        +FOURTH*ss(5,i,j,k)*(cc(i,j,k+1)-cc(i,j,k)) 

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
          !$OMP END PARALLEL DO
          !
          !   Lo/Hi j side
          !
       else if (side == -2 .or. side == 2) then

          j     = lo(2)
          jfine = ir(2)*j

          if (side == -2) then
             joff   = j+1
             sig_j  = 3
          else
             joff   = j-1
             sig_j  = 4
          end if

          !$OMP PARALLEL DO PRIVATE(i,k,lo_i_not,hi_i_not,lo_k_not,hi_k_not) &
          !$OMP PRIVATE(cell_mm,cell_pm,cell_mp,cell_pp,crse_flux,ifine,kfine)
          do k = lo(3),hi(3)
             do i = lo(1),hi(1)

                kfine = ir(3)*k
                ifine = ir(1)*i

                if (bc_dirichlet(mm_fine(ifine,jfine,kfine),1,0) .and. (.not.bc_dirichlet(mm_crse(i,j,k),1,0))) then

                   lo_i_not = .false.; hi_i_not = .false.; lo_k_not = .false.; hi_k_not = .false.

                   if (i == loflux(1)) then
                      if (.not. bc_neumann(mm_fine(ifine,jfine,kfine),1,-1)) lo_i_not = .true.
                   end if
                   if (i == hiflux(1)) then
                      if (.not. bc_neumann(mm_fine(ifine,jfine,kfine),1,+1)) hi_i_not = .true.
                   end if
                   if (k == loflux(3)) then
                      if (.not. bc_neumann(mm_fine(ifine,jfine,kfine),3,-1)) lo_k_not = .true.
                   end if
                   if (k == hiflux(3)) then
                      if (.not. bc_neumann(mm_fine(ifine,jfine,kfine),3,+1)) hi_k_not = .true.
                   end if

                   cell_mm = FOURTH*ss(sig_j,i,j,k)*(cc(i,joff,k)-cc(i,j,k)) &
                        +FOURTH*ss(2,i,j,k)*(cc(i-1,j,k)-cc(i,j,k)) &
                        +FOURTH*ss(6,i,j,k)*(cc(i,j,k-1)-cc(i,j,k)) 

                   cell_pm = FOURTH*ss(sig_j,i,j,k)*(cc(i,joff,k)-cc(i,j,k)) &
                        +FOURTH*ss(1,i,j,k)*(cc(i+1,j,k)-cc(i,j,k))& 
                        +FOURTH*ss(6,i,j,k)*(cc(i,j,k-1)-cc(i,j,k)) 

                   cell_mp = FOURTH*ss(sig_j,i,j,k)*(cc(i,joff,k)-cc(i,j,k)) &
                        +FOURTH*ss(2,i,j,k)*(cc(i-1,j,k)-cc(i,j,k)) &
                        +FOURTH*ss(5,i,j,k)*(cc(i,j,k+1)-cc(i,j,k)) 

                   cell_pp = FOURTH*ss(sig_j,i,j,k)*(cc(i,joff,k)-cc(i,j,k)) &
                        +FOURTH*ss(1,i,j,k)*(cc(i+1,j,k)-cc(i,j,k)) &
                        +FOURTH*ss(5,i,j,k)*(cc(i,j,k+1)-cc(i,j,k)) 

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
          !$OMP END PARALLEL DO
          !
          !   Lo/Hi k side
          !
       else if (side == -3 .or. side == 3) then

          k     = lo(3)
          kfine = ir(3)*k

          if (side == -3) then
             koff   = k+1
             sig_k  = 5
          else
             koff   = k-1
             sig_k  = 6
          end if

          !$OMP PARALLEL DO PRIVATE(i,j,lo_i_not,hi_i_not,lo_j_not,hi_j_not) &
          !$OMP PRIVATE(cell_mm,cell_pm,cell_mp,cell_pp,crse_flux,ifine,jfine)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                jfine = ir(2)*j
                ifine = ir(1)*i

                if (bc_dirichlet(mm_fine(ifine,jfine,kfine),1,0) .and. (.not.bc_dirichlet(mm_crse(i,j,k),1,0))) then

                   lo_i_not = .false.; hi_i_not = .false.; lo_j_not = .false.; hi_j_not = .false.

                   if (i == loflux(1)) then
                      if (.not. bc_neumann(mm_fine(ifine,jfine,kfine),1,-1)) lo_i_not = .true.
                   end if
                   if (i == hiflux(1)) then
                      if (.not. bc_neumann(mm_fine(ifine,jfine,kfine),1,+1)) hi_i_not = .true.
                   end if
                   if (j == loflux(2)) then
                      if (.not. bc_neumann(mm_fine(ifine,jfine,kfine),2,-1)) lo_j_not = .true.
                   end if
                   if (j == hiflux(2)) then
                      if (.not. bc_neumann(mm_fine(ifine,jfine,kfine),2,+1)) hi_j_not = .true.
                   end if

                   cell_mm = FOURTH*ss(sig_k,i,j,k)*(cc(i,j,koff)-cc(i,j,k)) &
                        +FOURTH*ss(2,i,j,k)*(cc(i-1,j,k)-cc(i,j,k)) &
                        +FOURTH*ss(4,i,j,k)*(cc(i,j-1,k)-cc(i,j,k)) 

                   cell_pm = FOURTH*ss(sig_k,i,j,k)*(cc(i,j,koff)-cc(i,j,k)) &
                        +FOURTH*ss(1,i,j,k)*(cc(i+1,j,k)-cc(i,j,k))& 
                        +FOURTH*ss(4,i,j,k)*(cc(i,j-1,k)-cc(i,j,k)) 

                   cell_mp = FOURTH*ss(sig_k,i,j,k)*(cc(i,j,koff)-cc(i,j,k)) &
                        +FOURTH*ss(2,i,j,k)*(cc(i-1,j,k)-cc(i,j,k)) &
                        +FOURTH*ss(3,i,j,k)*(cc(i,j+1,k)-cc(i,j,k)) 

                   cell_pp = FOURTH*ss(sig_k,i,j,k)*(cc(i,j,koff)-cc(i,j,k)) &
                        +FOURTH*ss(1,i,j,k)*(cc(i+1,j,k)-cc(i,j,k)) &
                        +FOURTH*ss(3,i,j,k)*(cc(i,j+1,k)-cc(i,j,k)) 

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
          !$OMP END PARALLEL DO
       end if

    end if

  end subroutine ml_interface_3d_nodal

  subroutine ml_fine_contrib(flux, res, mm, ratio, crse_domain, side)
    use bl_prof_module
    type(multifab), intent(inout) :: flux
    type(multifab), intent(inout) :: res
    type(imultifab), intent(in) :: mm
    type(box) :: crse_domain
    type(box) :: fbox
    integer :: side
    integer :: ratio(:)
    integer :: lof(get_dim(flux)), dm
    integer :: lo_dom(get_dim(flux)), hi_dom(get_dim(flux))
    integer :: i, n, dir
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: nc
    logical :: pmask(get_dim(res))
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_fine_contrib")

    nc = ncomp(res)

    if ( ncomp(res) /= ncomp(flux) ) then
       call bl_error("ML_FILL_FLUXES: res%nc /= flux%nc")
    end if

    lo_dom = lwb(crse_domain)
    hi_dom = upb(crse_domain)

    if ( nodal_q(res) ) hi_dom = hi_dom + 1

    dir   = iabs(side)
    pmask = get_pmask(get_layout(res))
    dm    = get_dim(flux)

    !$OMP PARALLEL DO PRIVATE(i,n,fbox,lof,fp,rp,mp)
    do i = 1, nfabs(flux)
       fbox   = get_ibox(flux,i)
       lof = lwb(fbox)
       fp => dataptr(flux, i)
       rp => dataptr(res, i)
       mp => dataptr(mm, i)
       do n = 1, nc
          if ( pmask(dir) .or. &
               (lof(dir) /= lo_dom(dir) .and. lof(dir) /= hi_dom(dir)) ) then
             select case(dm)
             case (1)
                call fine_edge_resid_1d(fp(:,1,1,n), rp(:,1,1,1), mp(:,1,1,1), ratio, side, lof)
             case (2)
                call fine_edge_resid_2d(fp(:,:,1,n), rp(:,:,1,1), mp(:,:,1,1), ratio, side, lof)
             case (3)
                call fine_edge_resid_3d(fp(:,:,:,n), rp(:,:,:,1), mp(:,:,:,1), ratio, side, lof)
             end select
          end if
       end do
    end do
    !$OMP END PARALLEL DO
    call destroy(bpt)
  end subroutine ml_fine_contrib

  subroutine fine_edge_resid_1d(dd, res, mm, ratio, side, lod)

    integer           , intent(in   ) :: lod(1)
    real (kind = dp_t), intent(inout) :: dd(lod(1):)
    real (kind = dp_t), intent(in   ) :: res(-1:)
    integer           , intent(in   ) ::  mm( 0:)
    integer, intent(in) :: ratio(:), side

    integer            :: i,ic,m,isign,ioff,nx,nxc
    real (kind = dp_t) :: fac,fac0

    nx  = size(mm,dim=1)-1
    nxc = size(dd,dim=1)

!   Lo/Hi i side
    if (side == -1) then
       i  = 0
       isign =  1 
    else
       i  = nx
       isign = -1
    end if

    ic = lod(1)

    dd(ic) = res(i)

!   Average towards the interior of the fine grid
    fac0 = ONE / ratio(1)
    do m = 1,ratio(1)-1
       fac = (ratio(1)-m) * fac0
       ioff = i+isign*m
       dd(ic) = dd(ic) + fac * res(ioff)
    end do

    if (.not.bc_dirichlet(mm(i),1,0)) dd(ic) = ZERO

  end subroutine fine_edge_resid_1d

  subroutine fine_edge_resid_2d(dd, res, mm, ratio, side, lod)

    use impose_neumann_bcs_module 

    integer           , intent(in   ) :: lod(2)
    real (kind = dp_t), intent(inout) :: dd(lod(1):,lod(2):)
    real (kind = dp_t), intent(inout) :: res(-1:,-1:)
    integer           , intent(in   ) ::  mm( 0:, 0:)
    integer, intent(in) :: ratio(:), side
    integer :: nx, ny, nxc, nyc
    integer :: hid(2)
    integer :: lo(2),ng_res
    integer :: i,j,ic,jc,m,n,isign,ioff,joff
    integer :: ileft,irght,jbot,jtop
    real (kind = dp_t) :: fac, fac0, fac1
    logical llo,lhi

    nx = size(mm,dim=1)-1
    ny = size(mm,dim=2)-1

    nxc = size(dd,dim=1)
    nyc = size(dd,dim=2)

    hid(1) = lod(1) + nxc-1
    hid(2) = lod(2) + nyc-1

    lo(:) = 0
    ng_res = 1
    call impose_neumann_bcs_2d(res,mm,lo,ng_res)

!   Lo/Hi i side
    if (side == -1 .or. side == 1) then

      if (side == -1) then
         i  = 0
         isign =  1
      else
         i  = nx
         isign = -1
      end if

      ic = lod(1)
      fac0 = ONE / ratio(2)

!     First average along the coarse-fine edge
      do jc = lod(2),hid(2)
         n = 0
         fac = HALF*ratio(2)*fac0
         j = (jc-lod(2))*ratio(2)
         if (j >  0) dd(ic,jc) = dd(ic,jc) + fac * res(i,j)
         if (j < ny) dd(ic,jc) = dd(ic,jc) + fac * res(i,j)

         do n = 1,ratio(2)-1
            fac = (ratio(2)-n)*fac0

            j = (jc-lod(2))*ratio(2) + n

            if (j < ny) then
               if (jc==lod(2)) then
                  if (.not. bc_dirichlet(mm(i,j),1,0)) fac = HALF * fac
               end if
               dd(ic,jc) = dd(ic,jc) + fac * res(i,j)
            end if

            j = (jc-lod(2))*ratio(2) - n

            if (j > 0) then
               if (jc==hid(2)) then
                  if (.not. bc_dirichlet(mm(i,j),1,0)) fac = HALF * fac
               end if
               dd(ic,jc) = dd(ic,jc) + fac * res(i,j)
            end if

         end do

      end do

      j = 0
      if (bc_neumann(mm(i,j),2,-1)) dd(ic,lod(2)) = TWO*dd(ic,lod(2))

      j = (hid(2)-lod(2))*ratio(2)
      if (bc_neumann(mm(i,j),2, 1)) dd(ic,hid(2)) = TWO*dd(ic,hid(2))

!     Now average towards the interior of the fine grid
      fac0 = fac0 / ratio(1)
      do n = 0,ratio(2)-1
         fac1 = (ratio(2)-n) * fac0
         if (n == 0) fac1 = HALF * fac1
         do m = 1,ratio(1)-1
            fac = (ratio(1)-m) * fac1
            ioff = i+isign*m
            do jc = lod(2),hid(2)
               j = (jc-lod(2))*ratio(2)
               jbot = j-n
               jtop = j+n

               if (j==0) then
                  if (bc_neumann(mm(i,j),2,-1)) jbot = jtop
               end if
               if (j==ny) then
                  if (bc_neumann(mm(i,j),2,+1)) jtop = jbot
               end if

               llo = .false.
               lhi = .false.

               if (j==0) then
                  if (.not. bc_neumann(mm(i,j),2,-1)) llo = .true.
               end if
               if (j==ny) then
                  if (.not. bc_neumann(mm(i,j),2,+1)) lhi = .true.
               end if

               if (llo) then
                  if (n==0) then
                     if (.not. bc_dirichlet(mm(ioff,j),1,0)) &
                          dd(ic,jc) = dd(ic,jc) + fac * res(ioff,j)
                  else
                     dd(ic,jc) = dd(ic,jc) + HALF * fac * res(ioff,jtop)
                  end if
               else if (lhi) then
                  if (n==0) then
                     if (.not. bc_dirichlet(mm(ioff,j),1,0)) &
                          dd(ic,jc) = dd(ic,jc) + fac * res(ioff,j)
                  else
                     dd(ic,jc) = dd(ic,jc) + HALF * fac * res(ioff,jbot)
                  end if
               else
                  dd(ic,jc) = dd(ic,jc) + fac * ( res(ioff,jtop) + &
                                                  res(ioff,jbot) )
               end if 
            end do
         end do
      end do

      do jc = lod(2),hid(2)
         if (.not.bc_dirichlet(mm(i,(jc-lod(2))*ratio(2)),1,0)) dd(ic,jc) = ZERO
      end do

!   Lo/Hi j side
    else if (side == -2 .or. side == 2) then

      if (side == -2) then
         j  = 0
         isign =  1
      else
         j  = ny
         isign = -1
      end if

      jc = lod(2)
      fac0 = ONE / ratio(1) 

!     First average along the coarse-fine edge
      do ic = lod(1),hid(1)
         do n = 0,ratio(1)-1
            fac = (ratio(1)-n)*fac0
            if (n == 0) fac = HALF * fac

            i = (ic-lod(1))*ratio(1) + n

            if (i == 0) then
               dd(ic,jc) = dd(ic,jc) + fac * res(i,j)
            else if (i < nx) then
               if (ic==lod(1) .and. n>0) then
                  if (.not. bc_dirichlet(mm(i,j),1,0)) fac = HALF * fac
               end if
               dd(ic,jc) = dd(ic,jc) + fac * res(i,j)
            end if

            i = (ic-lod(1))*ratio(1) - n
            if (i == nx) then
              dd(ic,jc) = dd(ic,jc) + fac * res(i,j)
            else if (i > 0) then
               if (ic==hid(1) .and. n>0) then
                  if (.not. bc_dirichlet(mm(i,j),1,0)) fac = HALF * fac
               end if
               dd(ic,jc) = dd(ic,jc) + fac * res(i,j)
            end if
         end do
      end do

      i = 0
      if (bc_neumann(mm(i,j),1,-1)) dd(lod(1),jc) = TWO*dd(lod(1),jc)

      i = (hid(1)-lod(1))*ratio(1)
      if (bc_neumann(mm(i,j),1, 1)) dd(hid(1),jc) = TWO*dd(hid(1),jc)

!     Now average towards the interior of the fine grid
      fac0 = fac0 / ratio(2)
      do n = 0,ratio(1)-1
         fac1 = (ratio(1)-n) * fac0
         if (n == 0) fac1 = HALF * fac1
         do m = 1,ratio(2)-1
            joff = j + isign*m
            fac = (ratio(2)-m) * fac1
            do ic = lod(1),hid(1)
               i = (ic-lod(1))*ratio(1)
               ileft = i-n
               irght = i+n

               if (i==0) then
                  if (bc_neumann(mm(i,j),1,-1)) ileft = irght
               end if
               if (i==nx) then
                  if (bc_neumann(mm(i,j),1,+1)) irght = ileft
               end if

               llo = .false.
               lhi = .false.

               if (i==0) then
                  if (.not. bc_neumann(mm(i,j),1,-1)) llo = .true.
               end if
               if (i==nx) then
                  if (.not. bc_neumann(mm(i,j),1,+1)) lhi = .true.
               end if

               if (llo) then
                  if (n==0) then
                     if (.not. bc_dirichlet(mm(i,joff),1,0)) &
                          dd(ic,jc) = dd(ic,jc) + fac * res(i,joff)
                  else
                     dd(ic,jc) = dd(ic,jc) + HALF * fac * res(irght,joff)
                  end if
               else if (lhi) then
                  if (n==0) then
                     if (.not. bc_dirichlet(mm(i,joff),1,0)) &
                          dd(ic,jc) = dd(ic,jc) + fac * res(i,joff)
                  else
                     dd(ic,jc) = dd(ic,jc) + HALF * fac * res(ileft,joff)
                  end if

               else
                  dd(ic,jc) = dd(ic,jc) + fac * ( res(irght,joff) + &
                                                  res(ileft,joff) )
               end if
            end do
         end do
      end do

      do ic = lod(1),hid(1)
         if (.not.bc_dirichlet(mm((ic-lod(1))*ratio(1),j),1,0)) dd(ic,jc) = ZERO
      end do

    end if

  end subroutine fine_edge_resid_2d

  subroutine fine_edge_resid_3d(dd, res, mm, ratio, side, lod)

    use impose_neumann_bcs_module 

    integer, intent(in) :: lod(:)
    real (kind = dp_t), intent(inout) :: dd(lod(1):,lod(2):,lod(3):)
    real (kind = dp_t), intent(inout) :: res(-1:,-1:,-1:)
    integer           , intent(in   ) ::  mm(0:,0:,0:)
    integer, intent(in) :: ratio(:),side
    integer :: nx, ny, nz, nxc, nyc, nzc
    integer :: hid(3),lo(3),ng_res
    integer :: i,j,k,l,ic,jc,kc,m,n
    integer :: isign,ioff,joff,koff
    integer :: ileft,irght,jbot,jtop,kdwn,kup
    real (kind = dp_t) :: fac, fac0, fac1, fac2, corner_fac
    logical ll1,ll2,ll3,lh1,lh2,lh3

    real (kind = dp_t), parameter :: VAL = (0.25_dp_t / 3.0_dp_t)

    nx = size(mm,dim=1)-1
    ny = size(mm,dim=2)-1
    nz = size(mm,dim=3)-1

    nxc = size(dd,dim=1)
    nyc = size(dd,dim=2)
    nzc = size(dd,dim=3)

    hid(1) = lod(1) + nxc-1
    hid(2) = lod(2) + nyc-1
    hid(3) = lod(3) + nzc-1

    lo     = 0
    ng_res = 1
    call impose_neumann_bcs_3d(res,mm,lo,ng_res)

    if (side == -1 .or. side == 1) then

      if (side == -1) then
         i     = 0
         isign =  1
      else
         i     = nx
         isign = -1
      end if

      ic   = lod(1)
      fac0 = 1.0_dp_t / (ratio(2)*ratio(3))
      !
      ! First average along the coarse-fine face.
      !
      do kc = lod(3),hid(3)
         do jc = lod(2),hid(2)
            do n = 0,ratio(2)-1
               fac2 = (ratio(2)-n)*fac0
               if (n == 0) fac2 = HALF * fac2

               j = (jc-lod(2))*ratio(2) + n
               if (j < ny) then
                  do l = 0,ratio(3)-1
                     fac = (ratio(3)-l)*fac2
                     if (l == 0) fac = HALF * fac

                     k = (kc-lod(3))*ratio(3) + l
                     if (k < nz) dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)

                     k = (kc-lod(3))*ratio(3) - l
                     if (k >  0) dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                  end do
               end if

               j = (jc-lod(2))*ratio(2) - n
               if (j > 0) then
                  do l = 0,ratio(3)-1
                     fac = (ratio(3)-l)*fac2
                     if (l == 0) fac = HALF * fac

                     k = (kc-lod(3))*ratio(3) + l
                     if (k < nz) dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)

                     k = (kc-lod(3))*ratio(3) - l
                     if (k >  0) dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                  end do
               end if

            end do
         end do
      end do

      jc = lod(2)
      kc = lod(3)
      j  = 0
      k  = 0
      if ( .not. bc_neumann(mm(i,j,k),2,-1) .and. (.not. bc_neumann(mm(i,j,k),3,-1)) ) then
         dd(ic,jc,kc) = dd(ic,jc,kc) + VAL * res(i,j,k)
      end if

      jc = hid(2)
      kc = lod(3)
      j  = ny
      k  = 0
      if ( .not. bc_neumann(mm(i,j,k),2,+1) .and. (.not. bc_neumann(mm(i,j,k),3,-1)) ) then
         dd(ic,jc,kc) = dd(ic,jc,kc) + VAL * res(i,j,k)
      end if

      jc = lod(2)
      kc = hid(3)
      j  = 0
      k  = nz
      if ( .not. bc_neumann(mm(i,j,k),2,-1) .and. (.not. bc_neumann(mm(i,j,k),3,+1)) ) then
         dd(ic,jc,kc) = dd(ic,jc,kc) + VAL * res(i,j,k)
      end if

      jc = hid(2)
      kc = hid(3)
      j  = ny
      k  = nz
      if ( .not. bc_neumann(mm(i,j,k),2,+1) .and. (.not. bc_neumann(mm(i,j,k),3,+1)) ) then
         dd(ic,jc,kc) = dd(ic,jc,kc) + VAL * res(i,j,k)
      end if

      j = 0
      do kc = lod(3),hid(3)
         k = (kc-lod(3))*ratio(3)
         if ( bc_neumann(mm(i,j,k),2,-1) ) dd(ic,lod(2),kc) = TWO*dd(ic,lod(2),kc)
      end do

      j = (hid(2)-lod(2))*ratio(2)
      do kc = lod(3),hid(3)
         k = (kc-lod(3))*ratio(3)
         if ( bc_neumann(mm(i,j,k),2, 1) ) dd(ic,hid(2),kc) = TWO*dd(ic,hid(2),kc)
      end do

      k = 0
      do jc = lod(2),hid(2)
         j = (jc-lod(2))*ratio(2)
         if ( bc_neumann(mm(i,j,k),3,-1) ) dd(ic,jc,lod(3)) = TWO*dd(ic,jc,lod(3))
      end do

      k = (hid(3)-lod(3))*ratio(3)
      do jc = lod(2),hid(2)
         j = (jc-lod(2))*ratio(2)
         if ( bc_neumann(mm(i,j,k),3, 1) ) dd(ic,jc,hid(3)) = TWO*dd(ic,jc,hid(3))
      end do
      !
      ! Now average towards the interior of the grid.
      !
      fac0 = fac0 / ratio(1)
      ic = lod(1)
      do l = 0, ratio(3)-1
        fac2 = (ratio(3)-l) * fac0
        if (l == 0) fac2 = HALF * fac2
        do n = 0, ratio(2)-1
          fac1 = (ratio(2)-n) * fac2
          if (n == 0) fac1 = HALF * fac1
          do m = 1, ratio(1)-1
            ioff = i+isign*m
            fac = (ratio(1)-m) * fac1
            if (m == 0) fac = HALF * fac

            do kc = lod(3),hid(3)
              k = (kc-lod(3))*ratio(3)
              do jc = lod(2),hid(2)
                j = (jc-lod(2))*ratio(2)
                jtop = j+n
                jbot = j-n
                kup  = k+l
                kdwn = k-l
                if (j==0  .and. bc_neumann(mm(i,j,k),2,-1)) jbot = jtop
                if (j==ny .and. bc_neumann(mm(i,j,k),2,+1)) jtop = jbot
                if (k==0  .and. bc_neumann(mm(i,j,k),3,-1)) kdwn = kup
                if (k==nz .and. bc_neumann(mm(i,j,k),3,+1)) kup  = kdwn

                ll2 = .false.; lh2 = .false.; ll3 = .false.; lh3 = .false.

                if (jc==lod(2) .and. (.not. bc_neumann(mm(i,j,k),2,-1))) ll2 = .true.
                if (jc==hid(2) .and. (.not. bc_neumann(mm(i,j,k),2,+1))) lh2 = .true.
                if (kc==lod(3) .and. (.not. bc_neumann(mm(i,j,k),3,-1))) ll3 = .true.
                if (kc==hid(3) .and. (.not. bc_neumann(mm(i,j,k),3,+1))) lh3 = .true.

                if ( ( ll2 .or. lh2 ) .and. ( ll3 .or. lh3 ) ) then
                   corner_fac = THIRD
                else if ( ( ll2 .or. lh2 ) .and. .not. ( ll3 .or. lh3 ) ) then
                   corner_fac = HALF
                else if ( .not. ( ll2 .or. lh2 ) .and. ( ll3 .or. lh3 ) ) then
                   corner_fac = HALF
                else
                   corner_fac = 1.0_dp_t
                end if

                ll2 = (j-n >  0); if (.not. ll2) ll2 = bc_neumann(mm(i,j,k),2,-1)
                lh2 = (j-n < ny); if (.not. lh2) lh2 = bc_neumann(mm(i,j,k),2,+1)
                ll3 = (k-l >  0); if (.not. ll3) ll3 = bc_neumann(mm(i,j,k),3,-1)
                lh3 = (k-l < nz); if (.not. lh3) lh3 = bc_neumann(mm(i,j,k),3,+1)
                ll1 = (k+l >  0); if (.not. ll1) ll1 = bc_neumann(mm(i,j,k),3,-1)
                lh1 = (k+l < nz); if (.not. lh1) lh1 = bc_neumann(mm(i,j,k),3,+1)

                if ( ll2 .and. lh2 ) then
                   if ( ll3 .and. lh3 ) dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * fac * res(ioff,jbot,kdwn) 
                   if ( ll1 .and. lh1 ) dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * fac * res(ioff,jbot,kup) 
                end if

                ll2 = (j+n >  0); if (.not. ll2) ll2 = bc_neumann(mm(i,j,k),2,-1)
                lh2 = (j+n < ny); if (.not. lh2) lh2 = bc_neumann(mm(i,j,k),2,+1)

                if ( ll2 .and. lh2 ) then
                   if ( ll3 .and. lh3 ) dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * fac * res(ioff,jtop,kdwn) 
                   if ( ll1 .and. lh1 ) dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * fac * res(ioff,jtop,kup) 
                end if

              end do
            end do

          end do
        end do
      end do

      do kc = lod(3),hid(3)
         do jc = lod(2),hid(2)
            j = (jc-lod(2))*ratio(2)
            k = (kc-lod(3))*ratio(3)
            if ( .not. bc_dirichlet(mm(i,j,k),1,0) ) dd(ic,jc,kc) = ZERO
         end do
      end do

    else if (side == -2 .or. side == 2) then

      if (side == -2) then
         j     = 0
         isign =  1
      else
         j     = ny
         isign = -1
      end if

      jc   = lod(2)
      fac0 = 1.0_dp_t / (ratio(1)*ratio(3))
      !
      ! First average along the coarse-fine face.
      !
      do kc = lod(3),hid(3)
         do ic = lod(1),hid(1)
            do n = 0,ratio(1)-1
               fac2 = (ratio(1)-n)*fac0
               if (n == 0) fac2 = HALF * fac2

               i = (ic-lod(1))*ratio(1) + n
               if (i < nx) then
                  do l = 0,ratio(3)-1
                     fac = (ratio(3)-l)*fac2
                     if (l == 0) fac = HALF * fac

                     k = (kc-lod(3))*ratio(3) + l
                     if (k < nz) dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)

                     k = (kc-lod(3))*ratio(3) - l
                     if (k >  0) dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                  end do
               end if

               i = (ic-lod(1))*ratio(1) - n
               if (i > 0) then
                  do l = 0,ratio(3)-1
                     fac = (ratio(3)-l)*fac2
                     if (l == 0) fac = HALF * fac

                     k = (kc-lod(3))*ratio(3) + l
                     if (k < nz) dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)

                     k = (kc-lod(3))*ratio(3) - l
                     if (k >  0) dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                  end do
               end if

            end do
         end do
      end do

      ic = lod(1)
      kc = lod(3)
      i  = 0
      k  = 0
      if ( .not. bc_neumann(mm(i,j,k),1,-1) .and. (.not. bc_neumann(mm(i,j,k),3,-1)) ) then
         dd(ic,jc,kc) = dd(ic,jc,kc) + VAL * res(i,j,k)
      end if

      ic = hid(1)
      kc = lod(3)
      i  = nx
      k  = 0
      if ( .not. bc_neumann(mm(i,j,k),1,+1) .and. (.not. bc_neumann(mm(i,j,k),3,-1)) ) then
         dd(ic,jc,kc) = dd(ic,jc,kc) + VAL * res(i,j,k)
      end if

      ic = lod(1)
      kc = hid(3)
      i  = 0
      k  = nz
      if ( .not. bc_neumann(mm(i,j,k),1,-1) .and. (.not. bc_neumann(mm(i,j,k),3,+1)) ) then
         dd(ic,jc,kc) = dd(ic,jc,kc) + VAL * res(i,j,k)
      end if

      ic = hid(1)
      kc = hid(3)
      i  = nx
      k  = nz
      if ( .not. bc_neumann(mm(i,j,k),1,+1) .and. (.not. bc_neumann(mm(i,j,k),3,+1)) ) then
         dd(ic,jc,kc) = dd(ic,jc,kc) + VAL * res(i,j,k)
      end if

      i = 0
      do kc = lod(3),hid(3)
         k = (kc-lod(3))*ratio(3)
         if ( bc_neumann(mm(i,j,k),1,-1) ) dd(lod(1),jc,kc) = TWO*dd(lod(1),jc,kc)
      end do

      i = (hid(1)-lod(1))*ratio(1)
      do kc = lod(3),hid(3)
         k = (kc-lod(3))*ratio(3)
         if ( bc_neumann(mm(i,j,k),1, 1) ) dd(hid(1),jc,kc) = TWO*dd(hid(1),jc,kc)
      end do

      k = 0
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         if ( bc_neumann(mm(i,j,k),3,-1) ) dd(ic,jc,lod(3)) = TWO*dd(ic,jc,lod(3))
      end do

      k = (hid(3)-lod(3))*ratio(3)
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         if ( bc_neumann(mm(i,j,k),3, 1) ) dd(ic,jc,hid(3)) = TWO*dd(ic,jc,hid(3))
      end do
      !
      ! Now average towards the interior of the grid.
      !
      fac0 = fac0 / ratio(2)
      jc = lod(2)
      do l = 0, ratio(3)-1
        fac2 = (ratio(3)-l) * fac0
        if (l == 0) fac2 = HALF * fac2
        do n = 0, ratio(1)-1
          fac1 = (ratio(1)-n) * fac2
          if (n == 0) fac1 = HALF * fac1
          do m = 1, ratio(2)-1
            joff = j+isign*m
            fac = (ratio(2)-m) * fac1
            if (m == 0) fac = HALF * fac

            do kc = lod(3),hid(3)
              k = (kc-lod(3))*ratio(3)
              do ic = lod(1),hid(1)
                i = (ic-lod(1))*ratio(1)
                irght = i+n
                ileft = i-n
                kup  = k+l
                kdwn = k-l
                if (i==0  .and. bc_neumann(mm(i,j,k),1,-1)) ileft = irght
                if (i==nx .and. bc_neumann(mm(i,j,k),1,+1)) irght = ileft
                if (k==0  .and. bc_neumann(mm(i,j,k),3,-1)) kdwn = kup
                if (k==nz .and. bc_neumann(mm(i,j,k),3,+1)) kup  = kdwn

                ll1 = .false.; lh1 = .false.; ll3 = .false.; lh3 = .false.

                if (ic==lod(1) .and. (.not. bc_neumann(mm(i,j,k),1,-1))) ll1 = .true.
                if (ic==hid(1) .and. (.not. bc_neumann(mm(i,j,k),1,+1))) lh1 = .true.
                if (kc==lod(3) .and. (.not. bc_neumann(mm(i,j,k),3,-1))) ll3 = .true.
                if (kc==hid(3) .and. (.not. bc_neumann(mm(i,j,k),3,+1))) lh3 = .true.

                if ( (  ll1 .or. lh1 ) .and. (  ll3 .or. lh3 ) ) then
                   corner_fac = THIRD
                else if ( ( ll1 .or. lh1 ) .and. .not. ( ll3 .or. lh3 ) ) then
                   corner_fac = HALF
                else if ( .not. &
                          ( ll1 .or. lh1 ) .and. ( ll3 .or. lh3 ) ) then
                   corner_fac = HALF
                else
                   corner_fac = 1.0_dp_t
                end if

                ll1 = (i-n >  0); if (.not. ll1) ll1 = bc_neumann(mm(i,j,k),1,-1)
                lh1 = (i-n < nx); if (.not. lh1) lh1 = bc_neumann(mm(i,j,k),1,+1)
                ll2 = (k-l >  0); if (.not. ll2) ll2 = bc_neumann(mm(i,j,k),3,-1)
                lh2 = (k-l < nz); if (.not. lh2) lh2 = bc_neumann(mm(i,j,k),3,+1)
                ll3 = (k+l >  0); if (.not. ll3) ll3 = bc_neumann(mm(i,j,k),3,-1)
                lh3 = (k+l < nz); if (.not. lh3) lh3 = bc_neumann(mm(i,j,k),3,+1)

                if ( ll1 .and. lh1 ) then
                   if ( ll2 .and. lh2 ) dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * fac * res(ileft,joff,kdwn) 
                   if ( ll3 .and. lh3 ) dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * fac * res(ileft,joff,kup) 
                end if

                ll1 = (i+n >  0); if (.not. ll1) ll1 = bc_neumann(mm(i,j,k),1,-1)
                lh1 = (i+n < nx); if (.not. lh1) lh1 = bc_neumann(mm(i,j,k),1,+1)

                if ( ll1 .and. lh1 ) then
                   if ( ll2 .and. lh2 ) dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * fac * res(irght,joff,kdwn) 
                   if ( ll3 .and. lh3 ) dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * fac * res(irght,joff,kup) 
                end if

              end do
            end do

          end do
        end do
      end do

      do kc = lod(3),hid(3)
         do ic = lod(1),hid(1)
            i = (ic-lod(1))*ratio(1)
            k = (kc-lod(3))*ratio(3)
            if ( .not. bc_dirichlet(mm(i,j,k),1,0) ) dd(ic,jc,kc) = ZERO
         end do
      end do

    else 

      if (side == -3) then
         k     = 0
         isign =  1
      else
         k     = nz
         isign = -1
      end if

      kc   = lod(3)
      fac0 = 1.0_dp_t / (ratio(1)*ratio(2))
      !
      ! First average along the coarse-fine face.
      !
      do jc = lod(2),hid(2)
         do ic = lod(1),hid(1)
            do n = 0,ratio(1)-1
               fac2 = (ratio(1)-n)*fac0
               if (n == 0) fac2 = HALF * fac2

               i = (ic-lod(1))*ratio(1) + n
               if (i < nx) then
                  do l = 0,ratio(2)-1
                     fac = (ratio(2)-l)*fac2
                     if (l == 0) fac = HALF * fac

                     j = (jc-lod(2))*ratio(2) + l
                     if (j < ny) dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)

                     j = (jc-lod(2))*ratio(2) - l
                     if (j >  0) dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                  end do
               end if

               i = (ic-lod(1))*ratio(1) - n
               if (i > 0) then
                  do l = 0,ratio(2)-1
                     fac = (ratio(2)-l)*fac2
                     if (l == 0) fac = HALF * fac

                     j = (jc-lod(2))*ratio(2) + l
                     if (j < ny) dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)

                     j = (jc-lod(2))*ratio(2) - l
                     if (j >  0) dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                  end do
               end if

            end do
         end do
      end do

      ic = lod(1)
      jc = lod(2)
      i  = 0
      j  = 0
      if ( .not. bc_neumann(mm(i,j,k),1,-1) .and. (.not. bc_neumann(mm(i,j,k),2,-1)) ) then
         dd(ic,jc,kc) = dd(ic,jc,kc) + VAL * res(i,j,k)
      end if

      ic = hid(1)
      jc = lod(2)
      i  = nx
      j  = 0
      if ( .not. bc_neumann(mm(i,j,k),1,+1) .and. (.not. bc_neumann(mm(i,j,k),2,-1)) ) then
         dd(ic,jc,kc) = dd(ic,jc,kc) + VAL * res(i,j,k)
      end if

      ic = lod(1)
      jc = hid(2)
      i  = 0
      j  = ny
      if ( .not. bc_neumann(mm(i,j,k),1,-1) .and. (.not. bc_neumann(mm(i,j,k),2,+1)) ) then
         dd(ic,jc,kc) = dd(ic,jc,kc) + VAL * res(i,j,k)
      end if

      ic = hid(1)
      jc = hid(2)
      i  = nx
      j  = ny
      if ( .not. bc_neumann(mm(i,j,k),1,+1) .and. (.not. bc_neumann(mm(i,j,k),2,+1)) ) then
         dd(ic,jc,kc) = dd(ic,jc,kc) + VAL * res(i,j,k)
      end if

      i = 0
      do jc = lod(2),hid(2)
         j = (jc-lod(2))*ratio(2)
         if ( bc_neumann(mm(i,j,k),1,-1) ) dd(lod(1),jc,kc) = TWO*dd(lod(1),jc,kc)
      end do

      i = (hid(1)-lod(1))*ratio(1)
      do jc = lod(2),hid(2)
         j = (jc-lod(2))*ratio(2)
         if ( bc_neumann(mm(i,j,k),1,+1) ) dd(hid(1),jc,kc) = TWO*dd(hid(1),jc,kc)
      end do

      j = 0
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         if ( bc_neumann(mm(i,j,k),2,-1) ) dd(ic,lod(2),kc) = TWO*dd(ic,lod(2),kc)
      end do

      j = (hid(2)-lod(2))*ratio(2)
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         if ( bc_neumann(mm(i,j,k),2,+1) ) dd(ic,hid(2),kc) = TWO*dd(ic,hid(2),kc)
      end do
      !
      ! Now average towards the interior of the grid.
      !
      fac0 = fac0 / ratio(3)
      kc = lod(3)
      do l = 0, ratio(2)-1
        fac2 = (ratio(2)-l) * fac0
        if (l == 0) fac2 = HALF * fac2
        do n = 0, ratio(1)-1
          fac1 = (ratio(1)-n) * fac2
          if (n == 0) fac1 = HALF * fac1
          do m = 1, ratio(3)-1
            koff = k+isign*m
            fac = (ratio(3)-m) * fac1
            if (m == 0) fac = HALF * fac

            do jc = lod(2),hid(2)
              j = (jc-lod(2))*ratio(2)
              do ic = lod(1),hid(1)
                i = (ic-lod(1))*ratio(1)
                irght = i+n
                ileft = i-n
                jtop  = j+l
                jbot  = j-l
                if (i==0  .and. bc_neumann(mm(i,j,k),1,-1)) ileft = irght
                if (i==nx .and. bc_neumann(mm(i,j,k),1,+1)) irght = ileft
                if (j==0  .and. bc_neumann(mm(i,j,k),2,-1)) jbot  = jtop
                if (j==ny .and. bc_neumann(mm(i,j,k),2,+1)) jtop  = jbot

                ll1 = .false.; lh1 = .false.; ll2 = .false.; lh2 = .false.

                if (ic==lod(1) .and. (.not. bc_neumann(mm(i,j,k),1,-1))) ll1 = .true.
                if (ic==hid(1) .and. (.not. bc_neumann(mm(i,j,k),1,+1))) lh1 = .true.
                if (jc==lod(2) .and. (.not. bc_neumann(mm(i,j,k),2,-1))) ll2 = .true.
                if (jc==hid(2) .and. (.not. bc_neumann(mm(i,j,k),2,+1))) lh2 = .true.

                if ( ( ll1 .or. lh1 ) .and. ( ll2 .or. lh2 ) ) then
                   corner_fac = THIRD
                else if ( ( ll1 .or. lh1 ) .and. .not. ( ll2 .or. lh2 ) ) then
                   corner_fac = HALF
                else if ( .not. &
                          ( ll1 .or. lh1 ) .and. ( ll2 .or. lh2 ) ) then
                   corner_fac = HALF
                else
                   corner_fac = 1.0_dp_t
                end if

                ll1 = (i-n >  0); if (.not. ll1) ll1 = bc_neumann(mm(i,j,k),1,-1) 
                lh1 = (i-n < nx); if (.not. lh1) lh1 = bc_neumann(mm(i,j,k),1,+1)
                ll2 = (j-l >  0); if (.not. ll2) ll2 = bc_neumann(mm(i,j,k),2,-1)
                lh2 = (j-l < ny); if (.not. lh2) lh2 = bc_neumann(mm(i,j,k),2,+1)
                ll3 = (j+l >  0); if (.not. ll3) ll3 = bc_neumann(mm(i,j,k),2,-1)
                lh3 = (j+l < ny); if (.not. lh3) lh3 = bc_neumann(mm(i,j,k),2,+1)

                if ( ll1 .and. lh1 ) then
                   if ( ll2 .and. lh2 ) dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * fac * res(ileft,jbot,koff) 
                   if ( ll3 .and. lh3 ) dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * fac * res(ileft,jtop,koff) 
                end if

                ll1 = (i+n >  0); if (.not. ll1) ll1 = bc_neumann(mm(i,j,k),1,-1) 
                lh1 = (i+n < nx); if (.not. lh1) lh1 = bc_neumann(mm(i,j,k),1,+1)

                if ( ll1 .and. lh1 ) then
                   if ( ll2 .and. lh2 ) dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * fac * res(irght,jbot,koff) 
                   if ( ll3 .and. lh3 ) dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * fac * res(irght,jtop,koff) 
                end if

              end do
            end do

          end do
        end do
      end do

      do jc = lod(2),hid(2)
         do ic = lod(1),hid(1)
            i = (ic-lod(1))*ratio(1)
            j = (jc-lod(2))*ratio(2)
            if ( .not. bc_dirichlet(mm(i,j,k),1,0) ) dd(ic,jc,kc) = ZERO
         end do
      end do

    end if

  end subroutine fine_edge_resid_3d

end module nodal_interface_stencil_module
