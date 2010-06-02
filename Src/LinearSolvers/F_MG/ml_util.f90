module ml_util_module

  use stencil_module
! use stencil_nodal_module

  implicit none

  private

  real (kind = dp_t), parameter :: ZERO = 0.0_dp_t
  real (kind = dp_t), parameter :: HALF = 0.5_dp_t
  real (kind = dp_t), parameter ::  ONE = 1.0_dp_t
  real (kind = dp_t), parameter ::  TWO = 2.0_dp_t

  public :: ml_fill_fluxes, ml_fill_fluxes_c, ml_fill_fine_fluxes, ml_fill_all_fluxes
  public :: ml_fine_contrib, ml_norm_inf, ml_norm_l2

contains

  subroutine ml_fill_fluxes(ss, flux, uu, mm, ratio, face, dim)
    use bl_prof_module
    type(multifab), intent(inout) :: flux
    type(multifab), intent(in) :: ss
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in) :: mm
    integer :: ratio
    integer :: face, dim
    integer :: i, n
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: ng
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_fill_fluxes")

    ng = uu%ng

    if ( uu%nc /= flux%nc ) then
       call bl_error("ML_FILL_FLUXES: uu%nc /= flux%nc")
    end if

    do i = 1, flux%nboxes
       if ( remote(flux, i) ) cycle
       fp => dataptr(flux, i)
       up => dataptr(uu, i)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       do n = 1, uu%nc
          select case(ss%dim)
          case (1)
             call stencil_flux_1d(sp(:,1,1,:), fp(:,1,1,n), up(:,1,1,n), &
                  mp(:,1,1,1), ng, ratio, face, dim)
          case (2)
             call stencil_flux_2d(sp(:,:,1,:), fp(:,:,1,n), up(:,:,1,n), &
                  mp(:,:,1,1), ng, ratio, face, dim)
          case (3)
             call stencil_flux_3d(sp(:,:,:,:), fp(:,:,:,n), up(:,:,:,n), &
                  mp(:,:,:,1), ng, ratio, face, dim)
          end select
       end do
    end do
    call destroy(bpt)
  end subroutine ml_fill_fluxes

  subroutine ml_fill_fluxes_c(ss, flux, cf, uu, cu, mm, ratio, face, dim)
    use bl_prof_module
    type(multifab), intent(inout) :: flux
    type(multifab), intent(in) :: ss
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in) :: mm
    integer, intent(in) :: cf, cu
    integer :: ratio
    integer :: face, dim
    integer :: i
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: ng
    logical :: lcross
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_fill_fluxes_c")

    ng = uu%ng

    lcross = ((ncomp(ss) == 5) .or. (ncomp(ss) == 7))

    call multifab_fill_boundary(uu, cross = lcross)
    do i = 1, flux%nboxes
       if ( remote(flux, i) ) cycle
       fp => dataptr(flux, i, cf)
       up => dataptr(uu, i, cu)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       select case(ss%dim)
       case (1)
          call stencil_flux_1d(sp(:,1,1,:), fp(:,1,1,1), up(:,1,1,1), &
               mp(:,1,1,1), ng, ratio, face, dim)
       case (2)
          call stencil_flux_2d(sp(:,:,1,:), fp(:,:,1,1), up(:,:,1,1), &
               mp(:,:,1,1), ng, ratio, face, dim)
       case (3)
          call stencil_flux_3d(sp(:,:,:,:), fp(:,:,:,1), up(:,:,:,1), &
               mp(:,:,:,1), ng, ratio, face, dim)
       end select
    end do
    call destroy(bpt)
  end subroutine ml_fill_fluxes_c

  subroutine ml_fill_fine_fluxes(ss, flux, uu, mm, face, dim)
    use bl_prof_module
    type(multifab), intent(inout) :: flux
    type(multifab), intent(in) :: ss
    type(multifab), intent(inout) :: uu
    type(imultifab), intent(in) :: mm
    integer :: face, dim
    integer :: i, n
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: ng
    logical :: lcross
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_fill_fine_fluxes")

    ng = uu%ng

    lcross = ((ncomp(ss) == 5) .or. (ncomp(ss) == 7))

    if ( uu%nc /= flux%nc ) then
       call bl_error("ML_FILL_FINE_FLUXES: uu%nc /= flux%nc")
    end if

    call multifab_fill_boundary(uu, cross = lcross)

    do i = 1, flux%nboxes
       if ( remote(flux, i) ) cycle
       fp => dataptr(flux, i)
       up => dataptr(uu, i)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       do n = 1, uu%nc
          select case(ss%dim)
          case (1)
             call stencil_fine_flux_1d(sp(:,1,1,:), fp(:,1,1,n), up(:,1,1,n), &
                  mp(:,1,1,1), ng, face, dim)
          case (2)
             call stencil_fine_flux_2d(sp(:,:,1,:), fp(:,:,1,n), up(:,:,1,n), &
                  mp(:,:,1,1), ng, face, dim)
          case (3)
             call stencil_fine_flux_3d(sp(:,:,:,:), fp(:,:,:,n), up(:,:,:,n), &
                  mp(:,:,:,1), ng, face, dim)
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine ml_fill_fine_fluxes

  subroutine ml_fill_all_fluxes(ss, flux, uu, mm)
    use bl_prof_module
    type( multifab), intent(in   ) :: ss
    type( multifab), intent(inout) :: flux(:)
    type( multifab), intent(inout) :: uu
    type(imultifab), intent(in   ) :: mm

    integer :: dim, i, ngu, ngf
    logical :: lcross

    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)

    type(bl_prof_timer), save :: bpt
    call build(bpt, "ml_fill_all_fluxes")

    ngu = uu%ng

    lcross = ((ncomp(ss) == 5) .or. (ncomp(ss) == 7))

    if ( uu%nc /= flux(1)%nc ) then
       call bl_error("ML_FILL_ALL_FLUXES: uu%nc /= flux%nc")
    end if

    call multifab_fill_boundary(uu, cross = lcross)

    do dim = 1, uu%dim
     do i = 1, flux(dim)%nboxes
       if ( remote(flux(dim), i) ) cycle
       ngf = flux(dim)%ng
       fp => dataptr(flux(dim), i)
       up => dataptr(uu, i)
       sp => dataptr(ss, i)
       mp => dataptr(mm, i)
       select case(ss%dim)
       case (1)
          call stencil_all_flux_1d(sp(:,1,1,:), fp(:,1,1,1), up(:,1,1,1), &
               mp(:,1,1,1), ngu, ngf)
       case (2)
          call stencil_all_flux_2d(sp(:,:,1,:), fp(:,:,1,1), up(:,:,1,1), &
               mp(:,:,1,1), ngu, ngf, dim)
       case (3)
          call stencil_all_flux_3d(sp(:,:,:,:), fp(:,:,:,1), up(:,:,:,1), &
               mp(:,:,:,1), ngu, ngf, dim)
       end select
     end do
    end do

    call destroy(bpt)

  end subroutine ml_fill_all_fluxes

  subroutine ml_fine_contrib(flux, res, mm, ratio, crse_domain, side)
    use bl_prof_module
    type(multifab), intent(inout) :: flux
    type(multifab), intent(inout) :: res
    type(imultifab), intent(in) :: mm
    type(box) :: crse_domain
    type(box) :: fbox
    integer :: side
    integer :: ratio(:)
    integer :: lof(flux%dim)
    integer :: lo_dom(flux%dim), hi_dom(flux%dim)
    integer :: i, n, dir
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer :: nc
    type(bl_prof_timer), save :: bpt

    call build(bpt, "ml_fine_contrib")

    nc = res%nc

    if ( res%nc /= flux%nc ) then
       call bl_error("ML_FILL_FLUXES: res%nc /= flux%nc")
    end if

    lo_dom = lwb(crse_domain)
    hi_dom = upb(crse_domain)
    if ( nodal_q(res) ) hi_dom = hi_dom + 1
    dir = iabs(side)

    do i = 1, flux%nboxes
       if ( remote(flux, i) ) cycle
       fbox   = get_ibox(flux,i)
       lof = lwb(fbox)
       fp => dataptr(flux, i)
       rp => dataptr(res, i)
       mp => dataptr(mm, i)
       do n = 1, nc
          if ( (res%la%lap%pmask(dir)) .or. &
               (lof(dir) /= lo_dom(dir) .and. lof(dir) /= hi_dom(dir)) ) then
             select case(flux%dim)
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
    real (kind = dp_t) :: fac, fac0, fac1, fac2
    real (kind = dp_t) :: corner_fac
    logical ll1,ll2,ll3,lh1,lh2,lh3

    nx = size(mm,dim=1)-1
    ny = size(mm,dim=2)-1
    nz = size(mm,dim=3)-1

    nxc = size(dd,dim=1)
    nyc = size(dd,dim=2)
    nzc = size(dd,dim=3)

    hid(1) = lod(1) + nxc-1
    hid(2) = lod(2) + nyc-1
    hid(3) = lod(3) + nzc-1

    lo(:) = 0
    ng_res = 1
    call impose_neumann_bcs_3d(res,mm,lo,ng_res)

    if (side == -1 .or. side == 1) then

      if (side == -1) then
         i  = 0
         isign =  1
      else
         i  = nx
         isign = -1
      end if

      ic = lod(1)
      fac0 = 1.0_dp_t / (ratio(2)*ratio(3))
      !
      ! First average along the coarse-fine face.
      !
      !$OMP PARALLEL DO PRIVATE(jc,kc,n,fac,fac2,j,l,k) IF((hid(3)-lod(3)).ge.3)
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
                     if (k < nz) then
                        dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                     end if

                     k = (kc-lod(3))*ratio(3) - l
                     if (k >  0) then
                        dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                     end if
                  end do
               end if

               j = (jc-lod(2))*ratio(2) - n
               if (j > 0) then
                  do l = 0,ratio(3)-1
                     fac = (ratio(3)-l)*fac2
                     if (l == 0) fac = HALF * fac

                     k = (kc-lod(3))*ratio(3) + l
                     if (k < nz) then
                        dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                     end if

                     k = (kc-lod(3))*ratio(3) - l
                     if (k >  0) then
                        dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                     end if
                  end do
               end if

            end do
         end do
      end do
      !$OMP END PARALLEL DO

      jc = lod(2)
      kc = lod(3)
      j = 0
      k = 0
      if (.not. bc_neumann(mm(i,j,k),2,-1)) then
          if (.not. bc_neumann(mm(i,j,k),3,-1)) &
               dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t
       end if

      jc = hid(2)
      kc = lod(3)
      j = ny
      k = 0
      if (.not. bc_neumann(mm(i,j,k),2,+1)) then
         if (.not. bc_neumann(mm(i,j,k),3,-1)) &
              dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t
      end if

      jc = lod(2)
      kc = hid(3)
      j = 0
      k = nz
      if (.not. bc_neumann(mm(i,j,k),2,-1)) then
         if (.not. bc_neumann(mm(i,j,k),3,+1)) &
              dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t
      end if

      jc = hid(2)
      kc = hid(3)
      j = ny
      k = nz
      if (.not. bc_neumann(mm(i,j,k),2,+1)) then
         if (.not. bc_neumann(mm(i,j,k),3,+1)) &
              dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t
      end if

      j = 0
      do kc = lod(3),hid(3)
         k = (kc-lod(3))*ratio(3)
         if (bc_neumann(mm(i,j,k),2,-1)) dd(ic,lod(2),kc) = TWO*dd(ic,lod(2),kc)
      end do

      j = (hid(2)-lod(2))*ratio(2)
      do kc = lod(3),hid(3)
         k = (kc-lod(3))*ratio(3)
         if (bc_neumann(mm(i,j,k),2, 1)) dd(ic,hid(2),kc) = TWO*dd(ic,hid(2),kc)
      end do

      k = 0
      do jc = lod(2),hid(2)
         j = (jc-lod(2))*ratio(2)
         if (bc_neumann(mm(i,j,k),3,-1)) dd(ic,jc,lod(3)) = TWO*dd(ic,jc,lod(3))
      end do

      k = (hid(3)-lod(3))*ratio(3)
      do jc = lod(2),hid(2)
         j = (jc-lod(2))*ratio(2)
         if (bc_neumann(mm(i,j,k),3, 1)) dd(ic,jc,hid(3)) = TWO*dd(ic,jc,hid(3))
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
            !$OMP PARALLEL DO PRIVATE(kc,k,jc,j,jtop,jbot,kup,kdwn) &
            !$OMP PRIVATE(ll1,lh1,ll2,lh2,ll3,lh3,corner_fac) IF((hid(3)-lod(3)).ge.3)
            do kc = lod(3),hid(3)
              k = (kc-lod(3))*ratio(3)
              do jc = lod(2),hid(2)
                j = (jc-lod(2))*ratio(2)
                jtop = j+n
                jbot = j-n
                kup  = k+l
                kdwn = k-l
                if (j==0) then
                   if (bc_neumann(mm(i,j,k),2,-1)) jbot = jtop
                end if
                if (j==ny) then
                   if (bc_neumann(mm(i,j,k),2,+1)) jtop = jbot
                end if
                if (k==0) then
                   if (bc_neumann(mm(i,j,k),3,-1)) kdwn = kup
                end if
                if (k==nz) then
                   if (bc_neumann(mm(i,j,k),3,+1)) kup  = kdwn
                end if

                ll2 = .false.
                lh2 = .false.
                ll3 = .false.
                lh3 = .false.

                if (jc==lod(2)) then
                   if (.not. bc_neumann(mm(i,j,k),2,-1)) ll2 = .true.
                end if
                if (jc==hid(2)) then
                   if (.not. bc_neumann(mm(i,j,k),2,+1)) lh2 = .true.
                end if
                if (kc==lod(3)) then
                   if (.not. bc_neumann(mm(i,j,k),3,-1)) ll3 = .true.
                end if
                if (kc==hid(3)) then
                   if (.not. bc_neumann(mm(i,j,k),3,+1)) lh3 = .true.
                end if

                if ( ( ll2 .or. lh2 ) .and. ( ll3 .or. lh3 ) ) then
                   corner_fac = 1.0_dp_t / 3.0_dp_t
                else if ( ( ll2 .or. lh2 ) .and. .not. ( ll3 .or. lh3 ) ) then
                   corner_fac = 1.0_dp_t  / 2.0_dp_t
                else if ( .not. ( ll2 .or. lh2 ) .and. ( ll3 .or. lh3 ) ) then
                   corner_fac = 1.0_dp_t / 2.0_dp_t
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
                   if ( ll3 .and. lh3 ) &
                        dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                        fac*res(ioff,jbot,kdwn) 
                   if ( ll1 .and. lh1 ) &
                        dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                        fac*res(ioff,jbot,kup) 
                end if

                ll2 = (j+n >  0); if (.not. ll2) ll2 = bc_neumann(mm(i,j,k),2,-1)
                lh2 = (j+n < ny); if (.not. lh2) lh2 = bc_neumann(mm(i,j,k),2,+1)

                if ( ll2 .and. lh2 ) then
                   if ( ll3 .and. lh3 ) &
                        dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                        fac*res(ioff,jtop,kdwn) 
                   if ( ll1 .and. lh1 ) &
                        dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                        fac*res(ioff,jtop,kup) 
                end if

              end do
            end do
            !$OMP END PARALLEL DO
          end do
        end do
      end do

      !$OMP PARALLEL DO PRIVATE(jc,kc,j,k) IF((hid(3)-lod(3)).ge.3)
      do kc = lod(3),hid(3)
         do jc = lod(2),hid(2)
            j = (jc-lod(2))*ratio(2)
            k = (kc-lod(3))*ratio(3)
            if (.not.bc_dirichlet(mm(i,j,k),1,0)) dd(ic,jc,kc) = ZERO
         end do
      end do
      !$OMP END PARALLEL DO

    else if (side == -2 .or. side == 2) then

      if (side == -2) then
         j  = 0
         isign =  1
      else
         j  = ny
         isign = -1
      end if
      jc = lod(2)
      fac0 = 1.0_dp_t / (ratio(1)*ratio(3))
      !
      ! First average along the coarse-fine face.
      !
      !$OMP PARALLEL DO PRIVATE(kc,ic,n,fac2,i,l,fac,k) IF((hid(3)-lod(3)).ge.3)
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
                     if (k < nz) then
                        dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                     end if

                     k = (kc-lod(3))*ratio(3) - l
                     if (k >  0) then
                        dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                     end if
                  end do
               end if

               i = (ic-lod(1))*ratio(1) - n
               if (i > 0) then
                  do l = 0,ratio(3)-1
                     fac = (ratio(3)-l)*fac2
                     if (l == 0) fac = HALF * fac

                     k = (kc-lod(3))*ratio(3) + l
                     if (k < nz) then
                        dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                     end if

                     k = (kc-lod(3))*ratio(3) - l
                     if (k >  0) then
                        dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                     end if
                  end do
               end if
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      ic = lod(1)
      kc = lod(3)
      i = 0
      k = 0
      if (.not. bc_neumann(mm(i,j,k),1,-1) ) then
         if (.not. bc_neumann(mm(i,j,k),3,-1) ) &
              dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t
      end if

      ic = hid(1)
      kc = lod(3)
      i = nx
      k = 0
      if (.not. bc_neumann(mm(i,j,k),1,+1) ) then
         if (.not. bc_neumann(mm(i,j,k),3,-1) ) &
              dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t
      end if

      ic = lod(1)
      kc = hid(3)
      i = 0
      k = nz
      if (.not. bc_neumann(mm(i,j,k),1,-1) ) then
         if (.not. bc_neumann(mm(i,j,k),3,+1) ) &
              dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t
      end if

      ic = hid(1)
      kc = hid(3)
      i = nx
      k = nz
      if (.not. bc_neumann(mm(i,j,k),1,+1) ) then
         if (.not. bc_neumann(mm(i,j,k),3,+1) ) &
              dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t
      end if

      i = 0
      do kc = lod(3),hid(3)
         k = (kc-lod(3))*ratio(3)
         if (bc_neumann(mm(i,j,k),1,-1)) dd(lod(1),jc,kc) = TWO*dd(lod(1),jc,kc)
      end do

      i = (hid(1)-lod(1))*ratio(1)
      do kc = lod(3),hid(3)
         k = (kc-lod(3))*ratio(3)
         if (bc_neumann(mm(i,j,k),1, 1)) dd(hid(1),jc,kc) = TWO*dd(hid(1),jc,kc)
      end do

      k = 0
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         if (bc_neumann(mm(i,j,k),3,-1)) dd(ic,jc,lod(3)) = TWO*dd(ic,jc,lod(3))
      end do

      k = (hid(3)-lod(3))*ratio(3)
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         if (bc_neumann(mm(i,j,k),3, 1)) dd(ic,jc,hid(3)) = TWO*dd(ic,jc,hid(3))
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
            !$OMP PARALLEL DO PRIVATE(kc,k,ic,i,irght,ileft,kup,kdwn) &
            !$OMP PRIVATE(ll1,lh1,ll2,lh2,ll3,lh3,corner_fac) IF((hid(3)-lod(3)).ge.3)
            do kc = lod(3),hid(3)
              k = (kc-lod(3))*ratio(3)
              do ic = lod(1),hid(1)
                i = (ic-lod(1))*ratio(1)
                irght = i+n
                ileft = i-n
                kup  = k+l
                kdwn = k-l
                if (i==0) then
                   if (bc_neumann(mm(i,j,k),1,-1)) ileft = irght
                end if
                if (i==nx) then
                   if (bc_neumann(mm(i,j,k),1,+1)) irght = ileft
                end if
                if (k==0) then
                   if (bc_neumann(mm(i,j,k),3,-1)) kdwn = kup
                end if
                if (k==nz) then
                   if (bc_neumann(mm(i,j,k),3,+1)) kup  = kdwn
                end if

                ll1 = .false.
                lh1 = .false.
                ll3 = .false.
                lh3 = .false.

                if (ic==lod(1)) then
                   if (.not. bc_neumann(mm(i,j,k),1,-1)) ll1 = .true.
                end if
                if (ic==hid(1)) then
                   if (.not. bc_neumann(mm(i,j,k),1,+1)) lh1 = .true.
                end if
                if (kc==lod(3)) then
                   if (.not. bc_neumann(mm(i,j,k),3,-1)) ll3 = .true.
                end if
                if (kc==hid(3)) then
                   if (.not. bc_neumann(mm(i,j,k),3,+1)) lh3 = .true.
                end if

                if ( (  ll1 .or. lh1 ) .and. (  ll3 .or. lh3 ) ) then
                   corner_fac = 1.0_dp_t / 3.0_dp_t
                else if ( ( ll1 .or. lh1 ) .and. .not. ( ll3 .or. lh3 ) ) then
                   corner_fac = 1.0_dp_t  / 2.0_dp_t
                else if ( .not. &
                          ( ll1 .or. lh1 ) .and. ( ll3 .or. lh3 ) ) then
                   corner_fac = 1.0_dp_t / 2.0_dp_t
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
                   if ( ll2 .and. lh2 ) &
                        dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                        fac*res(ileft,joff,kdwn) 
                   if ( ll3 .and. lh3 ) &
                        dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                        fac*res(ileft,joff,kup) 
                end if

                ll1 = (i+n >  0); if (.not. ll1) ll1 = bc_neumann(mm(i,j,k),1,-1)
                lh1 = (i+n < nx); if (.not. lh1) lh1 = bc_neumann(mm(i,j,k),1,+1)

                if ( ll1 .and. lh1 ) then
                   if ( ll2 .and. lh2 ) &
                        dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                        fac*res(irght,joff,kdwn) 
                   if ( ll3 .and. lh3 ) &
                        dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                        fac*res(irght,joff,kup) 
                end if
              end do
            end do
            !$OMP END PARALLEL DO
          end do
        end do
      end do

      do kc = lod(3),hid(3)
         do ic = lod(1),hid(1)
            i = (ic-lod(1))*ratio(1)
            k = (kc-lod(3))*ratio(3)
            if (.not.bc_dirichlet(mm(i,j,k),1,0)) dd(ic,jc,kc) = ZERO
         end do
      end do

    else 

      if (side == -3) then
         k  = 0
         isign =  1
      else
         k  = nz
         isign = -1
      end if
      kc = lod(3)
      fac0 = 1.0_dp_t / (ratio(1)*ratio(2))
      !
      ! First average along the coarse-fine face.
      !
      !$OMP PARALLEL DO PRIVATE(jc,ic,n,fac2,i,l,fac,j) IF((hid(2)-lod(2)).ge.3)
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
                     if (j < ny) then
                        dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                     end if

                     j = (jc-lod(2))*ratio(2) - l
                     if (j >  0) then
                        dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                     end if
                  end do
               end if

               i = (ic-lod(1))*ratio(1) - n
               if (i > 0) then
                  do l = 0,ratio(2)-1
                     fac = (ratio(2)-l)*fac2
                     if (l == 0) fac = HALF * fac

                     j = (jc-lod(2))*ratio(2) + l
                     if (j < ny) then
                        dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                     end if

                     j = (jc-lod(2))*ratio(2) - l
                     if (j >  0) then
                        dd(ic,jc,kc) = dd(ic,jc,kc) + fac * res(i,j,k)
                     end if
                  end do
               end if
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      ic = lod(1)
      jc = lod(2)
      i = 0
      j = 0
      if (.not. bc_neumann(mm(i,j,k),1,-1) ) then
         if (.not. bc_neumann(mm(i,j,k),2,-1) ) &
              dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t
      end if

      ic = hid(1)
      jc = lod(2)
      i = nx
      j = 0
      if (.not. bc_neumann(mm(i,j,k),1,+1) ) then
         if (.not. bc_neumann(mm(i,j,k),2,-1) ) &
              dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t
      end if

      ic = lod(1)
      jc = hid(2)
      i = 0
      j = ny
      if (.not. bc_neumann(mm(i,j,k),1,-1) ) then
         if (.not. bc_neumann(mm(i,j,k),2,+1) ) &
              dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t
      end if

      ic = hid(1)
      jc = hid(2)
      i = nx
      j = ny
      if (.not. bc_neumann(mm(i,j,k),1,+1) ) then
         if (.not. bc_neumann(mm(i,j,k),2,+1) ) &
              dd(ic,jc,kc) = dd(ic,jc,kc) + 0.25_dp_t * res(i,j,k) / 3.0_dp_t
      end if

      i = 0
      do jc = lod(2),hid(2)
         j = (jc-lod(2))*ratio(2)
         if (bc_neumann(mm(i,j,k),1,-1)) dd(lod(1),jc,kc) = TWO*dd(lod(1),jc,kc)
      end do

      i = (hid(1)-lod(1))*ratio(1)
      do jc = lod(2),hid(2)
         j = (jc-lod(2))*ratio(2)
         if (bc_neumann(mm(i,j,k),1,+1)) dd(hid(1),jc,kc) = TWO*dd(hid(1),jc,kc)
      end do

      j = 0
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         if (bc_neumann(mm(i,j,k),2,-1)) dd(ic,lod(2),kc) = TWO*dd(ic,lod(2),kc)
      end do

      j = (hid(2)-lod(2))*ratio(2)
      do ic = lod(1),hid(1)
         i = (ic-lod(1))*ratio(1)
         if (bc_neumann(mm(i,j,k),2,+1)) dd(ic,hid(2),kc) = TWO*dd(ic,hid(2),kc)
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
            !$OMP PARALLEL DO PRIVATE(jc,j,ic,i,irght,ileft,jtop,jbot) &
            !$OMP PRIVATE(ll1,lh1,ll2,lh2,ll3,lh3,corner_fac)  IF((hid(2)-lod(2)).ge.3)
            do jc = lod(2),hid(2)
              j = (jc-lod(2))*ratio(2)
              do ic = lod(1),hid(1)
                i = (ic-lod(1))*ratio(1)
                irght = i+n
                ileft = i-n
                jtop  = j+l
                jbot  = j-l
                if (i==0) then
                   if (bc_neumann(mm(i,j,k),1,-1)) ileft = irght
                end if
                if (i==nx) then
                   if (bc_neumann(mm(i,j,k),1,+1)) irght = ileft
                end if
                if (j==0) then
                   if (bc_neumann(mm(i,j,k),2,-1)) jbot  = jtop
                end if
                if (j==ny) then
                   if (bc_neumann(mm(i,j,k),2,+1)) jtop  = jbot
                end if

                ll1 = .false.
                lh1 = .false.
                ll2 = .false.
                lh2 = .false.

                if (ic==lod(1)) then
                   if (.not. bc_neumann(mm(i,j,k),1,-1)) ll1 = .true.
                end if
                if (ic==hid(1)) then
                   if (.not. bc_neumann(mm(i,j,k),1,+1)) lh1 = .true.
                end if
                if (jc==lod(2)) then
                   if (.not. bc_neumann(mm(i,j,k),2,-1)) ll2 = .true.
                end if
                if (jc==hid(2)) then
                   if (.not. bc_neumann(mm(i,j,k),2,+1)) lh2 = .true.
                end if

                if ( ( ll1 .or. lh1 ) .and. ( ll2 .or. lh2 ) ) then
                   corner_fac = 1.0_dp_t / 3.0_dp_t
                else if ( ( ll1 .or. lh1 ) .and. .not. ( ll2 .or. lh2 ) ) then
                   corner_fac = 1.0_dp_t  / 2.0_dp_t
                else if ( .not. &
                          ( ll1 .or. lh1 ) .and. ( ll2 .or. lh2 ) ) then
                   corner_fac = 1.0_dp_t / 2.0_dp_t
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
                   if ( ll2 .and. lh2 ) &
                        dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                        fac*res(ileft,jbot,koff) 
                   if ( ll3 .and. lh3 ) &
                        dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                        fac*res(ileft,jtop,koff) 
                end if

                ll1 = (i+n >  0); if (.not. ll1) ll1 = bc_neumann(mm(i,j,k),1,-1) 
                lh1 = (i+n < nx); if (.not. lh1) lh1 = bc_neumann(mm(i,j,k),1,+1)

                if ( ll1 .and. lh1 ) then
                   if ( ll2 .and. lh2 ) &
                        dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                        fac*res(irght,jbot,koff) 
                   if ( ll3 .and. lh3 ) &
                        dd(ic,jc,kc) = dd(ic,jc,kc) + corner_fac * &
                        fac*res(irght,jtop,koff) 
                end if

              end do
            end do
            !$OMP END PARALLEL DO
          end do
        end do
      end do

      do jc = lod(2),hid(2)
         do ic = lod(1),hid(1)
            i = (ic-lod(1))*ratio(1)
            j = (jc-lod(2))*ratio(2)
            if (.not. bc_dirichlet(mm(i,j,k),1,0)) dd(ic,jc,kc) = ZERO
         end do
      end do

    end if

  end subroutine fine_edge_resid_3d

  function ml_norm_inf(mf, mask) result(r)
    type( multifab), intent(in) :: mf(:)
    type(lmultifab), intent(in) :: mask(:)
    real(dp_t)                  :: r, r1
    integer                     :: n,nlevs
    nlevs = size(mf)
    r = norm_inf(mf(nlevs))
    do n = nlevs-1, 1, -1
       r1 = norm_inf(mf(n), mask(n))
       r = max(r1, r)
    end do
  end function ml_norm_inf

  function ml_norm_l2(mf, rr, mask) result(r)
    type( multifab), intent(in) :: mf(:)
    integer                     :: rr(:,:)
    type(lmultifab), intent(in) :: mask(:)
    real(dp_t)                  :: r
    integer                     :: n,nlevs
    nlevs = size(mf)
    r = norm_l2(mf(nlevs))**2
    do n = nlevs-1, 1, -1
       r =  r / product(rr(n,:)) &
          + norm_l2(mf(n), mask = mask(n))**2
    end do
    r = sqrt(r)
  end function ml_norm_l2

end module ml_util_module
