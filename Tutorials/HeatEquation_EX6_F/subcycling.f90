module subcycling_module

  use ml_layout_module
  use multifab_module
  use ml_restriction_module
  use multifab_fill_ghost_module
  use bndry_reg_module

  implicit none

  private

  public :: update_isrefined, apply_fluxfix_single_level, &
            flux_bndry_reg_copy_c,    &
            fine_flx_restriction_c,   &
            filter_assign, &
            filter_mean, &
            filter_quarter_of_diff, &
            filter_add, &
            filter_subtract_half, &
            filter_subtract_quarter

contains
  subroutine update_isrefined(mla,isrefined,nlevs)

    type(ml_layout), intent(inout) :: mla
    type(lmultifab) , intent(inout) :: isrefined(:)
    integer        , intent(inout) :: nlevs
    integer :: dm,n,ifab
    type(lmultifab),allocatable :: ischild(:)
    dm = mla%dim

    !prepare help data
    allocate(ischild(nlevs))
    do n=1,nlevs
     call lmultifab_build(ischild(n),mla%la(n),1,0)
    end do

    !#1 set isrefined=.false. to all cells
    do n=1,nlevs
       call setval(isrefined(n),.FALSE.,all=.TRUE.)
    enddo

    !#2 set isrefined=.false. to all ghost elements
    ! is done by all=.true. above

    !#3 identify each refined element and set isrefined=.true. where necessary
    ! elements from levels 2:nlevs are child elements
    call setval(ischild(1),.FALSE.,all=.TRUE.)

    do n=2,nlevs
      call setval(ischild(n),.TRUE.,all=.TRUE.)
    enddo

    ! all parent elements are refined elements
    do n=nlevs,2,-1
      call lml_cc_restriction(isrefined(n-1),ischild(n),mla%mba%rr(n-1,:))
    enddo
    !#3 end

    !#4 fill ghost cells only at the same level of refinement
    do n=1,nlevs
      ! fill ghost cells for two adjacent grids at the same level
      ! let elements know if neighbor is refined or not
      call lmultifab_fill_boundary(isrefined(n))
    enddo

    !destroy help data
    do n=1,nlevs
     call destroy(ischild(n))
    end do

  end subroutine update_isrefined

  subroutine apply_fluxfix_single_level(n,mla,phi,flux,isrefined,dx,dt)

    integer        , intent(in   ) :: n
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi(:)
    type(multifab) , intent(in   ) :: flux(:,:)
    type(lmultifab), intent(in   ) :: isrefined(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, ng_f, i, nlevs

    real(kind=dp_t), pointer ::  pp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)
    logical        , pointer :: isrefp(:,:,:,:)

    dm    = mla%dim
    nlevs = mla%nlevel

    ng_p = phi(1)%ng
    ng_f = flux(1,1)%ng

       do i=1,nfabs(phi(n))
          pp  => dataptr(phi(n),i)
          fxp => dataptr(flux(n,1),i)
          fyp => dataptr(flux(n,2),i)
          isrefp => dataptr(isrefined(n),i)
          lo = lwb(get_box(phi(n),i))
          hi = upb(get_box(phi(n),i))
          select case(dm)
          case (2)
             call apply_fluxfix_2d(pp(:,:,1,1), ng_p, &
                                fxp(:,:,1,1), fyp(:,:,1,1),&
                                isrefp(:,:,1,1), ng_f, &
                                lo, hi, dx(n), dt)
          case (3)
             fzp => dataptr(flux(n,3),i)
             call apply_fluxfix_3d(pp(:,:,:,1), ng_p, &
                                fxp(:,:,:,1), fyp(:,:,:,1), fzp(:,:,:,1),&
                                isrefp(:,:,:,1), ng_f, &
                                lo, hi, dx(n), dt)
          end select
       end do

  end subroutine apply_fluxfix_single_level

  subroutine apply_fluxfix_2d(phi, ng_p, fluxx, fluxy, isref, ng_f, lo, hi, dx, dt)

    integer          :: lo(2), hi(2), ng_p, ng_f
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:)
    logical          :: isref(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision :: dx, dt
    ! local variables
    integer i,j,ifab
    double precision :: ffixterm

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
         ! applies only to NOT refined elements
         if (.not.(isref(i,j))) then
            ffixterm = 0.0d0
            ! left  neighbor along x refined
            if (isref(i-1,j)) then
              ffixterm = ffixterm  + fluxx(i,j)
            endif
            ! right neighbor along x refined
            if (isref(i+1,j)) then
              ffixterm = ffixterm - fluxx(i+1,j)
            endif
            ! left  neighbor along y refined
            if (isref(i,j-1)) then
              ffixterm = ffixterm + fluxy(i,j)
            endif
            ! right neighbor along y refined
            if (isref(i,j+1)) then
              ffixterm = ffixterm - fluxy(i,j+1)
            endif
            if (ffixterm.ne.0) then
              phi(i,j) = phi(i,j) + dt * ffixterm / dx
            endif
         endif
       end do
    end do
  end subroutine apply_fluxfix_2d

  subroutine apply_fluxfix_3d(phi, ng_p, fluxx, fluxy, fluxz, isref, ng_f, lo, hi, dx, dt)

    integer          :: lo(3), hi(3), ng_p, ng_f
    double precision ::   phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    logical          :: isref(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision :: dx, dt
    ! local variables
    integer i,j,k,ifab
    double precision :: ffixterm
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
           ! applies only to NOT refined elements
           if (.not.(isref(i,j,k))) then
              ffixterm = 0.0d0
              ! left  neighbor along x refined
              if (isref(i-1,j,k)) then
                ffixterm = ffixterm  + fluxx(i,j,k)
              endif
              ! right neighbor along x refined
              if (isref(i+1,j,k)) then
                ffixterm = ffixterm - fluxx(i+1,j,k)
              endif
              ! left  neighbor along y refined
              if (isref(i,j-1,k)) then
                ffixterm = ffixterm + fluxy(i,j,k)
              endif
              ! right neighbor along y refined
              if (isref(i,j+1,k)) then
                ffixterm = ffixterm - fluxy(i,j+1,k)
              endif
              ! left  neighbor along z refined
              if (isref(i,j,k-1)) then
                ffixterm = ffixterm + fluxz(i,j,k)
              endif
              ! right neighbor along z refined
              if (isref(i,j,k+1)) then
                ffixterm = ffixterm - fluxz(i,j,k+1)
              endif
              if (ffixterm.ne.0) then
                phi(i,j,k) = phi(i,j,k) + dt * ffixterm / dx
              endif
           endif
         end do
      end do
    end do
    !$omp end parallel do
  end subroutine apply_fluxfix_3d

  subroutine flux_bndry_reg_copy_c(br, cb, mf, cm, nc, filter)
    type(multifab) , intent(in   ) :: mf(:)
    type(bndry_reg), intent(inout) :: br
    integer        , intent(in   ) :: cb, cm
    integer        , intent(in   ), optional :: nc
    integer :: i, f
    type(bl_prof_timer), save :: bpt
    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter

    call build(bpt, "fx_br_copy_c")
    do i = 1, br%dim
       do f = 0, 1
          call multifab_copy_c(br%bmf(i,f), cb, mf(i), cm, nc=nc, filter=filter)
       end do
    end do
    call destroy(bpt)

  end subroutine flux_bndry_reg_copy_c

  subroutine fine_flx_restriction_c(flux, cc, fine_flx, cf, ir, nc, filter)
    type(multifab) , intent(inout) :: flux(:)
    type(bndry_reg), intent(inout) :: fine_flx
    integer        , intent(in   ) :: cc, cf, ir(:)
    integer        , intent(in   ), optional  :: nc
    integer :: i, f
    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter

    do i = 1, fine_flx%dim
       do f = 0, 1
          call ml_br_edge_restriction_c(flux(i),cc,fine_flx%bmf(i,f),cf,ir,i,nc,&
                                        filter=filter)
       end do
    end do

  end subroutine fine_flx_restriction_c

  subroutine ml_br_edge_restriction_c(crse, cc, fine, cf, ir, face, nc, filter)
    use edge_restriction_module
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer       , intent(in   ) :: cc, cf, ir(:)
    integer       , intent(in   ) :: face
    integer, intent(in), optional :: nc
    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter

    integer             :: i, n, lnc, dm, len
    integer             :: lo(get_dim(fine)), hi(get_dim(fine)), loc(get_dim(fine)), lof(get_dim(fine))
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    type(layout)        :: lacfine
    type(layout)        :: laf
    type(multifab)      :: cfine
    logical             :: nodal(get_dim(fine))

    dm = get_dim(crse)

    lnc = 1; if ( present(nc) ) lnc = nc

    laf = get_layout(fine)

    call layout_build_coarse(lacfine, laf, ir)

    ! note that all multifabs in bndry_reg are NOT nodal, while coarse multiFab is.
    call multifab_build(cfine, lacfine, nc = ncomp(crse), ng = 0, nodal = nodal_flags(fine))

    !$OMP PARALLEL DO PRIVATE(i,n,lo,hi,loc,lof,fp,cp)
    do i = 1, nfabs(fine)
       lo  = lwb(get_ibox(cfine,i))
       hi  = upb(get_ibox(cfine,i))
       loc = lwb(get_pbox(cfine,i))
       lof = lwb(get_pbox(fine, i))
       do n = 1, lnc
          fp  => dataptr(fine,  i, n+cf-1, 1)
          cp  => dataptr(cfine, i, n,      1)
          select case (dm)
          case (1)
             call edge_restriction_1d(cp(:,1,1,1), loc, fp(:,1,1,1), lof, lo, hi, ir)
          case (2)
             call edge_restriction_2d(cp(:,:,1,1), loc, fp(:,:,1,1), lof, lo, hi, ir, face)
          case (3)
             call edge_restriction_3d(cp(:,:,:,1), loc, fp(:,:,:,1), lof, lo, hi, ir, face)
          end select
       enddo
    end do
    !$OMP END PARALLEL DO

    call copy(crse, cc, cfine, 1, lnc, filter=filter)

    call destroy(cfine)

  end subroutine ml_br_edge_restriction_c

  !---------------------
  ! filters
  !---------------------
  subroutine filter_assign(out, in)
     use bl_types
     real(dp_t), intent(inout) :: out(:,:,:,:)
     real(dp_t), intent(in   ) ::  in(:,:,:,:)
     integer                   :: i, j, k, n
     !
     ! out = out + in
     !
     do n = 1, size(out,4)
        !$omp parallel do private(i,j,k)
        do k = 1, size(out,3)
           do j = 1, size(out,2)
              do i = 1, size(out,1)
                 out(i,j,k,n) = in(i,j,k,n)
              end do
           end do
        end do
        !$omp end parallel do
     end do
  end subroutine filter_assign

  subroutine filter_mean(out, in)
     use bl_types
     real(dp_t), intent(inout) :: out(:,:,:,:)
     real(dp_t), intent(in   ) ::  in(:,:,:,:)
     integer                   :: i, j, k, n
     !
     ! out = out + in
     !
     do n = 1, size(out,4)
        !$omp parallel do private(i,j,k)
        do k = 1, size(out,3)
           do j = 1, size(out,2)
              do i = 1, size(out,1)
                 out(i,j,k,n) = out(i,j,k,n) + in(i,j,k,n)
                 out(i,j,k,n) = 0.5d0 * out(i,j,k,n)
              end do
           end do
        end do
        !$omp end parallel do
     end do
  end subroutine filter_mean

  subroutine filter_quarter_of_diff(out, in)
     use bl_types
     real(dp_t), intent(inout) :: out(:,:,:,:)
     real(dp_t), intent(in   ) ::  in(:,:,:,:)
     integer                   :: i, j, k, n
     !
     ! out = out + in
     !
     do n = 1, size(out,4)
        !$omp parallel do private(i,j,k)
        do k = 1, size(out,3)
           do j = 1, size(out,2)
              do i = 1, size(out,1)
                 out(i,j,k,n) = out(i,j,k,n) - in(i,j,k,n)
                 out(i,j,k,n) = 0.25d0 * out(i,j,k,n)
              end do
           end do
        end do
        !$omp end parallel do
     end do
  end subroutine filter_quarter_of_diff

  subroutine filter_add(out, in)
     use bl_types
     real(dp_t), intent(inout) :: out(:,:,:,:)
     real(dp_t), intent(in   ) ::  in(:,:,:,:)
     integer                   :: i, j, k, n
     !
     ! out = out + in
     !
     do n = 1, size(out,4)
        !$omp parallel do private(i,j,k)
        do k = 1, size(out,3)
           do j = 1, size(out,2)
              do i = 1, size(out,1)
                 out(i,j,k,n) = out(i,j,k,n) + in(i,j,k,n)
              end do
           end do
        end do
        !$omp end parallel do
     end do
  end subroutine filter_add

  subroutine filter_subtract_half(out, in)
     use bl_types
     real(dp_t), intent(inout) :: out(:,:,:,:)
     real(dp_t), intent(in   ) ::  in(:,:,:,:)
     integer                   :: i, j, k, n
     !
     ! out = out - in
     !
     do n = 1, size(out,4)
        !$omp parallel do private(i,j,k)
        do k = 1, size(out,3)
           do j = 1, size(out,2)
              do i = 1, size(out,1)
                 out(i,j,k,n) = out(i,j,k,n) - 0.5d0* in(i,j,k,n)
              end do
           end do
        end do
        !$omp end parallel do
     end do
  end subroutine filter_subtract_half

  subroutine filter_subtract_quarter(out, in)
     use bl_types
     real(dp_t), intent(inout) :: out(:,:,:,:)
     real(dp_t), intent(in   ) ::  in(:,:,:,:)
     integer                   :: i, j, k, n
     !
     ! out = out - in
     !
     do n = 1, size(out,4)
        !$omp parallel do private(i,j,k)
        do k = 1, size(out,3)
           do j = 1, size(out,2)
              do i = 1, size(out,1)
                 out(i,j,k,n) = out(i,j,k,n) - 0.25d0* in(i,j,k,n)
              end do
           end do
        end do
        !$omp end parallel do
     end do
  end subroutine filter_subtract_quarter

end module subcycling_module
