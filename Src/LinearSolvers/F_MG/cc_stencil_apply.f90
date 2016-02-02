module cc_stencil_apply_module

  use bl_constants_module
  use bl_types
  use bc_module
  use bc_functions_module
  use stencil_types_module
  use multifab_module

  implicit none

  private
  
  public  :: ml_fill_all_fluxes, &
       stencil_apply_1d, stencil_apply_2d, stencil_apply_3d, &
       stencil_flux_1d, stencil_flux_2d, stencil_flux_3d, &
       stencil_fine_flux_1d, stencil_fine_flux_2d, stencil_fine_flux_3d, &
       stencil_apply_ibc_2d, stencil_apply_ibc_3d

contains

  subroutine ml_fill_all_fluxes(ss, flux, uu, mm)

    use bl_prof_module
    use multifab_module
    use stencil_util_module, only : is_ibc_stencil

    type( multifab), intent(in   ) :: ss
    type( multifab), intent(inout) :: flux(:)
    type( multifab), intent(inout) :: uu
    type(imultifab), intent(in   ) :: mm

    integer :: dim, i, ngu, ngf, ndims

    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)

    type(bl_prof_timer), save :: bpt
    call build(bpt, "ml_fill_all_fluxes")

    ngu = nghost(uu)
    ndims = get_dim(uu)

    if ( ncomp(uu) /= ncomp(flux(1)) ) then
       call bl_error("ML_FILL_ALL_FLUXES: uu%nc /= flux%nc")
    end if

    call multifab_fill_boundary(uu)

    do dim = 1, ndims
       do i = 1, nfabs(flux(dim))
          ngf = nghost(flux(dim))
          fp => dataptr(flux(dim), i)
          up => dataptr(uu, i)
          sp => dataptr(ss, i)

          if (is_ibc_stencil(ss,i)) then
             select case (ndims)
             case (2)
                call stencil_all_flux_ibc_2d(sp(dim+1,1,1,1), fp(:,:,1,1), ngf, up(:,:,1,1), ngu, dim)
             case (3)
                call stencil_all_flux_ibc_3d(sp(dim+1,1,1,1), fp(:,:,:,1), ngf, up(:,:,:,1), ngu, dim)
             end select
          else
             mp => dataptr(mm, i)
             select case(ndims)
             case (1)
                call stencil_all_flux_1d(sp(:,:,1,1), fp(:,1,1,1), up(:,1,1,1), &
                                         mp(:,1,1,1), ngu, ngf)
             case (2)
                call stencil_all_flux_2d(sp(:,:,:,1), fp(:,:,1,1), up(:,:,1,1), &
                                         mp(:,:,1,1), ngu, ngf, dim)
             case (3)
                call stencil_all_flux_3d(sp(:,:,:,:), fp(:,:,:,1), up(:,:,:,1), &
                                         mp(:,:,:,1), ngu, ngf, dim)
             end select
          end if
       end do
    end do

    call destroy(bpt)

  end subroutine ml_fill_all_fluxes


  subroutine stencil_apply_1d(ss, dd, ng_d, uu, ng_u, mm, lo, hi, skwd)

    integer, intent(in) :: ng_d, ng_u, lo(:), hi(:)
    real (kind = dp_t), intent(in)  :: ss(0:,lo(1) :)
    real (kind = dp_t), intent(out) :: dd(lo(1)-ng_d:)
    real (kind = dp_t), intent(in)  :: uu(lo(1)-ng_u:)
    integer           , intent(in)  :: mm(lo(1):)
    logical           , intent(in), optional   :: skwd

    integer, parameter :: XBC = 3
    logical :: lskwd
    integer :: i
   
    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    do i = lo(1),hi(1)
       dd(i) = ss(0,i)*uu(i) + ss(1,i)*uu(i+1) + ss(2,i)*uu(i-1)
    end do

    if ( lskwd ) then
       if (hi(1) > lo(1)) then
          i = lo(1)
          if (bc_skewed(mm(i),1,+1)) then
             dd(i) = dd(i) + ss(XBC,i)*uu(i+2)
          end if
  
          i = hi(1)
          if (bc_skewed(mm(i),1,-1)) then
             dd(i) = dd(i) + ss(XBC,i)*uu(i-2)
          end if
       end if
    end if

  end subroutine stencil_apply_1d

  subroutine stencil_flux_1d(ss, flux, uu, mm, ng, ratio, face, dim, skwd)

    integer, intent(in) :: ng
    real (kind = dp_t), intent(in)  :: ss(0:,:)
    real (kind = dp_t), intent(out) :: flux(:)
    real (kind = dp_t), intent(in)  :: uu(1-ng:)
    integer           , intent(in)  :: mm(:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: ratio, face, dim
    integer nx
    integer i
    integer, parameter :: XBC = 3

    real (kind = dp_t) :: fac

    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=2)

    !   This factor is dx^fine / dx^crse
    fac = ONE / real(ratio, kind=dp_t)

    if ( dim == 1 ) then
       if ( face == -1 ) then
          i = 1
          if (bc_dirichlet(mm(1),1,-1)) then
             flux(1) = ss(1,i)*(uu(i+1)-uu(i)) + ss(2,i)*(uu(i-1)-uu(i)) &
                  - ss(2,i+1)*(uu(i+1)-uu(i))
             if (bc_skewed(mm(i),1,+1)) then
                flux(1) =  flux(1) + ss(XBC,i)*uu(i+2)
             end if
          else 
             flux(1) = Huge(flux)
          end if
          flux(1) = fac*flux(1)
       else if ( face == 1 ) then
          i = nx
          if (bc_dirichlet(mm(i),1,+1)) then
             flux(1) = ss(1,i)*(uu(i+1)-uu(i)) + ss(2,i)*(uu(i-1)-uu(i)) &
                  - ss(1,i-1)*(uu(i-1)-uu(i))
             if (bc_skewed(mm(i),1,-1)) then
                flux(1) =  flux(1) + ss(XBC,i)*uu(i-2)
             end if
          else 
             flux(1) = Huge(flux)
          end if
          flux(1) = fac*flux(1)
       end if
    end if

  end subroutine stencil_flux_1d

  subroutine stencil_apply_2d(ss, dd, ng_d, uu, ng_u, mm, lo, hi, skwd)

    integer           , intent(in   ) :: ng_d, ng_u, lo(:), hi(:)
    real (kind = dp_t), intent(in   ) :: ss(0:,lo(1):,lo(2):)
    real (kind = dp_t), intent(  out) :: dd(lo(1)-ng_d:,lo(2)-ng_d:)
    real (kind = dp_t), intent(in   ) :: uu(lo(1)-ng_u:,lo(2)-ng_u:)
    integer           , intent(in   ) :: mm(lo(1):,lo(2):)
    logical           , intent(in   ), optional :: skwd

    integer i,j

    integer, parameter :: XBC = 5, YBC = 6

    logical :: lskwd

    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    ! This is the Minion 4th order cross stencil.
    if (size(ss,dim=1) .eq. 9) then

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            dd(i,j) = &
                   ss(0,i,j) * uu(i,j) &
                 + ss(1,i,j) * uu(i-2,j) + ss(2,i,j) * uu(i-1,j) &
                 + ss(3,i,j) * uu(i+1,j) + ss(4,i,j) * uu(i+2,j) &
                 + ss(5,i,j) * uu(i,j-2) + ss(6,i,j) * uu(i,j-1) &
                 + ss(7,i,j) * uu(i,j+1) + ss(8,i,j) * uu(i,j+2)
          end do
       end do

    ! This is the Minion 4th order full stencil.
    else if (size(ss,dim=1) .eq. 25) then

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            dd(i,j) = ss( 0,i,j) * uu(i,j) &
                    + ss( 1,i,j) * uu(i-2,j-2) + ss( 2,i,j) * uu(i-1,j-2) & ! AT J-2
                    + ss( 3,i,j) * uu(i  ,j-2) + ss( 4,i,j) * uu(i+1,j-2) & ! AT J-2
                    + ss( 5,i,j) * uu(i+2,j-2)                            & ! AT J-2
                    + ss( 6,i,j) * uu(i-2,j-1) + ss( 7,i,j) * uu(i-1,j-1) & ! AT J-1
                    + ss( 8,i,j) * uu(i  ,j-1) + ss( 9,i,j) * uu(i+1,j-1) & ! AT J-1
                    + ss(10,i,j) * uu(i+2,j-1)                            & ! AT J-1
                    + ss(11,i,j) * uu(i-2,j  ) + ss(12,i,j) * uu(i-1,j  ) & ! AT J
                    + ss(13,i,j) * uu(i+1,j  ) + ss(14,i,j) * uu(i+2,j  ) & ! AT J
                    + ss(15,i,j) * uu(i-2,j+1) + ss(16,i,j) * uu(i-1,j+1) & ! AT J+1
                    + ss(17,i,j) * uu(i  ,j+1) + ss(18,i,j) * uu(i+1,j+1) & ! AT J+1
                    + ss(19,i,j) * uu(i+2,j+1)                            & ! AT J+1
                    + ss(20,i,j) * uu(i-2,j+2) + ss(21,i,j) * uu(i-1,j+2) & ! AT J+2
                    + ss(22,i,j) * uu(i  ,j+2) + ss(23,i,j) * uu(i+1,j+2) & ! AT J+2
                    + ss(24,i,j) * uu(i+2,j+2)                              ! AT J+2
          end do
       end do

    ! This is our standard 5-point Laplacian with a possible correction at boundaries
    else 

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             dd(i,j) = ss(0,i,j)*uu(i,j) &
                  + ss(1,i,j)*uu(i+1,j  ) + ss(2,i,j)*uu(i-1,j  ) &
                  + ss(3,i,j)*uu(i  ,j+1) + ss(4,i,j)*uu(i  ,j-1)
          end do
       end do

       if ( lskwd ) then
       ! Corrections for skewed stencils
       if (hi(1) > lo(1)) then
          do j = lo(2),hi(2)

             i = lo(1)
             if (bc_skewed(mm(i,j),1,+1)) then
                dd(i,j) = dd(i,j) + ss(XBC,i,j)*uu(i+2,j)
             end if

             i = hi(1)
             if (bc_skewed(mm(i,j),1,-1)) then
                dd(i,j) = dd(i,j) + ss(XBC,i,j)*uu(i-2,j)
             end if
          end do
       end if

       if (hi(2) > lo(2)) then
          do i = lo(1),hi(1)

             j = lo(2)
             if (bc_skewed(mm(i,j),2,+1)) then
                dd(i,j) = dd(i,j) + ss(YBC,i,j)*uu(i,j+2)
             end if

             j = hi(2)
             if (bc_skewed(mm(i,j),2,-1)) then
                dd(i,j) = dd(i,j) + ss(YBC,i,j)*uu(i,j-2)
             end if

          end do
       end if
       end if
    end if

  end subroutine stencil_apply_2d

  subroutine stencil_apply_ibc_2d(ss, dd, ng_d, uu, ng_u, lo, hi)

    integer           , intent(in   ) :: ng_d, ng_u, lo(:), hi(:)
    real (kind = dp_t), intent(in   ) :: ss(0:)
    real (kind = dp_t), intent(  out) :: dd(lo(1)-ng_d:,lo(2)-ng_d:)
    real (kind = dp_t), intent(in   ) :: uu(lo(1)-ng_u:,lo(2)-ng_u:)

    integer i,j

    ! This is our standard 5-point Laplacian without correction at boundaries
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          dd(i,j) = ss(0)*uu(i,j) &
               +    ss(1)*(uu(i-1,j) + uu(i+1,j)) &
               +    ss(2)*(uu(i,j-1) + uu(i,j+1))
       end do
    end do

  end subroutine stencil_apply_ibc_2d


  subroutine stencil_flux_2d(ss, flux, uu, mm, ng, ratio, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:)
    real (kind = dp_t), intent(out) :: flux(:,:)
    real (kind = dp_t), intent(in ) :: ss(0:,:,:)
    integer           , intent(in)  :: mm(:,:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: ratio, face, dim
    integer nx,ny
    integer i,j,ic,jc
    real (kind = dp_t) :: fac
    integer, parameter :: XBC = 5, YBC = 6
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=2)
    ny = size(ss,dim=3)

    !   Note that one factor of ratio is the tangential averaging, while the
    !     other is the normal factor
    fac = ONE/real(ratio*ratio, kind=dp_t)

!   Lo i face
    if ( dim == 1 ) then
       if (face == -1) then

          i = 1
          flux(1,:) = ZERO
          do j = 1,ny
             jc = (j-1)/ratio+1
             if (bc_dirichlet(mm(i,j),1,-1)) then
                flux(1,jc) = flux(1,jc)  &
                     + ss(1,i,j)*(uu(i+1,j)-uu(i,j)) &
                     + ss(2,i,j)*(uu(i-1,j)-uu(i,j)) - ss(2,i+1,j)*(uu(i+1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,+1)) &
                     flux(1,jc) = flux(1,jc) + ss(XBC,i,j)*(uu(i+2,j)-uu(i,j)) 
             else   
                flux(1,jc) = Huge(flux)
             end if
          end do
          flux(1,:) = fac * flux(1,:)

!      Hi i face
       else if (face == 1) then

          i = nx
          flux(1,:) = ZERO
          do j = 1,ny
             jc = (j-1)/ratio+1
             if (bc_dirichlet(mm(i,j),1,+1)) then

                flux(1,jc) = flux(1,jc) &
                     + ss(1,i,j)*(uu(i+1,j)-uu(i,j)) &
                     + ss(2,i,j)*(uu(i-1,j)-uu(i,j)) - ss(1,i-1,j)*(uu(i-1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,-1)) &
                     flux(1,jc) = flux(1,jc) + ss(XBC,i,j)*(uu(i-2,j)-uu(i,j))
             else 
                flux(1,jc) = Huge(flux)
             end if
          end do
          flux(1,:) = fac * flux(1,:)

       end if

!   Lo j face
    else if ( dim == 2 ) then
       if (face == -1) then

          j = 1
          flux(:,1) = ZERO
          do i = 1,nx
             ic = (i-1)/ratio+1
             if (bc_dirichlet(mm(i,j),2,-1)) then
                flux(ic,1) = flux(ic,1)  &
                     + ss(3,i,j)*(uu(i,j+1)-uu(i,j)) &
                     + ss(4,i,j)*(uu(i,j-1)-uu(i,j)) - ss(4,i,j+1)*(uu(i,j+1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,+1)) &
                     flux(ic,1) =  flux(ic,1) + ss(YBC,i,j)*(uu(i,j+2)-uu(i,j))
             else 
                flux(ic,1) = Huge(flux)
             end if
          end do
          flux(:,1) = fac * flux(:,1)


!      Hi j face
       else if (face == 1) then

          j = ny
          flux(:,1) = ZERO
          do i = 1,nx
             ic = (i-1)/ratio+1
             if (bc_dirichlet(mm(i,j),2,+1)) then
                flux(ic,1) = flux(ic,1)  &
                     + ss(3,i,j)*(uu(i,j+1)-uu(i,j)) &
                     + ss(4,i,j)*(uu(i,j-1)-uu(i,j)) - ss(3,i,j-1)*(uu(i,j-1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,-1)) &
                     flux(ic,1) = flux(ic,1) + ss(YBC,i,j)*(uu(i,j-2)-uu(i,j))
             else
                flux(ic,1) = Huge(flux)
             end if
          end do
          flux(:,1) = fac * flux(:,1)

       end if
    end if

  end subroutine stencil_flux_2d


  subroutine stencil_apply_3d(ss, dd, ng_d, uu, ng_u, mm, skwd, bottom_solver)

    integer           , intent(in ) :: ng_d,ng_u
    real (kind = dp_t), intent(in ) :: ss(0:,:,:,:)
    real (kind = dp_t), intent(out) :: dd(1-ng_d:,1-ng_d:,1-ng_d:)
    real (kind = dp_t), intent(in ) :: uu(1-ng_u:,1-ng_u:,1-ng_u:)
    integer           , intent(in ) :: mm(:,:,:)
    logical           , intent(in ), optional :: skwd, bottom_solver

    integer            :: nx,ny,nz,i,j,k
    integer, parameter :: XBC = 7, YBC = 8, ZBC = 9
    logical            :: lskwd, lbottom_solver

    lskwd          = .true.;   if ( present(skwd)          ) lskwd          = skwd
    lbottom_solver = .false. ; if ( present(bottom_solver) ) lbottom_solver = bottom_solver

    nx = size(ss,dim=2)
    ny = size(ss,dim=3)
    nz = size(ss,dim=4)

    !$OMP PARALLEL PRIVATE(i,j,k) IF(.not.lbottom_solver)

    if ( size(ss,dim=1) .eq. 13 ) then
       !
       ! This is the Minion 4th order cross stencil.
       ! 
       !$OMP DO
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                dd(i,j,k) = ss(0,i,j,k) * uu(i,j,k) &
                     + ss( 1,i,j,k) * uu(i-2,j,k) + ss( 2,i,j,k) * uu(i-1,j,k) &
                     + ss( 3,i,j,k) * uu(i+1,j,k) + ss( 4,i,j,k) * uu(i+2,j,k) &
                     + ss( 5,i,j,k) * uu(i,j-2,k) + ss( 6,i,j,k) * uu(i,j-1,k) &
                     + ss( 7,i,j,k) * uu(i,j+1,k) + ss( 8,i,j,k) * uu(i,j+2,k) &
                     + ss( 9,i,j,k) * uu(i,j,k-2) + ss(10,i,j,k) * uu(i,j,k-1) &
                     + ss(11,i,j,k) * uu(i,j,k+1) + ss(12,i,j,k) * uu(i,j,k+2)
             end do
          end do
       end do
       !$OMP END DO

    else if ( size(ss,dim=1) .eq. 61 ) then
       !
       ! This is the 4th order cross stencil for variable coefficients.
       !
       !$OMP DO
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                dd(i,j,k) = &
                       ss( 0,i,j,k) * uu(i  ,j  ,k  ) &
                       ! Contributions from k-2
                     + ss( 1,i,j,k) * uu(i  ,j-2,k-2) + ss( 2,i,j,k) * uu(i  ,j-1,k-2) &
                     + ss( 3,i,j,k) * uu(i-2,j  ,k-2) + ss( 4,i,j,k) * uu(i-1,j  ,k-2) &
                     + ss( 5,i,j,k) * uu(i  ,j  ,k-2) + ss( 6,i,j,k) * uu(i+1,j  ,k-2) &
                     + ss( 7,i,j,k) * uu(i+2,j  ,k-2) + ss( 8,i,j,k) * uu(i  ,j+1,k-2) &
                     + ss( 9,i,j,k) * uu(i  ,j+2,k-2)                                  &
                       ! Contributions from k-1
                     + ss(10,i,j,k) * uu(i  ,j-2,k-1) + ss(11,i,j,k) * uu(i  ,j-1,k-1) &
                     + ss(12,i,j,k) * uu(i-2,j  ,k-1) + ss(13,i,j,k) * uu(i-1,j  ,k-1) &
                     + ss(14,i,j,k) * uu(i  ,j  ,k-1) + ss(15,i,j,k) * uu(i+1,j  ,k-1) &
                     + ss(16,i,j,k) * uu(i+2,j  ,k-1) + ss(17,i,j,k) * uu(i  ,j+1,k-1) &
                     + ss(18,i,j,k) * uu(i  ,j+2,k-1)                                  &
                       ! Contributions from j-2,k
                     + ss(19,i,j,k) * uu(i-2,j-2,k  ) + ss(20,i,j,k) * uu(i-1,j-2,k  ) &
                     + ss(21,i,j,k) * uu(i  ,j-2,k  ) + ss(22,i,j,k) * uu(i+1,j-2,k  ) &
                     + ss(23,i,j,k) * uu(i+2,j-2,k  )                                  &
                       ! Contributions from j-1,k
                     + ss(24,i,j,k) * uu(i-2,j-1,k  ) + ss(25,i,j,k) * uu(i-1,j-1,k  ) &
                     + ss(26,i,j,k) * uu(i  ,j-1,k  ) + ss(27,i,j,k) * uu(i+1,j-1,k  ) &
                     + ss(28,i,j,k) * uu(i+2,j-1,k  )                                  &
                       ! Contributions from j  ,k
                     + ss(29,i,j,k) * uu(i-2,j  ,k  ) + ss(30,i,j,k) * uu(i-1,j  ,k  ) &
                                                      + ss(31,i,j,k) * uu(i+1,j  ,k  ) &
                     + ss(32,i,j,k) * uu(i+2,j  ,k  )                                  &
                       ! Contributions from j+1,k
                     + ss(33,i,j,k) * uu(i-2,j+1,k  ) + ss(34,i,j,k) * uu(i-1,j+1,k  ) &
                     + ss(35,i,j,k) * uu(i  ,j+1,k  ) + ss(36,i,j,k) * uu(i+1,j+1,k  ) &
                     + ss(37,i,j,k) * uu(i+2,j+1,k  )                                  &
                       ! Contributions from j+2,k
                     + ss(38,i,j,k) * uu(i-2,j+2,k  ) + ss(39,i,j,k) * uu(i-1,j+2,k  ) &
                     + ss(40,i,j,k) * uu(i  ,j+2,k  ) + ss(41,i,j,k) * uu(i+1,j+2,k  ) &
                     + ss(42,i,j,k) * uu(i+2,j+2,k  )                                  &
                       ! Contributions from k+1
                     + ss(43,i,j,k) * uu(i  ,j-2,k+1) + ss(44,i,j,k) * uu(i  ,j-1,k+1) &
                     + ss(45,i,j,k) * uu(i-2,j  ,k+1) + ss(46,i,j,k) * uu(i-1,j  ,k+1) &
                     + ss(47,i,j,k) * uu(i  ,j  ,k+1) + ss(48,i,j,k) * uu(i+1,j  ,k+1) &
                     + ss(49,i,j,k) * uu(i+2,j  ,k+1) + ss(50,i,j,k) * uu(i  ,j+1,k+1) &
                     + ss(51,i,j,k) * uu(i  ,j+2,k+1)                                  &
                       ! Contributions from k+2
                     + ss(52,i,j,k) * uu(i  ,j-2,k+2) + ss(53,i,j,k) * uu(i  ,j-1,k+2) &
                     + ss(54,i,j,k) * uu(i-2,j  ,k+2) + ss(55,i,j,k) * uu(i-1,j  ,k+2) &
                     + ss(56,i,j,k) * uu(i  ,j  ,k+2) + ss(57,i,j,k) * uu(i+1,j  ,k+2) &
                     + ss(58,i,j,k) * uu(i+2,j  ,k+2) + ss(59,i,j,k) * uu(i  ,j+1,k+2) &
                     + ss(60,i,j,k) * uu(i  ,j+2,k+2)
             end do
          end do
       end do
       !$OMP END DO

    else
       !
       ! This is the 2nd order cross stencil.
       !
       !$OMP DO
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                dd(i,j,k) = &
                     ss(0,i,j,k)*uu(i,j,k)       + &
                     ss(1,i,j,k)*uu(i+1,j  ,k  ) + &
                     ss(2,i,j,k)*uu(i-1,j  ,k  ) + &
                     ss(3,i,j,k)*uu(i  ,j+1,k  ) + &
                     ss(4,i,j,k)*uu(i  ,j-1,k  ) + &
                     ss(5,i,j,k)*uu(i  ,j  ,k+1) + &
                     ss(6,i,j,k)*uu(i  ,j  ,k-1)
             end do
          end do
       end do
       !$OMP END DO

    end if

    if ( lskwd ) then
       !
       ! Corrections for skewed stencils
       !
       if (nx > 1) then
          !$OMP DO
          do k = 1, nz
             do j = 1, ny
                i = 1
                if (bc_skewed(mm(i,j,k),1,+1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(XBC,i,j,k)*uu(i+2,j,k)
                end if

                i = nx
                if (bc_skewed(mm(i,j,k),1,-1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(XBC,i,j,k)*uu(i-2,j,k)
                end if
             end do
          end do
          !$OMP END DO
       end if

       if (ny > 1) then
          !$OMP DO
          do k = 1,nz
             do i = 1,nx
                j = 1
                if (bc_skewed(mm(i,j,k),2,+1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(YBC,i,j,k)*uu(i,j+2,k)
                end if

                j = ny
                if (bc_skewed(mm(i,j,k),2,-1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(YBC,i,j,k)*uu(i,j-2,k)
                end if
             end do
          end do
          !$OMP END DO
       end if

       if (nz > 1) then
          !$OMP DO
          do j = 1,ny
             do i = 1,nx
                k = 1
                if (bc_skewed(mm(i,j,k),3,+1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(ZBC,i,j,k)*uu(i,j,k+2)
                end if

                k = nz
                if (bc_skewed(mm(i,j,k),3,-1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(ZBC,i,j,k)*uu(i,j,k-2)
                end if
             end do
          end do
          !$OMP END DO
       end if
    end if

    !$OMP END PARALLEL

  end subroutine stencil_apply_3d

  subroutine stencil_apply_ibc_3d(ss, dd, ng_d, uu, ng_u, lo, hi)

    integer           , intent(in ) :: ng_d,ng_u, lo(:), hi(:)
    real (kind = dp_t), intent(in   ) :: ss(0:)
    real (kind = dp_t), intent(  out) :: dd(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:)
    real (kind = dp_t), intent(in   ) :: uu(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)

    integer i,j,k

    ! This is our standard 7-point Laplacian without correction at boundaries
    !$omp parallel do private(i,j,k) collapse(2)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             dd(i,j,k) = ss(0)*uu(i,j,k) &
                  +      ss(1)*(uu(i-1,j,k) + uu(i+1,j,k)) &
                  +      ss(2)*(uu(i,j-1,k) + uu(i,j+1,k)) &
                  +      ss(3)*(uu(i,j,k-1) + uu(i,j,k+1))
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine stencil_apply_ibc_3d

  subroutine stencil_flux_3d(ss, flux, uu, mm, ng, ratio, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), intent(out) :: flux(:,:,:)
    real (kind = dp_t), intent(in ) :: ss(0:,:,:,:)
    integer           , intent(in)  :: mm(:,:,:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: ratio, face, dim
    integer nx, ny, nz
    integer i,j,k,ic,jc,kc
    integer, parameter :: XBC = 7, YBC = 8, ZBC = 9
    real (kind = dp_t) :: fac
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=2)
    ny = size(ss,dim=3)
    nz = size(ss,dim=4)

    ! Note that two factors of ratio is from the tangential averaging, while the
    ! other is the normal factor
    fac = ONE/real(ratio*ratio*ratio, kind=dp_t)

    ! Note: Do not try to add OMP calls to this subroutine.  For example,
    !       in the first k loop below, kc may end up having the same value 
    !       on multiple threads, and then you try to update the same flux(1,jc,kc)
    !       memory simultaneously on different threads.

    !   Lo i face
    if ( dim ==  1 ) then
       if (face == -1) then

          !
          !   Lo i face
          !
          i = 1
          flux(1,:,:) = ZERO
          do k = 1,nz
             do j = 1,ny
                jc = (j-1)/ratio + 1
                kc = (k-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),1,-1)) then
                   flux(1,jc,kc) =  flux(1,jc,kc) &
                        + ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(2,i+1,j,k)*(uu(i+1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,+1)) &
                        flux(1,jc,kc) =  flux(1,jc,kc) + ss(XBC,i,j,k)*(uu(i+2,j,k)-uu(i,j,k))
                else 
                   flux(1,jc,kc) = Huge(flux)
                end if
             end do
          end do
          flux(1,:,:) = flux(1,:,:) * fac

          !   Hi i face
       else if (face ==  1) then

          !
          !   Hi i face
          !
          i = nx
          flux(1,:,:) = ZERO
          do k = 1,nz
             do j = 1,ny
                jc = (j-1)/ratio + 1
                kc = (k-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),1,+1)) then
                   flux(1,jc,kc) =  flux(1,jc,kc) &
                        + ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(1,i-1,j,k)*(uu(i-1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,-1)) &
                        flux(1,jc,kc) =  flux(1,jc,kc) + ss(XBC,i,j,k)*(uu(i-2,j,k)-uu(i,j,k))
                else 
                   flux(1,jc,kc) = Huge(flux)
                end if
             end do
          end do
          flux(1,:,:) = flux(1,:,:) * fac

       end if
       !   Lo j face
    else if ( dim == 2 ) then
       if (face == -1) then
          !
          !   Lo j face
          !
          j = 1
          flux(:,1,:) = ZERO
          do k = 1,nz
             do i = 1,nx
                ic = (i-1)/ratio + 1
                kc = (k-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),2,-1)) then
                   flux(ic,1,kc) =  flux(ic,1,kc) &
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(4,i,j+1,k)*(uu(i,j+1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,+1)) &
                        flux(ic,1,kc) =  flux(ic,1,kc) + ss(YBC,i,j,k)*(uu(i,j+2,k)-uu(i,j,k))
                else 
                   flux(ic,1,kc) = Huge(flux)
                end if
             end do
          end do
          flux(:,1,:) = flux(:,1,:) * fac

       else if (face ==  1) then
          !
          !   Hi j face
          !
          j = ny
          flux(:,1,:) = ZERO
          do k = 1,nz
             do i = 1,nx
                ic = (i-1)/ratio + 1
                kc = (k-1)/ratio + 1

                if (bc_dirichlet(mm(i,j,k),2,+1)) then
                   flux(ic,1,kc) =  flux(ic,1,kc) &
                        + ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(3,i,j-1,k)*(uu(i,j-1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,-1)) &
                        flux(ic,1,kc) =  flux(ic,1,kc) + ss(YBC,i,j,k)*(uu(i,j-2,k)-uu(i,j,k))
                else
                   flux(ic,1,kc) = Huge(flux)
                end if
             end do
          end do
          flux(:,1,:) = flux(:,1,:) * fac

       end if
    else if ( dim == 3 ) then
       if (face == -1) then

          !
          !   Lo k face
          !
          k = 1
          flux(:,:,1) = ZERO
          do j = 1,ny
             do i = 1,nx
                ic = (i-1)/ratio + 1
                jc = (j-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),3,-1)) then
                   flux(ic,jc,1) =  flux(ic,jc,1) &
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(6,i,j,k+1)*(uu(i,j,k+1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,+1)) &
                        flux(ic,jc,1) =  flux(ic,jc,1) + ss(ZBC,i,j,k)*(uu(i,j,k+2)-uu(i,j,k)) 
                else 
                   flux(ic,jc,1) = Huge(flux)
                end if
             end do
          end do
          flux(:,:,1) = flux(:,:,1) * fac

       else if (face ==  1) then

          !
          !   Hi k face
          !
          k = nz
          flux(:,:,1) = ZERO
          do j = 1,ny
             do i = 1,nx
                ic = (i-1)/ratio + 1
                jc = (j-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),3,+1)) then
                   flux(ic,jc,1) =  flux(ic,jc,1) &
                        + ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(5,i,j,k-1)*(uu(i,j,k-1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,-1)) &
                        flux(ic,jc,1) =  flux(ic,jc,1) + ss(ZBC,i,j,k)*(uu(i,j,k-2)-uu(i,j,k))
                else
                   flux(ic,jc,1) = Huge(flux)
                end if
             end do
          end do
          flux(:,:,1) = flux(:,:,1) * fac

       end if
    end if

  end subroutine stencil_flux_3d

  ! subroutine stencil_dense_apply_1d(ss, dd, ng_d, uu, ng_u)
  !   integer, intent(in) :: ng_d, ng_u
  !   real (kind = dp_t), intent(in   ) :: ss(0:,:)
  !   real (kind = dp_t), intent(  out) :: dd(1-ng_d:)
  !   real (kind = dp_t), intent(in   ) :: uu(1-ng_u:)
  !   integer i, nx
   
  !   nx = size(ss,dim=2)
  !   do i = 1, nx
  !     dd(i) = ss(1,i)*uu(i-1) + ss(0,i)*uu(i) + ss(2,i)*uu(i+1)
  !   end do

  ! end subroutine stencil_dense_apply_1d

  ! subroutine stencil_dense_apply_2d(ss, dd, ng_d, uu, ng_u)
  !   integer, intent(in) :: ng_d, ng_u
  !   real (kind = dp_t), intent(in   ) :: ss(0:,:,:)
  !   real (kind = dp_t), intent(  out) :: dd(1-ng_d:,1-ng_d:)
  !   real (kind = dp_t), intent(in   ) :: uu(1-ng_u:,1-ng_u:)
  !   integer i, j, nx, ny

  !   nx = size(ss,dim=2)
  !   ny = size(ss,dim=3)

  !   do j = 1, ny
  !      do i = 1, nx
  !         dd(i,j) = &
  !              + ss(1,i,j)*uu(i-1,j-1) + ss(2,i,j)*uu(i  ,j-1) + ss(3,i,j)*uu(i+1,j-1) &
  !              + ss(4,i,j)*uu(i-1,j  ) + ss(0,i,j)*uu(i  ,j  ) + ss(5,i,j)*uu(i+1,j  ) &
  !              + ss(6,i,j)*uu(i-1,j+1) + ss(7,i,j)*uu(i  ,j+1) + ss(8,i,j)*uu(i+1,j+1)
  !      end do
  !   end do

  ! end subroutine stencil_dense_apply_2d

  ! subroutine stencil_dense_apply_3d(ss, dd, ng_d, uu, ng_u)
  !   integer, intent(in) :: ng_d, ng_u
  !   real (kind = dp_t), intent(in   ) :: ss(0:,:,:,:)
  !   real (kind = dp_t), intent(in   ) :: uu(1-ng_u:,1-ng_u:,1-ng_u:)
  !   real (kind = dp_t), intent(  out) :: dd(1-ng_d:,1-ng_d:,1-ng_d:)
  !   integer i, j, k, nx, ny, nz

  !   nx = size(ss,dim=2)
  !   ny = size(ss,dim=3)
  !   nz = size(ss,dim=4)

  !   do k = 1, nz
  !      do j = 1, ny
  !         do i = 1, nx
  !            dd(i,j,k) = &
  !                 + ss( 1,i,j,k)*uu(i-1,j-1,k-1) &
  !                 + ss( 2,i,j,k)*uu(i  ,j-1,k-1) &
  !                 + ss( 3,i,j,k)*uu(i+1,j-1,k-1) &
  !                 + ss( 4,i,j,k)*uu(i-1,j  ,k-1) &
  !                 + ss( 5,i,j,k)*uu(i  ,j  ,k-1) &
  !                 + ss( 6,i,j,k)*uu(i+1,j  ,k-1) &
  !                 + ss( 7,i,j,k)*uu(i-1,j+1,k-1) &
  !                 + ss( 8,i,j,k)*uu(i  ,j+1,k-1) &
  !                 + ss( 9,i,j,k)*uu(i+1,j+1,k-1) &

  !                 + ss(10,i,j,k)*uu(i-1,j-1,k  ) &
  !                 + ss(11,i,j,k)*uu(i  ,j-1,k  ) &
  !                 + ss(12,i,j,k)*uu(i+1,j-1,k  ) &
  !                 + ss(13,i,j,k)*uu(i-1,j  ,k  ) &
  !                 + ss( 0,i,j,k)*uu(i  ,j  ,k  ) &
  !                 + ss(14,i,j,k)*uu(i+1,j  ,k  ) &
  !                 + ss(15,i,j,k)*uu(i-1,j+1,k  ) &
  !                 + ss(16,i,j,k)*uu(i  ,j+1,k  ) &
  !                 + ss(17,i,j,k)*uu(i+1,j+1,k  ) &

  !                 + ss(18,i,j,k)*uu(i-1,j-1,k+1) &
  !                 + ss(19,i,j,k)*uu(i  ,j-1,k+1) &
  !                 + ss(20,i,j,k)*uu(i+1,j-1,k+1) &
  !                 + ss(21,i,j,k)*uu(i-1,j  ,k+1) &
  !                 + ss(22,i,j,k)*uu(i  ,j  ,k+1) &
  !                 + ss(23,i,j,k)*uu(i+1,j  ,k+1) &
  !                 + ss(24,i,j,k)*uu(i-1,j+1,k+1) &
  !                 + ss(25,i,j,k)*uu(i  ,j+1,k+1) &
  !                 + ss(26,i,j,k)*uu(i+1,j+1,k+1)
  !         end do
  !      end do
  !   end do

  ! end subroutine stencil_dense_apply_3d

  subroutine stencil_fine_flux_1d(ss, flux, uu, mm, ng, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in)  :: ss(0:,:)
    real (kind = dp_t), intent(out) :: flux(:)
    real (kind = dp_t), intent(in)  :: uu(1-ng:)
    integer           , intent(in)  :: mm(:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: face, dim
    integer nx
    integer i
    integer, parameter :: XBC = 3
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=2)

    if ( dim == 1 ) then
       if ( face == -1 ) then
!         Lo i face
          i = 1
          if (bc_dirichlet(mm(1),1,-1)) then
             flux(1) = ss(1,i)*(uu(i+1)-uu(i)) + ss(2,i)*(uu(i-1)-uu(i)) &
                  - ss(2,i+1)*(uu(i+1)-uu(i))
             if (bc_skewed(mm(i),1,+1)) then
                flux(1) =  flux(1) + ss(XBC,i)*uu(i+2)
             end if
          else 
             flux(1) = ss(2,i)*(uu(i-1)-uu(i))
          end if
       else if ( face == 1 ) then

!         Hi i face
          i = nx
          if (bc_dirichlet(mm(i),1,+1)) then
             flux(1) = ss(1,i)*(uu(i+1)-uu(i)) + ss(2,i)*(uu(i-1)-uu(i)) &
                  - ss(1,i-1)*(uu(i-1)-uu(i))
             if (bc_skewed(mm(i),1,-1)) then
                flux(1) =  flux(1) + ss(XBC,i)*uu(i-2)
             end if
          else 
             flux(1) = ss(1,i)*(uu(i+1)-uu(i))
          end if
       end if
    end if

  end subroutine stencil_fine_flux_1d

  subroutine stencil_all_flux_1d(ss, flux, uu, mm, ngu, ngf, skwd)
    integer, intent(in) :: ngu, ngf
    real (kind = dp_t), intent(in ) ::   uu(-ngu:)
    real (kind = dp_t), intent(out) :: flux(-ngf:)
    real (kind = dp_t), intent(in ) :: ss(0:,0:)
    integer           , intent(in)  :: mm(0:)
    logical, intent(in), optional :: skwd
    integer nx
    integer i
    integer, parameter :: XBC = 3, YBC = 4
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=2)

    do i = 1,nx-2
      flux(i) = ss(2,i) * (uu(i)-uu(i-1)) 
    end do

    ! Must make sure we use stencil from interior fine cell, not fine cell
    ! next to c/f boundary.  Because we use ss(2,i,j) which looks "down", we
    ! only modify at the high side
    flux(nx-1) = ss(1,nx-2) * (uu(nx-1)-uu(nx-2)) 

    ! Lo i face
     i = 0
     if (bc_dirichlet(mm(i),1,-1)) then
        flux(0) = &
               ss(1,i)*(uu(i+1)-uu(i)) + ss(2,i  )*(uu(i-1)-uu(i)) &
                                       - ss(2,i+1)*(uu(i+1)-uu(i))
        if (bc_skewed(mm(i),1,+1)) &
             flux(0) = flux(0) + ss(XBC,i)*(uu(i+2)-uu(i)) 
        flux(0) = -flux(0)
     else if (bc_neumann(mm(i),1,-1)) then
        flux(0) = -ss(2,i)*uu(i-1)
        else   
        flux(0) = ss(2,i)*(uu(i)-uu(i-1))
     end if

    ! Hi i face
     i = nx-1
     if (bc_dirichlet(mm(i),1,+1)) then
        flux(nx) = &
               ss(1,i  )*(uu(i+1)-uu(i)) + ss(2,i)*(uu(i-1)-uu(i)) &
             - ss(1,i-1)*(uu(i-1)-uu(i))
        if (bc_skewed(mm(i),1,-1)) &
             flux(nx) = flux(nx) + ss(XBC,i)*(uu(i-2)-uu(i))
     else if (bc_neumann(mm(i),1,+1)) then
        flux(nx) = ss(1,i)*uu(i+1)
     else 
        flux(nx) = ss(1,i)*(uu(i+1)-uu(i))
     end if

  end subroutine stencil_all_flux_1d

  subroutine stencil_fine_flux_2d(ss, flux, uu, mm, ng, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:)
    real (kind = dp_t), intent(out) :: flux(:,:)
    real (kind = dp_t), intent(in ) :: ss(0:,:,:)
    integer           , intent(in)  :: mm(:,:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: face, dim
    integer nx,ny
    integer i,j
    integer, parameter :: XBC = 5, YBC = 6
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=2)
    ny = size(ss,dim=3)

!   Lo i face
    if ( dim == 1 ) then
       if (face == -1) then

          i = 1
          flux(1,:) = ZERO
          do j = 1,ny
             if (bc_dirichlet(mm(i,j),1,-1)) then
                flux(1,j) = &
                       ss(1,i,j)*(uu(i+1,j)-uu(i,j)) &
                     + ss(2,i,j)*(uu(i-1,j)-uu(i,j)) - ss(2,i+1,j)*(uu(i+1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,+1)) &
                     flux(1,j) = flux(1,j) + ss(XBC,i,j)*(uu(i+2,j)-uu(i,j)) 
             else if (bc_neumann(mm(i,j),1,-1)) then
                flux(1,j) = ss(2,i,j)*uu(i-1,j)
             else   
                flux(1,j) = ss(2,i,j)*(uu(i-1,j)-uu(i,j))
             end if
          end do

!      Hi i face
       else if (face == 1) then

          i = nx
          flux(1,:) = ZERO
          do j = 1,ny
             if (bc_dirichlet(mm(i,j),1,+1)) then
                flux(1,j) = &
                       ss(1,i,j)*(uu(i+1,j)-uu(i,j)) &
                     + ss(2,i,j)*(uu(i-1,j)-uu(i,j)) - ss(1,i-1,j)*(uu(i-1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,-1)) &
                     flux(1,j) = flux(1,j) + ss(XBC,i,j)*(uu(i-2,j)-uu(i,j))
             else if (bc_neumann(mm(i,j),1,+1)) then
                flux(1,j) = ss(1,i,j)*uu(i+1,j)
             else 
                flux(1,j) = ss(1,i,j)*(uu(i+1,j)-uu(i,j))
             end if
          end do

       end if

!   Lo j face
    else if ( dim == 2 ) then
       if (face == -1) then

          j = 1
          flux(:,1) = ZERO
          do i = 1,nx
             if (bc_dirichlet(mm(i,j),2,-1)) then
                flux(i,1) = &
                       ss(3,i,j)*(uu(i,j+1)-uu(i,j)) &
                     + ss(4,i,j)*(uu(i,j-1)-uu(i,j)) - ss(4,i,j+1)*(uu(i,j+1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,+1)) &
                     flux(i,1) =  flux(i,1) + ss(YBC,i,j)*(uu(i,j+2)-uu(i,j))
             else if (bc_neumann(mm(i,j),2,-1)) then
                flux(i,1) = ss(4,i,j)*uu(i,j-1)
             else 
                flux(i,1) = ss(4,i,j)*(uu(i,j-1)-uu(i,j))
             end if
          end do


!      Hi j face
       else if (face == 1) then

          j = ny
          flux(:,1) = ZERO
          do i = 1,nx
             if (bc_dirichlet(mm(i,j),2,+1)) then
                flux(i,1) = &
                       ss(3,i,j)*(uu(i,j+1)-uu(i,j)) &
                     + ss(4,i,j)*(uu(i,j-1)-uu(i,j)) - ss(3,i,j-1)*(uu(i,j-1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,-1)) &
                     flux(i,1) = flux(i,1) + ss(YBC,i,j)*(uu(i,j-2)-uu(i,j))
             else if (bc_neumann(mm(i,j),2,+1)) then
                flux(i,1) = ss(3,i,j)*uu(i,j+1)
             else
                flux(i,1) = ss(3,i,j)*(uu(i,j+1)-uu(i,j))
             end if
          end do

       end if
    end if

  end subroutine stencil_fine_flux_2d

  subroutine stencil_all_flux_2d(ss, flux, uu, mm, ngu, ngf, dim, skwd)
    integer, intent(in) :: ngu, ngf
    real (kind = dp_t), intent(in ) ::   uu(-ngu:,-ngu:)
    real (kind = dp_t), intent(out) :: flux(-ngf:,-ngf:)
    real (kind = dp_t), intent(in ) :: ss(0:,0:,0:)
    integer           , intent(in)  :: mm(0:,0:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: dim
    integer nx,ny
    integer i,j
    integer, parameter :: XBC = 5, YBC = 6
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=2)
    ny = size(ss,dim=3)

    if ( dim == 1 ) then
       do j = 0,ny-1
       do i = 1,nx-2
         flux(i,j) = ss(2,i,j) * (uu(i,j)-uu(i-1,j)) 
       end do
       end do

       ! Must make sure we use stencil from interior fine cell, not fine cell next to c/f boundary
       ! Because we use ss(2,i,j) which looks "down", we only modify at the high side
       do j = 0, ny-1
         flux(nx-1,j) = ss(1,nx-2,j) * (uu(nx-1,j)-uu(nx-2,j)) 
       end do

       ! Lo i face
        i = 0
        do j = 0,ny-1
             if (bc_dirichlet(mm(i,j),1,-1)) then
                flux(0,j) = ss(2,i,j)             *(uu(i-1,j)-uu(i,j)) &
                         + (ss(1,i,j)-ss(2,i+1,j))*(uu(i+1,j)-uu(i,j)) 
                if (bc_skewed(mm(i,j),1,+1)) &
                     flux(0,j) = flux(0,j) + ss(XBC,i,j)*(uu(i+2,j)-uu(i,j)) 
                flux(0,j) = -flux(0,j)
             else if (bc_neumann(mm(i,j),1,-1)) then
                flux(0,j) = -ss(2,i,j)*uu(i-1,j)
             else   
                flux(0,j) = ss(2,i,j)*(uu(i,j)-uu(i-1,j))
             end if
        end do

       ! Hi i face
        i = nx-1
        do j = 0,ny-1
             if (bc_dirichlet(mm(i,j),1,+1)) then
                flux(nx,j) = ss(1,i,j)             *(uu(i+1,j)-uu(i,j)) &
                          + (ss(2,i,j)-ss(1,i-1,j))*(uu(i-1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,-1)) &
                     flux(nx,j) = flux(nx,j) + ss(XBC,i,j)*(uu(i-2,j)-uu(i,j))
             else if (bc_neumann(mm(i,j),1,+1)) then
                flux(nx,j) = ss(1,i,j)*uu(i+1,j)
             else 
                flux(nx,j) = ss(1,i,j)*(uu(i+1,j)-uu(i,j))
             end if
        end do

    else if ( dim == 2 ) then

       do j = 1,ny-2
       do i = 0,nx-1
         flux(i,j) = ss(4,i,j) * (uu(i,j)-uu(i,j-1)) 
       end do
       end do

       ! Must make sure we use stencil from interior cell, not cell next to c/f boundary
       ! Because we use ss(4,i,j) which looks "down", we only modify at the high side
       do i = 0,nx-1
         flux(i,ny-1) = ss(3,i,ny-2) * (uu(i,ny-1)-uu(i,ny-2)) 
       end do

       ! Lo j face
       j = 0
       do i = 0,nx-1
             if (bc_dirichlet(mm(i,j),2,-1)) then
                flux(i,0) = &
                       ss(3,i,j)*(uu(i,j+1)-uu(i,j)) + ss(4,i,j  )*(uu(i,j-1)-uu(i,j)) & 
                                                     - ss(4,i,j+1)*(uu(i,j+1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,+1)) &
                     flux(i,0) =  flux(i,0) + ss(YBC,i,j)*(uu(i,j+2)-uu(i,j))
                flux(i,0) = -flux(i,0)
             else if (bc_neumann(mm(i,j),2,-1)) then
                flux(i,0) = -ss(4,i,j)*uu(i,j-1)
             else 
                flux(i,0) = ss(4,i,j)*(uu(i,j)-uu(i,j-1))
             end if
       end do

       ! Hi j face
       j = ny-1
       do i = 0,nx-1
             if (bc_dirichlet(mm(i,j),2,+1)) then
                flux(i,ny) = &
                       ss(3,i,j  )*(uu(i,j+1)-uu(i,j)) + ss(4,i,j)*(uu(i,j-1)-uu(i,j)) & 
                     - ss(3,i,j-1)*(uu(i,j-1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,-1)) &
                     flux(i,ny) = flux(i,ny) + ss(YBC,i,j)*(uu(i,j-2)-uu(i,j))
             else if (bc_neumann(mm(i,j),2,+1)) then
                flux(i,ny) = ss(3,i,j)*uu(i,j+1)
             else
                flux(i,ny) = ss(3,i,j)*(uu(i,j+1)-uu(i,j))
             end if
             flux(i,ny-1) = ss(4,i,ny-2) * (uu(i,j)-uu(i,j-1)) 
       end do

    end if

  end subroutine stencil_all_flux_2d

  subroutine stencil_fine_flux_3d(ss, flux, uu, mm, ng, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: ss(0:,:,:,:)
    real (kind = dp_t), intent(out) :: flux(:,:,:)
    real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:,1-ng:)
    integer           , intent(in)  :: mm(:,:,:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: face, dim
    integer nx, ny, nz
    integer i,j,k
    integer, parameter :: XBC = 7, YBC = 8, ZBC = 9
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=2)
    ny = size(ss,dim=3)
    nz = size(ss,dim=4)

    if ( dim ==  1 ) then
       !   Lo i face
       if (face == -1) then

          i = 1
          flux(1,:,:) = ZERO

          do k = 1,nz
             do j = 1,ny
                if (bc_dirichlet(mm(i,j,k),1,-1)) then
                   flux(1,j,k) =  &
                          ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(2,i+1,j,k)*(uu(i+1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,+1)) &
                        flux(1,j,k) =  flux(1,j,k) + ss(XBC,i,j,k)*(uu(i+2,j,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),1,-1)) then
                   flux(1,j,k) = ss(2,i,j,k)*uu(i-1,j,k)
                else 
                   flux(1,j,k) = ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k))
                end if
             end do
          end do

       !   Hi i face
       else if (face ==  1) then

          i = nx
          flux(1,:,:) = ZERO
          do k = 1,nz
             do j = 1,ny
                if (bc_dirichlet(mm(i,j,k),1,+1)) then
                   flux(1,j,k) = &
                          ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(1,i-1,j,k)*(uu(i-1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,-1)) &
                        flux(1,j,k) =  flux(1,j,k) + ss(XBC,i,j,k)*(uu(i-2,j,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),1,+1)) then
                   flux(1,j,k) = ss(1,i,j,k)*uu(i+1,j,k)
                else 
                   flux(1,j,k) = ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k))
                end if
             end do
          end do
       end if

    else if ( dim == 2 ) then

       !   Lo j face
       if (face == -1) then
          j = 1
          flux(:,1,:) = ZERO
          do k = 1,nz
             do i = 1,nx
                if (bc_dirichlet(mm(i,j,k),2,-1)) then
                   flux(i,1,k) = &
                          ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(4,i,j+1,k)*(uu(i,j+1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,+1)) &
                        flux(i,1,k) =  flux(i,1,k) + ss(YBC,i,j,k)*(uu(i,j+2,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),2,-1)) then
                   flux(i,1,k) = ss(4,i,j,k)*uu(i,j-1,k)
                else 
                   flux(i,1,k) = ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k))
                end if
             end do
          end do

       !   Hi j face
       else if (face ==  1) then

          j = ny
          flux(:,1,:) = ZERO
          do k = 1,nz
             do i = 1,nx
                if (bc_dirichlet(mm(i,j,k),2,+1)) then
                   flux(i,1,k) =  &
                          ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(3,i,j-1,k)*(uu(i,j-1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,-1)) &
                        flux(i,1,k) =  flux(i,1,k) + ss(YBC,i,j,k)*(uu(i,j-2,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),2,+1)) then
                   flux(i,1,k) = ss(3,i,j,k)*uu(i,j+1,k)
                else
                   flux(i,1,k) = ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k))
                end if
             end do
          end do
       end if

    else if ( dim == 3 ) then

       !   Lo k face
       if (face == -1) then

          k = 1
          flux(:,:,1) = ZERO
          do j = 1,ny
             do i = 1,nx
                if (bc_dirichlet(mm(i,j,k),3,-1)) then
                   flux(i,j,1) =  &
                          ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(6,i,j,k+1)*(uu(i,j,k+1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,+1)) &
                        flux(i,j,1) =  flux(i,j,1) + ss(ZBC,i,j,k)*(uu(i,j,k+2)-uu(i,j,k)) 
                else if (bc_neumann(mm(i,j,k),3,-1)) then
                   flux(i,j,1) = ss(6,i,j,k)*uu(i,j,k-1)
                else 
                   flux(i,j,1) = ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k))
                end if
             end do
          end do

       !   Hi k face
       else if (face ==  1) then

          k = nz
          flux(:,:,1) = ZERO
          do j = 1,ny
             do i = 1,nx
                if (bc_dirichlet(mm(i,j,k),3,+1)) then
                   flux(i,j,1) =  &
                          ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(5,i,j,k-1)*(uu(i,j,k-1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,-1)) &
                        flux(i,j,1) =  flux(i,j,1) + ss(ZBC,i,j,k)*(uu(i,j,k-2)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),3,+1)) then
                   flux(i,j,1) = ss(5,i,j,k)*uu(i,j,k+1)
                else
                   flux(i,j,1) = ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k))
                end if
             end do
          end do

       end if
    end if

  end subroutine stencil_fine_flux_3d

  subroutine stencil_all_flux_3d(ss, flux, uu, mm, ngu, ngf, dim, skwd)
    integer, intent(in) :: ngu,ngf
    real (kind = dp_t), intent(in ) ::   uu(-ngu:,-ngu:,-ngu:)
    real (kind = dp_t), intent(out) :: flux(-ngf:,-ngf:,-ngf:)
    real (kind = dp_t), intent(in ) :: ss(0:,0:,0:,0:)
    integer           , intent(in)  :: mm(0:,0:,0:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: dim
    integer nx, ny, nz
    integer i,j,k
    integer, parameter :: XBC = 7, YBC = 8, ZBC = 9
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=2)
    ny = size(ss,dim=3)
    nz = size(ss,dim=4)

    if ( dim ==  1 ) then

       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 0,nz-1
          do j = 0,ny-1
             do i = 1,nx-2
                flux(i,j,k) = ss(2,i,j,k) * (uu(i,j,k)-uu(i-1,j,k))
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       !
       ! Must make sure we use stencil from interior fine cell, not fine cell
       ! next to c/f boundary. Because we use ss(2,i,j) which looks "down", we
       ! only modify at the high side
       !
       !$OMP PARALLEL DO PRIVATE(j,k)
       do k = 0,nz-1
          do j = 0,ny-1
             flux(nx-1,j,k) = ss(1,nx-2,j,k) * (uu(nx-1,j,k)-uu(nx-2,j,k)) 
          end do
       end do
       !$OMP END PARALLEL DO

       !   Lo i face
       i = 0
       !$OMP PARALLEL DO PRIVATE(j,k)
       do k = 0,nz-1
             do j = 0,ny-1
                if (bc_dirichlet(mm(i,j,k),1,-1)) then
                   flux(0,j,k) =  &
                          ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(2,i+1,j,k)*(uu(i+1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,+1)) &
                        flux(0,j,k) =  flux(0,j,k) + ss(XBC,i,j,k)*(uu(i+2,j,k)-uu(i,j,k))
                   flux(0,j,k) = -flux(0,j,k)
                else if (bc_neumann(mm(i,j,k),1,-1)) then
                   flux(0,j,k) = -ss(2,i,j,k)*uu(i-1,j,k)
                else 
                   flux(0,j,k) = ss(2,i,j,k)*(uu(i,j,k)-uu(i-1,j,k))
                end if
             end do
       end do
       !$OMP END PARALLEL DO

       !   Hi i face
       i = nx-1
       !$OMP PARALLEL DO PRIVATE(j,k)
       do k = 0,nz-1
             do j = 0,ny-1
                if (bc_dirichlet(mm(i,j,k),1,+1)) then
                   flux(nx,j,k) = &
                          ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(2,i,j,k)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(1,i-1,j,k)*(uu(i-1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,-1)) &
                        flux(nx,j,k) =  flux(nx,j,k) + ss(XBC,i,j,k)*(uu(i-2,j,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),1,+1)) then
                   flux(nx,j,k) = ss(1,i,j,k)*uu(i+1,j,k)
                else 
                   flux(nx,j,k) = ss(1,i,j,k)*(uu(i+1,j,k)-uu(i,j,k))
                end if
             end do
       end do
       !$OMP END PARALLEL DO

    else if ( dim == 2 ) then

       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 0,nz-1
          do j = 1,ny-2
             do i = 0,nx-1
                flux(i,j,k) = ss(4,i,j,k) * (uu(i,j,k)-uu(i,j-1,k))
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       !
       ! Must make sure we use stencil from interior fine cell, not fine cell
       ! next to c/f boundary. Because we use ss(2,i,j) which looks "down", we
       ! only modify at the high side
       !
       !$OMP PARALLEL DO PRIVATE(i,k)
       do k = 0,nz-1
          do i = 0,nx-1
             flux(i,ny-1,k) = ss(3,i,ny-2,k) * (uu(i,ny-1,k)-uu(i,ny-2,k))
          end do
       end do
       !$OMP END PARALLEL DO

       !   Lo j face
       j = 0
       !$OMP PARALLEL DO PRIVATE(i,k)
       do k = 0,nz-1
             do i = 0,nx-1
                if (bc_dirichlet(mm(i,j,k),2,-1)) then
                   flux(i,0,k) = &
                          ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(4,i,j+1,k)*(uu(i,j+1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,+1)) &
                        flux(i,0,k) =  flux(i,0,k) + ss(YBC,i,j,k)*(uu(i,j+2,k)-uu(i,j,k))
                   flux(i,0,k) = -flux(i,0,k)
                else if (bc_neumann(mm(i,j,k),2,-1)) then
                   flux(i,0,k) = -ss(4,i,j,k)*uu(i,j-1,k)
                else 
                   flux(i,0,k) = ss(4,i,j,k)*(uu(i,j,k)-uu(i,j-1,k))
                end if
             end do
       end do
       !$OMP END PARALLEL DO

       !   Hi j face
       j = ny-1
       !$OMP PARALLEL DO PRIVATE(i,k)
       do k = 0,nz-1
             do i = 0,nx-1
                if (bc_dirichlet(mm(i,j,k),2,+1)) then
                   flux(i,ny,k) =  &
                          ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(4,i,j,k)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(3,i,j-1,k)*(uu(i,j-1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,-1)) &
                        flux(i,ny,k) =  flux(i,ny,k) + ss(YBC,i,j,k)*(uu(i,j-2,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),2,+1)) then
                   flux(i,ny,k) = ss(3,i,j,k)*uu(i,j+1,k)
                else
                   flux(i,ny,k) = ss(3,i,j,k)*(uu(i,j+1,k)-uu(i,j,k))
                end if
             end do
       end do
       !$OMP END PARALLEL DO

    else if ( dim == 3 ) then

       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 1,nz-2
          do j = 0,ny-1
             do i = 0,nx-1
                flux(i,j,k) = ss(6,i,j,k) * (uu(i,j,k)-uu(i,j,k-1))
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       !
       ! Must make sure we use stencil from interior fine cell, not fine cell
       ! next to c/f boundary. Because we use ss(2,i,j) which looks "down", we
       ! only modify at the high side
       !
       !$OMP PARALLEL DO PRIVATE(i,j)
       do j = 0,ny-1
          do i = 0,nx-1
             flux(i,j,nz-1) = ss(5,i,j,nz-2) * (uu(i,j,nz-1)-uu(i,j,nz-2))
          end do
       end do
       !$OMP END PARALLEL DO

       !   Lo k face
       k = 0
       !$OMP PARALLEL DO PRIVATE(i,j)
       do j = 0,ny-1
             do i = 0,nx-1
                if (bc_dirichlet(mm(i,j,k),3,-1)) then
                   flux(i,j,0) =  &
                          ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(6,i,j,k+1)*(uu(i,j,k+1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,+1)) &
                        flux(i,j,0) =  flux(i,j,0) + ss(ZBC,i,j,k)*(uu(i,j,k+2)-uu(i,j,k)) 
                   flux(i,j,0) = -flux(i,j,0)
                else if (bc_neumann(mm(i,j,k),3,-1)) then
                   flux(i,j,0) = -ss(6,i,j,k)*uu(i,j,k-1)
                else 
                   flux(i,j,0) = ss(6,i,j,k)*(uu(i,j,k)-uu(i,j,k-1))
                end if
             end do
       end do
       !$OMP END PARALLEL DO

       !   Hi k face
       k = nz-1
       !$OMP PARALLEL DO PRIVATE(i,j)
       do j = 0,ny-1
             do i = 0,nx-1
                if (bc_dirichlet(mm(i,j,k),3,+1)) then
                   flux(i,j,nz) =  &
                          ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(6,i,j,k)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(5,i,j,k-1)*(uu(i,j,k-1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,-1)) &
                        flux(i,j,nz) =  flux(i,j,nz) + ss(ZBC,i,j,k)*(uu(i,j,k-2)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),3,+1)) then
                   flux(i,j,nz) = ss(5,i,j,k)*uu(i,j,k+1)
                else
                   flux(i,j,nz) = ss(5,i,j,k)*(uu(i,j,k+1)-uu(i,j,k))
                end if
             end do
       end do
       !$OMP END PARALLEL DO

    end if

  end subroutine stencil_all_flux_3d

  subroutine stencil_all_flux_ibc_2d(ss, flux, ngf, uu, ngu, dim)
    integer, intent(in) :: ngu, ngf
    real (kind = dp_t), intent(in ) ::   uu(-ngu:,-ngu:)
    real (kind = dp_t), intent(out) :: flux(-ngf:,-ngf:)
    real (kind = dp_t), intent(in ) :: ss
    integer, intent(in) :: dim

    integer :: i, j, nx, ny

    nx = size(flux,dim=1) - 2*ngf
    ny = size(flux,dim=2) - 2*ngf
    
    if ( dim .eq. 1 ) then
       do j = 0,ny-1
          do i = 0,nx-1
             flux(i,j) = ss * (uu(i,j)-uu(i-1,j)) 
          end do
       end do

    else
       do j = 0,ny-1
          do i = 0,nx-1
             flux(i,j) = ss * (uu(i,j)-uu(i,j-1)) 
          end do
       end do
    end if
  end subroutine stencil_all_flux_ibc_2d

  subroutine stencil_all_flux_ibc_3d(ss, flux, ngf, uu, ngu, dim)
    integer, intent(in) :: ngu, ngf
    real (kind = dp_t), intent(in ) ::   uu(-ngu:,-ngu:,-ngu:)
    real (kind = dp_t), intent(out) :: flux(-ngf:,-ngf:,-ngf:)
    real (kind = dp_t), intent(in ) :: ss
    integer, intent(in) :: dim

    integer :: i, j, k, nx, ny, nz

    nx = size(flux,dim=1) - 2*ngf
    ny = size(flux,dim=2) - 2*ngf
    nz = size(flux,dim=3) - 2*ngf
    
    if ( dim .eq. 1 ) then
       !$omp parallel do private(i,j,k) collapse(2)
       do k = 0,nz-1
          do j = 0,ny-1
             do i = 0,nx-1
                flux(i,j,k) = ss * (uu(i,j,k)-uu(i-1,j,k)) 
             end do
          end do
       end do
       !$omp end parallel do
    else if ( dim .eq. 2 ) then
       !$omp parallel do private(i,j,k) collapse(2)
       do k = 0,nz-1
          do j = 0,ny-1
             do i = 0,nx-1
                flux(i,j,k) = ss * (uu(i,j,k)-uu(i,j-1,k)) 
             end do
          end do
       end do
       !$omp end parallel do
    else
       !$omp parallel do private(i,j,k) collapse(2)
       do k = 0,nz-1
          do j = 0,ny-1
             do i = 0,nx-1
                flux(i,j,k) = ss * (uu(i,j,k)-uu(i,j,k-1))                 
             end do
          end do
       end do
       !$omp end parallel do
    end if
  end subroutine stencil_all_flux_ibc_3d


  ! subroutine stencil_apply_n_2d(ss, dd, ng_d, uu, ng_u, mm, lo, hi, skwd)

  !   integer           , intent(in   ) :: ng_d, ng_u, lo(:), hi(:)
  !   real (kind = dp_t), intent(in   ) :: ss(0:,lo(1):,lo(2):)
  !   real (kind = dp_t), intent(  out) :: dd(lo(1)-ng_d:,lo(2)-ng_d:)
  !   real (kind = dp_t), intent(in   ) :: uu(lo(1)-ng_u:,lo(2)-ng_u:)
  !   integer           , intent(in   )  :: mm(lo(1):,lo(2):)
  !   logical           , intent(in   ), optional :: skwd

  !   integer i,j,n,nc,dm,nm1,nedge,nset
    
  !   integer, parameter :: XBC = 6, YBC = 7

  !   logical :: lskwd

  !   lskwd = .true.; if ( present(skwd) ) lskwd = skwd
    
  !   dm    = 2
  !   nset  = 1+3*dm
  !   nc    = (size(ss,dim=1)-1)/(nset+1)
  !   nedge = nc*nset

  !   do j = lo(2),hi(2)
  !      do i = lo(1),hi(1)
  !         dd(i,j) = ss(0,i,j)*uu(i,j)
  !      end do
  !   end do

  !   do n = 1,nc
  !      nm1 = (n-1)*nset
  !      do j = lo(2),hi(2)
  !         do i = lo(1),hi(1)
  !            dd(i,j) = dd(i,j) + &
  !                 (ss(1+nm1,i,j)*uu(i,j) &
  !                 + ss(2+nm1,i,j)*uu(i+1,j  ) + ss(3+nm1,i,j)*uu(i-1,j  ) &
  !                 + ss(4+nm1,i,j)*uu(i  ,j+1) + ss(5+nm1,i,j)*uu(i  ,j-1) &
  !                 )/ss(nedge+n,i,j)
  !         end do
  !      end do

  !      if ( lskwd ) then
  !      ! Corrections for skewed stencils
  !      if (hi(1) > lo(1)) then
  !         do j = lo(2),hi(2)

  !            i = lo(1)
  !            if (bc_skewed(mm(i,j),1,+1)) then
  !               dd(i,j) = dd(i,j) + ss(XBC+nm1,i,j)*uu(i+2,j)
  !            end if

  !            i = hi(1)
  !            if (bc_skewed(mm(i,j),1,-1)) then
  !               dd(i,j) = dd(i,j) + ss(XBC+nm1,i,j)*uu(i-2,j)
  !            end if
  !         end do
  !      end if

  !      if (hi(2) > lo(2)) then
  !         do i = lo(1),hi(1)

  !            j = lo(2)
  !            if (bc_skewed(mm(i,j),2,+1)) then
  !               dd(i,j) = dd(i,j) + ss(YBC+nm1,i,j)*uu(i,j+2)
  !            end if

  !            j = hi(2)
  !            if (bc_skewed(mm(i,j),2,-1)) then
  !               dd(i,j) = dd(i,j) + ss(YBC+nm1,i,j)*uu(i,j-2)
  !            end if
  !         end do
  !      end if
  !      end if
  !   end do

  ! end subroutine stencil_apply_n_2d


!   subroutine stencil_flux_n_2d(ss, flux, uu, mm, ng, ratio, face, dim, skwd)
!     integer, intent(in) :: ng
!     real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:)
!     real (kind = dp_t), intent(out) :: flux(:,:,1:)
!     real (kind = dp_t), intent(in ) :: ss(0:,:,:)
!     integer           , intent(in)  :: mm(:,:)
!     logical, intent(in), optional :: skwd
!     integer, intent(in) :: ratio, face, dim
!     integer nx,ny,dm,nc,nedge,nm1,nset
!     integer i,j,ic,jc,n
!     real (kind = dp_t) :: fac
!     integer, parameter :: XBC = 6, YBC = 7
!     logical :: lskwd

!     lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

!     nx = size(ss,dim=2)
!     ny = size(ss,dim=3)

!     dm    = 2
!     nset  = 1+3*dm
!     nc    = (size(ss,dim=1)-1)/(nset+1)
!     nedge = nc*nset

!     !   Note that one factor of ratio is the tangential averaging, while the
!     !     other is the normal factor
!     fac = ONE/real(ratio*ratio, kind=dp_t)

! !   Lo i face
!     if ( dim == 1 ) then
!        if (face == -1) then

!           i = 1
!           flux(1,:,:) = ZERO
!           do n = 1,nc
!              nm1  = (n-1)*nset
!              do j = 1,ny
!                 jc = (j-1)/ratio+1
!                 if (bc_dirichlet(mm(i,j),1,-1)) then
!                    flux(1,jc,n) = flux(1,jc,n)  &
!                      + ss(2+nm1,i,j)*(uu(i+1,j)-uu(i,j)) &
!                      + ss(3+nm1,i,j)*(uu(i-1,j)-uu(i,j)) - ss(3+nm1,i+1,j)*(uu(i+1,j)-uu(i,j))
!                    if (bc_skewed(mm(i,j),1,+1)) &
!                         flux(1,jc,n) = flux(1,jc,n) + ss(XBC+nm1,i,j)*(uu(i+2,j)-uu(i,j)) 
!                 else   
!                    flux(1,jc,n) = Huge(flux(:,:,n))
!                 end if
!              end do
!              flux(1,:,n) = fac * flux(1,:,n)
!           end do

! !      Hi i face
!        else if (face == 1) then

!           i = nx
!           flux(1,:,:) = ZERO
!           do n = 1,nc
!              nm1 = (n-1)*nset
!              do j = 1,ny
!                 jc = (j-1)/ratio+1
!                 if (bc_dirichlet(mm(i,j),1,+1)) then
!                    flux(1,jc,n) = flux(1,jc,n) &
!                      + ss(2+nm1,i,j)*(uu(i+1,j)-uu(i,j)) &
!                      + ss(3+nm1,i,j)*(uu(i-1,j)-uu(i,j)) - ss(2+nm1,i-1,j)*(uu(i-1,j)-uu(i,j))
!                    if (bc_skewed(mm(i,j),1,-1)) &
!                      flux(1,jc,n) = flux(1,jc,n) + ss(XBC+nm1,i,j)*(uu(i-2,j)-uu(i,j))
!                 else 
!                    flux(1,jc,n) = Huge(flux(:,:,n))
!                 end if
!              end do
!              flux(1,:,n) = fac * flux(1,:,n)
!           end do

!        end if

! !   Lo j face
!     else if ( dim == 2 ) then
!        if (face == -1) then

!           j = 1
!           flux(:,1,:) = ZERO
!           do n = 1,nc
!              nm1 = (n-1)*nset
!              do i = 1,nx
!                 ic = (i-1)/ratio+1
!                 if (bc_dirichlet(mm(i,j),2,-1)) then
!                    flux(ic,1,n) = flux(ic,1,n)  &
!                         + ss(4+nm1,i,j)*(uu(i,j+1)-uu(i,j)) &
!                         + ss(5+nm1,i,j)*(uu(i,j-1)-uu(i,j)) - ss(5+nm1,i,j+1)*(uu(i,j+1)-uu(i,j))
!                    if (bc_skewed(mm(i,j),2,+1)) &
!                         flux(ic,1,n) =  flux(ic,1,n) + ss(YBC+nm1,i,j)*(uu(i,j+2)-uu(i,j))
!                 else 
!                    flux(ic,1,n) = Huge(flux(:,:,n))
!                 end if
!              end do
!              flux(:,1,n) = fac * flux(:,1,n)
!           end do


! !      Hi j face
!        else if (face == 1) then

!           j = ny
!           flux(:,1,:) = ZERO
!           do n = 1,nc
!              nm1 = (n-1)*nset
!              do i = 1,nx
!                 ic = (i-1)/ratio+1
!                 if (bc_dirichlet(mm(i,j),2,+1)) then
!                    flux(ic,1,n) = flux(ic,1,n)  &
!                      + ss(4+nm1,i,j)*(uu(i,j+1)-uu(i,j)) &
!                      + ss(5+nm1,i,j)*(uu(i,j-1)-uu(i,j)) - ss(4+nm1,i,j-1)*(uu(i,j-1)-uu(i,j))
!                    if (bc_skewed(mm(i,j),2,-1)) &
!                      flux(ic,1,n) = flux(ic,1,n) + ss(YBC+nm1,i,j)*(uu(i,j-2)-uu(i,j))
!                 else
!                    flux(ic,1,n) = Huge(flux(:,:,n))
!                 end if
!              end do
!              flux(:,1,n) = fac * flux(:,1,n)
!           end do

!        end if
!     end if

!   end subroutine stencil_flux_n_2d

end module cc_stencil_apply_module
