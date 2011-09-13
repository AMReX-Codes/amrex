module cc_stencil_apply_module

  use bl_types
  use bc_module
  use bc_functions_module
  use multifab_module

  implicit none

  real(kind=dp_t), parameter, private :: ZERO = 0.0_dp_t
  real(kind=dp_t), parameter, private :: ONE  = 1.0_dp_t

  public  :: stencil_apply_1d, stencil_apply_2d, stencil_apply_3d
  private :: stencil_dense_apply_1d, stencil_dense_apply_2d, stencil_dense_apply_3d
  private :: stencil_all_flux_1d, stencil_all_flux_2d, stencil_all_flux_3d

contains

  subroutine stencil_apply_1d(ss, dd, ng_d, uu, ng_u, mm, lo, hi, skwd)

    integer, intent(in) :: ng_d, ng_u, lo(:), hi(:)
    real (kind = dp_t), intent(in)  :: ss(lo(1)     :,0:)
    real (kind = dp_t), intent(out) :: dd(lo(1)-ng_d:)
    real (kind = dp_t), intent(in)  :: uu(lo(1)-ng_u:)
    integer           , intent(in)  :: mm(lo(1):)
    logical, intent(in), optional   :: skwd

    integer, parameter :: XBC = 3
    logical :: lskwd
    integer :: i
   
    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    do i = lo(1),hi(1)
       dd(i) = ss(i,0)*uu(i) + ss(i,1)*uu(i+1) + ss(i,2)*uu(i-1)
    end do

    if ( lskwd ) then
       if (hi(1) > lo(1)) then
          i = lo(1)
          if (bc_skewed(mm(i),1,+1)) then
             dd(i) = dd(i) + ss(i,XBC)*uu(i+2)
          end if
  
          i = hi(1)
          if (bc_skewed(mm(i),1,-1)) then
             dd(i) = dd(i) + ss(i,XBC)*uu(i-2)
          end if
       end if
    end if

  end subroutine stencil_apply_1d

  subroutine stencil_flux_1d(ss, flux, uu, mm, ng, ratio, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in)  :: ss(:,0:)
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

    nx = size(ss,dim=1)

    !   This factor is dx^fine / dx^crse
    fac = ONE / real(ratio, kind=dp_t)

    if ( dim == 1 ) then
       if ( face == -1 ) then
          i = 1
          if (bc_dirichlet(mm(1),1,-1)) then
             flux(1) = ss(i,1)*(uu(i+1)-uu(i)) + ss(i,2)*(uu(i-1)-uu(i)) &
                  - ss(i+1,2)*(uu(i+1)-uu(i))
             if (bc_skewed(mm(i),1,+1)) then
                flux(1) =  flux(1) + ss(i,XBC)*uu(i+2)
             end if
          else 
             flux(1) = Huge(flux)
          end if
          flux(1) = fac*flux(1)
       else if ( face == 1 ) then
          i = nx
          if (bc_dirichlet(mm(i),1,+1)) then
             flux(1) = ss(i,1)*(uu(i+1)-uu(i)) + ss(i,2)*(uu(i-1)-uu(i)) &
                  - ss(i-1,1)*(uu(i-1)-uu(i))
             if (bc_skewed(mm(i),1,-1)) then
                flux(1) =  flux(1) + ss(i,XBC)*uu(i-2)
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
    real (kind = dp_t), intent(in   ) :: ss(lo(1):,lo(2):,0:)
    real (kind = dp_t), intent(  out) :: dd(lo(1)-ng_d:,lo(2)-ng_d:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng_u:,lo(2)-ng_u:)
    integer           , intent(in   )  :: mm(lo(1):,lo(2):)
    logical           , intent(in   ), optional :: skwd

    integer i,j

    integer, parameter :: XBC = 5, YBC = 6

    logical :: lskwd

    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    ! This is the Minion 4th order cross stencil.
    if (size(ss,dim=3) .eq. 9) then

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            dd(i,j) = &
                   ss(i,j,0) * uu(i,j) &
                 + ss(i,j,1) * uu(i-2,j) + ss(i,j,2) * uu(i-1,j) &
                 + ss(i,j,3) * uu(i+1,j) + ss(i,j,4) * uu(i+2,j) &
                 + ss(i,j,5) * uu(i,j-2) + ss(i,j,6) * uu(i,j-1) &
                 + ss(i,j,7) * uu(i,j+1) + ss(i,j,8) * uu(i,j+2)
          end do
       end do

    ! This is the Minion 4th order full stencil.
    else if (size(ss,dim=3) .eq. 25) then

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            dd(i,j) = ss(i,j, 0) * uu(i,j) &
                    + ss(i,j, 1) * uu(i-2,j-2) + ss(i,j, 2) * uu(i-1,j-2) & ! AT J-2
                    + ss(i,j, 3) * uu(i  ,j-2) + ss(i,j, 4) * uu(i+1,j-2) & ! AT J-2
                    + ss(i,j, 5) * uu(i+2,j-2)                            & ! AT J-2
                    + ss(i,j, 6) * uu(i-2,j-1) + ss(i,j, 7) * uu(i-1,j-1) & ! AT J-1
                    + ss(i,j, 8) * uu(i  ,j-1) + ss(i,j, 9) * uu(i+1,j-1) & ! AT J-1
                    + ss(i,j,10) * uu(i+2,j-1)                            & ! AT J-1
                    + ss(i,j,11) * uu(i-2,j  ) + ss(i,j,12) * uu(i-1,j  ) & ! AT J
                    + ss(i,j,13) * uu(i+1,j  ) + ss(i,j,14) * uu(i+2,j  ) & ! AT J
                    + ss(i,j,15) * uu(i-2,j+1) + ss(i,j,16) * uu(i-1,j+1) & ! AT J+1
                    + ss(i,j,17) * uu(i  ,j+1) + ss(i,j,18) * uu(i+1,j+1) & ! AT J+1
                    + ss(i,j,19) * uu(i+2,j+1)                            & ! AT J+1
                    + ss(i,j,20) * uu(i-2,j+2) + ss(i,j,21) * uu(i-1,j+2) & ! AT J+2
                    + ss(i,j,22) * uu(i  ,j+2) + ss(i,j,23) * uu(i+1,j+2) & ! AT J+2
                    + ss(i,j,24) * uu(i+2,j+2)                              ! AT J+2
          end do
       end do

    ! This is our standard 5-point Laplacian with a possible correction at boundaries
    else 

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             dd(i,j) = ss(i,j,0)*uu(i,j) &
                  + ss(i,j,1)*uu(i+1,j  ) + ss(i,j,2)*uu(i-1,j  ) &
                  + ss(i,j,3)*uu(i  ,j+1) + ss(i,j,4)*uu(i  ,j-1)
          end do
       end do

       if ( lskwd ) then
       ! Corrections for skewed stencils
       if (hi(1) > lo(1)) then
          do j = lo(2),hi(2)

             i = lo(1)
             if (bc_skewed(mm(i,j),1,+1)) then
                dd(i,j) = dd(i,j) + ss(i,j,XBC)*uu(i+2,j)
             end if

             i = hi(1)
             if (bc_skewed(mm(i,j),1,-1)) then
                dd(i,j) = dd(i,j) + ss(i,j,XBC)*uu(i-2,j)
             end if
          end do
       end if

       if (hi(2) > lo(2)) then
          do i = lo(1),hi(1)

             j = lo(2)
             if (bc_skewed(mm(i,j),2,+1)) then
                dd(i,j) = dd(i,j) + ss(i,j,YBC)*uu(i,j+2)
             end if

             j = hi(2)
             if (bc_skewed(mm(i,j),2,-1)) then
                dd(i,j) = dd(i,j) + ss(i,j,YBC)*uu(i,j-2)
             end if

          end do
       end if
       end if
    end if

  end subroutine stencil_apply_2d

subroutine stencil_apply_n_2d(ss, dd, ng_d, uu, ng_u, mm, lo, hi, skwd)
    integer           , intent(in   ) :: ng_d, ng_u, lo(:), hi(:)
    real (kind = dp_t), intent(in   ) :: ss(lo(1):,lo(2):,0:)
    real (kind = dp_t), intent(  out) :: dd(lo(1)-ng_d:,lo(2)-ng_d:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng_u:,lo(2)-ng_u:)
    integer           , intent(in   )  :: mm(lo(1):,lo(2):)
    logical           , intent(in   ), optional :: skwd

    integer i,j,n,nc,dm,nm1,nedge,nset
    
    integer, parameter :: XBC = 6, YBC = 7

    logical :: lskwd

    lskwd = .true.; if ( present(skwd) ) lskwd = skwd
    
    dm    = 2
    nset  = 1+3*dm
    nc    = (size(ss,dim=3)-1)/(nset+1)
    nedge = nc*nset

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          dd(i,j) = ss(i,j,0)*uu(i,j)
       end do
    end do

    do n = 1,nc
       nm1 = (n-1)*nset
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             dd(i,j) = dd(i,j) + &
                  (ss(i,j,1+nm1)*uu(i,j) &
                  + ss(i,j,2+nm1)*uu(i+1,j  ) + ss(i,j,3+nm1)*uu(i-1,j  ) &
                  + ss(i,j,4+nm1)*uu(i  ,j+1) + ss(i,j,5+nm1)*uu(i  ,j-1) &
                  )/ss(i,j,nedge+n)
          end do
       end do

       if ( lskwd ) then
       ! Corrections for skewed stencils
       if (hi(1) > lo(1)) then
          do j = lo(2),hi(2)

             i = lo(1)
             if (bc_skewed(mm(i,j),1,+1)) then
                dd(i,j) = dd(i,j) + ss(i,j,XBC+nm1)*uu(i+2,j)
             end if

             i = hi(1)
             if (bc_skewed(mm(i,j),1,-1)) then
                dd(i,j) = dd(i,j) + ss(i,j,XBC+nm1)*uu(i-2,j)
             end if
          end do
       end if

       if (hi(2) > lo(2)) then
          do i = lo(1),hi(1)

             j = lo(2)
             if (bc_skewed(mm(i,j),2,+1)) then
                dd(i,j) = dd(i,j) + ss(i,j,YBC+nm1)*uu(i,j+2)
             end if

             j = hi(2)
             if (bc_skewed(mm(i,j),2,-1)) then
                dd(i,j) = dd(i,j) + ss(i,j,YBC+nm1)*uu(i,j-2)
             end if
          end do
       end if
       end if
    end do

  end subroutine stencil_apply_n_2d

  subroutine stencil_flux_2d(ss, flux, uu, mm, ng, ratio, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:)
    real (kind = dp_t), intent(out) :: flux(:,:)
    real (kind = dp_t), intent(in ) :: ss(:,:,0:)
    integer           , intent(in)  :: mm(:,:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: ratio, face, dim
    integer nx,ny
    integer i,j,ic,jc
    real (kind = dp_t) :: fac
    integer, parameter :: XBC = 5, YBC = 6
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

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
                     + ss(i,j,1)*(uu(i+1,j)-uu(i,j)) &
                     + ss(i,j,2)*(uu(i-1,j)-uu(i,j)) - ss(i+1,j,2)*(uu(i+1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,+1)) &
                     flux(1,jc) = flux(1,jc) + ss(i,j,XBC)*(uu(i+2,j)-uu(i,j)) 
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
                     + ss(i,j,1)*(uu(i+1,j)-uu(i,j)) &
                     + ss(i,j,2)*(uu(i-1,j)-uu(i,j)) - ss(i-1,j,1)*(uu(i-1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,-1)) &
                     flux(1,jc) = flux(1,jc) + ss(i,j,XBC)*(uu(i-2,j)-uu(i,j))
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
                     + ss(i,j,3)*(uu(i,j+1)-uu(i,j)) &
                     + ss(i,j,4)*(uu(i,j-1)-uu(i,j)) - ss(i,j+1,4)*(uu(i,j+1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,+1)) &
                     flux(ic,1) =  flux(ic,1) + ss(i,j,YBC)*(uu(i,j+2)-uu(i,j))
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
                     + ss(i,j,3)*(uu(i,j+1)-uu(i,j)) &
                     + ss(i,j,4)*(uu(i,j-1)-uu(i,j)) - ss(i,j-1,3)*(uu(i,j-1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,-1)) &
                     flux(ic,1) = flux(ic,1) + ss(i,j,YBC)*(uu(i,j-2)-uu(i,j))
             else
                flux(ic,1) = Huge(flux)
             end if
          end do
          flux(:,1) = fac * flux(:,1)

       end if
    end if

  end subroutine stencil_flux_2d

  subroutine stencil_flux_n_2d(ss, flux, uu, mm, ng, ratio, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:)
    real (kind = dp_t), intent(out) :: flux(:,:,1:)
    real (kind = dp_t), intent(in ) :: ss(:,:,0:)
    integer           , intent(in)  :: mm(:,:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: ratio, face, dim
    integer nx,ny,dm,nc,nedge,nm1,nset
    integer i,j,ic,jc,n
    real (kind = dp_t) :: fac
    integer, parameter :: XBC = 6, YBC = 7
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

    dm    = 2
    nset  = 1+3*dm
    nc    = (size(ss,dim=3)-1)/(nset+1)
    nedge = nc*nset

    !   Note that one factor of ratio is the tangential averaging, while the
    !     other is the normal factor
    fac = ONE/real(ratio*ratio, kind=dp_t)

!   Lo i face
    if ( dim == 1 ) then
       if (face == -1) then

          i = 1
          flux(1,:,:) = ZERO
          do n = 1,nc
             nm1  = (n-1)*nset
             do j = 1,ny
                jc = (j-1)/ratio+1
                if (bc_dirichlet(mm(i,j),1,-1)) then
                   flux(1,jc,n) = flux(1,jc,n)  &
                     + ss(i,j,2+nm1)*(uu(i+1,j)-uu(i,j)) &
                     + ss(i,j,3+nm1)*(uu(i-1,j)-uu(i,j)) - ss(i+1,j,3+nm1)*(uu(i+1,j)-uu(i,j))
                   if (bc_skewed(mm(i,j),1,+1)) &
                        flux(1,jc,n) = flux(1,jc,n) + ss(i,j,XBC+nm1)*(uu(i+2,j)-uu(i,j)) 
                else   
                   flux(1,jc,n) = Huge(flux(:,:,n))
                end if
             end do
             flux(1,:,n) = fac * flux(1,:,n)
          end do

!      Hi i face
       else if (face == 1) then

          i = nx
          flux(1,:,:) = ZERO
          do n = 1,nc
             nm1 = (n-1)*nset
             do j = 1,ny
                jc = (j-1)/ratio+1
                if (bc_dirichlet(mm(i,j),1,+1)) then
                   flux(1,jc,n) = flux(1,jc,n) &
                     + ss(i,j,2+nm1)*(uu(i+1,j)-uu(i,j)) &
                     + ss(i,j,3+nm1)*(uu(i-1,j)-uu(i,j)) - ss(i-1,j,2+nm1)*(uu(i-1,j)-uu(i,j))
                   if (bc_skewed(mm(i,j),1,-1)) &
                     flux(1,jc,n) = flux(1,jc,n) + ss(i,j,XBC+nm1)*(uu(i-2,j)-uu(i,j))
                else 
                   flux(1,jc,n) = Huge(flux(:,:,n))
                end if
             end do
             flux(1,:,n) = fac * flux(1,:,n)
          end do

       end if

!   Lo j face
    else if ( dim == 2 ) then
       if (face == -1) then

          j = 1
          flux(:,1,:) = ZERO
          do n = 1,nc
             nm1 = (n-1)*nset
             do i = 1,nx
                ic = (i-1)/ratio+1
                if (bc_dirichlet(mm(i,j),2,-1)) then
                   flux(ic,1,n) = flux(ic,1,n)  &
                        + ss(i,j,4+nm1)*(uu(i,j+1)-uu(i,j)) &
                        + ss(i,j,5+nm1)*(uu(i,j-1)-uu(i,j)) - ss(i,j+1,5+nm1)*(uu(i,j+1)-uu(i,j))
                   if (bc_skewed(mm(i,j),2,+1)) &
                        flux(ic,1,n) =  flux(ic,1,n) + ss(i,j,YBC+nm1)*(uu(i,j+2)-uu(i,j))
                else 
                   flux(ic,1,n) = Huge(flux(:,:,n))
                end if
             end do
             flux(:,1,n) = fac * flux(:,1,n)
          end do


!      Hi j face
       else if (face == 1) then

          j = ny
          flux(:,1,:) = ZERO
          do n = 1,nc
             nm1 = (n-1)*nset
             do i = 1,nx
                ic = (i-1)/ratio+1
                if (bc_dirichlet(mm(i,j),2,+1)) then
                   flux(ic,1,n) = flux(ic,1,n)  &
                     + ss(i,j,4+nm1)*(uu(i,j+1)-uu(i,j)) &
                     + ss(i,j,5+nm1)*(uu(i,j-1)-uu(i,j)) - ss(i,j-1,4+nm1)*(uu(i,j-1)-uu(i,j))
                   if (bc_skewed(mm(i,j),2,-1)) &
                     flux(ic,1,n) = flux(ic,1,n) + ss(i,j,YBC+nm1)*(uu(i,j-2)-uu(i,j))
                else
                   flux(ic,1,n) = Huge(flux(:,:,n))
                end if
             end do
             flux(:,1,n) = fac * flux(:,1,n)
          end do

       end if
    end if

  end subroutine stencil_flux_n_2d

  subroutine stencil_apply_3d(ss, dd, ng_d, uu, ng_u, mm, skwd)

    integer           , intent(in ) :: ng_d,ng_u
    real (kind = dp_t), intent(in ) :: ss(:,:,:,0:)
    real (kind = dp_t), intent(out) :: dd(1-ng_d:,1-ng_d:,1-ng_d:)
    real (kind = dp_t), intent(in ) :: uu(1-ng_u:,1-ng_u:,1-ng_u:)
    integer           , intent(in ) :: mm(:,:,:)
    logical           , intent(in ), optional :: skwd

    integer nx,ny,nz,i,j,k
    integer, parameter :: XBC = 7, YBC = 8, ZBC = 9
    logical :: lskwd

    lskwd = .true.; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)

    ! This is the 4th order cross stencil for constant coefficients.
    if (size(ss,dim=4) .eq. 13) then
 
       !$OMP PARALLEL DO PRIVATE(i,j,k) IF(nz.ge.4)
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                dd(i,j,k) = ss(i,j,k,0) * uu(i,j,k) &
                     + ss(i,j,k, 1) * uu(i-2,j,k) + ss(i,j,k, 2) * uu(i-1,j,k) &
                     + ss(i,j,k, 3) * uu(i+1,j,k) + ss(i,j,k, 4) * uu(i+2,j,k) &
                     + ss(i,j,k, 5) * uu(i,j-2,k) + ss(i,j,k, 6) * uu(i,j-1,k) &
                     + ss(i,j,k, 7) * uu(i,j+1,k) + ss(i,j,k, 8) * uu(i,j+2,k) &
                     + ss(i,j,k, 9) * uu(i,j,k-2) + ss(i,j,k,10) * uu(i,j,k-1) &
                     + ss(i,j,k,11) * uu(i,j,k+1) + ss(i,j,k,12) * uu(i,j,k+2)
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    ! This is the 4th order cross stencil for variable coefficients.
    else if (size(ss,dim=4) .eq. 61) then

       !$OMP PARALLEL DO PRIVATE(i,j,k) IF(nz.ge.4)
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                dd(i,j,k) = &
                       ss(i,j,k, 0) * uu(i  ,j  ,k  ) &
                       ! Contributions from k-2
                     + ss(i,j,k, 1) * uu(i  ,j-2,k-2) + ss(i,j,k, 2) * uu(i  ,j-1,k-2) &
                     + ss(i,j,k, 3) * uu(i-2,j  ,k-2) + ss(i,j,k, 4) * uu(i-1,j  ,k-2) &
                     + ss(i,j,k, 5) * uu(i  ,j  ,k-2) + ss(i,j,k, 6) * uu(i+1,j  ,k-2) &
                     + ss(i,j,k, 7) * uu(i+2,j  ,k-2) + ss(i,j,k, 8) * uu(i  ,j+1,k-2) &
                     + ss(i,j,k, 9) * uu(i  ,j+2,k-2)                                  &
                       ! Contributions from k-1
                     + ss(i,j,k,10) * uu(i  ,j-2,k-1) + ss(i,j,k,11) * uu(i  ,j-1,k-1) &
                     + ss(i,j,k,12) * uu(i-2,j  ,k-1) + ss(i,j,k,13) * uu(i-1,j  ,k-1) &
                     + ss(i,j,k,14) * uu(i  ,j  ,k-1) + ss(i,j,k,15) * uu(i+1,j  ,k-1) &
                     + ss(i,j,k,16) * uu(i+2,j  ,k-1) + ss(i,j,k,17) * uu(i  ,j+1,k-1) &
                     + ss(i,j,k,18) * uu(i  ,j+2,k-1)                                  &
                       ! Contributions from j-2,k
                     + ss(i,j,k,19) * uu(i-2,j-2,k  ) + ss(i,j,k,20) * uu(i-1,j-2,k  ) &
                     + ss(i,j,k,21) * uu(i  ,j-2,k  ) + ss(i,j,k,22) * uu(i+1,j-2,k  ) &
                     + ss(i,j,k,23) * uu(i+2,j-2,k  ) 

                dd(i,j,k) = dd(i,j,k) &
                       ! Contributions from j-1,k
                     + ss(i,j,k,24) * uu(i-2,j-1,k  ) + ss(i,j,k,25) * uu(i-1,j-1,k  ) &
                     + ss(i,j,k,26) * uu(i  ,j-1,k  ) + ss(i,j,k,27) * uu(i+1,j-1,k  ) &
                     + ss(i,j,k,28) * uu(i+2,j-1,k  )                                  &
                       ! Contributions from j  ,k
                     + ss(i,j,k,29) * uu(i-2,j  ,k  ) + ss(i,j,k,30) * uu(i-1,j  ,k  ) &
                                                      + ss(i,j,k,31) * uu(i+1,j  ,k  ) &
                     + ss(i,j,k,32) * uu(i+2,j  ,k  )                                  &
                       ! Contributions from j+1,k
                     + ss(i,j,k,33) * uu(i-2,j+1,k  ) + ss(i,j,k,34) * uu(i-1,j+1,k  ) &
                     + ss(i,j,k,35) * uu(i  ,j+1,k  ) + ss(i,j,k,36) * uu(i+1,j+1,k  ) &
                     + ss(i,j,k,37) * uu(i+2,j+1,k  )                                  &
                       ! Contributions from j+2,k
                     + ss(i,j,k,38) * uu(i-2,j+2,k  ) + ss(i,j,k,39) * uu(i-1,j+2,k  ) &
                     + ss(i,j,k,40) * uu(i  ,j+2,k  ) + ss(i,j,k,41) * uu(i+1,j+2,k  ) &
                     + ss(i,j,k,42) * uu(i+2,j+2,k  )                                  

                dd(i,j,k) = dd(i,j,k) &
                       ! Contributions from k+1
                     + ss(i,j,k,43) * uu(i  ,j-2,k+1) + ss(i,j,k,44) * uu(i  ,j-1,k+1) &
                     + ss(i,j,k,45) * uu(i-2,j  ,k+1) + ss(i,j,k,46) * uu(i-1,j  ,k+1) &
                     + ss(i,j,k,47) * uu(i  ,j  ,k+1) + ss(i,j,k,48) * uu(i+1,j  ,k+1) &
                     + ss(i,j,k,49) * uu(i+2,j  ,k+1) + ss(i,j,k,50) * uu(i  ,j+1,k+1) &
                     + ss(i,j,k,51) * uu(i  ,j+2,k+1)                                  &
                       ! Contributions from k+2
                     + ss(i,j,k,52) * uu(i  ,j-2,k+2) + ss(i,j,k,53) * uu(i  ,j-1,k+2) &
                     + ss(i,j,k,54) * uu(i-2,j  ,k+2) + ss(i,j,k,55) * uu(i-1,j  ,k+2) &
                     + ss(i,j,k,56) * uu(i  ,j  ,k+2) + ss(i,j,k,57) * uu(i+1,j  ,k+2) &
                     + ss(i,j,k,58) * uu(i+2,j  ,k+2) + ss(i,j,k,59) * uu(i  ,j+1,k+2) &
                     + ss(i,j,k,60) * uu(i  ,j+2,k+2)                                   

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    ! This is the 2nd order cross stencil.
    else 

       !$OMP PARALLEL DO PRIVATE(i,j,k) IF(nz.ge.4)
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                dd(i,j,k) = &
                     ss(i,j,k,0)*uu(i,j,k)       + &
                     ss(i,j,k,1)*uu(i+1,j  ,k  ) + &
                     ss(i,j,k,2)*uu(i-1,j  ,k  ) + &
                     ss(i,j,k,3)*uu(i  ,j+1,k  ) + &
                     ss(i,j,k,4)*uu(i  ,j-1,k  ) + &
                     ss(i,j,k,5)*uu(i  ,j  ,k+1) + &
                     ss(i,j,k,6)*uu(i  ,j  ,k-1)
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    end if

    if ( lskwd ) then
       !
       ! Corrections for skewed stencils
       !
       if (nx > 1) then
          do k = 1, nz
             do j = 1, ny
                i = 1
                if (bc_skewed(mm(i,j,k),1,+1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(i,j,k,XBC)*uu(i+2,j,k)
                end if

                i = nx
                if (bc_skewed(mm(i,j,k),1,-1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(i,j,k,XBC)*uu(i-2,j,k)
                end if
             end do
          end do
       end if

       if (ny > 1) then
          do k = 1,nz
             do i = 1,nx
                j = 1
                if (bc_skewed(mm(i,j,k),2,+1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(i,j,k,YBC)*uu(i,j+2,k)
                end if

                j = ny
                if (bc_skewed(mm(i,j,k),2,-1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(i,j,k,YBC)*uu(i,j-2,k)
                end if
             end do
          end do
       end if

       if (nz > 1) then
          do j = 1,ny
             do i = 1,nx
                k = 1
                if (bc_skewed(mm(i,j,k),3,+1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(i,j,k,ZBC)*uu(i,j,k+2)
                end if

                k = nz
                if (bc_skewed(mm(i,j,k),3,-1)) then
                   dd(i,j,k) = dd(i,j,k) + ss(i,j,k,ZBC)*uu(i,j,k-2)
                end if
             end do
          end do
       end if
    end if
  end subroutine stencil_apply_3d

  subroutine stencil_flux_3d(ss, flux, uu, mm, ng, ratio, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), intent(out) :: flux(:,:,:)
    real (kind = dp_t), intent(in ) :: ss(:,:,:,0:)
    integer           , intent(in)  :: mm(:,:,:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: ratio, face, dim
    integer nx, ny, nz
    integer i,j,k,ic,jc,kc
    integer, parameter :: XBC = 7, YBC = 8, ZBC = 9
    real (kind = dp_t) :: fac
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)

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

          i = 1
          flux(1,:,:) = ZERO
          do k = 1,nz
             do j = 1,ny
                jc = (j-1)/ratio + 1
                kc = (k-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),1,-1)) then
                   flux(1,jc,kc) =  flux(1,jc,kc) &
                        + ss(i,j,k,1)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(i,j,k,2)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(i+1,j,k,2)*(uu(i+1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,+1)) &
                        flux(1,jc,kc) =  flux(1,jc,kc) + ss(i,j,k,XBC)*(uu(i+2,j,k)-uu(i,j,k))
                else 
                   flux(1,jc,kc) = Huge(flux)
                end if
             end do
          end do
          flux(1,:,:) = flux(1,:,:) * fac

          !   Hi i face
       else if (face ==  1) then

          i = nx
          flux(1,:,:) = ZERO
          do k = 1,nz
             do j = 1,ny
                jc = (j-1)/ratio + 1
                kc = (k-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),1,+1)) then
                   flux(1,jc,kc) =  flux(1,jc,kc) &
                        + ss(i,j,k,1)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(i,j,k,2)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(i-1,j,k,1)*(uu(i-1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,-1)) &
                        flux(1,jc,kc) =  flux(1,jc,kc) + ss(i,j,k,XBC)*(uu(i-2,j,k)-uu(i,j,k))
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
          j = 1
          flux(:,1,:) = ZERO
          do k = 1,nz
             do i = 1,nx
                ic = (i-1)/ratio + 1
                kc = (k-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),2,-1)) then
                   flux(ic,1,kc) =  flux(ic,1,kc) &
                        + ss(i,j,k,3)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(i,j,k,4)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(i,j+1,k,4)*(uu(i,j+1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,+1)) &
                        flux(ic,1,kc) =  flux(ic,1,kc) + ss(i,j,k,YBC)*(uu(i,j+2,k)-uu(i,j,k))
                else 
                   flux(ic,1,kc) = Huge(flux)
                end if
             end do
          end do
          flux(:,1,:) = flux(:,1,:) * fac

          !   Hi j face
       else if (face ==  1) then
          j = ny
          flux(:,1,:) = ZERO
          do k = 1,nz
             do i = 1,nx
                ic = (i-1)/ratio + 1
                kc = (k-1)/ratio + 1

                if (bc_dirichlet(mm(i,j,k),2,+1)) then
                   flux(ic,1,kc) =  flux(ic,1,kc) &
                        + ss(i,j,k,3)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(i,j,k,4)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(i,j-1,k,3)*(uu(i,j-1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,-1)) &
                        flux(ic,1,kc) =  flux(ic,1,kc) + ss(i,j,k,YBC)*(uu(i,j-2,k)-uu(i,j,k))
                else
                   flux(ic,1,kc) = Huge(flux)
                end if
             end do
          end do
          flux(:,1,:) = flux(:,1,:) * fac

          !   Lo k face
       end if
    else if ( dim == 3 ) then
       if (face == -1) then

          k = 1
          flux(:,:,1) = ZERO
          do j = 1,ny
             do i = 1,nx
                ic = (i-1)/ratio + 1
                jc = (j-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),3,-1)) then
                   flux(ic,jc,1) =  flux(ic,jc,1) &
                        + ss(i,j,k,5)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(i,j,k,6)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(i,j,k+1,6)*(uu(i,j,k+1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,+1)) &
                        flux(ic,jc,1) =  flux(ic,jc,1) + ss(i,j,k,ZBC)*(uu(i,j,k+2)-uu(i,j,k)) 
                else 
                   flux(ic,jc,1) = Huge(flux)
                end if
             end do
          end do
          flux(:,:,1) = flux(:,:,1) * fac

          !   Hi k face
       else if (face ==  1) then

          k = nz
          flux(:,:,1) = ZERO
          do j = 1,ny
             do i = 1,nx
                ic = (i-1)/ratio + 1
                jc = (j-1)/ratio + 1
                if (bc_dirichlet(mm(i,j,k),3,+1)) then
                   flux(ic,jc,1) =  flux(ic,jc,1) &
                        + ss(i,j,k,5)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(i,j,k,6)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(i,j,k-1,5)*(uu(i,j,k-1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,-1)) &
                        flux(ic,jc,1) =  flux(ic,jc,1) + ss(i,j,k,ZBC)*(uu(i,j,k-2)-uu(i,j,k))
                else
                   flux(ic,jc,1) = Huge(flux)
                end if
             end do
          end do
          flux(:,:,1) = flux(:,:,1) * fac

       end if
    end if

  end subroutine stencil_flux_3d

  subroutine stencil_dense_apply_1d(ss, dd, ng_d, uu, ng_u)
    integer, intent(in) :: ng_d, ng_u
    real (kind = dp_t), intent(in   ) :: ss(:,0:)
    real (kind = dp_t), intent(  out) :: dd(1-ng_d:)
    real (kind = dp_t), intent(in   ) :: uu(1-ng_u:)
    integer i, nx
   
    nx = size(ss,dim=1)
    do i = 1, nx
      dd(i) = ss(i,1)*uu(i-1) + ss(i,0)*uu(i) + ss(i,2)*uu(i+1)
    end do

  end subroutine stencil_dense_apply_1d

  subroutine stencil_dense_apply_2d(ss, dd, ng_d, uu, ng_u)
    integer, intent(in) :: ng_d, ng_u
    real (kind = dp_t), intent(in   ) :: ss(:,:,0:)
    real (kind = dp_t), intent(  out) :: dd(1-ng_d:,1-ng_d:)
    real (kind = dp_t), intent(in   ) :: uu(1-ng_u:,1-ng_u:)
    integer i, j, nx, ny

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

    do j = 1, ny
       do i = 1, nx
          dd(i,j) = &
               + ss(i,j,1)*uu(i-1,j-1) + ss(i,j,2)*uu(i  ,j-1) + ss(i,j,3)*uu(i+1,j-1) &
               + ss(i,j,4)*uu(i-1,j  ) + ss(i,j,0)*uu(i  ,j  ) + ss(i,j,5)*uu(i+1,j  ) &
               + ss(i,j,6)*uu(i-1,j+1) + ss(i,j,7)*uu(i  ,j+1) + ss(i,j,8)*uu(i+1,j+1)
       end do
    end do

  end subroutine stencil_dense_apply_2d

  subroutine stencil_dense_apply_3d(ss, dd, ng_d, uu, ng_u)
    integer, intent(in) :: ng_d, ng_u
    real (kind = dp_t), intent(in   ) :: ss(:,:,:,0:)
    real (kind = dp_t), intent(in   ) :: uu(1-ng_u:,1-ng_u:,1-ng_u:)
    real (kind = dp_t), intent(  out) :: dd(1-ng_d:,1-ng_d:,1-ng_d:)
    integer i, j, k, nx, ny, nz

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)

    !$OMP PARALLEL DO PRIVATE(i,j,k) IF(nz.ge.4)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             dd(i,j,k) = &
                  + ss(i,j,k, 1)*uu(i-1,j-1,k-1) &
                  + ss(i,j,k, 2)*uu(i  ,j-1,k-1) &
                  + ss(i,j,k, 3)*uu(i+1,j-1,k-1) &
                  + ss(i,j,k, 4)*uu(i-1,j  ,k-1) &
                  + ss(i,j,k, 5)*uu(i  ,j  ,k-1) &
                  + ss(i,j,k, 6)*uu(i+1,j  ,k-1) &
                  + ss(i,j,k, 7)*uu(i-1,j+1,k-1) &
                  + ss(i,j,k, 8)*uu(i  ,j+1,k-1) &
                  + ss(i,j,k, 9)*uu(i+1,j+1,k-1) &

                  + ss(i,j,k,10)*uu(i-1,j-1,k  ) &
                  + ss(i,j,k,11)*uu(i  ,j-1,k  ) &
                  + ss(i,j,k,12)*uu(i+1,j-1,k  ) &
                  + ss(i,j,k,13)*uu(i-1,j  ,k  ) &
                  + ss(i,j,k, 0)*uu(i  ,j  ,k  ) &
                  + ss(i,j,k,14)*uu(i+1,j  ,k  ) &
                  + ss(i,j,k,15)*uu(i-1,j+1,k  ) &
                  + ss(i,j,k,16)*uu(i  ,j+1,k  ) &
                  + ss(i,j,k,17)*uu(i+1,j+1,k  ) &

                  + ss(i,j,k,18)*uu(i-1,j-1,k+1) &
                  + ss(i,j,k,19)*uu(i  ,j-1,k+1) &
                  + ss(i,j,k,20)*uu(i+1,j-1,k+1) &
                  + ss(i,j,k,21)*uu(i-1,j  ,k+1) &
                  + ss(i,j,k,22)*uu(i  ,j  ,k+1) &
                  + ss(i,j,k,23)*uu(i+1,j  ,k+1) &
                  + ss(i,j,k,24)*uu(i-1,j+1,k+1) &
                  + ss(i,j,k,25)*uu(i  ,j+1,k+1) &
                  + ss(i,j,k,26)*uu(i+1,j+1,k+1)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine stencil_dense_apply_3d

  subroutine stencil_fine_flux_1d(ss, flux, uu, mm, ng, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in)  :: ss(:,0:)
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

    nx = size(ss,dim=1)

    if ( dim == 1 ) then
       if ( face == -1 ) then
!         Lo i face
          i = 1
          if (bc_dirichlet(mm(1),1,-1)) then
             flux(1) = ss(i,1)*(uu(i+1)-uu(i)) + ss(i,2)*(uu(i-1)-uu(i)) &
                  - ss(i+1,2)*(uu(i+1)-uu(i))
             if (bc_skewed(mm(i),1,+1)) then
                flux(1) =  flux(1) + ss(i,XBC)*uu(i+2)
             end if
          else 
             flux(1) = ss(i,2)*(uu(i-1)-uu(i))
          end if
       else if ( face == 1 ) then

!         Hi i face
          i = nx
          if (bc_dirichlet(mm(i),1,+1)) then
             flux(1) = ss(i,1)*(uu(i+1)-uu(i)) + ss(i,2)*(uu(i-1)-uu(i)) &
                  - ss(i-1,1)*(uu(i-1)-uu(i))
             if (bc_skewed(mm(i),1,-1)) then
                flux(1) =  flux(1) + ss(i,XBC)*uu(i-2)
             end if
          else 
             flux(1) = ss(i,1)*(uu(i+1)-uu(i))
          end if
       end if
    end if

  end subroutine stencil_fine_flux_1d

  subroutine ml_fill_all_fluxes(ss, flux, uu, mm)

    use bl_prof_module
    use multifab_module

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

    ngu = nghost(uu)

    lcross = ((ncomp(ss) == 5) .or. (ncomp(ss) == 7))

    if ( ncomp(uu) /= ncomp(flux(1)) ) then
       call bl_error("ML_FILL_ALL_FLUXES: uu%nc /= flux%nc")
    end if

    call multifab_fill_boundary(uu, cross = lcross)

    do dim = 1, get_dim(uu)
       do i = 1, nboxes(flux(dim))
          if ( remote(flux(dim), i) ) cycle
          ngf = nghost(flux(dim))
          fp => dataptr(flux(dim), i)
          up => dataptr(uu, i)
          sp => dataptr(ss, i)
          mp => dataptr(mm, i)
          select case(get_dim(ss))
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

    nx = size(ss,dim=1)

    do i = 1,nx-1
      flux(i) = ss(i,2) * (uu(i)-uu(i-1)) 
    end do

    ! Lo i face
     i = 0
     if (bc_dirichlet(mm(i),1,-1)) then
        flux(0) = &
               ss(i,1)*(uu(i+1)-uu(i)) + ss(i  ,2)*(uu(i-1)-uu(i)) &
                                             - ss(i+1,2)*(uu(i+1)-uu(i))
        if (bc_skewed(mm(i),1,+1)) &
             flux(0) = flux(0) + ss(i,XBC)*(uu(i+2)-uu(i)) 
        flux(0) = -flux(0)
     else if (bc_neumann(mm(i),1,-1)) then
        flux(0) = -ss(i,2)*uu(i-1)
        else   
        flux(0) = ss(i,2)*(uu(i)-uu(i-1))
     end if

    ! Hi i face
     i = nx-1
     if (bc_dirichlet(mm(i),1,+1)) then
        flux(nx) = &
               ss(i  ,1)*(uu(i+1)-uu(i)) + ss(i,2)*(uu(i-1)-uu(i)) &
             - ss(i-1,1)*(uu(i-1)-uu(i))
        if (bc_skewed(mm(i),1,-1)) &
             flux(nx) = flux(nx) + ss(i,XBC)*(uu(i-2)-uu(i))
     else if (bc_neumann(mm(i),1,+1)) then
        flux(nx) = ss(i,1)*uu(i+1)
     else 
        flux(nx) = ss(i,1)*(uu(i+1)-uu(i))
     end if

  end subroutine stencil_all_flux_1d

  subroutine stencil_fine_flux_2d(ss, flux, uu, mm, ng, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: uu(1-ng:,1-ng:)
    real (kind = dp_t), intent(out) :: flux(:,:)
    real (kind = dp_t), intent(in ) :: ss(:,:,0:)
    integer           , intent(in)  :: mm(:,:)
    logical, intent(in), optional :: skwd
    integer, intent(in) :: face, dim
    integer nx,ny
    integer i,j
    integer, parameter :: XBC = 5, YBC = 6
    logical :: lskwd

    lskwd = .true. ; if ( present(skwd) ) lskwd = skwd

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

!   Lo i face
    if ( dim == 1 ) then
       if (face == -1) then

          i = 1
          flux(1,:) = ZERO
          do j = 1,ny
             if (bc_dirichlet(mm(i,j),1,-1)) then
                flux(1,j) = &
                       ss(i,j,1)*(uu(i+1,j)-uu(i,j)) &
                     + ss(i,j,2)*(uu(i-1,j)-uu(i,j)) - ss(i+1,j,2)*(uu(i+1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,+1)) &
                     flux(1,j) = flux(1,j) + ss(i,j,XBC)*(uu(i+2,j)-uu(i,j)) 
             else if (bc_neumann(mm(i,j),1,-1)) then
                flux(1,j) = ss(i,j,2)*uu(i-1,j)
             else   
                flux(1,j) = ss(i,j,2)*(uu(i-1,j)-uu(i,j))
             end if
          end do

!      Hi i face
       else if (face == 1) then

          i = nx
          flux(1,:) = ZERO
          do j = 1,ny
             if (bc_dirichlet(mm(i,j),1,+1)) then
                flux(1,j) = &
                       ss(i,j,1)*(uu(i+1,j)-uu(i,j)) &
                     + ss(i,j,2)*(uu(i-1,j)-uu(i,j)) - ss(i-1,j,1)*(uu(i-1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,-1)) &
                     flux(1,j) = flux(1,j) + ss(i,j,XBC)*(uu(i-2,j)-uu(i,j))
             else if (bc_neumann(mm(i,j),1,+1)) then
                flux(1,j) = ss(i,j,1)*uu(i+1,j)
             else 
                flux(1,j) = ss(i,j,1)*(uu(i+1,j)-uu(i,j))
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
                       ss(i,j,3)*(uu(i,j+1)-uu(i,j)) &
                     + ss(i,j,4)*(uu(i,j-1)-uu(i,j)) - ss(i,j+1,4)*(uu(i,j+1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,+1)) &
                     flux(i,1) =  flux(i,1) + ss(i,j,YBC)*(uu(i,j+2)-uu(i,j))
             else if (bc_neumann(mm(i,j),2,-1)) then
                flux(i,1) = ss(i,j,4)*uu(i,j-1)
             else 
                flux(i,1) = ss(i,j,4)*(uu(i,j-1)-uu(i,j))
             end if
          end do


!      Hi j face
       else if (face == 1) then

          j = ny
          flux(:,1) = ZERO
          do i = 1,nx
             if (bc_dirichlet(mm(i,j),2,+1)) then
                flux(i,1) = &
                       ss(i,j,3)*(uu(i,j+1)-uu(i,j)) &
                     + ss(i,j,4)*(uu(i,j-1)-uu(i,j)) - ss(i,j-1,3)*(uu(i,j-1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,-1)) &
                     flux(i,1) = flux(i,1) + ss(i,j,YBC)*(uu(i,j-2)-uu(i,j))
             else if (bc_neumann(mm(i,j),2,+1)) then
                flux(i,1) = ss(i,j,3)*uu(i,j+1)
             else
                flux(i,1) = ss(i,j,3)*(uu(i,j+1)-uu(i,j))
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

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)

    if ( dim == 1 ) then
       do j = 0,ny-1
       do i = 1,nx-1
         flux(i,j) = ss(i,j,2) * (uu(i,j)-uu(i-1,j)) 
       end do
       end do

       ! Lo i face
        i = 0
        do j = 0,ny-1
             if (bc_dirichlet(mm(i,j),1,-1)) then
                flux(0,j) = &
                       ss(i,j,1)*(uu(i+1,j)-uu(i,j)) + ss(i  ,j,2)*(uu(i-1,j)-uu(i,j)) &
                                                     - ss(i+1,j,2)*(uu(i+1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,+1)) &
                     flux(0,j) = flux(0,j) + ss(i,j,XBC)*(uu(i+2,j)-uu(i,j)) 
                flux(0,j) = -flux(0,j)
             else if (bc_neumann(mm(i,j),1,-1)) then
                flux(0,j) = -ss(i,j,2)*uu(i-1,j)
             else   
                flux(0,j) = ss(i,j,2)*(uu(i,j)-uu(i-1,j))
             end if
        end do

       ! Hi i face
        i = nx-1
        do j = 0,ny-1
             if (bc_dirichlet(mm(i,j),1,+1)) then
                flux(nx,j) = &
                       ss(i  ,j,1)*(uu(i+1,j)-uu(i,j)) + ss(i,j,2)*(uu(i-1,j)-uu(i,j)) &
                     - ss(i-1,j,1)*(uu(i-1,j)-uu(i,j))
                if (bc_skewed(mm(i,j),1,-1)) &
                     flux(nx,j) = flux(nx,j) + ss(i,j,XBC)*(uu(i-2,j)-uu(i,j))
             else if (bc_neumann(mm(i,j),1,+1)) then
                flux(nx,j) = ss(i,j,1)*uu(i+1,j)
             else 
                flux(nx,j) = ss(i,j,1)*(uu(i+1,j)-uu(i,j))
             end if
        end do

    else if ( dim == 2 ) then
       do j = 1,ny-1
       do i = 0,nx-1
         flux(i,j) = ss(i,j,4) * (uu(i,j)-uu(i,j-1)) 
       end do
       end do

       ! Lo j face
       j = 0
       do i = 0,nx-1
             if (bc_dirichlet(mm(i,j),2,-1)) then
                flux(i,0) = &
                       ss(i,j,3)*(uu(i,j+1)-uu(i,j)) + ss(i,j  ,4)*(uu(i,j-1)-uu(i,j)) & 
                                                     - ss(i,j+1,4)*(uu(i,j+1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,+1)) &
                     flux(i,0) =  flux(i,0) + ss(i,j,YBC)*(uu(i,j+2)-uu(i,j))
                flux(i,0) = -flux(i,0)
             else if (bc_neumann(mm(i,j),2,-1)) then
                flux(i,0) = -ss(i,j,4)*uu(i,j-1)
             else 
                flux(i,0) = ss(i,j,4)*(uu(i,j)-uu(i,j-1))
             end if
       end do

       ! Hi j face
       j = ny-1
       do i = 0,nx-1
             if (bc_dirichlet(mm(i,j),2,+1)) then
                flux(i,ny) = &
                       ss(i,j  ,3)*(uu(i,j+1)-uu(i,j)) + ss(i,j,4)*(uu(i,j-1)-uu(i,j)) & 
                     - ss(i,j-1,3)*(uu(i,j-1)-uu(i,j))
                if (bc_skewed(mm(i,j),2,-1)) &
                     flux(i,ny) = flux(i,ny) + ss(i,j,YBC)*(uu(i,j-2)-uu(i,j))
             else if (bc_neumann(mm(i,j),2,+1)) then
                flux(i,ny) = ss(i,j,3)*uu(i,j+1)
             else
                flux(i,ny) = ss(i,j,3)*(uu(i,j+1)-uu(i,j))
             end if
       end do

    end if

  end subroutine stencil_all_flux_2d

  subroutine stencil_fine_flux_3d(ss, flux, uu, mm, ng, face, dim, skwd)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in ) :: ss(:,:,:,0:)
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

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)

    if ( dim ==  1 ) then
       !   Lo i face
       if (face == -1) then

          i = 1
          flux(1,:,:) = ZERO

          do k = 1,nz
             do j = 1,ny
                if (bc_dirichlet(mm(i,j,k),1,-1)) then
                   flux(1,j,k) =  &
                          ss(i,j,k,1)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(i,j,k,2)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(i+1,j,k,2)*(uu(i+1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,+1)) &
                        flux(1,j,k) =  flux(1,j,k) + ss(i,j,k,XBC)*(uu(i+2,j,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),1,-1)) then
                   flux(1,j,k) = ss(i,j,k,2)*uu(i-1,j,k)
                else 
                   flux(1,j,k) = ss(i,j,k,2)*(uu(i-1,j,k)-uu(i,j,k))
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
                          ss(i,j,k,1)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(i,j,k,2)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(i-1,j,k,1)*(uu(i-1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,-1)) &
                        flux(1,j,k) =  flux(1,j,k) + ss(i,j,k,XBC)*(uu(i-2,j,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),1,+1)) then
                   flux(1,j,k) = ss(i,j,k,1)*uu(i+1,j,k)
                else 
                   flux(1,j,k) = ss(i,j,k,1)*(uu(i+1,j,k)-uu(i,j,k))
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
                          ss(i,j,k,3)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(i,j,k,4)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(i,j+1,k,4)*(uu(i,j+1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,+1)) &
                        flux(i,1,k) =  flux(i,1,k) + ss(i,j,k,YBC)*(uu(i,j+2,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),2,-1)) then
                   flux(i,1,k) = ss(i,j,k,4)*uu(i,j-1,k)
                else 
                   flux(i,1,k) = ss(i,j,k,4)*(uu(i,j-1,k)-uu(i,j,k))
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
                          ss(i,j,k,3)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(i,j,k,4)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(i,j-1,k,3)*(uu(i,j-1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,-1)) &
                        flux(i,1,k) =  flux(i,1,k) + ss(i,j,k,YBC)*(uu(i,j-2,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),2,+1)) then
                   flux(i,1,k) = ss(i,j,k,3)*uu(i,j+1,k)
                else
                   flux(i,1,k) = ss(i,j,k,3)*(uu(i,j+1,k)-uu(i,j,k))
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
                          ss(i,j,k,5)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(i,j,k,6)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(i,j,k+1,6)*(uu(i,j,k+1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,+1)) &
                        flux(i,j,1) =  flux(i,j,1) + ss(i,j,k,ZBC)*(uu(i,j,k+2)-uu(i,j,k)) 
                else if (bc_neumann(mm(i,j,k),3,-1)) then
                   flux(i,j,1) = ss(i,j,k,6)*uu(i,j,k-1)
                else 
                   flux(i,j,1) = ss(i,j,k,6)*(uu(i,j,k-1)-uu(i,j,k))
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
                          ss(i,j,k,5)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(i,j,k,6)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(i,j,k-1,5)*(uu(i,j,k-1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,-1)) &
                        flux(i,j,1) =  flux(i,j,1) + ss(i,j,k,ZBC)*(uu(i,j,k-2)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),3,+1)) then
                   flux(i,j,1) = ss(i,j,k,5)*uu(i,j,k+1)
                else
                   flux(i,j,1) = ss(i,j,k,5)*(uu(i,j,k+1)-uu(i,j,k))
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

    nx = size(ss,dim=1)
    ny = size(ss,dim=2)
    nz = size(ss,dim=3)

    if ( dim ==  1 ) then

       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 0,nz-1
          do j = 0,ny-1
             do i = 0,nx-1
                flux(i,j,k) = ss(i,j,k,2) * (uu(i,j,k)-uu(i-1,j,k))
             end do
          end do
       end do
       !$OMP END PARALLEL DO

       !   Lo i face
       i = 0
       do k = 0,nz-1
             do j = 0,ny-1
                if (bc_dirichlet(mm(i,j,k),1,-1)) then
                   flux(0,j,k) =  &
                          ss(i,j,k,1)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(i,j,k,2)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(i+1,j,k,2)*(uu(i+1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,+1)) &
                        flux(0,j,k) =  flux(0,j,k) + ss(i,j,k,XBC)*(uu(i+2,j,k)-uu(i,j,k))
                   flux(0,j,k) = -flux(0,j,k)
                else if (bc_neumann(mm(i,j,k),1,-1)) then
                   flux(0,j,k) = -ss(i,j,k,2)*uu(i-1,j,k)
                else 
                   flux(0,j,k) = ss(i,j,k,2)*(uu(i,j,k)-uu(i-1,j,k))
                end if
             end do
       end do

       !   Hi i face
       i = nx-1
       do k = 0,nz-1
             do j = 0,ny-1
                if (bc_dirichlet(mm(i,j,k),1,+1)) then
                   flux(nx,j,k) = &
                          ss(i,j,k,1)*(uu(i+1,j,k)-uu(i,j,k)) &
                        + ss(i,j,k,2)*(uu(i-1,j,k)-uu(i,j,k)) &
                        - ss(i-1,j,k,1)*(uu(i-1,j,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),1,-1)) &
                        flux(nx,j,k) =  flux(nx,j,k) + ss(i,j,k,XBC)*(uu(i-2,j,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),1,+1)) then
                   flux(nx,j,k) = ss(i,j,k,1)*uu(i+1,j,k)
                else 
                   flux(nx,j,k) = ss(i,j,k,1)*(uu(i+1,j,k)-uu(i,j,k))
                end if
             end do
       end do

    else if ( dim == 2 ) then

       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 0,nz-1
          do j = 0,ny-1
             do i = 0,nx-1
                flux(i,j,k) = ss(i,j,k,4) * (uu(i,j,k)-uu(i,j-1,k))
             end do
          end do
       end do
       !$OMP END PARALLEL DO

       !   Lo j face
       j = 0
       do k = 0,nz-1
             do i = 0,nx-1
                if (bc_dirichlet(mm(i,j,k),2,-1)) then
                   flux(i,0,k) = &
                          ss(i,j,k,3)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(i,j,k,4)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(i,j+1,k,4)*(uu(i,j+1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,+1)) &
                        flux(i,0,k) =  flux(i,0,k) + ss(i,j,k,YBC)*(uu(i,j+2,k)-uu(i,j,k))
                   flux(i,0,k) = -flux(i,0,k)
                else if (bc_neumann(mm(i,j,k),2,-1)) then
                   flux(i,0,k) = -ss(i,j,k,4)*uu(i,j-1,k)
                else 
                   flux(i,0,k) = ss(i,j,k,4)*(uu(i,j,k)-uu(i,j-1,k))
                end if
             end do
       end do

       !   Hi j face
       j = ny-1
       do k = 0,nz-1
             do i = 0,nx-1
                if (bc_dirichlet(mm(i,j,k),2,+1)) then
                   flux(i,ny,k) =  &
                          ss(i,j,k,3)*(uu(i,j+1,k)-uu(i,j,k)) &
                        + ss(i,j,k,4)*(uu(i,j-1,k)-uu(i,j,k)) &
                        - ss(i,j-1,k,3)*(uu(i,j-1,k)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),2,-1)) &
                        flux(i,ny,k) =  flux(i,1,ny) + ss(i,j,k,YBC)*(uu(i,j-2,k)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),2,+1)) then
                   flux(i,ny,k) = ss(i,j,k,3)*uu(i,j+1,k)
                else
                   flux(i,ny,k) = ss(i,j,k,3)*(uu(i,j+1,k)-uu(i,j,k))
                end if
             end do
       end do

    else if ( dim == 3 ) then

       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k = 0,nz-1
          do j = 0,ny-1
             do i = 0,nx-1
                flux(i,j,k) = ss(i,j,k,6) * (uu(i,j,k)-uu(i,j,k-1))
             end do
          end do
       end do
       !$OMP END PARALLEL DO

       !   Lo k face
       k = 0
       do j = 0,ny-1
             do i = 0,nx-1
                if (bc_dirichlet(mm(i,j,k),3,-1)) then
                   flux(i,j,0) =  &
                          ss(i,j,k,5)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(i,j,k,6)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(i,j,k+1,6)*(uu(i,j,k+1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,+1)) &
                        flux(i,j,0) =  flux(i,j,0) + ss(i,j,k,ZBC)*(uu(i,j,k+2)-uu(i,j,k)) 
                   flux(i,j,0) = -flux(i,j,0)
                else if (bc_neumann(mm(i,j,k),3,-1)) then
                   flux(i,j,0) = -ss(i,j,k,6)*uu(i,j,k-1)
                else 
                   flux(i,j,0) = ss(i,j,k,6)*(uu(i,j,k)-uu(i,j,k-1))
                end if
             end do
       end do

       !   Hi k face
       k = nz-1
       do j = 0,ny-1
             do i = 0,nx-1
                if (bc_dirichlet(mm(i,j,k),3,+1)) then
                   flux(i,j,nz) =  &
                          ss(i,j,k,5)*(uu(i,j,k+1)-uu(i,j,k)) &
                        + ss(i,j,k,6)*(uu(i,j,k-1)-uu(i,j,k)) &
                        - ss(i,j,k-1,5)*(uu(i,j,k-1)-uu(i,j,k))
                   if (bc_skewed(mm(i,j,k),3,-1)) &
                        flux(i,j,nz) =  flux(i,j,nz) + ss(i,j,k,ZBC)*(uu(i,j,k-2)-uu(i,j,k))
                else if (bc_neumann(mm(i,j,k),3,+1)) then
                   flux(i,j,nz) = ss(i,j,k,5)*uu(i,j,k+1)
                else
                   flux(i,j,nz) = ss(i,j,k,5)*(uu(i,j,k+1)-uu(i,j,k))
                end if
             end do
       end do

    end if

  end subroutine stencil_all_flux_3d

end module cc_stencil_apply_module
