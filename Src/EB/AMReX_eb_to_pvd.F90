module amrex_eb_to_vtk

  use iso_c_binding    , only: c_int
  use amrex_fort_module, only: amrex_real

  implicit none

  integer :: np = 0
  integer :: nc = 0

  real(amrex_real), allocatable :: points(:,:)
  integer,          allocatable :: connectivity(:,:)

  ! parameter to exclude small cells
  real(amrex_real), parameter :: stol = 1.0d-12

contains

  subroutine amrex_eb_to_polygon (problo, dx, lo, hi, flag, fglo, fghi, bcent, blo, bhi, &
       apx, axlo, axhi,  apy, aylo, ayhi, apz, azlo, azhi) &
       bind(C, name="amrex_eb_to_polygon")

  use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell, is_single_valued_cell

  implicit none

  integer, intent(in   ) :: lo(3),   hi(3)
  integer, intent(in   ) :: axlo(3), axhi(3)
  integer, intent(in   ) :: aylo(3), ayhi(3)
  integer, intent(in   ) :: azlo(3), azhi(3)
  integer, intent(in   ) :: fglo(3), fghi(3)
  integer, intent(in   ) :: blo(3),  bhi(3)

  integer,          intent(in   ) :: &
       flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

  real(amrex_real), intent(in   ) :: problo(3), dx(3), &
       bcent(blo(1):bhi(1),blo(2):bhi(2),blo(3):bhi(3),3),   &
       apx(axlo(1):axhi(1),axlo(2):axhi(2),axlo(3):axhi(3)), &
       apy(aylo(1):ayhi(1),aylo(2):ayhi(2),aylo(3):ayhi(3)), &
       apz(azlo(1):azhi(1),azlo(2):azhi(2),azlo(3):azhi(3))

  real(amrex_real) :: axm, axp, aym, ayp, azm, azp
  real(amrex_real) :: apnorm, apnorminv

  real(amrex_real) :: centroid(3), normal(3)
  real(amrex_real) :: distance
  real(amrex_real) :: n0(3), p, n0_d(3), p_d, tol
  real(amrex_real) :: vertex(8,3), alpha(12), alpha_d(12), apoints(12,3)
  logical          :: alpha_intersect(12), alpha_d_intersect(12)

  integer :: i, j, k, lc1, lc2, count, count_d

  real(amrex_real), parameter :: ihat(3) = (/1.0d0, 0.0d0, 0.0d0/)
  real(amrex_real), parameter :: jhat(3) = (/0.0d0, 1.0d0, 0.0d0/)
  real(amrex_real), parameter :: khat(3) = (/0.0d0, 0.0d0, 1.0d0/)

  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

           ! NOTE: do not skip fully enclosed cells (is_covered_cell), as this seems
           ! to skip thin walls in the domain:
           !if(.not.is_regular_cell(flag(i,j,k)) .and. &
           !   .not.is_covered_cell(flag(i,j,k))) then

           ! Instead only look for EBs
           ! if( .not.is_regular_cell(flag(i,j,k))) then

           ! If covered cells are accounted for in this loop, a FPE arises
           ! since apnorm is zero.
           if (is_single_valued_cell(flag(i,j,k))) then

              ! Calculate unit normal
              axm = apx(i,  j  , k  )
              axp = apx(i+1,j  , k  )
              aym = apy(i,  j  , k  )
              ayp = apy(i,  j+1, k  )
              azm = apz(i,  j  , k  )
              azp = apz(i,  j  , k+1)

              apnorm = sqrt((axm-axp)**2 + (aym-ayp)**2 + (azm-azp)**2)


              apnorminv = 1.d0 / apnorm
              normal(1) = (axp-axm) * apnorminv
              normal(2) = (ayp-aym) * apnorminv
              normal(3) = (azp-azm) * apnorminv

              ! convert bcent to global coordinate system centered at plo
              centroid(1) = problo(1) + bcent(i,j,k,1)*dx(1) + (dble(i) + 0.5d0)*dx(1)
              centroid(2) = problo(2) + bcent(i,j,k,2)*dx(2) + (dble(j) + 0.5d0)*dx(2)
              centroid(3) = problo(3) + bcent(i,j,k,3)*dx(3) + (dble(k) + 0.5d0)*dx(3)

              ! vertices of bounding cell (i,j,k)
              vertex(1,:) = (/problo(1)+dble(i  )*dx(1), problo(2)+dble(j  )*dx(2), problo(3)+dble(k  )*dx(3)/)
              vertex(2,:) = (/problo(1)+dble(i+1)*dx(1), problo(2)+dble(j  )*dx(2), problo(3)+dble(k  )*dx(3)/)
              vertex(3,:) = (/problo(1)+dble(i  )*dx(1), problo(2)+dble(j+1)*dx(2), problo(3)+dble(k  )*dx(3)/)
              vertex(4,:) = (/problo(1)+dble(i+1)*dx(1), problo(2)+dble(j+1)*dx(2), problo(3)+dble(k  )*dx(3)/)
              vertex(5,:) = (/problo(1)+dble(i  )*dx(1), problo(2)+dble(j  )*dx(2), problo(3)+dble(k+1)*dx(3)/)
              vertex(6,:) = (/problo(1)+dble(i+1)*dx(1), problo(2)+dble(j  )*dx(2), problo(3)+dble(k+1)*dx(3)/)
              vertex(7,:) = (/problo(1)+dble(i  )*dx(1), problo(2)+dble(j+1)*dx(2), problo(3)+dble(k+1)*dx(3)/)
              vertex(8,:) = (/problo(1)+dble(i+1)*dx(1), problo(2)+dble(j+1)*dx(2), problo(3)+dble(k+1)*dx(3)/)

              ! NOTE: this seems to be unncessary:
              ! skip cells that have a tiny intersection and cells that have
              ! the centroid on a face/edge/corner
              !if(apnorm > stol    .and. &
              !   vertex(1,1) < centroid(1) .and. centroid(1) < vertex(8,1) .and. &
              !   vertex(1,2) < centroid(2) .and. centroid(2) < vertex(8,2) .and. &
              !   vertex(1,3) < centroid(3) .and. centroid(3) < vertex(8,3)) then

                 ! Compute EB facets for current cell.
                 call calc_hesse(distance, n0, p, normal, centroid)
                 call calc_alpha(alpha, distance, n0, p, vertex, dx)
                 call calc_intercects(count, alpha_intersect, alpha)

                 ! If the number of facet "contained" in does not describe a facet:
                 ! ... I.e. there's less than 3 (not even a triangle) or more than 6
                 ! ... (I have no idea what that is):
                 !   => Move the centroid a little back and forth along the normal
                 !      to see if that makes a difference:
                 if((count < 3) .or. (count > 6)) then
                    tol = min(dx(1), dx(2), dx(3)) / 100  ! bit of a fudge factor

                    call calc_hesse(distance, n0_d, p_d, normal, centroid + tol*normal)
                    call calc_alpha(alpha_d, distance, n0_d, p_d, vertex, dx)
                    call calc_intercects(count_d, alpha_d_intersect, alpha_d)
                    if(count_d >= 3 .and. count_d <= 6) then
                        count = count_d
                        alpha_intersect = alpha_d_intersect
                    endif

                    call calc_hesse(distance, n0_d, p_d, normal, centroid - tol*normal)
                    call calc_alpha(alpha_d, distance, n0_d, p_d, vertex, dx)
                    call calc_intercects(count_d, alpha_d_intersect, alpha_d)
                    if((count_d >= 3) .and. (count_d <= 6)) then
                        count = count_d
                        alpha_intersect = alpha_d_intersect
                    endif
                 endif

                 ! I know this was a bit of a hack, but it's the only way I prevent
                 ! missing facets...

                 if((count >= 3) .and. (count <= 6)) then

                    call grow_connectivity(nc, connectivity)
                    connectivity(nc,0) = 0

                    ! calculate intersection points.
                    apoints( 1,:) = vertex(1,:) + ihat*dx(1)*alpha( 1)
                    apoints( 2,:) = vertex(2,:) + jhat*dx(2)*alpha( 2)
                    apoints( 3,:) = vertex(3,:) + ihat*dx(1)*alpha( 3)
                    apoints( 4,:) = vertex(1,:) + jhat*dx(2)*alpha( 4)
                    apoints( 5,:) = vertex(1,:) + khat*dx(3)*alpha( 5)
                    apoints( 6,:) = vertex(2,:) + khat*dx(3)*alpha( 6)
                    apoints( 7,:) = vertex(4,:) + khat*dx(3)*alpha( 7)
                    apoints( 8,:) = vertex(3,:) + khat*dx(3)*alpha( 8)
                    apoints( 9,:) = vertex(5,:) + ihat*dx(1)*alpha( 9)
                    apoints(10,:) = vertex(6,:) + jhat*dx(2)*alpha(10)
                    apoints(11,:) = vertex(7,:) + ihat*dx(1)*alpha(11)
                    apoints(12,:) = vertex(5,:) + jhat*dx(2)*alpha(12)

                    ! store intersections with grid cell alpha in [0,1]
                    do lc1=1,12
                       if(alpha_intersect(lc1)) then
                          call grow_points(np, points)
                          points(np,:) = apoints(lc1,:)
                          lc2 = connectivity(nc,0) + 1
                          connectivity(nc,0) = lc2
                          connectivity(nc,lc2) = np
                       endif
                    enddo

                    ! reorder points for a convex polygon
                    call reorder_polygon(points, connectivity(nc,:), n0, centroid)

                 endif

              !endif

           end if ! if .not. regular
        end do
     end do
  end do

end subroutine amrex_eb_to_polygon

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
subroutine amrex_write_eb_vtp(myID) bind(C, name="amrex_write_eb_vtp")

  implicit none

  integer(c_int), intent(in   ) :: myID

  character(len=8) :: cID
  integer :: lc1, lc2

  write(cID,"(I8.8)") myID

  open(unit=100, file='eb_'//cID//'.vtp', status='unknown')

! Write the necessary header information for a PolyData file type
  write(100,"(A)")'<?xml version="1.0"?>'
  write(100,"(2A)") '<VTKFile type="PolyData"',&
       ' version="0.1" byte_order="LittleEndian">'
  write(100,"(3x,A)") '<PolyData>'

! Write Piece tag and identify the number of particles in the system.
  write(100,"(6x,a,i10.10,a,a,i10.10,a)") &
       '<Piece NumberOfPoints="',np, '" NumberOfVerts="0" ', &
       'NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="',nc,'">'

  write(100,"(9x,a)") '<Points>'
  write(100,"(12x,a,a)") '<DataArray type="Float32" ',&
       'NumberOfComponents="3" format="ascii">'
  do lc1 = 1,np
     write (100,"(15x,3(es13.6,3x))") real(points(lc1,1:3))
  enddo
  write(100,"(12x,a,/9x,a)")'</DataArray>','</Points>'

  write(100,"(9x,a)") '<Polys>'
  write(100,"(12x,a,a)") '<DataArray type="Int32" ',&
       'Name="connectivity" format="ascii">'
  do lc1 = 1,nc
     do lc2 = 1, connectivity(lc1,0)
        write(100,"(1x,i6)",advance='no') connectivity(lc1,lc2)-1
     enddo
     write(100,"(' ')")
  enddo
  write(100,"(12x,a)")'</DataArray>'

  write(100,"(12x,a,a)") '<DataArray type="Int32" ',&
       'Name="offsets" format="ascii">'
  write(100,"(15x)",advance='no')
  lc2 = 0
  do lc1 = 1,nc
     lc2 = lc2 + connectivity(lc1,0)
     write(100,"(1x,i8)",advance='no') lc2
  enddo
  write(100,"(' ')")
  write(100,"(12x,a)")'</DataArray>'
  write(100,"( 9x,a)")'</Polys>'

  ! Write tags for data not included (vtp format style)
  write(100,"(9x,a)")'<PointData></PointData>'
  write(100,"(9x,a)")'<CellData></CellData>'
  write(100,"(6x,a,/3x,a,/a)")&
       '</Piece>','</PolyData>','</VTKFile>'

  close(100)

  if(allocated(points      )) deallocate(points      )
  if(allocated(connectivity)) deallocate(connectivity)

end subroutine amrex_write_eb_vtp



!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
subroutine amrex_write_pvtp(nProcs) bind(C, name="amrex_write_pvtp")

  implicit none

  integer(c_int), intent(in   ) :: nProcs

  integer :: lc1
  character(len=8) :: clc1

  open(unit=100, file='eb.pvtp', status='unknown')

! Write the necessary header information for a PolyData file type
  write(100,"(A)")'<?xml version="1.0"?>'
  write(100,"(2a)") '<VTKFile type="PPolyData"',&
       ' version="0.1" byte_order="LittleEndian">'

  write(100,"(2x,a)") '<PPolyData GhostLevel="0">'
  write(100,"(4x,a)") '<PPointData/>'
  write(100,"(4x,a)") '<PCellData/>'
  write(100,"(4x,a)") '<PPoints>'
  write(100,"(6x,a,a)") '<PDataArray type="Float32" ',&
       'NumberOfComponents="3"/>'
  write(100,"(4x,a)")'</PPoints>'

  do lc1=0,nProcs-1
     write(clc1,"(I8.8)") lc1
     write(100,"(4x,a)") '<Piece Source="eb_'//clc1//'.vtp"/>'
  enddo
  write(100,"(2x,a)") '</PPolyData>'
  write(100,"(a)") '</VTKFile>'

  close(100)

  if(allocated(points      )) deallocate(points      )
  if(allocated(connectivity)) deallocate(connectivity)

end subroutine amrex_write_pvtp

!.......................................................................!
! SUBROUTINE CALC_HESSE(distance, n0, p, normal, centroid)              !
!   Computes the Hesse Normal Form corresponding to normal and centroid !
!.......................................................................!
  subroutine calc_hesse(distance, n0, p, normal, centroid)
    real(amrex_real), intent(  out) :: distance, n0(3), p
    real(amrex_real), intent(in   ) :: normal(3), centroid(3)

    real(amrex_real) :: sign_of_dist

    ! General equation of a plane: Ax + By + Cz + D = 0
    ! here D := distance
    distance = -dot_product(normal, centroid)

    ! Get the sign of the distance
    sign_of_dist = -distance/dabs(distance)

    ! Get Hessian form
    n0 = sign_of_dist*normal/dot_product(normal, normal)
    p  = sign_of_dist*(-distance)

  end subroutine calc_hesse

!.......................................................................!
! SOUBROUTINE CALC_ALPHA(alpha, distance, n0, p, vertex, dx)            !
!   Fills the alpha array                                               !
!.......................................................................!
  subroutine calc_alpha(alpha, distance, n0, p, vertex, dx)
    real(amrex_real), intent(  out) :: alpha(12)
    real(amrex_real), intent(in   ) :: distance, n0(3), p, vertex(8,3), dx(3)


    ! default (large) value
    alpha = 10.0

    ! Ray-xAxis intersection
    if(abs(n0(1)) > epsilon(0.0d0)) then
       alpha( 1) = (p - dot_product(n0,vertex(1,:)))/(n0(1)*dx(1))
       alpha( 3) = (p - dot_product(n0,vertex(3,:)))/(n0(1)*dx(1))
       alpha( 9) = (p - dot_product(n0,vertex(5,:)))/(n0(1)*dx(1))
       alpha(11) = (p - dot_product(n0,vertex(7,:)))/(n0(1)*dx(1))
    endif

    ! Ray-yAxis intersection
    if(abs(n0(2)) > epsilon(0.0d0)) then
       alpha( 2) = (p - dot_product(n0,vertex(2,:)))/(n0(2)*dx(2))
       alpha( 4) = (p - dot_product(n0,vertex(1,:)))/(n0(2)*dx(2))
       alpha(10) = (p - dot_product(n0,vertex(6,:)))/(n0(2)*dx(2))
       alpha(12) = (p - dot_product(n0,vertex(5,:)))/(n0(2)*dx(2))
    endif

    ! Ray-zAxis intersection
    if(abs(n0(3)) > epsilon(0.0d0)) then
       alpha( 5) = (p - dot_product(n0,vertex(1,:)))/(n0(3)*dx(3))
       alpha( 6) = (p - dot_product(n0,vertex(2,:)))/(n0(3)*dx(3))
       alpha( 7) = (p - dot_product(n0,vertex(4,:)))/(n0(3)*dx(3))
       alpha( 8) = (p - dot_product(n0,vertex(3,:)))/(n0(3)*dx(3))
    endif

  end subroutine calc_alpha


!.......................................................................!
! SUBROUTINE CALC_INTERSECTS(int_count, intersect_flags, alpha)         !
!   Fills count and flags selecting the alphas which are in (0,1)       !
!.......................................................................!
  subroutine calc_intercects(int_count, intersect_flags, alpha)
    logical,          intent(  out) :: intersect_flags(12)
    integer,          intent(  out) :: int_count
    real(amrex_real), intent(in   ) :: alpha(12)

    integer :: lc1

    intersect_flags(:) = .false.
    int_count          = 0

    do lc1=1,12
       if(intersects(alpha(lc1))) then
          int_count = int_count + 1
          intersect_flags(lc1) = .true.
       endif
    enddo
  end subroutine


!.......................................................................!
! FUNCTION INTERSECTS(val)                                              !
!   Returns .true. iff val in [0,1]                                     !
!.......................................................................!
  logical function intersects(val)
    real(amrex_real), intent(in) :: val
    intersects = ((val > 0.0d0) .and. (val < 1.0d0))
  end function intersects


!.......................................................................!
!                                                                       !
!                                                                       !
!.......................................................................!
  subroutine grow_points(nPoints, lPoints)

    integer,      intent(inout)              :: nPoints
    real(amrex_real), intent(inout), allocatable :: lPoints(:,:)

    integer :: csize
    real(amrex_real),   allocatable :: rtmp(:,:)

    nPoints = nPoints + 1

    ! Increase real data
    if(.not.(allocated(lPoints))) then
       allocate(lPoints(1024,3))
    else
       csize = size(lPoints,1)
       if(nPoints >= csize) then
          allocate(rtmp(2*csize,3))
          rtmp(1:csize,:) = lPoints(1:csize,:)
          call move_alloc(rtmp,lPoints)
       endif
    endif

    return
  end subroutine grow_points


!.......................................................................!
!                                                                       !
!                                                                       !
!.......................................................................!
  subroutine grow_connectivity(nConnections, lConnectivity)

    integer, intent(inout)              :: nConnections
    integer, intent(inout), allocatable :: lConnectivity(:,:)

    integer :: csize
    integer, allocatable :: itmp(:,:)

    nConnections = nConnections + 1

    ! Increase real data
    if(.not.(allocated(lConnectivity))) then
       allocate(lConnectivity(1024,0:6))
       lConnectivity(:,0) = 0
    else
       csize = size(lConnectivity,1)
       if(nConnections >= csize) then
          allocate(itmp(2*csize,0:6))
          itmp(1:csize,:) = lConnectivity(1:csize,:)
          call move_alloc(itmp,lConnectivity)
          lConnectivity(nConnections:,0) = 0
       endif
    endif

    return
  end subroutine grow_connectivity


!.......................................................................!
!                                                                       !
!                                                                       !
!.......................................................................!
  subroutine reorder_polygon(lpoints, lconnect, lnormal, lcentroid)

    implicit none

    real(amrex_real), intent(in   ) :: lpoints(:,:)
    integer,          intent(inout) :: lconnect(0:6)
    real(amrex_real), intent(in   ) :: lnormal(3)
    real(amrex_real), intent(in   ) :: lcentroid(3)

    integer :: i, k, pi, pk, longest

    real(amrex_real) :: ref_angle, angle
    real(amrex_real) :: center(3)

    longest = 3
    if(abs(lnormal(1)) > abs(lnormal(2))) then
       if(abs(lnormal(1)) > abs(lnormal(3))) longest = 1
    else
       if(abs(lnormal(2)) > abs(lnormal(3))) longest = 2
    endif

    center = 0.0d0
    do i=1, lconnect(0)
       center = center + points(lconnect(i),:)
    enddo
    center = center/dble(lconnect(0))

    if(longest == 1) then
       do i=1, lconnect(0)-1
          pi = lconnect(i)
          ref_angle = atan2(lpoints(pi,3)-center(3), &
               lpoints(pi,2)-center(2))
          do k = i+1, lconnect(0)
             pk = lconnect(k)
             angle = atan2(lpoints(pk,3)-center(3),  &
                           lpoints(pk,2)-center(2))
             if(angle < ref_angle) then
                ref_angle = angle
                lconnect(k) = pi
                lconnect(i) = pk
                pi = pk
             endif
          enddo
       enddo

    else if(longest == 2) then
       do i=1, lconnect(0)-1
          pi = lconnect(i)
          ref_angle = atan2(lpoints(pi,1)-center(1), &
                            lpoints(pi,3)-center(3))
          do k = i+1, lconnect(0)
             pk = lconnect(k)
             angle = atan2(lpoints(pk,1)-center(1), &
                           lpoints(pk,3)-center(3))
             if(angle < ref_angle) then
                ref_angle = angle
                lconnect(k) = pi
                lconnect(i) = pk
                pi = pk
             endif
          enddo
       enddo

    else if(longest == 3) then
       do i=1, lconnect(0)-1
          pi = lconnect(i)
          ref_angle = atan2(lpoints(pi,2)-center(2), &
                            lpoints(pi,1)-center(1))
          do k = i+1, lconnect(0)
             pk = lconnect(k)
             angle = atan2(lpoints(pk,2)-center(2), &
                           lpoints(pk,1)-center(1))
             if(angle < ref_angle) then
                ref_angle = angle
                lconnect(k) = pi
                lconnect(i) = pk
                pi = pk
             endif
          enddo
       enddo
    endif

  end subroutine reorder_polygon

!.......................................................................!
!                                                                       !
!                                                                       !
!.......................................................................!
  subroutine amrex_eb_grid_coverage (myID, problo, dx, lo, hi, flag, fglo, fghi)&
       bind(C, name="amrex_eb_grid_coverage")

  use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell

  implicit none

  integer, intent(in   ) :: myID
  integer, intent(in   ) :: lo(3),   hi(3)
  integer, intent(in   ) :: fglo(3), fghi(3)

  integer,          intent(in   ) :: &
       flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

  real(amrex_real), intent(in   ) :: problo(3), dx(3)

  integer :: i, j, k, lc1

  integer, save :: grid = 0

  character(len=4) :: cgrid, cID

  integer :: nodes(3)

  lc1 = 0
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)

          !if(is_regular_cell(flag(i,j,k)) .or. &
          !     .not.is_covered_cell(flag(i,j,k))) lc1 = lc1 + 1
          if(is_regular_cell(flag(i,j,k))) lc1 = lc1 + 1

        end do
     end do
  end do

  grid = grid + 1

  if(lc1 == 0) return

  write(cID,  "(i4.4)") myID
  write(cgrid,"(i4.4)") grid

  nodes = hi-lo + 1

  open(unit=100, file='eb_grid_'//cID//'_'//cgrid//'.vtr', status='unknown')

  write(100,'(A)') '<?xml version="1.0"?>'
  write(100,'(A)') '<VTKFile type="RectilinearGrid" &
       &version="0.1" byte_order="LittleEndian">'

  write(100,'(A,6I6,A)') '<RectilinearGrid &
       &WholeExtent="',0,nodes(1),0,nodes(2),0,nodes(3),'">'

  write(100,'(A,6I6,A)') '<Piece Extent="',0,&
     nodes(1),0,nodes(2),0,nodes(3),'">'
  write(100,'(A)') '<Coordinates>'

  call data_array(problo(1), lo(1), nodes(1), dx(1))
  call data_array(problo(2), lo(2), nodes(2), dx(2))
  call data_array(problo(3), lo(3), nodes(3), dx(3))

  write(100,'("</Coordinates>")')
  write(100,'("</Piece>")')
  write(100,'("</RectilinearGrid>")')
  write(100,'("</VTKFile>")')

  close(100)

contains

  subroutine data_array(problo, lo, lnodes, ldx)

    implicit none

    integer,          intent(in) :: lo, lnodes
    real(amrex_real), intent(in) :: problo
    real(amrex_real), intent(in) :: ldx

    integer :: llc
    real(amrex_real), allocatable :: lines(:)
    real(amrex_real)              :: grid_start 

    if(allocated(lines)) deallocate(lines)
    allocate(lines(0:lnodes))

    grid_start = problo + real(lo,amrex_real)*ldx

    do llc = 0, lnodes
       lines(llc) = grid_start + real(llc,amrex_real)*ldx
    enddo

    write(100,'(A,F14.8,A,F14.8,A)') '<DataArray &
         &type="Float32" format="ascii" &
         &RangeMin="',lines(0),'" RangeMax="',lines(lnodes),'">'
    write(100,'(10F14.8)') lines
    write(100,'(A)') '</DataArray>'

    deallocate(lines)

  end subroutine data_array

end subroutine amrex_eb_grid_coverage

end module amrex_eb_to_vtk
