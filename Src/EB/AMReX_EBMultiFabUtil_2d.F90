module amrex_eb_util_module

  use amrex_fort_module, only : amrex_real
  use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell, is_single_valued_cell
  use amrex_constants_module, only : half, zero, one
  implicit none

  private
  public :: amrex_eb_avgdown_sv, amrex_eb_avgdown, amrex_eb_avgdown_faces, &
       amrex_eb_avgdown_boundaries, amrex_compute_eb_divergence, &
       amrex_eb_avg_fc_to_cc, amrex_eb_set_covered_nodes, &
       amrex_eb_interpolate_to_face_centroid, &
       amrex_eb_interpolate_to_face_centroid_per_cell

contains

  subroutine amrex_eb_avgdown_sv (lo, hi, fine, flo, fhi, crse, clo, chi, &
       fv, fvlo, fvhi, vfrc, vflo, vfhi, lrat, ncomp) bind(c,name='amrex_eb_avgdown_sv')
    implicit none
    integer, intent(in) :: lo(2), hi(2), flo(2), fhi(2), clo(2), chi(2), &
         fvlo(2), fvhi(2), vflo(2), vfhi(2), lrat(2), ncomp
    real(amrex_real), intent(in   ) :: fine( flo(1): fhi(1), flo(2): fhi(2),ncomp)
    real(amrex_real), intent(inout) :: crse( clo(1): chi(1), clo(2): chi(2),ncomp)
    real(amrex_real), intent(in   ) :: fv  (fvlo(1):fvhi(1),fvlo(2):fvhi(2))
    real(amrex_real), intent(in   ) :: vfrc(vflo(1):vfhi(1),vflo(2):vfhi(2))
    
    integer :: i, j, ii, jj, n, iref, jref
    real(amrex_real) :: cv
    
    do n = 1, ncomp
       do j     = lo(2), hi(2)
          jj    = j * lrat(2)
          do i  = lo(1), hi(1)
             ii = i * lrat(1)
             crse(i,j,n) = 0.d0
             cv          = 0.d0
             do    jref = 0, lrat(2)-1
                do iref = 0, lrat(1)-1
                   cv          = cv          +                         (fv(ii+iref,jj+jref) &
                        * vfrc(ii+iref,jj+jref))
                   crse(i,j,n) = crse(i,j,n) + fine(ii+iref,jj+jref,n)*(fv(ii+iref,jj+jref) &
                        * vfrc(ii+iref,jj+jref))
                end do
             end do
             if (cv .gt. 1.d-30) then
                crse(i,j,n) = crse(i,j,n) / cv
             else
                crse(i,j,n) = fine(ii,jj,n)  ! covered cell
             end if
          end do
       end do
    end do
  end subroutine amrex_eb_avgdown_sv


  subroutine amrex_eb_avgdown (lo, hi, fine, flo, fhi, crse, clo, chi, &
       vfrc, vflo, vfhi, lrat, ncomp) bind(c,name='amrex_eb_avgdown')
    implicit none
    integer, intent(in) :: lo(2), hi(2), flo(2), fhi(2), clo(2), chi(2), &
         vflo(2), vfhi(2), lrat(2), ncomp
    real(amrex_real), intent(in   ) :: fine( flo(1): fhi(1), flo(2): fhi(2),ncomp)
    real(amrex_real), intent(inout) :: crse( clo(1): chi(1), clo(2): chi(2),ncomp)
    real(amrex_real), intent(in   ) :: vfrc(vflo(1):vfhi(1),vflo(2):vfhi(2))
    
    integer :: i, j, ii, jj, n, iref, jref
    real(amrex_real) :: cv
    
    do n = 1, ncomp
       do j     = lo(2), hi(2)
          jj    = j * lrat(2)
          do i  = lo(1), hi(1)
             ii = i * lrat(1)
             crse(i,j,n) = 0.d0
             cv          = 0.d0
             do    jref = 0, lrat(2)-1
                do iref = 0, lrat(1)-1
                   cv  = cv + vfrc(ii+iref,jj+jref)
                   crse(i,j,n) = crse(i,j,n) + fine(ii+iref,jj+jref,n) * vfrc(ii+iref,jj+jref)
                end do
             end do
             if (cv .gt. 1.d-30) then
                crse(i,j,n) = crse(i,j,n) / cv
             else
                crse(i,j,n) = fine(ii,jj,n)  ! covered cell
             end if
          end do
       end do
    end do
  end subroutine amrex_eb_avgdown

  subroutine amrex_eb_avgdown_faces (lo, hi, fine, flo, fhi, crse, clo, chi, &
       ap, aplo, aphi, lrat, idir, ncomp) bind(c,name='amrex_eb_avgdown_faces')
    implicit none
    integer, intent(in) :: lo(2), hi(2), flo(2), fhi(2), clo(2), chi(2), aplo(2), aphi(2), lrat(2), idir ,ncomp 
    real(amrex_real), intent(in   )   :: fine( flo(1): fhi(1), flo(2): fhi(2), ncomp)
    real(amrex_real), intent(inout)   :: crse( clo(1): chi(1), clo(2): chi(2), ncomp) 
    real(amrex_real), intent(in   )   ::  ap(aplo(1):aphi(1),aplo(2):aphi(2)) 

    integer :: i, j, ii, jj, n, iref, jref, facx, facy
    real(amrex_real) :: fa 
  
    facx = lrat(1) 
    facy = lrat(2)
    if(idir == 0) then
      do n = 1, ncomp 
         do j       = lo(2), hi(2) 
            jj      = j*facy
            do i    = lo(1), hi(1)
               ii   = i*facx
               crse(i,j,n) = 0.d0 
               fa          = 0.d0 
               do    jref  = 0, facy-1
                     fa          = fa + ap(ii, jj+jref) 
                     crse(i,j,n) = crse(i,j,n) + ap(ii,jj+jref)*fine(ii,jj+jref,n)
               enddo
               if(fa.gt.1.d-30) then 
                 crse(i,j,n) = crse(i,j,n)/fa
               else 
                 crse(i,j,n) = fine(ii,jj,n) !covered face
               endif 
            enddo
         enddo 
      enddo
    else
      do n = 1, ncomp
         do j       = lo(2), hi(2)
            jj      = j*facy
            do i    = lo(1), hi(1) 
               ii   = i*facx
               crse(i,j,n) = 0.d0
               fa          = 0.d0 
               do    iref  = 0, facx-1
                     fa          = fa + ap(ii+iref, jj) 
                     crse(i,j,n) = crse(i,j,n) + ap(ii+iref,jj)*fine(ii+iref,jj,n)
               enddo
               if(fa.gt.1.d-30) then 
                 crse(i,j,n) = crse(i,j,n)/fa
               else
                 crse(i,j,n) = fine(ii,jj,n) !covered face 
               endif
             enddo
         enddo
      enddo 
   endif
  end subroutine amrex_eb_avgdown_faces 


  subroutine amrex_eb_avgdown_boundaries (lo, hi, fine, flo, fhi, crse, clo, chi, &
       ba, blo, bhi, lrat, ncomp) bind(c,name='amrex_eb_avgdown_boundaries')
    implicit none
    integer, dimension(2), intent(in) :: lo, hi, flo, fhi, clo, chi, blo, bhi, lrat
    integer, intent(in) :: ncomp
    real(amrex_real), intent(in   )   :: fine( flo(1): fhi(1), flo(2): fhi(2), ncomp)
    real(amrex_real), intent(inout)   :: crse( clo(1): chi(1), clo(2): chi(2), ncomp) 
    real(amrex_real), intent(in   )   :: ba  ( blo(1): bhi(1), blo(2): bhi(2)) 

    integer :: i, j, ii, jj, n, iref, jref, facx, facy
    real(amrex_real) :: fa 
  
    facx = lrat(1) 
    facy = lrat(2)
    do n = 1, ncomp 
       do j       = lo(2), hi(2) 
          jj      = j*facy
          do i    = lo(1), hi(1)
             ii   = i*facx
             crse(i,j,n) = 0.d0 
             fa          = 0.d0 
             do    jref  = 0, facy-1
                do iref  = 0, facx-1
                   fa          = fa          + ba(ii+iref,jj+jref)
                   crse(i,j,n) = crse(i,j,n) + ba(ii+iref,jj+jref)*fine(ii+iref,jj+jref,n)
                end do
             enddo
             if(fa.gt.1.d-30) then 
                crse(i,j,n) = crse(i,j,n)/fa
             else 
                crse(i,j,n) = zero
             endif
          enddo
       enddo
    enddo
  end subroutine amrex_eb_avgdown_boundaries


  subroutine amrex_compute_eb_divergence (lo, hi, divu, dlo, dhi, u, ulo, uhi, v, vlo, vhi, &
       ccm, cmlo, cmhi, flag, flo, fhi, vfrc, klo, khi, &
       apx, axlo, axhi, apy, aylo, ayhi, fcx, cxlo, cxhi, fcy, cylo, cyhi, &
       dxinv) bind(c, name='amrex_compute_eb_divergence')
    implicit none
    integer, dimension(2), intent(in) :: lo, hi, dlo, dhi, ulo, uhi, vlo, vhi, klo, khi, &
         cmlo, cmhi, flo, fhi, axlo, axhi, aylo, ayhi, cxlo, cxhi, cylo, cyhi
    real(amrex_real), intent(inout) :: divu(dlo(1):dhi(1),dlo(2):dhi(2))
    real(amrex_real), intent(in   ) ::    u(ulo(1):uhi(1),ulo(2):uhi(2))
    real(amrex_real), intent(in   ) ::    v(vlo(1):vhi(1),vlo(2):vhi(2))
    integer         , intent(in   ) ::  ccm(cmlo(1):cmhi(1),cmlo(2):cmhi(2))
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2))
    real(amrex_real), intent(in   ) :: vfrc( klo(1): khi(1), klo(2): khi(2))
    real(amrex_real), intent(in   ) ::  apx(axlo(1):axhi(1),axlo(2):axhi(2))
    real(amrex_real), intent(in   ) ::  apy(aylo(1):ayhi(1),aylo(2):ayhi(2))
    real(amrex_real), intent(in   ) ::  fcx(cxlo(1):cxhi(1),cxlo(2):cxhi(2))
    real(amrex_real), intent(in   ) ::  fcy(cylo(1):cyhi(1),cylo(2):cyhi(2))
    real(amrex_real), intent(in) :: dxinv(2)

    integer :: i,j,ii,jj
    real(amrex_real) :: fxm, fxp, fym, fyp, fracx, fracy

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (is_covered_cell(flag(i,j))) then
             divu(i,j) = zero
          else if (is_regular_cell(flag(i,j))) then
             divu(i,j) = dxinv(1) * (u(i+1,j)-u(i,j)) + dxinv(2) * (v(i,j+1)-v(i,j))
          else
             fxm = u(i,j)
             if (apx(i,j).ne.zero .and. apx(i,j).ne.one) then
                jj = j + int(sign(one,fcx(i,j)))
                fracy = abs(fcx(i,j))*real(ior(ccm(i-1,jj),ccm(i,jj)),amrex_real)
                fxm = (one-fracy)*fxm + fracy*u(i,jj)
             end if

             fxp = u(i+1,j)
             if (apx(i+1,j).ne.zero .and. apx(i+1,j).ne.one) then
                jj = j + int(sign(one,fcx(i+1,j)))
                fracy = abs(fcx(i+1,j))*real(ior(ccm(i,jj),ccm(i+1,jj)),amrex_real)
                fxp = (one-fracy)*fxp + fracy*u(i+1,jj)
             end if

             fym = v(i,j)
             if (apy(i,j).ne.zero .and. apy(i,j).ne.one) then
                ii = i + int(sign(one,fcy(i,j)))
                fracx = abs(fcy(i,j))*real(ior(ccm(ii,j-1),ccm(ii,j)),amrex_real)
                fym = (one-fracx)*fym + fracx*v(ii,j)
             end if

             fyp = v(i,j+1)
             if (apy(i,j+1).ne.zero .and. apy(i,j+1).ne.one) then
                ii = i + int(sign(one,fcy(i,j+1)))
                fracx = abs(fcy(i,j+1))*real(ior(ccm(ii,j),ccm(ii,j+1)),amrex_real)
                fyp = (one-fracx)*fyp + fracx*v(ii,j+1)
             end if

             divu(i,j) = (one/vfrc(i,j)) * &
                  ( dxinv(1) * (apx(i+1,j)*fxp-apx(i,j)*fxm) &
                  + dxinv(2) * (apy(i,j+1)*fyp-apy(i,j)*fym) ) 

          end if
       end do
    end do
  end subroutine amrex_compute_eb_divergence


  subroutine amrex_eb_avg_fc_to_cc (lo, hi, cc, cclo, cchi, fx, fxlo, fxhi, fy, fylo, fyhi, &
       ax, axlo, axhi, ay, aylo, ayhi, flag, flo, fhi) bind(c,name='amrex_eb_avg_fc_to_cc')
    integer, dimension(2), intent(in) :: lo, hi, cclo, cchi, fxlo, fxhi, fylo, fyhi, &
         axlo, axhi, aylo, ayhi, flo, fhi
    real(amrex_real), intent(inout) :: cc  (cclo(1):cchi(1),cclo(2):cchi(2),2)
    real(amrex_real), intent(in   ) :: fx  (fxlo(1):fxhi(1),fxlo(2):fxhi(2))
    real(amrex_real), intent(in   ) :: fy  (fylo(1):fyhi(1),fylo(2):fyhi(2))
    real(amrex_real), intent(in   ) :: ax  (axlo(1):axhi(1),axlo(2):axhi(2))
    real(amrex_real), intent(in   ) :: ay  (aylo(1):ayhi(1),aylo(2):ayhi(2))
    integer         , intent(in   ) :: flag( flo(1): fhi(1), flo(2): fhi(2))
    
    integer :: i,j

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (is_covered_cell(flag(i,j))) then
             cc(i,j,1) = zero
             cc(i,j,2) = zero
          else
             if (ax(i,j) .eq. zero) then
                cc(i,j,1) = fx(i+1,j)
             else if (ax(i+1,j) .eq. zero) then
                cc(i,j,1) = fx(i,j)
             else
                cc(i,j,1) = half * (fx(i,j) + fx(i+1,j))
             end if

             if (ay(i,j) .eq. zero) then
                cc(i,j,2) = fy(i,j+1)
             else if (ay(i,j+1) .eq. zero) then
                cc(i,j,2) = fy(i,j)
             else
                cc(i,j,2) = half * (fy(i,j) + fy(i,j+1))
             end if
          end if
       end do
    end do
  end subroutine amrex_eb_avg_fc_to_cc

  subroutine amrex_eb_set_covered_nodes (lo, hi, d, dlo, dhi, f, flo, fhi, v, nc) &
       bind(c,name='amrex_eb_set_covered_nodes')
    integer, intent(in) :: lo(2), hi(2), dlo(2), dhi(2), flo(2), fhi(2), nc
    real(amrex_real), intent(inout) :: d(dlo(1):dhi(1),dlo(2):dhi(2),nc)
    integer, intent(in) :: f(flo(1):fhi(1),flo(2):fhi(2))
    real(amrex_real), intent(in) :: v(nc)

    integer :: i, j, n

    do n = 1, nc
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (is_covered_cell(f(i-1,j-1)) .and. &
                 is_covered_cell(f(i  ,j-1)) .and. &
                 is_covered_cell(f(i-1,j  )) .and. &
                 is_covered_cell(f(i  ,j  ))) then
                d(i,j,n) = v(n)
             end if
          end do
       end do
    end do

  end subroutine amrex_eb_set_covered_nodes

  ! Interpolate face-based variable from face center to face centroid -- this version 
  !   does one face on all grids
  subroutine amrex_eb_interpolate_to_face_centroid ( lo, hi, ivar, var, vlo, vhi, ncomp, &
        areafrac, alo, ahi, cent, clo, chi, flags, flo, fhi, face_type  ) &
       bind(c,name='amrex_eb_interpolate_to_face_centroid')

      use amrex_ebcellflag_module, only: is_covered_cell, get_neighbor_cells

      ! Tile bounds ( face centered )
      integer,  intent(in   ) :: lo(2),  hi(2)

      ! Array Bounds
      integer,  intent(in   ) :: vlo(2), vhi(2)
      integer,  intent(in   ) :: alo(2), ahi(2)
      integer,  intent(in   ) :: clo(2), chi(2)
      integer,  intent(in   ) :: flo(2), fhi(2)
      integer,  intent(in   ) :: ncomp

      ! Type of face (1=x, 2=y)
      integer,  intent(in   ) :: face_type

      ! Arrays
      real(amrex_real), intent(inout) ::                            &
           & ivar(vlo(1):vhi(1),vlo(2):vhi(2),ncomp)  ! Interpolated Variable

      real(amrex_real), intent(inout) ::                            &
           &      var(vlo(1):vhi(1),vlo(2):vhi(2),ncomp), &
           & areafrac(alo(1):ahi(1),alo(2):ahi(2)      ),  &
           &     cent(clo(1):chi(1),clo(2):chi(2),    2)

      integer, intent(in   ) :: &
           & flags(flo(1):fhi(1),flo(2):fhi(2))

      ! Local variables
      integer          :: i, j, n, nbr(-1:1,-1:1)
      real(amrex_real) :: fracx, fracy

      select case ( face_type )
      case(1) ! >>>>>>>>>>>>>>>>>>>>>>  X-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !
         do n = 1, ncomp
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if ( ( areafrac(i,j) > zero ) .and. ( areafrac(i,j) < one ) ) then

                        call get_neighbor_cells( flags(i,j), nbr )

                        if ( cent(i,j,1) < zero ) then
                           fracy = - cent(i,j,1) * nbr(0,-1)
                           ivar(i,j,n) = fracy * var(i,j-1,n) + (one-fracy) * var(i,j,n)
                        else
                           fracy = cent(i,j,1) * nbr(0,1)
                           ivar(i,j,n) = fracy * var(i,j+1,n) + (one-fracy) * var(i,j,n)
                        end if
                     else
                        ivar(i,j,n) = var(i,j,n)
                     end if
               end do
            end do
         end do

      case(2)  ! >>>>>>>>>>>>>>>>>>>>>>  Y-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

         do n = 1, ncomp
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if ( ( areafrac(i,j) > zero ) .and. ( areafrac(i,j) < one ) ) then

                        call get_neighbor_cells( flags(i,j), nbr )

                        if ( cent(i,j,1) < zero ) then
                           fracx = - cent(i,j,1) * nbr(-1,0)
                           ivar(i,j,n) = fracx * var(i-1,j,n) + (one-fracx) * var(i,j,n)
                        else
                           fracx = cent(i,j,1) * nbr(1,0)
                           ivar(i,j,n) = fracx * var(i+1,j,n) + (one-fracx) * var(i,j,n)
                        end if
                     else
                        ivar(i,j,n) = var(i,j,n)
                     end if
               end do
            end do
         end do


      case default

         write(*,*) "amrex_eb_interpolate_to_face_centroid(): face_type = ", face_type, " but valid values are 1,2"
         stop

      end select

  end subroutine amrex_eb_interpolate_to_face_centroid

   !
   ! Returns flux at face centroid in direction dir for just cell (i,j) -- 
   !         note nbr is passed in 
   !
   function amrex_eb_interpolate_to_face_centroid_per_cell ( i, j, dir, var, vlo,  n,  &
        afrac, alo, cent, clo, nbr )  result(ivar)

      use amrex_ebcellflag_module, only: is_covered_cell
      use amrex_error_module,      only: amrex_abort

      ! Face indices: these must be consistent with a staggered indexing
      ! and therefore consistent with the value of dir
      integer,  intent(in   ) :: i, j

      ! Direction of staggering (1=x, 2=y): this specify how (i,j) must
      ! be interpreted, i.e. which staggered numbering the indexing refer to
      integer,  intent(in   ) :: dir

      ! The component to interpolate
      integer,  intent(in   ) :: n

      ! Array Bounds ( only start index )
      integer,  intent(in   ) :: vlo(2), alo(2), clo(2)

      ! Arrays
      real(amrex_real), intent(in   ) ::  &
           &   var(vlo(1):, vlo(2):, 1:), &
           & afrac(alo(1):, alo(2):    ), &
           &  cent(clo(1):, clo(2):, 1:)

      ! Neighbors information
      integer,  intent(in   ) :: nbr(-1:1,-1:1)

      ! Output: the interpolated value
      real(amrex_real)               :: ivar

      ! Local variables
      real(amrex_real)               :: fracx, fracy

      select case ( dir )
      case(1) ! >>>>>>>>>>>>>>>>>>>>>>  X-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

         if ( afrac(i,j) == zero ) then
            ivar = zero
         else if ( afrac(i,j) == one ) then
            ivar = var(i,j,n)
         else
            if ( cent(i,j,1) < zero ) then
               fracy = - cent(i,j,1) * nbr(0,-1)
               ivar = fracy * var(i,j-1,n) + (one-fracy) * var(i,j,n)
            else
               fracy = cent(i,j,1) * nbr(0,1)
               ivar = fracy * var(i,j+1,n) + (one-fracy) * var(i,j,n)
            end if
         end if


      case(2)  ! >>>>>>>>>>>>>>>>>>>>>>  Y-face <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

         if ( afrac(i,j) == zero ) then
            ivar = zero
         else if ( afrac(i,j) == one ) then
            ivar = var(i,j,n)
         else
            if ( cent(i,j,1) < zero ) then
               fracx = - cent(i,j,1) * nbr(-1,0)
               ivar = fracx * var(i-1,j,n) + (one-fracx) * var(i,j,n)
            else
               fracx = cent(i,j,1) * nbr(1,0)
               ivar = fracx * var(i+1,j,n) + (one-fracx) * var(i,j,n)
            end if
         end if

      case default

         call amrex_abort( "interpolate_to_face_centroid(): value of 'dir'"&
              //" is invalid. Must be either 1 or 2")

      end select

   end function amrex_eb_interpolate_to_face_centroid_per_cell

end module amrex_eb_util_module
