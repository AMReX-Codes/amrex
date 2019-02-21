module amrex_eb_levelset_module
    use amrex_fort_module, only: c_real => amrex_real
    use iso_c_binding,     only: c_int

    use amrex_ebcellflag_module, only: is_single_valued_cell

    implicit none

contains

    !-----------------------------------------------------------------------------------------------------------------
    !!
    !>   pure subroutine INIT_LEVELSET
    !!
    !!   Purpose: Initializes level-set array to the fortran huge(real(c_real)) value. This way these values of the
    !!   level-set function will be overwritten by the min() function (used in these levelset_update function).
    !!
    !!   Comments: If you want to "clear" the whole level-set array (phi), make sure that lo and hi match the grown
    !!   tile box (i.e. mfi.growntilebox()).
    !!
    !-----------------------------------------------------------------------------------------------------------------

    pure subroutine amrex_eb_init_levelset(lo,  hi,          &
                                           phi, phlo, phhi ) &
                    bind(C, name="amrex_eb_init_levelset")

        implicit none

        ! ** define I/O dummy variables
        integer,      dimension(3), intent(in   ) :: lo, hi, phlo, phhi
        real(c_real),               intent(  out) :: phi ( phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )

        ! ** define internal variables
        !    i, j, k: loop index variables
        integer :: i, j, k

        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    phi(i, j, k) = huge(phi)
                end do
            end do
        end do

      end subroutine amrex_eb_init_levelset



    !----------------------------------------------------------------------------------------------------------------
    !!
    !>   pure subroutine FILL_LEVELSET
    !!
    !!   Purpose: given a list of EB-facets, fill the level-set MultiFab between `lo` and `hi` with the closets
    !!   distance to the EB-facests. Also fill a iMultiFab with 0's and 1's. 0 Indicating that the closest distance
    !!   was not the result of a projection onto a facet's plane, but instead onto an edge/corner.
    !!
    !!   Comments: Distances are **signed** ( < 0 => inside an EB). Note that at this point the algorithm assumes an
    !!   edge case lies inside an EB (i.e. its distance is negative). These points are marked using valid = 0. We
    !!   recommend that the EB's implicit function is used to check these cases (i.e. validate the level-set).
    !!
    !----------------------------------------------------------------------------------------------------------------

    pure subroutine amrex_eb_fill_levelset(lo,      hi,          &
                                           eb_list, l_eb,        &
                                           valid,   vlo,  vhi,   &
                                           phi,     phlo, phhi,  &
                                           dx,      dx_eb      ) &
                     bind(C, name="amrex_eb_fill_levelset")

        implicit none

        ! ** define I/O dummy variables
        integer,                       intent(in   ) :: l_eb
        integer,      dimension(3),    intent(in   ) :: lo, hi, vlo, vhi, phlo, phhi
        real(c_real), dimension(l_eb), intent(in   ) :: eb_list
        real(c_real),                  intent(  out) :: phi     (phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3))
        integer,                       intent(  out) :: valid   ( vlo(1):vhi(1),   vlo(2):vhi(2),   vlo(3):vhi(3) )
        real(c_real), dimension(3),    intent(in   ) :: dx, dx_eb


        ! ** define internal variables
        !    pos_node:      position of the level-set MultiFab node (where level-set is evaluated)
        !    levelset_node: value of the signed-distance function at pos_node
        real(c_real), dimension(3) :: pos_node
        real(c_real)               :: levelset_node
        !    ii, jj, kk: loop index variables
        !    valid_cell: .true. iff levelset_node is signed (if .false., levelset_node needs to be validated by IF)
        integer :: ii, jj, kk
        logical :: valid_cell

        do kk = lo(3), hi(3)
            do jj = lo(2), hi(2)
                do ii = lo(1), hi(1)
                    pos_node      = (/ ii*dx(1), jj*dx(2), kk*dx(3) /)
                    call closest_dist ( levelset_node, valid_cell, eb_list, l_eb, dx_eb, pos_node)

                    phi(ii, jj, kk) = levelset_node;

                    if ( valid_cell ) then
                        valid(ii, jj, kk) = 1
                    else
                        valid(ii, jj, kk) = 0
                    end if
                end do
            end do
        end do

    end subroutine amrex_eb_fill_levelset


    !---------------------------------------------------------------------------
    !!
    !>   pure subroutine FILL_LEVELSET_LOC
    !!
    !!   Purpose: given a list of EB-facets, fill the level-set MultiFab between
    !!   `lo` and `hi` with the closets distance to the EB-facests. For all
    !!   `abs(ls_guess) > ls_thres `, level-set filling is aborted and the
    !!   level-set is set to +/- `ls_thres`. This avoid unnecessary filling far
    !!   from EB boundaries.
    !!
    !!   Also fill a iMultiFab with 0's and 1's. 0 Indicating that the closest
    !!   distance was not the result of a projection onto a facet's plane, but
    !!   instead onto an edge/corner.
    !!
    !!   Comments: Distances are **signed** ( < 0 => inside an EB). Note that at
    !!   this point the algorithm assumes an edge case lies inside an EB (i.e.
    !!   its distance is negative). These points are marked using valid = 0. We
    !!   recommend that the EB's implicit function is used to check these cases
    !!   (i.e. validate the level-set).
    !!
    !---------------------------------------------------------------------------

    pure subroutine amrex_eb_fill_levelset_loc(lo,       hi,              &
                                               eb_list,  l_eb,            &
                                               valid,    v_lo,   v_hi,    &
                                               phi,      p_lo,   p_hi,    &
                                               ls_guess, lsg_lo, lsg_hi,  &
                                               ls_thres, dx,     dx_eb  ) &
                    bind(C, name="amrex_eb_fill_levelset_loc")

      implicit none

      ! ** define I/O dummy variables
      integer,                       intent(in   ) :: l_eb
      integer,      dimension(3),    intent(in   ) :: lo, hi, v_lo, v_hi, p_lo, p_hi, lsg_lo, lsg_hi
      real(c_real), dimension(l_eb), intent(in   ) :: eb_list
      real(c_real), dimension(3),    intent(in   ) :: dx, dx_eb
      real(c_real),                  intent(in   ) :: ls_thres

      real(c_real),  intent(  out) :: phi     (  p_lo(1):p_hi(1),     p_lo(2):p_hi(2),     p_lo(3):p_hi(3)  )
      integer,       intent(  out) :: valid   (  v_lo(1):v_hi(1),     v_lo(2):v_hi(2),     v_lo(3):v_hi(3)  )
      real(c_real),  intent(in   ) :: ls_guess(lsg_lo(1):lsg_hi(1), lsg_lo(2):lsg_hi(2), lsg_lo(3):lsg_hi(3))

      ! ** define internal variables
      !    pos_node:      position of the level-set MultiFab node (where level-set is evaluated)
      !    levelset_node: value of the signed-distance function at pos_node
      real(c_real), dimension(3) :: pos_node
      real(c_real)               :: levelset_node
      !    ii, jj, kk: loop index variables
      !    valid_cell: .true. iff levelset_node is signed (if .false., levelset_node needs to be validated by IF)
      integer      :: ii, jj, kk
      logical      :: valid_cell
      real(c_real) :: phi_th

      phi_th = ls_thres
      if (phi_th < 0) phi_th = huge(phi_th)

      do kk = lo(3), hi(3)
         do jj = lo(2), hi(2)
            do ii = lo(1), hi(1)
               if ( abs(ls_guess(ii, jj, kk)) .lt.  phi_th ) then

                  pos_node = (/ ii*dx(1), jj*dx(2), kk*dx(3) /)
                  call closest_dist ( levelset_node, valid_cell, eb_list, l_eb, dx_eb, pos_node)

                  phi(ii, jj, kk) = levelset_node;

                  if ( valid_cell ) then
                     valid(ii, jj, kk) = 1
                  else
                     valid(ii, jj, kk) = 0
                  end if

               else if ( ls_guess(ii, jj, kk) .le. -phi_th ) then

                  phi(ii, jj, kk) = - phi_th
                  valid(ii, jj, kk) = 0

               else if ( ls_guess(ii, jj, kk) .ge.  phi_th ) then

                  phi(ii, jj, kk) = phi_th
                  valid(ii, jj, kk) = 0

               end if
            end do
         end do
      end do

    end subroutine amrex_eb_fill_levelset_loc


    pure subroutine amrex_eb_fill_levelset_bcs( phi,      philo, phihi, &
                                                valid,    vlo,   vhi,   &
                                                periodic, domlo, domhi, &
                                                eb_list, l_eb,          &
                                                dx, dx_eb              )&
                    bind(C, name="amrex_eb_fill_levelset_bcs")

        implicit none

        integer,      dimension(3), intent(in   ) :: philo, phihi, vlo, vhi, periodic, domlo, domhi
        real(c_real),               intent(  out) :: phi(philo(1):phihi(1), philo(2):phihi(2), philo(3):phihi(3))
        integer,                    intent(  out) :: valid(vlo(1):vhi(1),     vlo(2):vhi(2),     vlo(3):vhi(3))

        ! ** extra data used by fill levelset operation
        real(c_real), dimension(3),    intent(in   ) :: dx, dx_eb
        integer,                       intent(in   ) :: l_eb
        real(c_real), dimension(l_eb), intent(in   ) :: eb_list


        integer, dimension(3) :: lo, hi

        !-------------------------------------------------------------------------------------------------------------
        ! Iterate over each of the 6 "faces" of the rectangular domain
        !-------------------------------------------------------------------------------------------------------------

        ! 2 i-j faces => k is in [philo(3), domlo(3)) U (domhi(3), phihi(3)]

        if ( (periodic(3).eq.0) .and. (philo(3).lt.domlo(3)) ) then
            lo(:) = philo(:)      ! k = philo(3), domlo(3) - 1
            hi(:) = phihi(:)      ! j = philo(2), phihi(2)
            hi(3) = domlo(3) - 1  ! i = philo(1), phihi(1)

            call amrex_eb_fill_levelset ( &
                lo, hi,                   &
                eb_list, l_eb,            &
                valid, vlo,   vhi,        &
                phi,   philo, phihi,      &
                dx, dx_eb                 &
            )
        end if

        if ( (periodic(3).eq.0) .and. (phihi(3).gt.domhi(3)) ) then
            lo(:) = philo(:)      ! k = domhi(3) + 1, phihi(3)
            hi(:) = phihi(:)      ! j = philo(2), phihi(2)
            lo(3) = domhi(3) + 1  ! i = philo(1), phihi(1)

            call amrex_eb_fill_levelset ( &
                lo, hi,                   &
                eb_list, l_eb,            &
                valid, vlo,   vhi,        &
                phi,   philo, phihi,      &
                dx, dx_eb                 &
            )
        end if


        ! 2 i-k faces => j is in [philo(2), domlo(2)) U (domhi(2), phihi(2)]

        if ( (periodic(2).eq.0) .and. (philo(2).lt.domlo(2)) ) then
            lo(:) = philo(:)      ! k = philo(3), phihi(3)
            hi(:) = phihi(:)      ! j = philo(2), domlo(2) - 1
            hi(2) = domlo(2) - 1  ! i = philo(1), phihi(1)

            call amrex_eb_fill_levelset ( &
                lo, hi,                   &
                eb_list, l_eb,            &
                valid, vlo,   vhi,        &
                phi,   philo, phihi,      &
                dx, dx_eb                 &
            )
        end if

        if ( (periodic(2).eq.0) .and. (phihi(2).gt.domhi(2)) ) then
            lo(:) = philo(:)      ! k = philo(3), phihi(3)
            hi(:) = phihi(:)      ! j = domhi(2) + 1, phihi(2)
            lo(2) = domhi(2) + 1  ! i = philo(1), phihi(1)

            call amrex_eb_fill_levelset ( &
                lo, hi,                   &
                eb_list, l_eb,            &
                valid, vlo,   vhi,        &
                phi,   philo, phihi,      &
                dx, dx_eb                 &
            )
        end if


        ! 2 j-k faces => i is in [philo(1), domlo(1)) U (domhi(1), phihi(1)]

        if ( (periodic(1).eq.0) .and. (philo(1).lt.domlo(1)) ) then
            lo(:) = philo(:)      ! k = philo(3), phihi(3)
            hi(:) = phihi(:)      ! j = philo(2), phihi(2)
            hi(1) = domlo(1) - 1  ! i = philo(1), domlo(1) - 1

            call amrex_eb_fill_levelset ( &
                lo, hi,                   &
                eb_list, l_eb,            &
                valid, vlo,   vhi,        &
                phi,   philo, phihi,      &
                dx, dx_eb                 &
            )
        end if

        if ( (periodic(1).eq.0) .and. (phihi(1).gt.domhi(1)) ) then
            lo(:) = philo(:)      ! k = philo(3), phihi(3)
            hi(:) = phihi(:)      ! j = philo(2), phihi(2)
            lo(1) = domhi(1) + 1  ! i = domhi(1) + 1, phihi(1)

            call amrex_eb_fill_levelset ( &
                lo, hi,                   &
                eb_list, l_eb,            &
                valid, vlo,   vhi,        &
                phi,   philo, phihi,      &
                dx, dx_eb                 &
            )
        end if

    end subroutine amrex_eb_fill_levelset_bcs


    !------------------------------------------------------------------------------------------------------------
    !!
    !>   pure subroutine CLOSEST_DIST
    !!
    !!   Purpose: Find the distance to the closets point on the surface defined by the EB-facet list (from the
    !!   point `pos`). Note that this distance is **signed** => if the vector `eb_center - pos` points
    !!   **towards** the surface.
    !!
    !!   Comments: sometimes the closes point is on an EB-facet edge. In this case, the surface normal is not
    !!   trivial, and this algorithm defaults to a negative distance. However this point is given an `valid`
    !!   flag of false. It is recommended that the EB's implicit function is used to determine the weather the
    !!   lies in the EB interior.
    !!
    !-----------------------------------------------------------------------------------------------------------

    pure subroutine closest_dist(min_dist, proj_valid,  &
                                 eb_data,  l_eb, dx_eb, &
                                 pos                   )

      use amrex_eb_geometry_module, only: facets_nearest_pt

      implicit none

      ! ** define I/O dummy variables
      integer,                       intent(in   ) :: l_eb
      logical,                       intent(  out) :: proj_valid
      real(c_real),                  intent(  out) :: min_dist
      real(c_real), dimension(3),    intent(in   ) :: pos, dx_eb
      real(c_real), dimension(l_eb), intent(in   ) :: eb_data


      ! ** define internal variables
      !    i:         loop index variable
      !    i_nearest: index of facet nearest to ps
      integer                    :: i, i_nearest
      !    vi_pt, vi_cent: vector indices (in MultiFab index-space) of:
      !       +------|---> the projection point on the nearest EB facet
      !              +---> the center of the nearest EB facet
      integer,      dimension(3) :: vi_pt, vi_cent
      !    dist_proj:        projected (minimal) distance to the nearest EB facet
      !    dist2, min_dist2: squred distance to the EB facet centre, and square distance to the nearest EB facet
      real(c_real)               :: dist_proj, dist2, min_dist2, min_edge_dist2
      !    ind_dx:           inverse of dx_eb (used to allocate MultiFab indices to position vector)
      !    eb_norm, eb_cent: EB normal and center (LATER: of the nearest EB facet)
      !    eb_min_pt, c_vec: projected point on EB facet (c_vec: onto facet edge)
      real(c_real), dimension(3) :: inv_dx, eb_norm, eb_cent, eb_min_pt, c_vec

      inv_dx(:)  = 1.d0 / dx_eb(:)

      min_dist   = huge(min_dist)
      min_dist2  = huge(min_dist2)
      i_nearest  = 0
      proj_valid = .false.

      ! Find nearest EB facet
      do i = 1, l_eb, 6
         eb_cent(:)   = eb_data(i     : i + 2)
         eb_norm(:)   = eb_data(i + 3 : i + 5)

         dist2        = dot_product( pos(:) - eb_cent(:), pos(:) - eb_cent(:) )

         if ( dist2 < min_dist2 ) then
            min_dist2 = dist2
            i_nearest = i
         end if
      end do


      ! Test if pos "projects onto" the nearest EB facet's interior
      eb_cent(:)   = eb_data(i_nearest     : i_nearest + 2)
      eb_norm(:)   = eb_data(i_nearest + 3 : i_nearest + 5)

      dist_proj = dot_product( pos(:) - eb_cent(:), -eb_norm(:) )
      eb_min_pt(:) = pos(:) + eb_norm(:) * dist_proj

      vi_cent(:) = floor( eb_cent(:) * inv_dx)
      vi_pt(:)   = floor( eb_min_pt(:) * inv_dx);

      ! If projects onto nearest EB facet, then return projected distance
      ! Alternatively: find the nearest point on the EB edge
      if ( all( vi_pt == vi_cent ) ) then
         ! this is a signed distance function
         min_dist   = dist_proj
         proj_valid = .true.
      else
         ! fallback: find the nearest point on the EB edge
         c_vec = facets_nearest_pt(vi_pt, vi_cent, pos, eb_norm, eb_cent, dx_eb)
         min_edge_dist2 = dot_product( c_vec(:) - pos(:), c_vec(:) - pos(:))
         min_dist       = -sqrt( min(min_dist2, min_edge_dist2) )
      end if

    end subroutine closest_dist

    !---------------------------------------------------------------------------
    !!
    !>   pure subroutine THRESHOLD_LEVELSET
    !!
    !!   PURPOSE: sets max/min threshold to the level-set function. This ensures
    !!   that the level-set function is a local description near boundaries.
    !!
    !!   COMMENTS: some applications (such as filling level-set from discrete)
    !!   EB facets become less accurate far from the boundary. Thresholding can
    !!   therefore eliminate spurious level-set values
    !!
    !---------------------------------------------------------------------------

    pure subroutine amrex_eb_threshold_levelset(lo,  hi,     threshold, &
                                                phi, phi_lo, phi_hi    )&
                    bind(C, name="amrex_eb_threshold_levelset")

      implicit none

      integer,      dimension(3), intent(in   ) :: lo, hi, phi_lo, phi_hi
      real(c_real),               intent(in   ) :: threshold
      real(c_real),               intent(inout) :: phi (phi_lo(1):phi_hi(1), &
                                                        phi_lo(2):phi_hi(2), &
                                                        phi_lo(3):phi_hi(3))

      integer      :: ii, jj, kk
      real(c_real) :: phi_th


      phi_th = threshold
      if (phi_th < 0) phi_th = huge(phi_th)


      do kk = lo(3), hi(3)
         do jj = lo(2), hi(2)
            do ii = lo(1), hi(1)

               if ( phi(ii, jj, kk) >  phi_th ) phi(ii, jj, kk) =  phi_th
               if ( phi(ii, jj, kk) < -phi_th ) phi(ii, jj, kk) = -phi_th

            end do
         end do
      end do


    end subroutine amrex_eb_threshold_levelset

    !----------------------------------------------------------------------------------------------------------------
    !!
    !>   pure subroutine VALIDATE_LEVELSET
    !!
    !!   Purpose: ensure that the sign of the level-set (`phi`) function is the same as the sign of the implicit
    !!   function (`impf`), if `valid == 0`.
    !!
    !!   Comments: Note the role of `valid` above, `valid` is 0 for points where the closest distance to the surface
    !!   is on an edge. In that case, fill_levelset_eb defaults to negative distances. Hence this function can be
    !!   used for checking this assumption, and flipping the level-set value's sign if necessary.
    !!
    !----------------------------------------------------------------------------------------------------------------

    pure subroutine amrex_eb_validate_levelset(lo,    hi,   n_pad, &
                                               impf,  imlo, imhi,  &
                                               valid, vlo,  vhi,   &
                                               phi,   phlo, phhi  )&
                    bind(C, name="amrex_eb_validate_levelset")

        implicit none

        integer,      dimension(3), intent(in   ) :: lo, hi, imlo, imhi, vlo, vhi, phlo, phhi
        integer,                    intent(in   ) :: n_pad
        real(c_real),               intent(in   ) :: impf  ( imlo(1):imhi(1), imlo(2):imhi(2), imlo(3):imhi(3) )
        integer,                    intent(in   ) :: valid (  vlo(1):vhi(1),   vlo(2):vhi(2),   vlo(3):vhi(3)  )
        real(c_real),               intent(inout) :: phi   ( phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )

        integer      :: ii, jj, kk
        real(c_real) :: levelset_node

        do kk = lo(3), hi(3)
            do jj = lo(2), hi(2)
                do ii = lo(1), hi(1)
                    if ( valid(ii, jj, kk)  == 0 ) then
                        levelset_node = abs( phi(ii, jj, kk) )
                        if ( impf(ii, jj, kk) <= 0 ) then
                            phi(ii, jj, kk) = levelset_node
                        else
                            phi(ii, jj, kk) = -levelset_node
                        end if
                     end if
                end do
            end do
        end do

    end subroutine amrex_eb_validate_levelset


    pure subroutine amrex_eb_validate_levelset_bcs( phi,      phlo,  phhi,  &
                                                    valid,    vlo,   vhi,   &
                                                    periodic, domlo, domhi, &
                                                    impf,     imlo,  imhi ) &
                    bind(C, name="amrex_eb_validate_levelset_bcs")

        implicit none

        integer, dimension(3), intent(in   ) :: phlo, phhi, vlo, vhi, periodic, domlo, domhi, imlo, imhi
        integer,               intent(in   ) :: valid (  vlo(1):vhi(1),   vlo(2):vhi(2),   vlo(3):vhi(3)  )
        real(c_real),          intent(inout) :: phi   ( phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )
        real(c_real),          intent(in   ) :: impf(imlo(1):imhi(1), imlo(2):imhi(2), imlo(3):imhi(3))

        integer, dimension(3) :: lo, hi
        !-------------------------------------------------------------------------------------------------------------
        ! Iterate over each of the 6 "faces" of the rectangular domain
        !-------------------------------------------------------------------------------------------------------------

        ! 2 i-j faces => k is in [vlo(3), domlo(3)) U (domhi(3), vhi(3)]

        if ( (periodic(3).eq.0) .and. (phlo(3).lt.domlo(3)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), domlo(3) - 1
            hi(:) = phhi(:)       ! j = phlo(2), phhi(2)
            hi(3) = domlo(3) - 1  ! i = phlo(1), phhi(1)

            call amrex_eb_validate_levelset ( &
                lo, hi, 0,                    &
                impf,  imlo, imhi,            &
                valid, vlo,  vhi,             &
                phi,   phlo, phhi             &
            )
        end if

        if ( (periodic(3).eq.0) .and. (phhi(3).gt.domhi(3)) ) then
            lo(:) = phlo(:)       ! k = domhi(3) + 1, phihi(3)
            hi(:) = phhi(:)       ! j = phlo(2), phhi(2)
            lo(3) = domhi(3) + 1  ! i = phlo(1), phhi(1)

            call amrex_eb_validate_levelset ( &
                lo, hi, 0,                    &
                impf,  imlo, imhi,            &
                valid, vlo,  vhi,             &
                phi,   phlo, phhi             &
            )
        end if


        ! 2 i-k faces => j is in [vlo(2), domlo(2)) U (domhi(2), vhi(2)]

        if ( (periodic(2).eq.0) .and. (phlo(2).lt.domlo(2)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), phhi(3)
            hi(:) = phhi(:)       ! j = phlo(2), domlo(2) - 1
            hi(2) = domlo(2) - 1  ! i = phlo(1), phhi(1)

            call amrex_eb_validate_levelset ( &
                lo, hi, 0,                    &
                impf,  imlo, imhi,            &
                valid, vlo,  vhi,             &
                phi,   phlo, phhi             &
            )
        end if

        if ( (periodic(2).eq.0) .and. (phhi(2).gt.domhi(2)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), phhi(3)
            hi(:) = phhi(:)       ! j = domhi(2) + 1, phhi(2)
            lo(2) = domhi(2) + 1  ! i = phlo(1), phhi(1)

            call amrex_eb_validate_levelset ( &
                lo, hi, 0,                    &
                impf,  imlo, imhi,            &
                valid, vlo,  vhi,             &
                phi,   phlo, phhi             &
            )
        end if


        ! 2 j-k faces => i is in [vlo(1), domlo(1)) U (domhi(1), vhi(1)]

        if ( (periodic(1).eq.0) .and. (phlo(1).lt.domlo(1)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), phhi(3)
            hi(:) = phhi(:)       ! j = phlo(2), phhi(2)
            hi(1) = domlo(1) - 1  ! i = phlo(1), domlo(1) - 1

            call amrex_eb_validate_levelset ( &
                lo, hi, 0,                    &
                impf,  imlo, imhi,            &
                valid, vlo,  vhi,             &
                phi,   phlo, phhi             &
            )
        end if

        if ( (periodic(1).eq.0) .and. (phhi(1).gt.domhi(1)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), phhi(3)
            hi(:) = phhi(:)       ! j = phlo(2), phhi(2)
            lo(1) = domhi(1) + 1  ! i = domhi(1) + 1, phhi(1)

            call amrex_eb_validate_levelset ( &
                lo, hi, 0,                    &
                impf,  imlo, imhi,            &
                valid, vlo,  vhi,             &
                phi,   phlo, phhi             &
            )
        end if

    end subroutine amrex_eb_validate_levelset_bcs


    !----------------------------------------------------------------------------------------------------------------
    !!
    !>   pure subroutine UPDATE_LEVELSET_INTERSECTION
    !!
    !!   Purpose: Update level set using the "intersection" selection rule: the minimum value between `phi` and
    !!   `ls_in` is stored in `phi`. This is the level-set equivalent of the `GeometryShop::IntersectionIF`.
    !!
    !!   Comments: The role of `valid` here is to flag cells near (or in) regions of negative level-set.
    !!
    !----------------------------------------------------------------------------------------------------------------

    pure subroutine amrex_eb_update_levelset_intersection(lo,    hi,           &
                                                          v_in,  vilo, vihi,   &
                                                          ls_in, lslo, lshi,   &
                                                          valid, vlo,  vhi,    &
                                                          phi,   phlo, phhi  ) &
                    bind(C, name="amrex_eb_update_levelset_intersection")

        implicit none

        integer,      dimension(3), intent(in   ) :: lo, hi, vilo, vihi, lslo, lshi, vlo, vhi, phlo, phhi
        integer,                    intent(in   ) :: v_in  (vilo(1):vihi(1),vilo(2):vihi(2),vilo(3):vihi(3))
        real(c_real),               intent(in   ) :: ls_in (lslo(1):lshi(1),lslo(2):lshi(2),lslo(3):lshi(3))
        integer,                    intent(inout) :: valid ( vlo(1):vhi(1),  vlo(2):vhi(2),  vlo(3):vhi(3) )
        real(c_real),               intent(inout) :: phi   (phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))

        real(c_real) :: ls_node, in_node
        integer      :: ii, jj, kk

        do kk = lo(3), hi(3)
            do jj = lo(2), hi(2)
                do ii = lo(1), hi(1)

                    if (v_in(ii, jj, kk) .eq. 1 ) then
                        in_node = ls_in(ii, jj, kk)
                        ls_node = phi(ii, jj, kk)

                        if ( in_node .lt. ls_node ) then
                            phi(ii, jj, kk) = in_node
                            if ( ls_node .le. 0 ) then
                                valid(ii, jj, kk) = 1
                            end if
                        end if
                    end if

                end do
            end do
        end do

    end subroutine amrex_eb_update_levelset_intersection



    pure subroutine amrex_eb_update_levelset_intersection_bcs( v_in,     vilo,  vihi,   &
                                                               ls_in,    lslo,  lshi,   &
                                                               valid,    vlo,   vhi,    &
                                                               phi,      phlo,  phhi,   &
                                                               periodic, domlo, domhi ) &
                    bind(C, name="amrex_eb_update_levelset_intersection_bcs")

        implicit none

        integer, dimension(3), intent(in   ) :: vilo, vihi, lslo, lshi, vlo, vhi, phlo, phhi, periodic, domlo, domhi
        integer,               intent(in   ) :: v_in  (vilo(1):vihi(1),vilo(2):vihi(2),vilo(3):vihi(3))
        real(c_real),          intent(in   ) :: ls_in (lslo(1):lshi(1),lslo(2):lshi(2),lslo(3):lshi(3))
        integer,               intent(  out) :: valid ( vlo(1):vhi(1),  vlo(2):vhi(2),  vlo(3):vhi(3) )
        real(c_real),          intent(  out) :: phi   (phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))

        integer, dimension(3) :: lo, hi

        !-------------------------------------------------------------------------------------------------------------
        ! Iterate over each of the 6 "faces" of the rectangular domain
        !-------------------------------------------------------------------------------------------------------------

        ! 2 i-j faces => k is in [vlo(3), domlo(3)) U (domhi(3), vhi(3)]

        if ( (periodic(3).eq.0) .and. (phlo(3).lt.domlo(3)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), domlo(3) - 1
            hi(:) = phhi(:)       ! j = phlo(2), phhi(2)
            hi(3) = domlo(3) - 1  ! i = phlo(1), phhi(1)

            call amrex_eb_update_levelset_intersection ( &
                lo,    hi,                               &
                v_in,  vilo, vihi,                       &
                ls_in, lslo, lshi,                       &
                valid, vlo,  vhi,                        &
                phi,   phlo, phhi                        &
            )
        end if

        if ( (periodic(3).eq.0) .and. (phhi(3).gt.domhi(3)) ) then
            lo(:) = phlo(:)       ! k = domhi(3) + 1, phihi(3)
            hi(:) = phhi(:)       ! j = phlo(2), phhi(2)
            lo(3) = domhi(3) + 1  ! i = phlo(1), phhi(1)

            call amrex_eb_update_levelset_intersection ( &
                lo,    hi,                               &
                v_in,  vilo, vihi,                       &
                ls_in, lslo, lshi,                       &
                valid, vlo,  vhi,                        &
                phi,   phlo, phhi                        &
            )
        end if


        ! 2 i-k faces => j is in [vlo(2), domlo(2)) U (domhi(2), vhi(2)]

        if ( (periodic(2).eq.0) .and. (phlo(2).lt.domlo(2)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), phhi(3)
            hi(:) = phhi(:)       ! j = phlo(2), domlo(2) - 1
            hi(2) = domlo(2) - 1  ! i = phlo(1), phhi(1)

            call amrex_eb_update_levelset_intersection ( &
                lo,    hi,                               &
                v_in,  vilo, vihi,                       &
                ls_in, lslo, lshi,                       &
                valid, vlo,  vhi,                        &
                phi,   phlo, phhi                        &
            )
        end if

        if ( (periodic(2).eq.0) .and. (phhi(2).gt.domhi(2)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), phhi(3)
            hi(:) = phhi(:)       ! j = domhi(2) + 1, phhi(2)
            lo(2) = domhi(2) + 1  ! i = phlo(1), phhi(1)

            call amrex_eb_update_levelset_intersection ( &
                lo,    hi,                               &
                v_in,  vilo, vihi,                       &
                ls_in, lslo, lshi,                       &
                valid, vlo,  vhi,                        &
                phi,   phlo, phhi                        &
            )
        end if


        ! 2 j-k faces => i is in [vlo(1), domlo(1)) U (domhi(1), vhi(1)]

        if ( (periodic(1).eq.0) .and. (phlo(1).lt.domlo(1)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), phhi(3)
            hi(:) = phhi(:)       ! j = phlo(2), phhi(2)
            hi(1) = domlo(1) - 1  ! i = phlo(1), domlo(1) - 1

            call amrex_eb_update_levelset_intersection ( &
                lo,    hi,                               &
                v_in,  vilo, vihi,                       &
                ls_in, lslo, lshi,                       &
                valid, vlo,  vhi,                        &
                phi,   phlo, phhi                        &
            )
        end if

        if ( (periodic(1).eq.0) .and. (phhi(1).gt.domhi(1)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), phhi(3)
            hi(:) = phhi(:)       ! j = phlo(2), phhi(2)
            lo(1) = domhi(1) + 1  ! i = domhi(1) + 1, phhi(1)

            call amrex_eb_update_levelset_intersection ( &
                lo,    hi,                               &
                v_in,  vilo, vihi,                       &
                ls_in, lslo, lshi,                       &
                valid, vlo,  vhi,                        &
                phi,   phlo, phhi                        &
            )
        end if

    end subroutine amrex_eb_update_levelset_intersection_bcs



    !----------------------------------------------------------------------------------------------------------------
    !!
    !>   pure subroutine UPDATE_LEVELSET_UNION
    !!
    !!   Purpose: Update level set using the "union" selection rule: the maximum value between `phi` and `ls_in` is
    !!   stored in `phi`. This is the level-set equivalent of the `GeometryShop::IntersectionIF`.
    !!
    !!   Comments: The role of `valid` here is to flag cells near (or in) regions of negative level-set.
    !!
    !----------------------------------------------------------------------------------------------------------------

    pure subroutine amrex_eb_update_levelset_union(lo,    hi,           &
                                                   v_in,  vilo, vihi,   &
                                                   ls_in, lslo, lshi,   &
                                                   valid, vlo,  vhi,    &
                                                   phi,   phlo, phhi  ) &
                    bind(C, name="amrex_eb_update_levelset_union")

        implicit none

        integer,      dimension(3), intent(in   ) :: lo, hi, vilo, vihi, lslo, lshi, vlo, vhi, phlo, phhi
        integer,                    intent(in   ) :: v_in  (vilo(1):vihi(1),vilo(2):vihi(2),vilo(3):vihi(3))
        real(c_real),               intent(in   ) :: ls_in (lslo(1):lshi(1),lslo(2):lshi(2),lslo(3):lshi(3))
        integer,                    intent(  out) :: valid ( vlo(1):vhi(1),  vlo(2):vhi(2),  vlo(3):vhi(3) )
        real(c_real),               intent(  out) :: phi   (phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))

        real(c_real) :: ls_node, in_node
        integer      :: ii, jj, kk

        do kk = lo(3), hi(3)
            do jj = lo(2), hi(2)
                do ii = lo(1), hi(1)

                    if (v_in(ii, jj, kk) .eq. 1 ) then
                        in_node = ls_in(ii, jj, kk)
                        ls_node = phi(ii, jj, kk)

                        if ( in_node .gt. ls_node ) then
                            phi(ii, jj, kk) = in_node
                            if ( ls_node .le. 0 ) then
                                valid(ii, jj, kk) = 1
                            end if
                        end if
                    end if

                end do
            end do
        end do

    end subroutine amrex_eb_update_levelset_union

    pure subroutine amrex_eb_update_levelset_union_bcs( v_in,     vilo,  vihi,   &
                                                        ls_in,    lslo,  lshi,   &
                                                        valid,    vlo,   vhi,    &
                                                        phi,      phlo,  phhi,   &
                                                        periodic, domlo, domhi ) &
                    bind(C, name="amrex_eb_update_levelset_union_bcs")

        implicit none

        integer, dimension(3), intent(in   ) :: vilo, vihi, lslo, lshi, vlo, vhi, phlo, phhi, periodic, domlo, domhi
        integer,               intent(in   ) :: v_in  (vilo(1):vihi(1),vilo(2):vihi(2),vilo(3):vihi(3))
        real(c_real),          intent(in   ) :: ls_in (lslo(1):lshi(1),lslo(2):lshi(2),lslo(3):lshi(3))
        integer,               intent(  out) :: valid ( vlo(1):vhi(1),  vlo(2):vhi(2),  vlo(3):vhi(3) )
        real(c_real),          intent(  out) :: phi   (phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))

        !integer      :: i, j, k
        !real(c_real) :: ls_node, in_node

        integer, dimension(3) :: lo, hi

        !-------------------------------------------------------------------------------------------------------------
        ! Iterate over each of the 6 "faces" of the rectangular domain
        !-------------------------------------------------------------------------------------------------------------

        ! 2 i-j faces => k is in [vlo(3), domlo(3)) U (domhi(3), vhi(3)]
        if ( (periodic(3).eq.0) .and. (phlo(3).lt.domlo(3)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), domlo(3) - 1
            hi(:) = phhi(:)       ! j = phlo(2), phhi(2)
            hi(3) = domlo(3) - 1  ! i = phlo(1), phhi(1)

            call amrex_eb_update_levelset_union ( &
                lo,    hi,                        &
                v_in,  vilo, vihi,                &
                ls_in, lslo, lshi,                &
                valid, vlo,  vhi,                 &
                phi,   phlo, phhi                 &
            )
        end if

        if ( (periodic(3).eq.0) .and. (phhi(3).gt.domhi(3)) ) then
            lo(:) = phlo(:)       ! k = domhi(3) + 1, phihi(3)
            hi(:) = phhi(:)       ! j = phlo(2), phhi(2)
            lo(3) = domhi(3) + 1  ! i = phlo(1), phhi(1)

            call amrex_eb_update_levelset_union ( &
                lo,    hi,                        &
                v_in,  vilo, vihi,                &
                ls_in, lslo, lshi,                &
                valid, vlo,  vhi,                 &
                phi,   phlo, phhi                 &
            )
        end if


        ! 2 i-k faces => j is in [vlo(2), domlo(2)) U (domhi(2), vhi(2)]

        if ( (periodic(2).eq.0) .and. (phlo(2).lt.domlo(2)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), phhi(3)
            hi(:) = phhi(:)       ! j = phlo(2), domlo(2) - 1
            hi(2) = domlo(2) - 1  ! i = phlo(1), phhi(1)

            call amrex_eb_update_levelset_union ( &
                lo,    hi,                        &
                v_in,  vilo, vihi,                &
                ls_in, lslo, lshi,                &
                valid, vlo,  vhi,                 &
                phi,   phlo, phhi                 &
            )
        end if

        if ( (periodic(2).eq.0) .and. (phhi(2).gt.domhi(2)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), phhi(3)
            hi(:) = phhi(:)       ! j = domhi(2) + 1, phhi(2)
            lo(2) = domhi(2) + 1  ! i = phlo(1), phhi(1)

            call amrex_eb_update_levelset_union ( &
                lo,    hi,                        &
                v_in,  vilo, vihi,                &
                ls_in, lslo, lshi,                &
                valid, vlo,  vhi,                 &
                phi,   phlo, phhi                 &
            )
        end if


        ! 2 j-k faces => i is in [vlo(1), domlo(1)) U (domhi(1), vhi(1)]

        if ( (periodic(1).eq.0) .and. (phlo(1).lt.domlo(1)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), phhi(3)
            hi(:) = phhi(:)       ! j = phlo(2), phhi(2)
            hi(1) = domlo(1) - 1  ! i = phlo(1), domlo(1) - 1

            call amrex_eb_update_levelset_union ( &
                lo,    hi,                        &
                v_in,  vilo, vihi,                &
                ls_in, lslo, lshi,                &
                valid, vlo,  vhi,                 &
                phi,   phlo, phhi                 &
            )
        end if

        if ( (periodic(1).eq.0) .and. (phhi(1).gt.domhi(1)) ) then
            lo(:) = phlo(:)       ! k = phlo(3), phhi(3)
            hi(:) = phhi(:)       ! j = phlo(2), phhi(2)
            lo(1) = domhi(1) + 1  ! i = domhi(1) + 1, phhi(1)

            call amrex_eb_update_levelset_union ( &
                lo,    hi,                        &
                v_in,  vilo, vihi,                &
                ls_in, lslo, lshi,                &
                valid, vlo,  vhi,                 &
                phi,   phlo, phhi                 &
            )
        end if

    end subroutine amrex_eb_update_levelset_union_bcs

    !----------------------------------------------------------------------------------------------------------------
    !!
    !>   pure subroutine FILL_VALID
    !!
    !!   Purpose: Fills elements of valid with 1 whenever it corresponds to a position with n_pad of the level set
    !!   value being negative (i.e. phi < 0), and 0 otherwise.
    !!
    !----------------------------------------------------------------------------------------------------------------

    pure subroutine amrex_eb_fill_valid(       lo,     hi, &
                                        valid, vlo,   vhi, &
                                        phi,   phlo, phhi, &
                                        n_pad            ) &
                    bind(C, name="amrex_eb_fill_valid")

        implicit none

        integer, dimension(3), intent(in   ) :: lo, hi, vlo, vhi, phlo, phhi
        integer,               intent(in   ) :: n_pad
        integer,               intent(  out) :: valid( vlo(1):vhi(1),   vlo(2):vhi(2),   vlo(3):vhi(3))
        real(c_real),          intent(in   ) :: phi  (phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3))

        integer :: i, j, k
        logical :: valid_cell


        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    valid_cell = neighbour_is_valid(phi, phlo, phhi, i, j, k, n_pad)
                    if ( valid_cell ) then
                        valid(i, j, k) = 1
                    end if
                end do
            end do
        end do

    end subroutine amrex_eb_fill_valid



    !----------------------------------------------------------------------------------------------------------------
    !!
    !>   pure subroutine FILL_VALID_BCS
    !!
    !!   Purpose: Fills elements of valid with 1 whenever it corresponds to a position with n_pad of the level set
    !!   value being negative (i.e. phi < 0), and 0 otherwise. For all ghost cells outside the domain
    !!
!----------------------------------------------------------------------------------------------------------------


    pure subroutine amrex_eb_fill_valid_bcs( valid, vlo, vhi, periodic, domlo, domhi ) &
                    bind(C, name="amrex_eb_fill_valid_bcs")

        implicit none

        integer, dimension(3), intent(in   ) :: vlo, vhi, periodic, domlo, domhi
        integer,               intent(  out) :: valid(vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3))

        integer :: i, j, k

        !-------------------------------------------------------------------------------------------------------------
        ! Iterate over each of the 6 "faces" of the rectangular domain
        !-------------------------------------------------------------------------------------------------------------

        ! 2 i-j faces => k is in [vlo(3), domlo(3)) U (domhi(3), vhi(3)]

        if ( periodic(3).eq.0 ) then
            do k = vlo(3), domlo(3) - 1
                do j = vlo(2), vhi(2)
                    do i = vlo(1), vhi(1)
                        valid(i, j, k) = 1
                    end do
                end do
            end do

            do k = domhi(3) + 1, vhi(3)
                do j = vlo(2), vhi(2)
                    do i = vlo(1), vhi(1)
                        valid(i, j, k) = 1
                    end do
                end do
            end do
        end if

        ! 2 i-k faces => j is in [vlo(2), domlo(2)) U (domhi(2), vhi(2)]

        if ( periodic(2).eq.0 ) then
            if ( vlo(2).lt.domlo(2) ) then
                do k = vlo(3), vhi(3)
                    do j = vlo(2), domlo(2) - 1
                        do i = vlo(1), vhi(1)
                            valid(i, j, k) = 1
                        end do
                    end do
                end do
            end if

            if (vhi(2).gt.domhi(2) ) then
                do k = vlo(3), vhi(3)
                    do j = domhi(2) + 1, vhi(2)
                        do i = vlo(1), vhi(1)
                            valid(i, j, k) = 1
                        end do
                    end do
                end do
            end if
        end if

        ! 2 j-k faces => i is in [vlo(1), domlo(1)) U (domhi(1), vhi(1)]

        if ( periodic(1).eq.0 ) then
            if ( vlo(1).lt.domlo(1) ) then
                do k = vlo(3), vhi(3)
                    do j = vlo(2), vhi(2)
                        do i = vlo(1), domlo(1) - 1
                            valid(i, j, k) = 1
                        end do
                    end do
                end do
            end if

            if ( vhi(1).gt.domhi(1) ) then
                do k = vlo(3), vhi(3)
                    do j = vlo(2), vhi(2)
                        do i = domhi(1) + 1, vhi(1)
                            valid(i, j, k) = 1
                        end do
                    end do
                end do
            end if
        end if

    end subroutine amrex_eb_fill_valid_bcs



    pure function neighbour_is_valid(phi, phlo, phhi, i, j, k, n_pad)
        implicit none

        ! ** output type
        logical :: neighbour_is_valid

        ! ** input types
        integer,      dimension(3), intent(in) :: phlo, phhi
        real(c_real),               intent(in) :: phi( phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )
        integer,                    intent(in) :: i, j, k, n_pad


        ! ** declare local variables
        ! ii, jj, kk : loop variables itterating over neighbour stencil
        ! klo ... ihi: boundaries of stencil which will be checked for valid cells
        !              a cell is valid if phi <= 0
        integer :: ii, jj, kk, klo, khi, jlo, jhi, ilo, ihi


        !----------------------------------------------------------------------------------------------------
        ! build neighbour stencil of size n_pad
        ! note: stencil could be out-of-bounds  => bounds-checking
        !----------------------------------------------------------------------------------------------------

        ilo = max(i-n_pad, phlo(1))
        ihi = min(i+n_pad, phhi(1))

        jlo = max(j-n_pad, phlo(2))
        jhi = min(j+n_pad, phhi(2))

        klo = max(k-n_pad, phlo(3))
        khi = min(k+n_pad, phhi(3))


        !---------------------------------------------------------------------------------------------------
        ! check members of neighbour stencil:
        !       cell is "valid" whenever at least one cell in the neighbour stencil has a level-set phi
        !       less than, or equal to, 0
        !---------------------------------------------------------------------------------------------------

        neighbour_is_valid = .false.

        do kk = klo, khi
            do jj = jlo, jhi
                do ii = ilo, ihi
                    if ( phi(ii, jj, kk) .le. 0 ) then
                        neighbour_is_valid = .true.
                        return
                    end if
                end do
            end do
        end do

    end function neighbour_is_valid



    pure subroutine amrex_eb_count_facets(lo, hi, flag, flo, fhi, n_facets) &
                    bind(C, name="amrex_eb_count_facets")

        implicit none

        integer, dimension(3), intent(in   ) :: lo, hi, flo, fhi
        integer,               intent(in   ) :: flag ( flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3) )
        integer,               intent(inout) :: n_facets

        integer :: i, j, k

        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    if ( is_single_valued_cell( flag(i, j, k) ) ) then
                        n_facets = n_facets + 1
                    end if
                end do
            end do
        end do

    end subroutine amrex_eb_count_facets



    pure subroutine amrex_eb_as_list(lo,       hi,   c_facets,  &
                               flag,     flo,  fhi,       &
                               norm,     nlo,  nhi,       &
                               bcent,    blo,  bhi,       &
                               list_out, lsize,           &
                               dx                       ) &
                    bind(C, name="amrex_eb_as_list")

        implicit none

        integer,                        intent(in   ) :: lsize
        integer, dimension(3),          intent(in   ) :: lo, hi, flo, fhi, nlo, nhi, blo, bhi
        real(c_real),                   intent(in   ) :: dx(3)
        integer,                        intent(in   ) :: flag  ( flo(1):fhi(1), flo(2):fhi(2), flo(3):fhi(3) )
        real(c_real),                   intent(in   ) :: norm  ( nlo(1):nhi(1), nlo(2):nhi(2), nlo(3):nhi(3), 3)
        real(c_real),                   intent(in   ) :: bcent ( blo(1):bhi(1), blo(2):bhi(2), blo(3):bhi(3), 3)
        real(c_real), dimension(lsize), intent(  out) :: list_out
        integer,                        intent(inout) :: c_facets

        integer                    :: i, j, k, i_facet
        real(c_real), dimension(3) :: eb_cent

        i_facet = 6 * c_facets + 1

        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    if ( is_single_valued_cell( flag(i, j, k) ) ) then

                        c_facets = c_facets + 1

                        eb_cent(:) = ( bcent(i, j, k, :)                     &
                                       + (/ dble(i), dble(j), dble(k) /)     &
                                       + (/ 0.5d0, 0.5d0, 0.5d0 /) ) * dx(:)

                        list_out( i_facet     : i_facet + 2) = eb_cent(:)
                        list_out( i_facet + 3 : i_facet + 5) = norm(i, j, k, :)

                        i_facet = i_facet + 6

                    end if
                end do
            end do
        end do

    end subroutine amrex_eb_as_list



    pure subroutine amrex_eb_interp_levelset(pos, plo,  n_refine, &
                                             phi, phlo, phhi,     &
                                             dx,  phi_interp    ) &
                    bind(C, name="amrex_eb_interp_levelset")


        implicit none

        real(c_real), dimension(3), intent(in   ) :: pos, plo
        integer,      dimension(3), intent(in   ) :: phlo, phhi
        integer,                    intent(in   ) :: n_refine
        real(c_real),               intent(in   ) :: phi( phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )
        real(c_real), dimension(3), intent(in   ) :: dx
        real(c_real),               intent(  out) :: phi_interp

        integer                    :: i, j, k
        real(c_real)               :: xp, yp, zp, lx, ly, lz, wx_lo, wx_hi, wy_lo, wy_hi, wz_lo, wz_hi
        real(c_real), dimension(3) :: inv_dx

        inv_dx(:) = n_refine / dx(:)

        xp = pos(1) - plo(1)
        yp = pos(2) - plo(2)
        zp = pos(3) - plo(3)

        lx = xp * inv_dx(1)
        ly = yp * inv_dx(2)
        lz = zp * inv_dx(3)

        i = floor(lx)
        j = floor(ly)
        k = floor(lz)

        wx_hi = lx - i
        wy_hi = ly - j
        wz_hi = lz - k

        wx_lo = 1.0d0 - wx_hi
        wy_lo = 1.0d0 - wy_hi
        wz_lo = 1.0d0 - wz_hi

        phi_interp = phi(i,   j,   k  ) * wx_lo * wy_lo * wz_lo &
                   + phi(i+1, j,   k  ) * wx_hi * wy_lo * wz_lo &
                   + phi(i,   j+1, k  ) * wx_lo * wy_hi * wz_lo &
                   + phi(i,   j,   k+1) * wx_lo * wy_lo * wz_hi &
                   + phi(i+1, j+1, k  ) * wx_hi * wy_hi * wz_lo &
                   + phi(i,   j+1, k+1) * wx_lo * wy_hi * wz_hi &
                   + phi(i+1, j,   k+1) * wx_hi * wy_lo * wz_hi &
                   + phi(i+1, j+1, k+1) * wx_hi * wy_hi * wz_hi

    end subroutine amrex_eb_interp_levelset



    pure subroutine amrex_eb_normal_levelset(pos, plo,   n_refine, &
                                             phi, phlo,  phhi,     &
                                             dx,  normal         ) &
                    bind(C, name="amrex_eb_normal_levelset")

        implicit none

        real(c_real), dimension(3), intent(in   ) :: pos, plo
        integer,      dimension(3), intent(in   ) :: phlo, phhi
        integer,                    intent(in   ) :: n_refine
        real(c_real),               intent(in   ) :: phi( phlo(1):phhi(1), phlo(2):phhi(2), phlo(3):phhi(3) )
        real(c_real), dimension(3), intent(in   ) :: dx
        real(c_real), dimension(3), intent(  out) :: normal

        integer                    :: i, j, k
        real(c_real)               :: xp, yp, zp, lx, ly, lz, wx_lo, wx_hi, wy_lo, wy_hi, wz_lo, wz_hi
        real(c_real), dimension(3) :: inv_dx

        real(c_real) :: inv_norm

        inv_dx = n_refine / dx

        xp = pos(1) - plo(1)
        yp = pos(2) - plo(2)
        zp = pos(3) - plo(3)

        lx = xp * inv_dx(1)
        ly = yp * inv_dx(2)
        lz = zp * inv_dx(3)

        i = floor(lx)
        j = floor(ly)
        k = floor(lz)

        wx_hi = lx - i
        wy_hi = ly - j
        wz_hi = lz - k

        wx_lo = 1.0d0 - wx_hi
        wy_lo = 1.0d0 - wy_hi
        wz_lo = 1.0d0 - wz_hi

        normal(1) = - phi(i,   j,   k  )*inv_dx(1) * wy_lo * wz_lo  &
                    + phi(i+1, j,   k  )*inv_dx(1) * wy_lo * wz_lo  &
                    - phi(i,   j+1, k  )*inv_dx(1) * wy_hi * wz_lo  &
                    + phi(i+1, j+1, k  )*inv_dx(1) * wy_hi * wz_lo  &
                    - phi(i,   j,   k+1)*inv_dx(1) * wy_lo * wz_hi  &
                    + phi(i+1, j,   k+1)*inv_dx(1) * wy_lo * wz_hi  &
                    - phi(i,   j+1, k+1)*inv_dx(1) * wy_hi * wz_hi  &
                    + phi(i+1, j+1, k+1)*inv_dx(1) * wy_hi * wz_hi

        normal(2) = - phi(i,   j,   k  )*inv_dx(2) * wx_lo * wz_lo  &
                    + phi(i,   j+1, k  )*inv_dx(2) * wx_lo * wz_lo  &
                    - phi(i+1, j,   k  )*inv_dx(2) * wx_hi * wz_lo  &
                    + phi(i+1, j+1, k  )*inv_dx(2) * wx_hi * wz_lo  &
                    - phi(i,   j,   k+1)*inv_dx(2) * wx_lo * wz_hi  &
                    + phi(i,   j+1, k+1)*inv_dx(2) * wx_lo * wz_hi  &
                    - phi(i+1, j,   k+1)*inv_dx(2) * wx_hi * wz_hi  &
                    + phi(i+1, j+1, k+1)*inv_dx(2) * wx_hi * wz_hi

        normal(3) = - phi(i,   j,   k  )*inv_dx(3) * wx_lo * wy_lo  &
                    + phi(i,   j,   k+1)*inv_dx(3) * wx_lo * wy_lo  &
                    - phi(i+1, j,   k  )*inv_dx(3) * wx_hi * wy_lo  &
                    + phi(i+1, j,   k+1)*inv_dx(3) * wx_hi * wy_lo  &
                    - phi(i,   j+1, k  )*inv_dx(3) * wx_lo * wy_hi  &
                    + phi(i,   j+1, k+1)*inv_dx(3) * wx_lo * wy_hi  &
                    - phi(i+1, j+1, k  )*inv_dx(3) * wx_hi * wy_hi  &
                    + phi(i+1, j+1, k+1)*inv_dx(3) * wx_hi * wy_hi

        ! this might not be necessary if the phi grid is dense enough...
        inv_norm = 1.0d0 / sqrt(normal(1)**2 + normal(2)**2 + normal(3)**2)
        normal(:) = normal(:) * inv_norm

    end subroutine amrex_eb_normal_levelset

end module amrex_eb_levelset_module
