module eb_levelset
    use amrex_fort_module, only: c_real => amrex_real
    use iso_c_binding,     only: c_int

    use amrex_ebcellflag_module, only: is_single_valued_cell

    implicit none

contains

    !-----------------------------------------------------------------------------------------------------------------!
    !                                                                                                                 !
    !   pure subroutine INIT_LEVELSET                                                                                 !
    !                                                                                                                 !
    !   Purpose: Initializes level-set array to the fortran huge(real(c_real)) value. This way these values of the    !
    !   level-set function will be overwritten by the min() function (used in these levelset_update function).        !
    !                                                                                                                 !
    !   Comments: If you want to "clear" the whole level-set array (phi), make sure that lo and hi match the grown    !
    !   tile box (i.e. mfi.growntilebox()).                                                                           !
    !                                                                                                                 !
    !-----------------------------------------------------------------------------------------------------------------!

    pure subroutine init_levelset(lo,  hi,          &
                                  phi, phlo, phhi ) &
                    bind(C, name="init_levelset")

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

    end subroutine init_levelset



    !----------------------------------------------------------------------------------------------------------------!
    !                                                                                                                !
    !   pure subroutine FILL_LEVELSET                                                                                !
    !                                                                                                                !
    !   Purpose: given a list of EB-facets, fill the level-set MultiFab between `lo` and `hi` with the closets       !
    !   distance to the EB-facests. Also fill a iMultiFab with 0's and 1's. 0 Indicating that the closest distance   !
    !   was not the result of a projection onto a facet's plane, but instead onto an edge/corner.                    !
    !                                                                                                                !
    !   Comments: Distances are **signed** ( < 0 => inside an EB). Note that at this point the algorithm assumes an  !
    !   edge case lies inside an EB (i.e. its distance is negative). These points are marked using valid = 0. We     !
    !   recommend that the EB's implicit function is used to check these cases (i.e. validate the level-set).        !
    !                                                                                                                !
    !----------------------------------------------------------------------------------------------------------------!

    pure subroutine fill_levelset_eb(lo,      hi,          &
                                     eb_list, l_eb,        &
                                     valid,   vlo,  vhi,   &
                                     phi,     phlo, phhi,  &
                                     dx,      dx_eb      ) &
                     bind(C, name="fill_levelset_eb")

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

    contains


        !------------------------------------------------------------------------------------------------------------!
        !                                                                                                            !
        !   pure subroutine CLOSEST_DIST                                                                             !
        !                                                                                                            !
        !   Purpose: Find the distance to the closets point on the surface defined by the EB-facet list (from the    !
        !   point `pos`). Note that this distance is **signed** => if the vector `eb_center - pos` points            !
        !   **towards** the surface.                                                                                 !
        !                                                                                                            !
        !   Comments: sometimes the closes point is on an EB-facet edge. In this case, the surface normal is not     !
        !   trivial, and this algorithm defaults to a negative distance. However this point is given an `valid`      !
        !   flag of false. It is recommended that the EB's implicit function is used to determine the weather the    !
        !   lies in the EB interior.                                                                                 !
        !                                                                                                            !
        !----------------------------------------------------------------------------------------------------------- !

        pure subroutine closest_dist(min_dist, proj_valid,  &
                                     eb_data,  l_eb, dx_eb, &
                                     pos                   )

            use eb_geometry, only: facets_nearest_pt

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

    end subroutine fill_levelset_eb



    !----------------------------------------------------------------------------------------------------------------!
    !                                                                                                                !
    !   pure subroutine VALIDATE_LEVELSET                                                                            !
    !                                                                                                                !
    !   Purpose: ensure that the sign of the level-set (`phi`) function is the same as the sign of the implicit      !
    !   function (`impf`), if `valid == 0`.                                                                          !
    !                                                                                                                !
    !   Comments: Note the role of `valid` above, `valid` is 0 for points where the closest distance to the surface  !
    !   is on an edge. In that case, fill_levelset_eb defaults to negative distances. Hence this function can be     !
    !   used for checking this assumption, and flipping the level-set value's sign if necessary.                     !
    !                                                                                                                !
    !----------------------------------------------------------------------------------------------------------------!

    pure subroutine validate_levelset(lo,    hi,   n_pad, &
                                      impf,  imlo, imhi,  &
                                      valid, vlo,  vhi,   &
                                      phi,   phlo, phhi  )&
                    bind(C, name="validate_levelset")

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
                        levelset_node = dabs( phi(ii, jj, kk) )
                        if ( impf(ii, jj, kk) <= 0 ) then
                            phi(ii, jj, kk) = levelset_node
                        else
                            phi(ii, jj, kk) = -levelset_node
                        end if
                    end if
                end do
            end do
        end do

    end subroutine validate_levelset



    !----------------------------------------------------------------------------------------------------------------!
    !                                                                                                                !
    !   pure subroutine UPDATE_LEVELSET_INTERSECTION                                                                 !
    !                                                                                                                !
    !   Purpose: Update level set using the "intersection" selection rule: the minimum value between `phi` and       !
    !   `ls_in` is stored in `phi`. This is the level-set equivalent of the `GeometryShop::IntersectionIF`.          !
    !                                                                                                                !
    !   Comments: The role of `valid` here is to flag cells near (or in) regions of negative level-set.              !
    !                                                                                                                !
    !----------------------------------------------------------------------------------------------------------------!

    pure subroutine update_levelset_intersection(lo,    hi,           &
                                                 v_in,  vilo, vihi,   &
                                                 ls_in, lslo, lshi,   &
                                                 valid, vlo,  vhi,    &
                                                 phi,   phlo, phhi,   &
                                                 dx,    n_pad       ) &
                    bind(C, name="update_levelset_intersection")

        implicit none

        integer,      dimension(3), intent(in   ) :: lo, hi, vilo, vihi, lslo, lshi, vlo, vhi, phlo, phhi
        integer,                    intent(in   ) :: v_in  (vilo(1):vihi(1),vilo(2):vihi(2),vilo(3):vihi(3))
        real(c_real),               intent(in   ) :: ls_in (lslo(1):lshi(1),lslo(2):lshi(2),lslo(3):lshi(3))
        integer,                    intent(  out) :: valid ( vlo(1):vhi(1),  vlo(2):vhi(2),  vlo(3):vhi(3) )
        real(c_real),               intent(  out) :: phi   (phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))
        real(c_real), dimension(3), intent(in   ) :: dx
        integer,                    intent(in   ) :: n_pad

        real(c_real) :: ls_node, in_node
        integer      :: i, j, k, ii, jj, kk
        logical      :: valid_cell

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

    end subroutine update_levelset_intersection



    !----------------------------------------------------------------------------------------------------------------!
    !                                                                                                                !
    !   pure subroutine UPDATE_LEVELSET_UNION                                                                        !
    !                                                                                                                !
    !   Purpose: Update level set using the "union" selection rule: the maximum value between `phi` and `ls_in` is   !
    !   stored in `phi`. This is the level-set equivalent of the `GeometryShop::IntersectionIF`.                     !
    !                                                                                                                !
    !   Comments: The role of `valid` here is to flag cells near (or in) regions of negative level-set.              !
    !                                                                                                                !
    !----------------------------------------------------------------------------------------------------------------!

    pure subroutine update_levelset_union(lo,    hi,           &
                                          v_in,  vilo, vihi,   &
                                          ls_in, lslo, lshi,   &
                                          valid, vlo,  vhi,    &
                                          phi,   phlo, phhi,   &
                                          dx,    n_pad       ) &
                    bind(C, name="update_levelset_union")

        implicit none

        integer,      dimension(3), intent(in   ) :: lo, hi, vilo, vihi, lslo, lshi, vlo, vhi, phlo, phhi
        integer,                    intent(in   ) :: v_in  (vilo(1):vihi(1),vilo(2):vihi(2),vilo(3):vihi(3))
        real(c_real),               intent(in   ) :: ls_in (lslo(1):lshi(1),lslo(2):lshi(2),lslo(3):lshi(3))
        integer,                    intent(  out) :: valid ( vlo(1):vhi(1),  vlo(2):vhi(2),  vlo(3):vhi(3) )
        real(c_real),               intent(  out) :: phi   (phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))
        real(c_real), dimension(3), intent(in   ) :: dx
        integer,                    intent(in   ) :: n_pad

        real(c_real) :: ls_node, in_node
        integer      :: i, j, k, ii, jj, kk
        logical      :: valid_cell

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

    end subroutine update_levelset_union



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


        !----------------------------------------------------------------------------------------------------!
        ! build neighbour stencil of size 2 * n_pad                                                          !
        !                                ^^^                                                                 !
        !                  *** fudge factor: finite-sized particle can overlap with neighbouring cells       !
        ! note: stencil could be out-of-bounds (due to fudge factor) => bounds-checking                      !
        !----------------------------------------------------------------------------------------------------!

        klo = k - 2 * n_pad
        if ( klo .lt. phlo(3) ) then
            klo = phlo(3)
        end if

        khi = k + 2 * n_pad
        if ( khi .gt. phhi(3) ) then
            khi = phhi(3)
        end if

        jlo = j - 2 * n_pad
        if ( jlo .lt. phlo(2) ) then
            jlo = phlo(2)
        end if

        jhi = j + 2 * n_pad
        if ( jhi .gt. phhi(2) ) then
            jhi = phhi(2)
        end if

        ilo = i - 2 * n_pad
        if ( ilo .lt. phlo(1) ) then
            ilo = phlo(1)
        end if

        ihi = i + 2 * n_pad
        if ( ihi .gt. phhi(1) ) then
            ihi = phhi(1)
        end if


        !---------------------------------------------------------------------------------------------------!
        ! check members of neighbour stencil:                                                               !
        !       cell is "valid" whenever at least one cell in the neighbour stencil has a level-set phi     !
        !       less than, or equal to, 0                                                                   !
        !---------------------------------------------------------------------------------------------------!

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



    pure subroutine count_eb_facets(lo, hi, flag, flo, fhi, n_facets) &
                    bind(C, name="count_eb_facets")

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

    end subroutine count_eb_facets



    pure subroutine eb_as_list(lo,       hi,   c_facets,  &
                               flag,     flo,  fhi,       &
                               norm,     nlo,  nhi,       &
                               bcent,    blo,  bhi,       &
                               list_out, lsize,           &
                               dx                       ) &
                    bind(C, name="eb_as_list")

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

    end subroutine eb_as_list



    pure subroutine interp_levelset(pos, plo,  n_refine, &
                                    phi, phlo, phhi,     &
                                    dx,  phi_interp    ) &
                    bind(C, name="interp_levelset")


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

    end subroutine interp_levelset



    pure subroutine normal_levelset(pos, plo,   n_refine, &
                                    phi, phlo,  phhi,     &
                                    dx,  normal         ) &
                    bind(C, name="normal_levelset")

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

    end subroutine normal_levelset

end module eb_levelset
