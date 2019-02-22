module amrex_eb_geometry_module

    use amrex_fort_module, only: c_real => amrex_real

    implicit none

contains

   !----------------------------------------------------------------------
   !!
   !>  Pure Function: DOT_3D_REAL
   !!
   !!  Purpose: Returns the cartesian dot product for two vectors in three
   !!  dimensions.
   !!
   !!  Comments: Vectors are represented as one-dimensional arrays of type
   !!  real(c_real) and of dimension(3) (i.e. indices range from 1..3).
   !----------------------------------------------------------------------
    pure function dot_3d_real (v1, v2)
        implicit none

        real(c_real)                           :: dot_3d_real
        real(c_real), dimension(3), intent(in) :: v1, v2

        ! really naive implementation
        dot_3d_real = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)

    end function dot_3d_real

   !----------------------------------------------------------------------
   !!
   !>  Pure Function: CROSS_3D_REAL
   !!
   !!  Purpose: Returns the cartesian cross product for two vectors in
   !!  three dimensions.
   !!
   !!  Comments: Vectors are represented as one-dimensional arrays of type
   !!  real(c_real) and of dimension(3) (i.e. indices range from 1..3).
   !----------------------------------------------------------------------
    pure function cross_3d_real (v1, v2)
        implicit none

        real(c_real), dimension(3)             :: cross_3d_real
        real(c_real), dimension(3), intent(in) :: v1, v2

        cross_3d_real(1) = v1(2)*v2(3) - v1(3)*v2(2)
        cross_3d_real(2) = v1(3)*v2(1) - v1(1)*v2(3)
        cross_3d_real(3) = v1(1)*v2(2) - v1(2)*v2(1)

    end function cross_3d_real

   !----------------------------------------------------------------------
   !!
   !>  Pure Function: PT_IN_BOX
   !!
   !!  Purpose: Returns true if the coordinate vector represents a point
   !!  inside a three-dimensional cube of size dx and origin specified by
   !!  integer vector. As a final input, the user specifies a dimension
   !!  (axis) to ignore. This allows IN-BOX checking on a box face.
   !!
   !!  Comments: Position vectors are represented as one-dimensional arrays
   !!  of type real(c_real) and of dimension(3) (i.e. indices range
   !!  from 1..3). Cells are enumerated using arrays of type integer and
   !!  dimension(3).
   !----------------------------------------------------------------------
   ! pure function pt_in_box (pt, id, id_ignore)
   !     implicit none

   !     logical                                :: pt_in_box
   !     real(c_real), dimension(3), intent(in) :: pt
   !     integer,      dimension(3), intent(in) :: id
   !     integer,                    intent(in) :: id_ignore

   !     real(c_real), dimension(3) :: box_max, box_min
   !     integer                    :: i

   !     ! Determine box boundaries
   !     box_min(:) = id(:) * dx(:)
   !     box_max(:) = (id(:) + 1.) * dx(:)

   !     pt_in_box = .true.

   !     ! Check each coordinate. Skip ignored coordinate.
   !     do i = 1, 3
   !         if (.not. (i .eq. id_ignore) ) then
   !             if ( pt(i) .le. box_min(i) ) then
   !                 pt_in_box = .false.
   !                 exit
   !             end if

   !             if ( pt(i) .ge. box_max(i) ) then
   !                 pt_in_box = .false.
   !                 exit
   !             end if
   !         end if
   !     end do

   ! end function pt_in_box

   !----------------------------------------------------------------------
   !!
   !>  Pure Subroutine: CALC_FACET_EDGE
   !!
   !!  Purpose: Calculates the line (represented by a position and a
   !!  direction vector) given by the intersection of two planes (defined
   !!  by two normal (n1, n2) and two positions (h1 = n1.p1, h2 = n2.p2).
   !!
   !!  When one plane is the EB surface, and the other is a face of the
   !!  cell. Then this line represents the edge of the EB facet.
   !!
   !!  Comments: Vectors are represented as one-dimensional arrays of type
   !!  real(c_real) and of dimension(3) (i.e. indices range from 1..3).
   !----------------------------------------------------------------------
    pure subroutine calc_facet_edge (p0, v, h1, h2, n1, n2)
        implicit none

        real(c_real), dimension(3), intent(  out) :: p0, v
        real(c_real), dimension(3), intent(in   ) :: n1, n2
        real(c_real),               intent(in   ) :: h1, h2

        real(c_real) :: c1, c2, c_dp, c_norm

        c_dp = dot_3d_real(n1, n2)
        c_norm = 1 - c_dp * c_dp

        c1 = ( h1 - h2 * c_dp ) / c_norm
        c2 = ( h2 - h1 * c_dp ) / c_norm

        p0(:) = c1 * n1(:) + c2 * n2(:)
        v = cross_3d_real(n1, n2)

    end subroutine calc_facet_edge

   !----------------------------------------------------------------------
   !!
   !>  Pure Subroutine: LINES_NEAREST_PT
   !!
   !!  Purpose: Given an a line an a point, this subroutine finds the point
   !!  one the line which minimizes the cartesian distance. It also finds
   !!  the corresponing distance along the line corresponding to this point
   !!
   !!  Comments: Vectors are represented as one-dimensional arrays of type
   !!  real(c_real) and of dimension(3) (i.e. indices range from 1..3).
   !----------------------------------------------------------------------
    pure subroutine lines_nearest_pt (lambda_min, nearest_pt, p0, v, pt)
        implicit none

        real(c_real),               intent(  out) :: lambda_min
        real(c_real), dimension(3), intent(  out) :: nearest_pt
        real(c_real), dimension(3), intent(in   ) :: p0, v, pt

        real(c_real), dimension(3) :: c

        c(:) = p0(:) - pt(:)
        lambda_min = - dot_3d_real(v, c) / dot_3d_real(v, v)

        nearest_pt(:) = p0(:) + lambda_min*v(:)

    end subroutine lines_nearest_pt

   !----------------------------------------------------------------------
   !!
   !>  Pure Subroutine: SWAP_REALS
   !!
   !!  Purpose: Stupid little subroutine which swaps the values of its
   !!  inputs.
   !!
   !!  Comments: Inputs are of type real(c_real)
   !!
   !----------------------------------------------------------------------
    pure subroutine swap_reals(a, b)
        implicit none

        real(c_real), intent(inout) :: a, b
        real(c_real)                :: bucket

        bucket = a
        a      = b
        b      = bucket

    end subroutine swap_reals

   !----------------------------------------------------------------------
   !!
   !>  Pure Subroutine: LAMBDA_BOUNDS
   !!
   !!  Purpose: Given a line which passes through a box in three dimensions
   !!  (it can pass through the edges). Let lambda be a real value
   !!  representing the coordinate along the line. This subroutine finds
   !!  teh min/max values of lambda, in order for the point described by
   !!  lambda to be contained within the box.
   !!
   !!  Comments: Vectors are represented as one-dimensional arrays of type
   !!  real(c_real) and of dimension(3) (i.e. indices range from 1..3).
   !----------------------------------------------------------------------
    pure subroutine lambda_bounds(lambda_min, lambda_max, id_cell, p0, v, dx)
        implicit none

        real(c_real),               intent(  out) :: lambda_min, lambda_max
        integer,      dimension(3), intent(in   ) :: id_cell
        real(c_real), dimension(3), intent(in   ) :: p0, v, dx

        ! c... are the preliminary boundaries
        real(c_real) :: cx_lo, cy_lo, cz_lo, cx_hi, cy_hi, cz_hi

        ! defaults such that if skipped, min/max will not choose these values anyway
        cx_lo = -huge(cx_lo)
        cy_lo = -huge(cy_lo)
        cz_lo = -huge(cz_lo)

        cx_hi = huge(cx_hi)
        cy_hi = huge(cy_hi)
        cz_hi = huge(cz_hi)

        ! if the line runs parrallel to any of these dimensions (which is true for
        ! EB edges), then skip -> the min/max functions at the end will skip them
        ! due to the +/-huge(c...) defaults (above).

        if ( abs(v(1)) .gt. epsilon(v) ) then
            cx_lo = -( p0(1) - dble(id_cell(1)) * dx(1) ) / v(1)
            cx_hi = -( p0(1) - ( dble(id_cell(1)) + 1. ) * dx(1) ) / v(1)

            if ( v(1) .lt. 0. ) then
                call swap_reals(cx_lo, cx_hi)
            end if
        end if

        if ( abs(v(2)) .gt. epsilon(v) ) then
            cy_lo = -( p0(2) - dble(id_cell(2)) * dx(2) ) / v(2)
            cy_hi = -( p0(2) - ( dble(id_cell(2)) + 1. ) * dx(2) ) / v(2)

            if ( v(2) .lt. 0. ) then
                call swap_reals(cy_lo, cy_hi)
            end if
        end if

        if ( abs(v(3)) .gt. epsilon(v) )  then
            cz_lo = -( p0(3) - dble(id_cell(3)) * dx(3) ) / v(3)
            cz_hi = -( p0(3) - ( dble(id_cell(3)) + 1. ) * dx(3) ) / v(3)

            if ( v(3) .lt. 0. ) then
                call swap_reals(cz_lo, cz_hi)
            end if
        endif

        lambda_min = max(cx_lo, cy_lo, cz_lo)
        lambda_max = min(cx_hi, cy_hi, cz_hi)

    end subroutine lambda_bounds

   !----------------------------------------------------------------------
   !!
   !>  Pure Function: FACETS_NEAREST_PT
   !!
   !!  Purpose: Given a collision between particle and EB surface, and
   !!  given that a neighbour cell owns the EB surface, a collision between
   !!  the particle and the EDGE of the EB facet might occur. This
   !!  function returns the coordinates of the closest point on the edge of
   !!  an EB facet. This function does not check of collisions.
   !!
   !!  Comments: Position and normal vectors are represented as
   !!  one-dimensional arrays of type real(c_real) and of dimension(3)
   !!  (i.e. indices range from 1..3). Cells are enumerated using arrays of
   !!  type integer and of dimension(3).
   !----------------------------------------------------------------------
    pure function facets_nearest_pt(ind_pt, ind_loop, r_vec, eb_normal, eb_p0, dx)
        implicit none

        real(c_real), dimension(3)             :: facets_nearest_pt
        integer,      dimension(3), intent(in) :: ind_pt, ind_loop
        real(c_real), dimension(3), intent(in) :: r_vec, eb_normal, eb_p0, dx

        integer,      dimension(3) :: ind_facets
        integer                    :: n_facets, i_facet, tmp_facet, ind_cell, ind_nb
        real(c_real), dimension(3) :: c_vec, c_vec_tmp, rc_vec
        real(c_real), dimension(3) :: facet_normal, facet_p0, edge_p0, edge_v
        real(c_real)               :: min_dist, min_dist_tmp, eb_h, facet_h

        ! variables keeping track of coordinates on EB edges
        ! lambda_tmp: current lambda-value being used in testing for edge collions
        ! lambda: minimum (closets to bcentre) lambda value satisfying potential collision
        real(c_real) :: f_c, lambda_tmp, lambda_max, lambda_min


        ! Enumerate the possible EB facet edges invovlved.
        n_facets = 0

        if ( .not. (ind_pt(1) .eq. ind_loop(1)) ) then
            n_facets = n_facets + 1
            ind_facets(n_facets) = 1
        end if

        if ( .not. (ind_pt(2) .eq. ind_loop(2)) ) then
            n_facets = n_facets + 1
            ind_facets(n_facets) = 2
        end if

        if ( .not. (ind_pt(3) .eq. ind_loop(3)) ) then
            n_facets = n_facets + 1
            ind_facets(n_facets) = 3
        end if

        ! scalar characterizing EB facet position
        eb_h = dot_3d_real(eb_normal, eb_p0)

        ! itterate over EB facet edges and find whichever has the closest nearest point
        min_dist = huge(min_dist)
        do i_facet = 1, n_facets
            tmp_facet = ind_facets(i_facet)

            ! determine the normal of the cell's facet (cube faces)
            facet_normal = (/ 0., 0., 0. /)
            facet_normal(tmp_facet) = 1.  ! whether facing inwards or outwards is not important here

            ind_cell = ind_loop(tmp_facet)
            ind_nb = ind_pt(tmp_facet)

            ! determine position of the cell's facet
            if (ind_cell .lt. ind_nb) then
                f_c = ( dble(ind_cell) + 1.0 ) * dx(tmp_facet)
            else ! if (ind_cell .gt. ind_nb) then
                f_c = dble(ind_cell) * dx(tmp_facet)
            end if

            facet_p0 = (/                             &
                ( dble(ind_loop(1)) + 0.5 ) * dx(1) , &
                ( dble(ind_loop(2)) + 0.5 ) * dx(2) , &
                ( dble(ind_loop(3)) + 0.5 ) * dx(3)   &
            /)
            facet_p0(tmp_facet) = f_c

            ! scalar characterizing cell facet position
            facet_h = dot_3d_real(facet_normal, facet_p0)

            ! compute EB facet edge by finding the intercept between EB surface (first plane)
            ! and the cell's facet (second plane)
            call calc_facet_edge (edge_p0, edge_v, eb_h, facet_h, eb_normal, facet_normal)
            ! this solution is a line representing the closest EB edge, now compute the point
            ! on the line which minimizes the distance to the particle
            call lines_nearest_pt (lambda_tmp, c_vec_tmp, edge_p0, edge_v, r_vec)

            ! IMPORTANT: this point might be outside the cell
            !  -> in that case, it will be one of the cell's corners

            ! but don't use the PT_IN_BOX function, as it can yield false positives / negatives
            !  -> but if you do want to use it, include test [1] below to avoid rounding errors
            !if (.not. pt_in_box(c_vec_tmp, ind_loop, tmp_facet)) then

            ! if closest point is outside cell, determine the furthest we can go along the
            ! EB edge line whilst staying within the cell.
            call lambda_bounds(lambda_min, lambda_max, ind_loop, edge_p0, edge_v, dx)

            if (lambda_tmp .lt. lambda_min) then
                lambda_tmp = lambda_min
            elseif ( lambda_tmp .gt. lambda_max) then
                lambda_tmp = lambda_max
            end if
            c_vec_tmp(:) = edge_p0(:) + lambda_tmp*edge_v(:)

            !end if

            ! determine new distance to particle
            rc_vec(:) = c_vec_tmp(:) - r_vec(:)
            min_dist_tmp = dot_3d_real(rc_vec, rc_vec)

            ! minimize distance
            if (min_dist_tmp .lt. min_dist) then
                min_dist = min_dist_tmp
                c_vec(:) = c_vec_tmp(:)
            end if
        end do

        facets_nearest_pt(:) = c_vec(:)

    end function facets_nearest_pt

end module amrex_eb_geometry_module
