#include "AMReX_EB_geom.H"

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BoxArray.H>



namespace amrex {
namespace geom {


inline Real dot_3d_real (const RealVect & v1, const RealVect & v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}



inline RealVect cross_3d_real (const RealVect & v1, const RealVect & v2) {
    RealVect cross;

    cross[0] = v1[1]*v2[2] - v1[2]*v2[1];
    cross[1] = v1[2]*v2[0] - v1[0]*v2[2];
    cross[2] = v1[0]*v2[1] - v1[1]*v2[0];

    return cross;
}



void calc_facet_edge (      RealVect & p0,       RealVect & v,
                            Real       h1,       Real       h2,
                      const RealVect & n1, const RealVect & n2) {

    Real c_dp   = dot_3d_real(n1, n2);
    Real c_norm = 1 - c_dp * c_dp;

    Real c1 = ( h1 - h2 * c_dp ) / c_norm;
    Real c2 = ( h2 - h1 * c_dp ) / c_norm;

    for (int d=0; d<AMREX_SPACEDIM; ++d)
        p0[d] = c1 * n1[d] + c2 * n2[d];
    v = cross_3d_real(n1, n2);
}



void lines_nearest_pt (Real & lambda_min, RealVect & nearest_pt,
                       const RealVect & p0, const RealVect & v, const RealVect & pt) {

    RealVect c;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        c[d] = p0[d] - pt[d];

    lambda_min = - dot_3d_real(v, c) / dot_3d_real (v, v);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        nearest_pt[d] = p0[d] + lambda_min*v[d];
}



inline void swap_reals (Real & a, Real & b) {

    Real c = a;
    a      = b;
    b      = c;

}



void lambda_bounds (Real & lambda_min, Real & lambda_max, const IntVect & id_cell,
                    const RealVect & p0, const RealVect & v, const RealVect & dx) {


    Real cx_lo = -std::numeric_limits<Real>::max();
    Real cy_lo = -std::numeric_limits<Real>::max();
    Real cz_lo = -std::numeric_limits<Real>::max();

    Real cx_hi = std::numeric_limits<Real>::max();
    Real cy_hi = std::numeric_limits<Real>::max();
    Real cz_hi = std::numeric_limits<Real>::max();

    Real eps = std::numeric_limits<Real>::epsilon();

    if ( amrex::Math::abs(v[0]) > eps ) {
        cx_lo = -( p0[0] -   id_cell[0]       * dx[0] ) / v[0];
        cx_hi = -( p0[0] - ( id_cell[0] + 1 ) * dx[0] ) / v[0];

        if ( v[0] < 0. ) swap_reals(cx_lo, cx_hi);
    }

    if ( amrex::Math::abs(v[1]) > eps ) {
        cy_lo = -( p0[1] -   id_cell[1]       * dx[1] ) / v[1];
        cy_hi = -( p0[1] - ( id_cell[1] + 1 ) * dx[1] ) / v[1];

        if ( v[1] < 0. ) swap_reals(cy_lo, cy_hi);
    }

    if ( amrex::Math::abs(v[2]) > eps ) {
        cz_lo = -( p0[2] -   id_cell[2]       * dx[2] ) / v[2];
        cz_hi = -( p0[2] - ( id_cell[2] + 1 ) * dx[2] ) / v[2];

        if ( v[2] < 0. ) swap_reals(cz_lo, cz_hi);
    }


    lambda_min = std::max(cx_lo, std::max(cy_lo, cz_lo));
    lambda_max = std::min(cx_hi, std::min(cy_hi, cz_hi));
}



bool norm_all_eq (const RealVect & a, const RealVect & b) {

    Real diff = 0, eps = std::numeric_limits<Real>::epsilon();
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        Real abs_a = amrex::Math::abs(a[d]);
        Real abs_b = amrex::Math::abs(b[d]);
        diff += amrex::Math::abs(abs_a - abs_b);
    }

    return diff <= AMREX_SPACEDIM*eps;
}



bool all_eq (const IntVect & a, const IntVect & b) {

    for (int d=0; d<AMREX_SPACEDIM; ++d)
        if (a[d] != b[d])
            return false;

    return true;
}



RealVect facets_nearest_pt (const IntVect  & ind_pt, const IntVect  & ind_loop,
                            const RealVect & r_vec,  const RealVect & eb_normal,
                            const RealVect & eb_p0,  const RealVect & dx ) {

    int n_facets = 0;
    IntVect ind_facets {AMREX_D_DECL(0, 0, 0)};

    if ( ! (ind_pt[0] == ind_loop[0]) ) {
        n_facets = n_facets + 1;
        ind_facets[n_facets - 1] = 1;
    }

    if ( ! (ind_pt[1] == ind_loop[1]) ) {
        n_facets = n_facets + 1;
        ind_facets[n_facets - 1] = 2;
    }

    if ( ! (ind_pt[2] == ind_loop[2]) ) {
        n_facets = n_facets + 1;
        ind_facets[n_facets - 1] = 3;
    }


    Real eb_h = dot_3d_real(eb_normal, eb_p0);


    Real min_dist = std::numeric_limits<Real>::max();
    RealVect c_vec;
    for (int i_facet=0; i_facet<n_facets; ++i_facet) {

        int tmp_facet = ind_facets[i_facet];

        RealVect facet_normal {AMREX_D_DECL(0., 0., 0.)};
        facet_normal[tmp_facet - 1] = 1.;


        if (norm_all_eq(eb_normal, facet_normal)) continue;

        int ind_cell = ind_loop[tmp_facet - 1];
        int ind_nb   = ind_pt[tmp_facet - 1];

        Real f_c;
        if (ind_cell < ind_nb) {
            f_c = ( ind_cell + 1 ) * dx[tmp_facet - 1];
        } else {
            f_c =   ind_cell       * dx[tmp_facet - 1];
        }

        RealVect facet_p0;
        for (int d=0; d<AMREX_SPACEDIM; ++d)
            facet_p0[d] = (ind_loop[d] + 0.5) * dx[d];

        facet_p0[tmp_facet - 1] = f_c;

        Real facet_h = dot_3d_real(facet_normal, facet_p0);

        RealVect edge_p0, edge_v;
        calc_facet_edge(edge_p0, edge_v, eb_h, facet_h, eb_normal, facet_normal);

        Real lambda_tmp;
        RealVect c_vec_tmp;
        lines_nearest_pt(lambda_tmp, c_vec_tmp, edge_p0, edge_v, r_vec);

        Real lambda_min, lambda_max;
        lambda_bounds(lambda_min, lambda_max, ind_loop, edge_p0, edge_v, dx);


        if (lambda_tmp < lambda_min) {
            lambda_tmp = lambda_min;
        } else if ( lambda_tmp > lambda_max) {
            lambda_tmp = lambda_max;
        }

        RealVect rc_vec;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            c_vec_tmp[d] = edge_p0[d] + lambda_tmp*edge_v[d];
            rc_vec[d] = c_vec_tmp[d] - r_vec[d];
        }

        Real min_dist_tmp = dot_3d_real(rc_vec, rc_vec);

        if (min_dist_tmp < min_dist) {
            min_dist = min_dist_tmp;
            for (int d=0; d<AMREX_SPACEDIM; ++d)
                c_vec[d] = c_vec_tmp[d];
        }
    }

    return c_vec;
}



void closest_dist (Real & min_dist, bool & proj_valid,
                   const Vector<Real> & eb_data,
                   const RealVect & dx_eb, const RealVect & pos) {

    int l_eb = eb_data.size();

    RealVect inv_dx;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        inv_dx[d] = 1. / dx_eb[d];

    min_dist   = std::numeric_limits<Real>::max();
    proj_valid = false;

    Real min_dist2 = std::numeric_limits<Real>::max();
    int i_nearest  = 0;


    for (int i=0; i<l_eb; i+=2*AMREX_SPACEDIM) {

        RealVect eb_cent, eb_norm;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            eb_cent[d] = eb_data[i + d];
            eb_norm[d] = eb_data[i + d + AMREX_SPACEDIM];
        }

        Real dist2 = dot_3d_real( pos - eb_cent, pos - eb_cent );

        if ( dist2 < min_dist2 ) {
            min_dist2 = dist2;
            i_nearest = i;
        }
    }

    RealVect eb_cent, eb_norm;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        eb_cent[d] = eb_data[i_nearest + d];
        eb_norm[d] = eb_data[i_nearest + d + AMREX_SPACEDIM];
    }

    Real dist_proj = dot_3d_real( pos - eb_cent, -eb_norm );

    RealVect eb_min_pt;
    IntVect  vi_cent, vi_pt;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        eb_min_pt[d] = pos[d] + eb_norm[d] * dist_proj;

        vi_cent[d] = amrex::Math::floor( eb_cent[d] * inv_dx[d]);
        vi_pt[d]   = amrex::Math::floor( eb_min_pt[d] * inv_dx[d]);
    }


    bool min_pt_valid = false;
    if ( all_eq( vi_pt, vi_cent ) ) {
        min_pt_valid = true;
    } else {
        for (int k_shift = -1; k_shift <= 1; ++k_shift) {
            for (int j_shift = -1; j_shift <= 1; ++j_shift) {
                for (int i_shift = -1; i_shift <= 1; ++i_shift) {

                    RealVect shift {AMREX_D_DECL(1.e-6*i_shift,
                                                 1.e-6*j_shift,
                                                 1.e-6*k_shift)};

                    for (int d=0; d<AMREX_SPACEDIM; ++d)
                        vi_pt[d] = amrex::Math::floor(
                                (eb_min_pt[d] + shift[d]*dx_eb[d]) * inv_dx[d]
                            );

                    if ( all_eq( vi_pt, vi_cent ) ) min_pt_valid = true;
                }
            }
        }
    }


    if ( min_pt_valid ) {
        min_dist   = dist_proj;
        proj_valid = true;
    } else {
        RealVect c_vec = facets_nearest_pt(vi_pt, vi_cent, pos, eb_norm, eb_cent, dx_eb);
        Real min_edge_dist2 = dot_3d_real( c_vec - pos, c_vec - pos);
        min_dist = -std::sqrt( std::min(min_dist2, min_edge_dist2) );
    }

}


}
}
