/**
 * \file AMReX_EB_StateRedistItracker.cpp
 * @{
 *
 */

#include <AMReX_EB_Redistribution.H>
#include <AMReX_EB_StateRedistTracker_K.H>
#include <AMReX_GpuMemory.H>

namespace amrex {

using state_redist_detail::imap;
using state_redist_detail::jmap;
#if (AMREX_SPACEDIM == 3)
using state_redist_detail::kmap;
#endif

namespace {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
initialize_state_data (int i, int j, int k,
                       Box const& bxg3, Box const& bxg4,
                       Array4<EBCellFlag const> const& flag,
                       Array4<Real const> const& ccent,
                       Array4<int> const& nrs,
                       Array4<Real> const& alpha,
                       Array4<Real> const& cent_hat) noexcept
{
    nrs(i,j,k) = 1;

    if (bxg4.contains(i,j,k)) {
        auto const cell_flag = flag(i,j,k);
        alpha(i,j,k,0) = cell_flag.isCovered() ? 0.0_rt : 1.0_rt;

        if (bxg3.contains(i,j,k)) {
            if (cell_flag.isRegular()) {
                AMREX_D_TERM(cent_hat(i,j,k,0) = 0.0_rt;,
                             cent_hat(i,j,k,1) = 0.0_rt;,
                             cent_hat(i,j,k,2) = 0.0_rt;);
            } else if (cell_flag.isCovered()) {
                AMREX_D_TERM(cent_hat(i,j,k,0) = eb_covered_val;,
                             cent_hat(i,j,k,1) = eb_covered_val;,
                             cent_hat(i,j,k,2) = eb_covered_val;);
            } else {
                AMREX_D_TERM(cent_hat(i,j,k,0) = ccent(i,j,k,0);,
                             cent_hat(i,j,k,1) = ccent(i,j,k,1);,
                             cent_hat(i,j,k,2) = ccent(i,j,k,2););
            }
        }
    }
}

int
find_active_cells (Box const& bx, Array4<Real const> const& vfrac,
                   Dim3* active_cells, Real target_volfrac)
{
    Gpu::DeviceScalar<int> active_count(0);
    int* const pactive_count = active_count.dataPtr();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real const vf = vfrac(i,j,k);
        if (vf > 0.0 && vf < target_volfrac) {
            int const active_index = Gpu::Atomic::Add(pactive_count, 1);
            active_cells[active_index] = Dim3{.x = i, .y = j, .z = k};
        }
    });

    return active_count.dataValue();
}

#if (AMREX_SPACEDIM == 2)

// Neighbor ID ordering; component 0 stores the neighbor count.
//
// ^  6 7 8
// |  4   5
// j  1 2 3
//   i --->

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
construct_itracker_cell (int i, int j, int k,
                         Array4<Real const> const& apx,
                         Array4<Real const> const& apy,
                         Array4<Real const> const& vfrac,
                         Array4<int> const& itracker,
                         Box const& domain,
                         bool is_periodic_x, bool is_periodic_y,
                         Real target_volfrac) noexcept
{
#ifdef AMREX_USE_FLOAT
    const Real small_norm_diff = Real(1.e-4);
#else
    const Real small_norm_diff = Real(1.e-8);
#endif

    Real apnorm, apnorm_inv;
    const Real dapx = apx(i+1,j  ,k  ) - apx(i,j,k);
    const Real dapy = apy(i  ,j+1,k  ) - apy(i,j,k);
    apnorm = std::sqrt(dapx*dapx+dapy*dapy);
    apnorm_inv = Real(1)/apnorm;
    Real nx = dapx * apnorm_inv;
    Real ny = dapy * apnorm_inv;

    bool nx_eq_ny = ( (std::abs(nx-ny) < small_norm_diff) ||
                      (std::abs(nx+ny) < small_norm_diff)  ) ? true : false;

    int n1;
    // As a first pass, choose just based on the normal
    if (std::abs(nx) > std::abs(ny))
    {
        if (nx > 0) {
            n1 = 5;
        } else {
            n1 = 4;
        }

    } else {
        if (ny > 0) {
            n1 = 7;
        } else {
            n1 = 2;
        }
    }

    bool xdir_mns_ok = (is_periodic_x || (i > domain.smallEnd(0)));
    bool xdir_pls_ok = (is_periodic_x || (i < domain.bigEnd(0)  ));
    bool ydir_mns_ok = (is_periodic_y || (j > domain.smallEnd(1)));
    bool ydir_pls_ok = (is_periodic_y || (j < domain.bigEnd(1)  ));

    // Override above logic if trying to reach outside a domain boundary (and non-periodic).
    // Note that the replacement must stay inside the domain as well -- we prefer the
    // direction the normal points in, but reverse the sign if that side is outside.
    // Otherwise, at a domain corner, the second override could undo the first one and
    // leave us pointing outside the domain again.
    if ( (!xdir_mns_ok && (n1 == 4)) ||
         (!xdir_pls_ok && (n1 == 5)) )
    {
        n1 = ((ny > 0 && ydir_pls_ok) || !ydir_mns_ok) ? 7 : 2;
    }
    if ( (!ydir_mns_ok && (n1 == 2)) ||
         (!ydir_pls_ok && (n1 == 7)) )
    {
        n1 = ((nx > 0 && xdir_pls_ok) || !xdir_mns_ok) ? 5 : 4;
    }

    // (i,j) merges with at least one cell now
    int num_neighbors = 1;

    // (i+ioff,j+joff) is in the nbhd of (i,j)
    int ioff = imap[n1];
    int joff = jmap[n1];

    // Variables for vfrac fallback
    bool select_by_vfrac = false;
    Real second_vfrac = Real(-1.0);
    int second_idx = -1;

    // Sanity check - if selected neighbor is covered, use vfrac fallback
    Real first_vfrac = vfrac(i+ioff,j+joff,k);
    if (first_vfrac == 0.)
    {
        Real vfrac_x = (xdir_pls_ok) ? vfrac(i+1,j,k) : Real(-1.0);
        int idx_x = 5;
        if (xdir_mns_ok && vfrac(i-1,j,k) > vfrac_x) {
            vfrac_x = vfrac(i-1,j,k);
            idx_x = 4;
        }

        Real vfrac_y = (ydir_pls_ok) ? vfrac(i,j+1,k) : Real(-1.0);
        int idx_y = 7;
        if (ydir_mns_ok && vfrac(i,j-1,k) > vfrac_y) {
            vfrac_y = vfrac(i,j-1,k);
            idx_y = 2;
        }

        // Select largest and second-largest vfrac
        int first_idx;
        if (vfrac_x >= vfrac_y) {
            first_idx = idx_x;
            first_vfrac = vfrac_x;
            second_idx = idx_y;
            second_vfrac = vfrac_y;
        } else {
            first_idx = idx_y;
            first_vfrac = vfrac_y;
            second_idx = idx_x;
            second_vfrac = vfrac_x;
        }

        if (first_vfrac <= 0.0) {
            amrex::Abort("Trying to merge with covered cell - all cardinal neighbors are also covered or outside domain");
        }

        n1 = first_idx;

        select_by_vfrac = true;
    }

    itracker(i,j,k,1) = n1;
    Real const vf = vfrac(i,j,k);
    Real sum_vol = vf + first_vfrac;

#if 0
    if (debug_verbose > 0)
        amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vf <<
                          " trying to merge with " << IntVect(i+imap[n1],j+jmap[n1]) <<
                          " with volfrac " << vfrac(i+imap[n1],j+jmap[n1],k) <<
                          " to get new sum_vol " <<  sum_vol << '\n';
#endif

    int n2 = 0;
    // If the merged cell isn't large enough, we try to merge in the other direction
    if (sum_vol < target_volfrac || nx_eq_ny)
    {
        // If first neighbor was selected by vfrac, use second-best vfrac candidate for second neighbor
        if (select_by_vfrac && second_vfrac > 0.0)
        {
            n2 = second_idx;
            itracker(i,j,k,2) = n2;
            num_neighbors += 1;
            sum_vol += second_vfrac;
        }
        else if (!select_by_vfrac)
        {
            // Original normal-based perpendicular selection logic
            // Original offset was in y-direction, so we will add to the x-direction
            // Note that if we can't because it would go outside the domain, we don't
            if (ioff == 0) {
                if (nx >= 0 && xdir_pls_ok)
                {
                    n2 = 5;
                    num_neighbors += 1;
                }
                else if (nx <= 0 && xdir_mns_ok)
                {
                    n2 = 4;
                    num_neighbors += 1;
                }

            // Original offset was in x-direction, so we will add to the y-direction
            // Note that if we can't because it would go outside the domain, we don't
            } else {
                if (ny >= 0 && ydir_pls_ok)
                {
                    n2 = 7;
                    num_neighbors += 1;
                }
                else if (ny <= 0 && ydir_mns_ok)
                {
                    n2 = 2;
                    num_neighbors += 1;
                }
            }

            if (num_neighbors > 1)
            {
                // (i+ioff2,j+joff2) is in the nbhd of (i,j)
                int ioff2 = imap[n2];
                int joff2 = jmap[n2];

                sum_vol += vfrac(i+ioff2,j+joff2,k);
                itracker(i,j,k,2) = n2;
#if 0
                if (debug_verbose > 0)
                    amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vf <<
                                      " trying to ALSO merge with " << IntVect(i+ioff2,j+joff2) <<
                                      " with volfrac " << vfrac(i+ioff2,j+joff2,k) <<
                                       " to get new sum_vol " <<  sum_vol << '\n';
#endif
            }
        }  // end else if (!select_by_vfrac)
    }

    // Now we merge in the corner direction if we have already claimed two
    if (num_neighbors == 2)
    {
        int n3;
        // We already have two offsets, and we know they are in different directions
        ioff = imap[n1] + imap[n2];
        joff = jmap[n1] + jmap[n2];

        if (ioff > 0 && joff > 0) {
            n3 = 8;
        } else if (ioff < 0 && joff > 0) {
            n3 = 6;
        } else if (ioff > 0 && joff < 0) {
            n3 = 3;
        } else {
            n3 = 1;
        }

        itracker(i,j,k,3) = n3;

        // (i,j) merges with at least three cells now
        num_neighbors += 1;

        sum_vol += vfrac(i+ioff,j+joff,k);
#if 0
        if (debug_verbose > 0)
            amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vf <<
                              " trying to ALSO merge with " << IntVect(i+ioff,j+joff) <<
                              " with volfrac " << vfrac(i+ioff,j+joff,k) <<
                              " to get new sum_vol " <<  sum_vol << '\n';
#endif
    }
    // A cell next to a non-periodic domain boundary can only merge with the cells
    // that are available inside the domain, so its neighborhood may fall short of
    // target_volfrac even though we merged with everything we were allowed to use.
    // We accept the largest neighborhood we could build in that case: reaching
    // outside the domain instead would silently break conservation, because
    // StateRedistribute leaves cells outside the domain out of Qhat even though
    // MakeStateRedistUtils counts their volume in nbhd_vol.
    bool nbhd_limited_by_domain = !(xdir_mns_ok && xdir_pls_ok &&
                                    ydir_mns_ok && ydir_pls_ok);

    if (sum_vol < target_volfrac && !nbhd_limited_by_domain)
    {
#if 0
      amrex::Print() << "Couldn't merge with enough cells to raise volume at " <<
                        IntVect(i,j) << " so stuck with sum_vol " << sum_vol << '\n';
#endif
      amrex::Abort("Couldn't merge with enough cells to raise volume greater than target_volfrac");
    }
    itracker(i,j,k,0) = num_neighbors;
}

#elif (AMREX_SPACEDIM == 3)

// Neighbor ID ordering; component 0 stores the neighbor count.
//
//    at k-1   |   at k   |   at k+1
//
// ^  15 16 17 |  6 7 8  |  24 25 26
// |  12 13 14 |  4   5  |  21 22 23
// j  9  10 11 |  1 2 3  |  18 19 20
//   i --->

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
construct_itracker_cell (int i, int j, int k,
                         Array4<Real const> const& apx,
                         Array4<Real const> const& apy,
                         Array4<Real const> const& apz,
                         Array4<Real const> const& vfrac,
                         Array4<int> const& itracker,
                         Box const& domain,
                         bool is_periodic_x, bool is_periodic_y,
                         bool is_periodic_z,
                         Real target_volfrac) noexcept
{
#ifdef AMREX_USE_FLOAT
    const Real small_norm_diff = Real(1.e-4);
#else
    const Real small_norm_diff = Real(1.e-8);
#endif

    Real apnorm, apnorm_inv;
    const Real dapx = apx(i+1,j  ,k  ) - apx(i,j,k);
    const Real dapy = apy(i  ,j+1,k  ) - apy(i,j,k);
    const Real dapz = apz(i  ,j  ,k+1) - apz(i,j,k);
    apnorm = std::sqrt(dapx*dapx+dapy*dapy+dapz*dapz);
    apnorm_inv = Real(1.0)/apnorm;
    Real nx = dapx * apnorm_inv;
    Real ny = dapy * apnorm_inv;
    Real nz = dapz * apnorm_inv;

    bool nx_eq_ny = ( (std::abs(nx-ny) < small_norm_diff) ||
                      (std::abs(nx+ny) < small_norm_diff)  ) ? true : false;
    bool nx_eq_nz = ( (std::abs(nx-nz) < small_norm_diff) ||
                      (std::abs(nx+nz) < small_norm_diff)  ) ? true : false;
    bool ny_eq_nz = ( (std::abs(ny-nz) < small_norm_diff) ||
                      (std::abs(ny+nz) < small_norm_diff)  ) ? true : false;

    bool xdir_mns_ok = (is_periodic_x || (i > domain.smallEnd(0)));
    bool xdir_pls_ok = (is_periodic_x || (i < domain.bigEnd(0)  ));
    bool ydir_mns_ok = (is_periodic_y || (j > domain.smallEnd(1)));
    bool ydir_pls_ok = (is_periodic_y || (j < domain.bigEnd(1)  ));
    bool zdir_mns_ok = (is_periodic_z || (k > domain.smallEnd(2)));
    bool zdir_pls_ok = (is_periodic_z || (k < domain.bigEnd(2)  ));

    int n1;
    // x-component of normal is greatest
    if ( (std::abs(nx) > std::abs(ny)) &&
         (std::abs(nx) > std::abs(nz)) )
    {
        if (nx > 0) {
            n1 = 5;
        } else {
            n1 = 4;
        }

    // y-component of normal is greatest
    } else if ( (std::abs(ny) >= std::abs(nx)) &&
                (std::abs(ny) > std::abs(nz)) )
    {
        if (ny > 0) {
            n1 = 7;
        } else {
            n1 = 2;
        }

    // z-component of normal is greatest
    } else {
        if (nz > 0) {
            n1 = 22;
        } else {
            n1 = 13;
        }
    }

    // Is the direction the normal points in open in each coordinate direction?
    bool xdir_nrm_ok = (nx > 0) ? xdir_pls_ok : xdir_mns_ok;
    bool ydir_nrm_ok = (ny > 0) ? ydir_pls_ok : ydir_mns_ok;
    bool zdir_nrm_ok = (nz > 0) ? zdir_pls_ok : zdir_mns_ok;

    // Override above logic if trying to reach outside a domain boundary (and non-periodic).
    // Note that the replacement must stay inside the domain as well: we prefer the other
    // coordinate direction along which the normal points inward, and only reverse a sign
    // if neither of the two remaining directions is open on the side the normal points to.
    // Otherwise, at a domain corner or edge, a later override could undo an earlier one and
    // leave us pointing outside the domain again.
    if ( (!xdir_mns_ok && (n1 == 4)) ||
         (!xdir_pls_ok && (n1 == 5)) )
    {
        if ( (std::abs(ny) > std::abs(nz) && ydir_nrm_ok) || !zdir_nrm_ok ) {
            n1 = ((ny > 0 && ydir_pls_ok) || !ydir_mns_ok) ? 7 : 2;
        } else {
            n1 = (nz > 0) ? 22 : 13;
        }
    }

    if ( (!ydir_mns_ok && (n1 == 2)) ||
         (!ydir_pls_ok && (n1 == 7)) )
    {
        if ( (std::abs(nx) > std::abs(nz) && xdir_nrm_ok) || !zdir_nrm_ok ) {
            n1 = ((nx > 0 && xdir_pls_ok) || !xdir_mns_ok) ? 5 : 4;
        } else {
            n1 = (nz > 0) ? 22 : 13;
        }
    }

    if ( (!zdir_mns_ok && (n1 == 13)) ||
         (!zdir_pls_ok && (n1 == 22)) )
    {
        if ( (std::abs(nx) > std::abs(ny) && xdir_nrm_ok) || !ydir_nrm_ok ) {
            n1 = ((nx > 0 && xdir_pls_ok) || !xdir_mns_ok) ? 5 : 4;
        } else {
            n1 = (ny > 0) ? 7 : 2;
        }
    }

    // (i,j,k) merges with at least one cell now
    int num_neighbors = 1;

    // (i+ioff,j+joff,k+koff) is now the first cell in the nbhd of (i,j,k)
    int ioff = imap[n1];
    int joff = jmap[n1];
    int koff = kmap[n1];

    // Variables for vfrac fallback
    bool select_by_vfrac = false;
    Real second_vfrac = Real(-1.0);
    int second_idx = -1;

    // Sanity check - if selected neighbor is covered, use vfrac fallback
    Real first_vfrac = vfrac(i+ioff,j+joff,k+koff);
    if (first_vfrac == 0.)
    {
        Real vfrac_x = (xdir_pls_ok) ? vfrac(i+1,j,k) : Real(-1.0);
        int idx_x = 5;
        if (xdir_mns_ok && vfrac(i-1,j,k) > vfrac_x) {  // Strict > means +x wins ties
            vfrac_x = vfrac(i-1,j,k);
            idx_x = 4;
        }

        Real vfrac_y = (ydir_pls_ok) ? vfrac(i,j+1,k) : Real(-1.0);
        int idx_y = 7;
        if (ydir_mns_ok && vfrac(i,j-1,k) > vfrac_y) {
            vfrac_y = vfrac(i,j-1,k);
            idx_y = 2;
        }

        Real vfrac_z = (zdir_pls_ok) ? vfrac(i,j,k+1) : Real(-1.0);
        int idx_z = 22;
        if (zdir_mns_ok && vfrac(i,j,k-1) > vfrac_z) {
            vfrac_z = vfrac(i,j,k-1);
            idx_z = 13;
        }

        // Find largest and second-largest vfrac
        int first_idx;

        if (vfrac_x >= vfrac_y && vfrac_x >= vfrac_z) {
            // X is largest
            first_idx = idx_x;
            first_vfrac = vfrac_x;
            if (vfrac_y >= vfrac_z) {
                second_idx = idx_y;
                second_vfrac = vfrac_y;
            } else {
                second_idx = idx_z;
                second_vfrac = vfrac_z;
            }
        } else if (vfrac_y >= vfrac_z) {
            // Y is largest
            first_idx = idx_y;
            first_vfrac = vfrac_y;
            if (vfrac_x >= vfrac_z) {
                second_idx = idx_x;
                second_vfrac = vfrac_x;
            } else {
                second_idx = idx_z;
                second_vfrac = vfrac_z;
            }
        } else {
            // Z is largest
            first_idx = idx_z;
            first_vfrac = vfrac_z;
            if (vfrac_x >= vfrac_y) {
                second_idx = idx_x;
                second_vfrac = vfrac_x;
            } else {
                second_idx = idx_y;
                second_vfrac = vfrac_y;
            }
        }

        if (first_vfrac <= 0.0) {
            amrex::Abort("Trying to merge with covered cell - all cardinal neighbors are also covered or outside domain");
        }

        n1 = first_idx;

        ioff = imap[first_idx];
        joff = jmap[first_idx];
        koff = kmap[first_idx];

        select_by_vfrac = true;
    }

    itracker(i,j,k,1) = n1;
    Real const vf = vfrac(i,j,k);
    Real sum_vol = vf + first_vfrac;

#if 0
    if (debug_print)
        amrex::Print() << "Cell " << IntVect(i,j,k) << " with volfrac " << vf <<
                          " trying to merge with " << IntVect(i+ioff,j+joff,k+koff) <<
                          " with volfrac " << vfrac(i+ioff,j+joff,k+koff) <<
                          " to get new sum_vol " <<  sum_vol << '\n';
#endif

    bool just_broke_symmetry = ( ( (joff == 0 && koff == 0) && (nx_eq_ny || nx_eq_nz) ) ||
                                 ( (ioff == 0 && koff == 0) && (nx_eq_ny || ny_eq_nz) ) ||
                                 ( (ioff == 0 && joff == 0) && (nx_eq_nz || ny_eq_nz) ) );

    int n2 = 0;
    // If the merged cell isn't large enough, or if we broke symmetry by the current merge,
    // we merge in one of the other directions.  Note that the direction of the next merge
    // is first set by a symmetry break, but if that isn't happening, we choose the next largest normal

    if ( (sum_vol < target_volfrac) || just_broke_symmetry )
    {
        // If first neighbor was selected by vfrac, use second-best vfrac candidate for second neighbor
        if (select_by_vfrac && second_vfrac > 0.0)
        {
            n2 = second_idx;
            itracker(i,j,k,2) = n2;
            num_neighbors += 1;

            sum_vol += second_vfrac;

        }
        else if (!select_by_vfrac)
        {
            // Original normal-based perpendicular selection logic
            // Original offset was in x-direction
            if (joff == 0 && koff == 0)
            {
                if (nx_eq_ny) {
                    n2 = (ny > 0) ? 7 : 2;
                } else if (nx_eq_nz) {
                    n2 = (nz > 0) ? 22 : 13;
                } else if ( (std::abs(ny) > std::abs(nz)) ) {
                    n2 = (ny > 0) ? 7 : 2;
                } else {
                    n2 = (nz > 0) ? 22 : 13;
                }

            // Original offset was in y-direction
            } else if (ioff == 0 && koff == 0)
            {
                if (nx_eq_ny) {
                    n2 = (nx > 0) ? 5 : 4;
                } else if (ny_eq_nz) {
                    n2 = (nz > 0) ? 22 : 13;
                } else if ( (std::abs(nx) > std::abs(nz)) ) {
                    n2 = (nx > 0) ? 5 : 4;
                } else {
                    n2 = (nz > 0) ? 22 : 13;
                }

            // Original offset was in z-direction
            } else if (ioff == 0 && joff == 0)
            {
                if (nx_eq_nz) {
                    n2 = (nx > 0) ? 5 : 4;
                } else if (ny_eq_nz) {
                    n2 = (ny > 0) ? 7 : 2;
                } else if ( (std::abs(nx) > std::abs(ny)) ) {
                    n2 = (nx > 0) ? 5 : 4;
                } else {
                    n2 = (ny > 0) ? 7 : 2;
                }
            }

            // (i+ioff2,j+joff2,k+koff2) would be in the nbhd of (i,j,k)
            int ioff2 = imap[n2];
            int joff2 = jmap[n2];
            int koff2 = kmap[n2];

            // As in 2D, we don't reach outside a domain boundary (and non-periodic).
            // The first merge used one coordinate direction and this one wants another,
            // so if the chosen direction is outside the domain we try the one remaining
            // direction instead, and only skip this merge if that is blocked as well.
            bool nbor2_ok = (ioff2 >= 0 || xdir_mns_ok) && (ioff2 <= 0 || xdir_pls_ok) &&
                            (joff2 >= 0 || ydir_mns_ok) && (joff2 <= 0 || ydir_pls_ok) &&
                            (koff2 >= 0 || zdir_mns_ok) && (koff2 <= 0 || zdir_pls_ok);

            if (!nbor2_ok)
            {
                if (ioff == 0 && ioff2 == 0) {
                    n2 = ((nx > 0 && xdir_pls_ok) || !xdir_mns_ok) ? 5 : 4;
                    nbor2_ok = (n2 == 5) ? xdir_pls_ok : xdir_mns_ok;
                } else if (joff == 0 && joff2 == 0) {
                    n2 = ((ny > 0 && ydir_pls_ok) || !ydir_mns_ok) ? 7 : 2;
                    nbor2_ok = (n2 == 7) ? ydir_pls_ok : ydir_mns_ok;
                } else {
                    n2 = ((nz > 0 && zdir_pls_ok) || !zdir_mns_ok) ? 22 : 13;
                    nbor2_ok = (n2 == 22) ? zdir_pls_ok : zdir_mns_ok;
                }

                ioff2 = imap[n2];
                joff2 = jmap[n2];
                koff2 = kmap[n2];
            }

            if (!nbor2_ok)
            {
                // Leave the neighbor count alone so that nothing downstream -- including
                // the corner merge and the 2x2x2 expansion below -- sees this neighbor
                n2 = 0;
            }
            else
            {
                // (i,j,k) merges with at least two cells now
                num_neighbors += 1;

                sum_vol += vfrac(i+ioff2,j+joff2,k+koff2);
            }
            itracker(i,j,k,2) = n2;
#if 0
            if (debug_print)
                amrex::Print() << "Cell " << IntVect(i,j,k) << " with volfrac " << vf <<
                                  " trying to ALSO merge with " << IntVect(i+ioff2,j+joff2,k+koff2) <<
                                  " with volfrac " << vfrac(i+ioff2,j+joff2,k+koff2) <<
                                   " to get new sum_vol " <<  sum_vol << '\n';
#endif
        }  // end else if (!select_by_vfrac)
    }

    // If the merged cell has merged in two directions, we now merge in the corner direction within the current plane
    if (num_neighbors >= 2)
    {
        int n3;
        // We already have two offsets, and we know they are in different directions
        ioff = imap[n1] + imap[n2];
        joff = jmap[n1] + jmap[n2];
        koff = kmap[n1] + kmap[n2];

        // Both nbors are in the koff=0 plane
        if (koff == 0)
        {
            if (ioff > 0 && joff > 0) {
                n3 = 8;
            } else if (ioff < 0 && joff > 0) {
                n3 = 6;
            } else if (ioff > 0 && joff < 0) {
                n3 = 3;
            } else {
                n3 = 1;
            }

        // Both nbors are in the joff=0 plane
        } else if (joff == 0) {
            if (ioff > 0 && koff > 0) {
                n3 = 23;
            } else if (ioff < 0 && koff > 0) {
                n3 = 21;
            } else if (ioff > 0 && koff < 0) {
                n3 = 14;
            } else {
                n3 = 12;
            }

        // Both nbors are in the ioff=0 plane
        } else {
            if (joff > 0 && koff > 0) {
                n3 = 25;
            } else if (joff < 0 && koff > 0) {
                n3 = 19;
            } else if (joff > 0 && koff < 0) {
                n3 = 16;
            } else {
                n3 = 10;
            }
        }

        itracker(i,j,k,3) = n3;

        // (i,j,k) merges with at least three cells now
        num_neighbors += 1;

        sum_vol += vfrac(i+ioff,j+joff,k+koff);

#if 0
        int ioff3 = imap[n3];
        int joff3 = jmap[n3];
        int koff3 = kmap[n3];
        if (debug_print)
            amrex::Print() << "Cell " << IntVect(i,j,k) << " with volfrac " << vf <<
                              " trying to ALSO merge with " << IntVect(i+ioff3,j+joff3,k+koff3) <<
                              " with volfrac " << vfrac(i+ioff3,j+joff3,k+koff3) << '\n';
#endif

        // All nbors are currently in one of three planes
        just_broke_symmetry = ( ( (koff == 0) && (nx_eq_nz || ny_eq_nz) ) ||
                                ( (joff == 0) && (nx_eq_ny || ny_eq_nz) ) ||
                                ( (ioff == 0) && (nx_eq_ny || nx_eq_nz) ) );

        // If with a nbhd of four cells we have still not reached vfrac > target_volfrac, we add another four
        //    cells to the nbhd to make a 2x2x2 block.  We use the direction of the remaining
        //    normal to know whether to go lo or hi in the new direction.
        if (sum_vol < target_volfrac || just_broke_symmetry)
        {
#if 0
            if (debug_print)
                if (just_broke_symmetry)
                    amrex::Print() << "Expanding neighborhood of " << IntVect(i,j,k) <<
                                      " from 4 to 8 since we just broke symmetry with the last merge " << '\n';
                else
                    amrex::Print() << "Expanding neighborhood of " << IntVect(i,j,k) <<
                                      " from 4 to 8 since sum_vol with 4 was only " << sum_vol << " " << '\n';
#endif

            for (int n = 4; n < 8; ++n) {
                int nbor;
                if (koff == 0) {
                    bool const expand_plus = (nz > 0 && zdir_pls_ok) || !zdir_mns_ok;
                    if (n == 4) {
                        nbor = expand_plus ? 22 : 13;
                    } else if (n == 5) {
                        nbor = expand_plus ? ((ioff > 0) ? 23 : 21)
                                           : ((ioff > 0) ? 14 : 12);
                    } else if (n == 6) {
                        nbor = expand_plus ? ((joff > 0) ? 25 : 19)
                                           : ((joff > 0) ? 16 : 10);
                    } else {
                        if (ioff > 0 && joff > 0) {
                            nbor = expand_plus ? 26 : 17;
                        } else if (ioff < 0 && joff > 0) {
                            nbor = expand_plus ? 24 : 15;
                        } else if (ioff > 0 && joff < 0) {
                            nbor = expand_plus ? 20 : 11;
                        } else {
                            nbor = expand_plus ? 18 : 9;
                        }
                    }
                } else if (joff == 0) {
                    bool const expand_plus = (ny > 0 && ydir_pls_ok) || !ydir_mns_ok;
                    if (n == 4) {
                        nbor = expand_plus ? 7 : 2;
                    } else if (n == 5) {
                        nbor = expand_plus ? ((ioff > 0) ? 8 : 6)
                                           : ((ioff > 0) ? 3 : 1);
                    } else if (n == 6) {
                        nbor = expand_plus ? ((koff > 0) ? 25 : 16)
                                           : ((koff > 0) ? 19 : 10);
                    } else {
                        if (ioff > 0 && koff > 0) {
                            nbor = expand_plus ? 26 : 20;
                        } else if (ioff < 0 && koff > 0) {
                            nbor = expand_plus ? 24 : 18;
                        } else if (ioff > 0 && koff < 0) {
                            nbor = expand_plus ? 17 : 11;
                        } else {
                            nbor = expand_plus ? 15 : 9;
                        }
                    }
                } else {
                    // The two cardinal neighbors leave ioff == 0.
                    bool const expand_plus = (nx > 0 && xdir_pls_ok) || !xdir_mns_ok;
                    if (n == 4) {
                        nbor = expand_plus ? 5 : 4;
                    } else if (n == 5) {
                        nbor = expand_plus ? ((joff > 0) ? 8 : 3)
                                           : ((joff > 0) ? 6 : 1);
                    } else if (n == 6) {
                        nbor = expand_plus ? ((koff > 0) ? 23 : 14)
                                           : ((koff > 0) ? 21 : 12);
                    } else {
                        if (joff > 0 && koff > 0) {
                            nbor = expand_plus ? 26 : 24;
                        } else if (joff < 0 && koff > 0) {
                            nbor = expand_plus ? 20 : 18;
                        } else if (joff > 0 && koff < 0) {
                            nbor = expand_plus ? 17 : 15;
                        } else {
                            nbor = expand_plus ? 11 : 9;
                        }
                    }
                }

                int const ioffn = imap[nbor];
                int const joffn = jmap[nbor];
                int const koffn = kmap[nbor];
                sum_vol += vfrac(i+ioffn,j+joffn,k+koffn);
                itracker(i,j,k,n) = nbor;
#if 0
                if (debug_print)
                    amrex::Print() << "Cell " << IntVect(i,j,k) << " with volfrac " << vf <<
                      " trying to ALSO merge with " << IntVect(i+ioffn,j+joffn,k+koffn) <<
                      " with volfrac " << vfrac(i+ioffn,j+joffn,k+koffn) <<
                      " to get new sum_vol " << sum_vol << '\n';
#endif
            }

            // (i,j,k) has a 2x2x2 neighborhood now
            num_neighbors += 4;
        }
    }
    // A cell next to a non-periodic domain boundary can only merge with the cells
    // that are available inside the domain, so its neighborhood may fall short of
    // target_volfrac even though we merged with everything we were allowed to use.
    // We accept the largest neighborhood we could build in that case: reaching
    // outside the domain instead would silently break conservation, because
    // StateRedistribute leaves cells outside the domain out of Qhat even though
    // MakeStateRedistUtils counts their volume in nbhd_vol.
    bool nbhd_limited_by_domain = !(xdir_mns_ok && xdir_pls_ok &&
                                    ydir_mns_ok && ydir_pls_ok &&
                                    zdir_mns_ok && zdir_pls_ok);

    if (sum_vol < target_volfrac && !nbhd_limited_by_domain)
    {
#if 0
      amrex::Print() << "Couldn't merge with enough cells to raise volume at " <<
                        IntVect(i,j,k) << " so stuck with sum_vol " << sum_vol << '\n';
#endif
      amrex::Abort("Couldn't merge with enough cells to raise volume greater than target_volfrac");
    }
    itracker(i,j,k,0) = num_neighbors;
}

#endif

template <bool initialize_state>
int
make_itracker ( Box const& bx,
               AMREX_D_DECL(Array4<Real const> const& apx,
                            Array4<Real const> const& apy,
                            Array4<Real const> const& apz),
               Array4<Real const> const& vfrac,
               Array4<int> const& itracker,
               Array4<EBCellFlag const> const& flag,
               Array4<Real const> const& ccent,
               Array4<int> const& nrs,
               Array4<Real> const& alpha,
               Array4<Real> const& cent_hat,
               Dim3* active_cells,
               Geometry const& lev_geom,
               Real target_volfrac)
{
    Box const domain = lev_geom.Domain();
    AMREX_D_TERM(bool const is_periodic_x = lev_geom.isPeriodic(0);,
                 bool const is_periodic_y = lev_geom.isPeriodic(1);,
                 bool const is_periodic_z = lev_geom.isPeriodic(2););

    Box const bxg3 = amrex::grow(bx,3);
    Box const bxg4 = amrex::grow(bx,4);

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::ignore_unused(bxg3, bxg4, flag, ccent, nrs, alpha, cent_hat);
        itracker(i,j,k,0) = 0;

        if constexpr (initialize_state) {
            initialize_state_data(i, j, k, bxg3, bxg4, flag, ccent,
                                  nrs, alpha, cent_hat);
        }
    });

    Box domain_per_grown = domain;
    AMREX_D_TERM(if (is_periodic_x) { domain_per_grown.grow(0,4); },
                 if (is_periodic_y) { domain_per_grown.grow(1,4); },
                 if (is_periodic_z) { domain_per_grown.grow(2,4); });

    Box const bx_per_g4 = domain_per_grown & bxg4;
    int num_active_cells = 0;

    if constexpr (initialize_state) {
        num_active_cells = find_active_cells(bx_per_g4, vfrac, active_cells,
                                             target_volfrac);
        amrex::ParallelFor(num_active_cells,
        [=] AMREX_GPU_DEVICE (int active_index) noexcept
        {
            Dim3 const cell = active_cells[active_index];
            construct_itracker_cell(
                cell.x, cell.y, cell.z, AMREX_D_DECL(apx, apy, apz),
                vfrac, itracker, domain,
                AMREX_D_DECL(is_periodic_x, is_periodic_y, is_periodic_z),
                target_volfrac);
        });
    } else {
        amrex::ParallelFor(bx_per_g4,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const vf = vfrac(i,j,k);
            if (vf > 0.0 && vf < target_volfrac) {
                construct_itracker_cell(
                    i, j, k, AMREX_D_DECL(apx, apy, apz), vfrac, itracker,
                    domain,
                    AMREX_D_DECL(is_periodic_x, is_periodic_y, is_periodic_z),
                    target_volfrac);
            }
        });
    }
    return num_active_cells;
}

}

void
MakeITracker ( Box const& bx,
               AMREX_D_DECL(Array4<Real const> const& apx,
                            Array4<Real const> const& apy,
                            Array4<Real const> const& apz),
               Array4<Real const> const& vfrac,
               Array4<int> const& itracker,
               Geometry const& lev_geom,
               Real target_volfrac)
{
    make_itracker<false>(bx, AMREX_D_DECL(apx, apy, apz), vfrac, itracker,
                         {}, {}, {}, {}, {}, nullptr, lev_geom,
                         target_volfrac);
}

int
MakeITracker ( Box const& bx,
               AMREX_D_DECL(Array4<Real const> const& apx,
                            Array4<Real const> const& apy,
                            Array4<Real const> const& apz),
               Array4<Real const> const& vfrac,
               Array4<int> const& itracker,
               Array4<EBCellFlag const> const& flag,
               Array4<Real const> const& ccent,
               Array4<int> const& nrs,
               Array4<Real> const& alpha,
               Array4<Real> const& cent_hat,
               Dim3* active_cells,
               Geometry const& lev_geom,
               Real target_volfrac)
{
    return make_itracker<true>(
        bx, AMREX_D_DECL(apx, apy, apz), vfrac, itracker, flag, ccent, nrs,
        alpha, cent_hat, active_cells, lev_geom, target_volfrac);
}

}

/** @} */
