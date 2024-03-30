/**
 * \file StateRedistItracer.cpp
 * @{
 *
 */

#include <AMReX_EB_Redistribution.H>

namespace amrex {

#if (AMREX_SPACEDIM == 2)

void
MakeITracker ( Box const& bx,
               Array4<Real const> const& apx,
               Array4<Real const> const& apy,
               Array4<Real const> const& vfrac,
               Array4<int> const& itracker,
               Geometry const& lev_geom,
               Real target_volfrac)
{
#if 0
    int debug_verbose = 0;
#endif

    const Real small_norm_diff = 1.e-8;

    const Box domain = lev_geom.Domain();

    // Note that itracker has 4 components -- the first component stores the
    //      number of neighbors and each additional component stores where that neighbor is
    //
    // We identify the neighbors with the following ordering
    //
    // ^  6 7 8
    // |  4   5
    // j  1 2 3
    //   i --->

    // Note the first component of imap/jmap should never be used
    Array<int,9> imap{0,-1,0,1,-1,1,-1,0,1};
    Array<int,9> jmap{0,-1,-1,-1,0,0,1,1,1};

    const auto& is_periodic_x = lev_geom.isPeriodic(0);
    const auto& is_periodic_y = lev_geom.isPeriodic(1);

//  if (debug_verbose > 0)
//      amrex::Print() << " IN MAKE_ITRACKER DOING BOX " << bx << '\n';

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        itracker(i,j,k,0) = 0;
    });

    Box domain_per_grown = domain;
    if (is_periodic_x) { domain_per_grown.grow(0,4); }
    if (is_periodic_y) { domain_per_grown.grow(1,4); }

    Box const& bxg4 = amrex::grow(bx,4);
    Box bx_per_g4= domain_per_grown & bxg4;

    amrex::ParallelFor(bx_per_g4,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       if (vfrac(i,j,k) > 0.0 && vfrac(i,j,k) < target_volfrac)
       {
           Real apnorm, apnorm_inv;
           const Real dapx = apx(i+1,j  ,k  ) - apx(i,j,k);
           const Real dapy = apy(i  ,j+1,k  ) - apy(i,j,k);
           apnorm = std::sqrt(dapx*dapx+dapy*dapy);
           apnorm_inv = 1.0/apnorm;
           Real nx = dapx * apnorm_inv;
           Real ny = dapy * apnorm_inv;

           bool nx_eq_ny = ( (std::abs(nx-ny) < small_norm_diff) ||
                             (std::abs(nx+ny) < small_norm_diff)  ) ? true : false;

           // As a first pass, choose just based on the normal
           if (std::abs(nx) > std::abs(ny))
           {
               if (nx > 0) {
                   itracker(i,j,k,1) = 5;
               } else {
                   itracker(i,j,k,1) = 4;
               }

           } else {
               if (ny > 0) {
                   itracker(i,j,k,1) = 7;
               } else {
                   itracker(i,j,k,1) = 2;
               }
           }

           bool xdir_mns_ok = (is_periodic_x || (i > domain.smallEnd(0)));
           bool xdir_pls_ok = (is_periodic_x || (i < domain.bigEnd(0)  ));
           bool ydir_mns_ok = (is_periodic_y || (j > domain.smallEnd(1)));
           bool ydir_pls_ok = (is_periodic_y || (j < domain.bigEnd(1)  ));

           // Override above logic if trying to reach outside a domain boundary (and non-periodic)
           if ( (!xdir_mns_ok && (itracker(i,j,k,1) == 4)) ||
                (!xdir_pls_ok && (itracker(i,j,k,1) == 5)) )
           {
               itracker(i,j,k,1) = (ny > 0) ? 7 : 2;
           }
           if ( (!ydir_mns_ok && (itracker(i,j,k,1) == 2)) ||
                (!ydir_pls_ok && (itracker(i,j,k,1) == 7)) )
           {
               itracker(i,j,k,1) = (nx > 0) ? 5 : 4;
           }

           // (i,j) merges with at least one cell now
           itracker(i,j,k,0) += 1;

           // (i+ioff,j+joff) is in the nbhd of (i,j)
           int ioff = imap[itracker(i,j,k,1)];
           int joff = jmap[itracker(i,j,k,1)];

           // Sanity check
           if (vfrac(i+ioff,j+joff,k) == 0.) {
               amrex::Abort(" Trying to merge with covered cell");
           }

           Real sum_vol = vfrac(i,j,k) + vfrac(i+ioff,j+joff,k);

#if 0
           if (debug_verbose > 0)
               amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vfrac(i,j,k) <<
                                 " trying to merge with " << IntVect(i+ioff,j+joff) <<
                                 " with volfrac " << vfrac(i+ioff,j+joff,k) <<
                                 " to get new sum_vol " <<  sum_vol << '\n';
#endif

           // If the merged cell isn't large enough, we try to merge in the other direction
           if (sum_vol < target_volfrac || nx_eq_ny)
           {
               // Original offset was in y-direction, so we will add to the x-direction
               // Note that if we can't because it would go outside the domain, we don't
               if (ioff == 0) {
                   if (nx >= 0 && xdir_pls_ok)
                   {
                       itracker(i,j,k,2) = 5;
                       itracker(i,j,k,0) += 1;
                   }
                   else if (nx <= 0 && xdir_mns_ok)
                   {
                       itracker(i,j,k,2) = 4;
                       itracker(i,j,k,0) += 1;
                   }

               // Original offset was in x-direction, so we will add to the y-direction
               // Note that if we can't because it would go outside the domain, we don't
               } else {
                   if (ny >= 0 && ydir_pls_ok)
                   {
                       itracker(i,j,k,2) = 7;
                       itracker(i,j,k,0) += 1;
                   }
                   else if (ny <= 0 && ydir_mns_ok)
                   {
                       itracker(i,j,k,2) = 2;
                       itracker(i,j,k,0) += 1;
                   }
               }

               if (itracker(i,j,k,0) > 1)
               {
                   // (i+ioff2,j+joff2) is in the nbhd of (i,j)
                   int ioff2 = imap[itracker(i,j,k,2)];
                   int joff2 = jmap[itracker(i,j,k,2)];

                   sum_vol += vfrac(i+ioff2,j+joff2,k);
#if 0
                   if (debug_verbose > 0)
                       amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vfrac(i,j,k) <<
                                         " trying to ALSO merge with " << IntVect(i+ioff2,j+joff2) <<
                                         " with volfrac " << vfrac(i+ioff2,j+joff2,k) <<
                                          " to get new sum_vol " <<  sum_vol << '\n';
#endif
               }
           }

           // Now we merge in the corner direction if we have already claimed two
           if (itracker(i,j,k,0) == 2)
           {
               // We already have two offsets, and we know they are in different directions
               ioff = imap[itracker(i,j,k,1)] + imap[itracker(i,j,k,2)];
               joff = jmap[itracker(i,j,k,1)] + jmap[itracker(i,j,k,2)];

               if (ioff > 0 && joff > 0) {
                   itracker(i,j,k,3) = 8;
               } else if (ioff < 0 && joff > 0) {
                   itracker(i,j,k,3) = 6;
               } else if (ioff > 0 && joff < 0) {
                   itracker(i,j,k,3) = 3;
               } else {
                   itracker(i,j,k,3) = 1;
               }

               // (i,j) merges with at least three cells now
               itracker(i,j,k,0) += 1;

               sum_vol += vfrac(i+ioff,j+joff,k);
#if 0
               if (debug_verbose > 0)
                   amrex::Print() << "Cell " << IntVect(i,j) << " with volfrac " << vfrac(i,j,k) <<
                                     " trying to ALSO merge with " << IntVect(i+ioff,j+joff) <<
                                     " with volfrac " << vfrac(i+ioff,j+joff,k) <<
                                     " to get new sum_vol " <<  sum_vol << '\n';
#endif
           }
           if (sum_vol < target_volfrac)
           {
#if 0
             amrex::Print() << "Couldnt merge with enough cells to raise volume at " <<
                               IntVect(i,j) << " so stuck with sum_vol " << sum_vol << '\n';
#endif
             amrex::Abort("Couldnt merge with enough cells to raise volume greater than target_volfrac");
           }
       }
    });
}

#elif (AMREX_SPACEDIM == 3)

void
MakeITracker ( Box const& bx,
               Array4<Real const> const& apx,
               Array4<Real const> const& apy,
               Array4<Real const> const& apz,
               Array4<Real const> const& vfrac,
               Array4<int> const& itracker,
               Geometry const& lev_geom,
               Real target_volfrac)
{
#if 0
     bool debug_print = false;
#endif

    const Real small_norm_diff = Real(1.e-8);

    const Box domain = lev_geom.Domain();

    // Note that itracker has 8 components -- the first component stores the
    //      number of neighbors and each additional component stores where that neighbor is
    //
    // We identify the neighbors with the following ordering
    //
    //    at k-1   |   at k  |   at k+1
    //
    // ^  15 16 17 |  6 7 8  |  24 25 26
    // |  12 13 14 |  4   5  |  21 22 23
    // j  9  10 11 |  1 2 3  |  18 19 20
    //   i --->
    //
    // Note the first component of imap/jmap/kmap should never be used
    Array<int,27>    imap{0,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    Array<int,27>    jmap{0,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    Array<int,27>    kmap{0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    const auto& is_periodic_x = lev_geom.isPeriodic(0);
    const auto& is_periodic_y = lev_geom.isPeriodic(1);
    const auto& is_periodic_z = lev_geom.isPeriodic(2);

#if 0
    if (debug_print)
        amrex::Print() << " IN MERGE_REDISTRIBUTE DOING BOX " << bx << '\n';
#endif

    amrex::ParallelFor(Box(itracker),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        itracker(i,j,k,0) = 0;
    });

    Box domain_per_grown = domain;
    if (is_periodic_x) { domain_per_grown.grow(0,4); }
    if (is_periodic_y) { domain_per_grown.grow(1,4); }
#if (AMREX_SPACEDIM == 3)
    if (is_periodic_z) { domain_per_grown.grow(2,4); }
#endif

    Box const& bxg4 = amrex::grow(bx,4);
    Box bx_per_g4= domain_per_grown & bxg4;

    amrex::ParallelFor(bx_per_g4,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
       if (vfrac(i,j,k) > 0.0 && vfrac(i,j,k) < target_volfrac)
       {
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

           // x-component of normal is greatest
           if ( (std::abs(nx) > std::abs(ny)) &&
                (std::abs(nx) > std::abs(nz)) )
           {
               if (nx > 0) {
                   itracker(i,j,k,1) = 5;
               } else {
                   itracker(i,j,k,1) = 4;
               }

           // y-component of normal is greatest
           } else if ( (std::abs(ny) >= std::abs(nx)) &&
                       (std::abs(ny) > std::abs(nz)) )
           {
               if (ny > 0) {
                   itracker(i,j,k,1) = 7;
               } else {
                   itracker(i,j,k,1) = 2;
               }

           // z-component of normal is greatest
           } else {
               if (nz > 0) {
                   itracker(i,j,k,1) = 22;
               } else {
                   itracker(i,j,k,1) = 13;
               }
           }

           // Override above logic if trying to reach outside a domain boundary (and non-periodic)
           if ( (!xdir_mns_ok && (itracker(i,j,k,1) == 4)) ||
                (!xdir_pls_ok && (itracker(i,j,k,1) == 5)) )
           {
               if ( (std::abs(ny) > std::abs(nz)) ) {
                   itracker(i,j,k,1) = (ny > 0) ? 7 : 2;
               } else {
                   itracker(i,j,k,1) = (nz > 0) ? 22 : 13;
               }
           }

           if ( (!ydir_mns_ok && (itracker(i,j,k,1) == 2)) ||
                (!ydir_pls_ok && (itracker(i,j,k,1) == 7)) )
           {
               if ( (std::abs(nx) > std::abs(nz)) ) {
                   itracker(i,j,k,1) = (nx > 0) ? 5 : 4;
               } else {
                   itracker(i,j,k,1) = (nz > 0) ? 22 : 13;
               }
           }

           if ( (!zdir_mns_ok && (itracker(i,j,k,1) == 13)) ||
                (!zdir_pls_ok && (itracker(i,j,k,1) == 22)) )
           {
               if ( (std::abs(nx) > std::abs(ny)) ) {
                   itracker(i,j,k,1) = (nx > 0) ? 5 : 4;
               } else {
                   itracker(i,j,k,1) = (ny > 0) ? 7 : 2;
               }
           }

           // (i,j,k) merges with at least one cell now
           itracker(i,j,k,0) += 1;

           // (i+ioff,j+joff,k+koff) is now the first cell in the nbhd of (i,j,k)
           int ioff = imap[itracker(i,j,k,1)];
           int joff = jmap[itracker(i,j,k,1)];
           int koff = kmap[itracker(i,j,k,1)];

           // Sanity check
           if (vfrac(i+ioff,j+joff,k+koff) == 0.)
           {
               // amrex::Print() << "Cell " << IntVect(i,j,k) << " is trying to merge with cell " << IntVect(i+ioff,j+joff,k+koff) << '\n';
               amrex::Abort(" Trying to merge with covered cell");
           }

           Real sum_vol = vfrac(i,j,k) + vfrac(i+ioff,j+joff,k+koff);

#if 0
           if (debug_print)
               amrex::Print() << "Cell " << IntVect(i,j,k) << " with volfrac " << vfrac(i,j,k) <<
                                 " trying to merge with " << IntVect(i+ioff,j+joff,k+koff) <<
                                 " with volfrac " << vfrac(i+ioff,j+joff,k+koff) <<
                                 " to get new sum_vol " <<  sum_vol << '\n';
#endif

           bool just_broke_symmetry = ( ( (joff == 0 && koff == 0) && (nx_eq_ny || nx_eq_nz) ) ||
                                        ( (ioff == 0 && koff == 0) && (nx_eq_ny || ny_eq_nz) ) ||
                                        ( (ioff == 0 && joff == 0) && (nx_eq_nz || ny_eq_nz) ) );

           // If the merged cell isn't large enough, or if we broke symmetry by the current merge,
           // we merge in one of the other directions.  Note that the direction of the next merge
           // is first set by a symmetry break, but if that isn't happening, we choose the next largest normal

           if ( (sum_vol < target_volfrac) || just_broke_symmetry )
           {
               // Original offset was in x-direction
               if (joff == 0 && koff == 0)
               {
                   if (nx_eq_ny) {
                       itracker(i,j,k,2) = (ny > 0) ? 7 : 2;
                   } else if (nx_eq_nz) {
                       itracker(i,j,k,2) = (nz > 0) ? 22 : 13;
                   } else if ( (std::abs(ny) > std::abs(nz)) ) {
                       itracker(i,j,k,2) = (ny > 0) ? 7 : 2;
                   } else {
                       itracker(i,j,k,2) = (nz > 0) ? 22 : 13;
                   }

               // Original offset was in y-direction
               } else if (ioff == 0 && koff == 0)
               {
                   if (nx_eq_ny) {
                       itracker(i,j,k,2) = (nx > 0) ? 5 : 4;
                   } else if (ny_eq_nz) {
                       itracker(i,j,k,2) = (nz > 0) ? 22 : 13;
                   } else if ( (std::abs(nx) > std::abs(nz)) ) {
                       itracker(i,j,k,2) = (nx > 0) ? 5 : 4;
                   } else {
                       itracker(i,j,k,2) = (nz > 0) ? 22 : 13;
                   }

               // Original offset was in z-direction
               } else if (ioff == 0 && joff == 0)
               {
                   if (nx_eq_nz) {
                       itracker(i,j,k,2) = (nx > 0) ? 5 : 4;
                   } else if (ny_eq_nz) {
                       itracker(i,j,k,2) = (ny > 0) ? 7 : 2;
                   } else if ( (std::abs(nx) > std::abs(ny)) ) {
                       itracker(i,j,k,2) = (nx > 0) ? 5 : 4;
                   } else {
                       itracker(i,j,k,2) = (ny > 0) ? 7 : 2;
                   }
               }

               // (i,j,k) merges with at least two cells now
               itracker(i,j,k,0) += 1;

               // (i+ioff2,j+joff2,k+koff2) is in the nbhd of (i,j,k)
               int ioff2 = imap[itracker(i,j,k,2)];
               int joff2 = jmap[itracker(i,j,k,2)];
               int koff2 = kmap[itracker(i,j,k,2)];

               sum_vol += vfrac(i+ioff2,j+joff2,k+koff2);
#if 0
               if (debug_print)
                   amrex::Print() << "Cell " << IntVect(i,j,k) << " with volfrac " << vfrac(i,j,k) <<
                                     " trying to ALSO merge with " << IntVect(i+ioff2,j+joff2,k+koff2) <<
                                     " with volfrac " << vfrac(i+ioff2,j+joff2,k+koff2) <<
                                      " to get new sum_vol " <<  sum_vol << '\n';
#endif
           }

           // If the merged cell has merged in two directions, we now merge in the corner direction within the current plane
           if (itracker(i,j,k,0) >= 2)
           {
               // We already have two offsets, and we know they are in different directions
               ioff = imap[itracker(i,j,k,1)] + imap[itracker(i,j,k,2)];
               joff = jmap[itracker(i,j,k,1)] + jmap[itracker(i,j,k,2)];
               koff = kmap[itracker(i,j,k,1)] + kmap[itracker(i,j,k,2)];

               // Both nbors are in the koff=0 plane
               if (koff == 0)
               {
                   if (ioff > 0 && joff > 0) {
                       itracker(i,j,k,3) = 8;
                   } else if (ioff < 0 && joff > 0) {
                       itracker(i,j,k,3) = 6;
                   } else if (ioff > 0 && joff < 0) {
                       itracker(i,j,k,3) = 3;
                   } else {
                       itracker(i,j,k,3) = 1;
                   }

               // Both nbors are in the joff=0 plane
               } else if (joff == 0) {
                   if (ioff > 0 && koff > 0) {
                       itracker(i,j,k,3) = 23;
                   } else if (ioff < 0 && koff > 0) {
                       itracker(i,j,k,3) = 21;
                   } else if (ioff > 0 && koff < 0) {
                       itracker(i,j,k,3) = 14;
                   } else {
                       itracker(i,j,k,3) = 12;
                   }

               // Both nbors are in the ioff=0 plane
               } else {
                   if (joff > 0 && koff > 0) {
                       itracker(i,j,k,3) = 25;
                   } else if (joff < 0 && koff > 0) {
                       itracker(i,j,k,3) = 19;
                   } else if (joff > 0 && koff < 0) {
                       itracker(i,j,k,3) = 16;
                   } else {
                       itracker(i,j,k,3) = 10;
                   }
               }

               // (i,j,k) merges with at least three cells now
               itracker(i,j,k,0) += 1;

               sum_vol += vfrac(i+ioff,j+joff,k+koff);

#if 0
               int ioff3 = imap[itracker(i,j,k,3)];
               int joff3 = jmap[itracker(i,j,k,3)];
               int koff3 = kmap[itracker(i,j,k,3)];
               if (debug_print)
                   amrex::Print() << "Cell " << IntVect(i,j,k) << " with volfrac " << vfrac(i,j,k) <<
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
                   // All nbors are currently in the koff=0 plane
                   if (koff == 0)
                   {
                       if (nz > 0)
                       {
                           itracker(i,j,k,4) = 22;

                           if (ioff > 0) {
                               itracker(i,j,k,5) =  23;
                           } else {
                               itracker(i,j,k,5) =  21;
                           }
                           if (joff > 0) {
                               itracker(i,j,k,6) =  25;
                           } else {
                               itracker(i,j,k,6) =  19;
                           }

                           if (ioff > 0 && joff > 0) {
                               itracker(i,j,k,7) = 26;
                           } else if (ioff < 0 && joff > 0) {
                               itracker(i,j,k,7) = 24;
                           } else if (ioff > 0 && joff < 0) {
                               itracker(i,j,k,7) = 20;
                           } else {
                               itracker(i,j,k,7) = 18;
                           }

                       } else { // nz <= 0

                           itracker(i,j,k,4) = 13;

                           if (ioff > 0) {
                               itracker(i,j,k,5) =  14;
                           } else {
                               itracker(i,j,k,5) =  12;
                           }
                           if (joff > 0) {
                               itracker(i,j,k,6) =  16;
                           } else {
                               itracker(i,j,k,6) =  10;
                           }

                           if (ioff > 0 && joff > 0) {
                               itracker(i,j,k,7) = 17;
                           } else if (ioff < 0 && joff > 0) {
                               itracker(i,j,k,7) = 15;
                           } else if (ioff > 0 && joff < 0) {
                               itracker(i,j,k,7) = 11;
                           } else {
                               itracker(i,j,k,7) =  9;
                           }
                       }
                   } else if (joff == 0) {
                       if (ny > 0)
                       {
                           itracker(i,j,k,4) = 7;

                           if (ioff > 0) {
                               itracker(i,j,k,5) =  8;
                           } else {
                               itracker(i,j,k,5) =  6;
                           }
                           if (koff > 0) {
                               itracker(i,j,k,6) =  25;
                           } else {
                               itracker(i,j,k,6) =  16;
                           }

                           if (ioff > 0 && koff > 0) {
                               itracker(i,j,k,7) = 26;
                           } else if (ioff < 0 && koff > 0) {
                               itracker(i,j,k,7) = 24;
                           } else if (ioff > 0 && koff < 0) {
                               itracker(i,j,k,7) = 17;
                           } else {
                               itracker(i,j,k,7) = 15;
                           }

                       } else { // ny <= 0

                           itracker(i,j,k,4) = 2;

                           if (ioff > 0) {
                               itracker(i,j,k,5) =  3;
                           } else {
                               itracker(i,j,k,5) =  1;
                           }
                           if (koff > 0) {
                               itracker(i,j,k,6) =  19;
                           } else {
                               itracker(i,j,k,6) =  10;
                           }

                           if (ioff > 0 && koff > 0) {
                               itracker(i,j,k,7) = 20;
                           } else if (ioff < 0 && koff > 0) {
                               itracker(i,j,k,7) = 18;
                           } else if (ioff > 0 && koff < 0) {
                               itracker(i,j,k,7) = 11;
                           } else {
                               itracker(i,j,k,7) =  9;
                           }
                       }
                   } else if (ioff == 0) {

                       if (nx > 0)
                       {
                           itracker(i,j,k,4) = 5;

                           if (joff > 0) {
                               itracker(i,j,k,5) =  8;
                           } else {
                               itracker(i,j,k,5) =  3;
                           }
                           if (koff > 0) {
                               itracker(i,j,k,6) =  23;
                           } else {
                               itracker(i,j,k,6) =  14;
                           }

                           if (joff > 0 && koff > 0) {
                               itracker(i,j,k,7) = 26;
                           } else if (joff < 0 && koff > 0) {
                               itracker(i,j,k,7) = 20;
                           } else if (joff > 0 && koff < 0) {
                               itracker(i,j,k,7) = 17;
                           } else {
                               itracker(i,j,k,7) = 11;
                           }
                       } else { // nx <= 0

                           itracker(i,j,k,4) = 4;

                           if (joff > 0) {
                               itracker(i,j,k,5) =  6;
                           } else {
                               itracker(i,j,k,5) =  1;
                           }
                           if (koff > 0) {
                               itracker(i,j,k,6) =  21;
                           } else {
                               itracker(i,j,k,6) =  12;
                           }

                           if (joff > 0 && koff > 0) {
                               itracker(i,j,k,7) = 24;
                           } else if (joff < 0 && koff > 0) {
                               itracker(i,j,k,7) = 18;
                           } else if (joff > 0 && koff < 0) {
                               itracker(i,j,k,7) = 15;
                           } else {
                               itracker(i,j,k,7) =  9;
                           }
                       }
                   }

                   for(int n(4); n<8; n++){
                     int ioffn = imap[itracker(i,j,k,n)];
                     int joffn = jmap[itracker(i,j,k,n)];
                     int koffn = kmap[itracker(i,j,k,n)];
                     sum_vol += vfrac(i+ioffn,j+joffn,k+koffn);
#if 0
                     if (debug_print)
                       amrex::Print() << "Cell " << IntVect(i,j,k) << " with volfrac " << vfrac(i,j,k) <<
                         " trying to ALSO merge with " << IntVect(i+ioffn,j+joffn,k+koffn) <<
                         " with volfrac " << vfrac(i+ioffn,j+joffn,k+koffn) <<
                         " to get new sum_vol " <<  sum_vol << '\n';
#endif
                   }

                   // (i,j,k) has a 2x2x2 neighborhood now
                   itracker(i,j,k,0) += 4;
               }
           }
           if (sum_vol < target_volfrac)
           {
#if 0
             amrex::Print() << "Couldnt merge with enough cells to raise volume at " <<
                               IntVect(i,j,k) << " so stuck with sum_vol " << sum_vol << '\n';
#endif
             amrex::Abort("Couldnt merge with enough cells to raise volume greater than target_volfrac");
           }
       }
    });
}

#endif

}
