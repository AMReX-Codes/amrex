#ifndef FIELDGATHER_H_
#define FIELDGATHER_H_

#include "ShapeFactors.H"

/* \brief Field gather for particles handled by thread thread_num
 * /param xp, yp, zp   : Pointer to arrays of particle positions.
 * \param Exp, Eyp, Ezp: Pointer to array of electric field on particles.
 * \param Bxp, Byp, Bzp: Pointer to array of magnetic field on particles.
 * \param ex_arr ey_arr: Array4 of current density, either full array or tile.
 * \param ez_arr bx_arr: Array4 of current density, either full array or tile.
 * \param by_arr bz_arr: Array4 of current density, either full array or tile.
 * \param np_to_gather : Number of particles for which field is gathered.
 * \param dx           : 3D cell size
 * \param xyzmin       : Physical lower bounds of domain.
 * \param lo           : Index lower bounds of domain.
 * \param stagger_shift: 0 if nodal, 0.5 if staggered.
 */
template <int depos_order, int lower_in_v>
void doGatherShapeN(const amrex::Real * const xp,
                    const amrex::Real * const yp,
                    const amrex::Real * const zp,
                    amrex::Real * const Exp, amrex::Real * const Eyp,
                    amrex::Real * const Ezp, amrex::Real * const Bxp,
                    amrex::Real * const Byp, amrex::Real * const Bzp,
                    const amrex::Array4<const amrex::Real>& ex_arr,
                    const amrex::Array4<const amrex::Real>& ey_arr,
                    const amrex::Array4<const amrex::Real>& ez_arr,
                    const amrex::Array4<const amrex::Real>& bx_arr,
                    const amrex::Array4<const amrex::Real>& by_arr,
                    const amrex::Array4<const amrex::Real>& bz_arr,
                    const long np_to_gather,
                    const std::array<amrex::Real, 3>& dx,
                    const std::array<amrex::Real, 3> xyzmin,
                    const amrex::Dim3 lo,
                    const amrex::Real stagger_shift)
{
    const amrex::Real dxi = 1.0/dx[0];
    const amrex::Real dzi = 1.0/dx[2];
#if (AMREX_SPACEDIM == 3)
    const amrex::Real dyi = 1.0/dx[1];
#endif

    const amrex::Real xmin = xyzmin[0];
    const amrex::Real ymin = xyzmin[1];
    const amrex::Real zmin = xyzmin[2];

    // Loop over particles and gather fields from
    // {e,b}{x,y,z}_arr to {E,B}{xyz}p.
    amrex::ParallelFor(
        np_to_gather,
        [=] AMREX_GPU_DEVICE (long ip) {
            // --- Compute shape factors
            // x direction
            // Get particle position
            const amrex::Real x = (xp[ip]-xmin)*dxi;
            // Compute shape factors for node-centered quantities
            amrex::Real AMREX_RESTRICT sx [depos_order + 1];
            // j: leftmost grid point (node-centered) that particle touches
            const int j  = compute_shape_factor<depos_order>(sx, x);
            // Compute shape factors for cell-centered quantities
            amrex::Real AMREX_RESTRICT sx0[depos_order + 1 - lower_in_v];
            // j0: leftmost grid point (cell-centered) that particle touches
            const int j0 = compute_shape_factor<depos_order - lower_in_v>(
                sx0, x-stagger_shift);
#if (AMREX_SPACEDIM == 3)
            // y direction
            const amrex::Real y = (yp[ip]-ymin)*dyi;
            amrex::Real AMREX_RESTRICT sy [depos_order + 1];
            const int k  = compute_shape_factor<depos_order>(sy, y);
            amrex::Real AMREX_RESTRICT sy0[depos_order + 1 - lower_in_v];
            const int k0 = compute_shape_factor<depos_order-lower_in_v>(
                sy0, y-stagger_shift);
#endif
            // z direction
            const amrex::Real z = (zp[ip]-zmin)*dzi;
            amrex::Real AMREX_RESTRICT sz [depos_order + 1];
            const int l  = compute_shape_factor<depos_order>(sz, z);
            amrex::Real AMREX_RESTRICT sz0[depos_order + 1 - lower_in_v];
            const int l0 = compute_shape_factor<depos_order - lower_in_v>(
                sz0, z-stagger_shift);

            // Set fields on particle to zero
            Exp[ip] = 0;
            Eyp[ip] = 0;
            Ezp[ip] = 0;
            Bxp[ip] = 0;
            Byp[ip] = 0;
            Bzp[ip] = 0;
            // Each field is gathered in a separate block of 
            // AMREX_SPACEDIM nested loops because the deposition
            // order can differ for each component of each field
            // when lower_in_v is set to 1
#if (AMREX_SPACEDIM == 2)
            // Gather field on particle Eyp[i] from field on grid ey_arr
            for (int iz=0; iz<=depos_order; iz++){
                for (int ix=0; ix<=depos_order; ix++){
                    Eyp[ip] += sx[ix]*sz[iz]*
                        ey_arr(lo.x+j+ix, lo.y+l+iz, 0);
                }
            }
            // Gather field on particle Exp[i] from field on grid ex_arr
            // Gather field on particle Bzp[i] from field on grid bz_arr
            for (int iz=0; iz<=depos_order; iz++){
                for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                    Exp[ip] += sx0[ix]*sz[iz]*
                        ex_arr(lo.x+j0+ix, lo.y+l +iz, 0);
                    Bzp[ip] += sx0[ix]*sz[iz]*
                        bz_arr(lo.x+j0+ix, lo.y+l +iz, 0);
                }
            }
            // Gather field on particle Ezp[i] from field on grid ez_arr
            // Gather field on particle Bxp[i] from field on grid bx_arr
            for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                for (int ix=0; ix<=depos_order; ix++){
                    Ezp[ip] += sx[ix]*sz0[iz]*
                        ez_arr(lo.x+j+ix, lo.y+l0 +iz, 0);
                    Bxp[ip] += sx[ix]*sz0[iz]*
                        bx_arr(lo.x+j+ix, lo.y+l0 +iz, 0);
                }
            }
            // Gather field on particle Byp[i] from field on grid by_arr
            for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                    Byp[ip] += sx0[ix]*sz0[iz]*
                        by_arr(lo.x+j0+ix, lo.y+l0+iz, 0);
                }
            }
#else // (AMREX_SPACEDIM == 3)
            // Gather field on particle Exp[i] from field on grid ex_arr
            for (int iz=0; iz<=depos_order; iz++){
                for (int iy=0; iy<=depos_order; iy++){
                    for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                        Exp[ip] += sx0[ix]*sy[iy]*sz[iz]*
                            ex_arr(lo.x+j0+ix, lo.y+k+iy, lo.z+l+iz);
                    }
                }
            }
            // Gather field on particle Eyp[i] from field on grid ey_arr
            for (int iz=0; iz<=depos_order; iz++){
                for (int iy=0; iy<=depos_order-lower_in_v; iy++){
                    for (int ix=0; ix<=depos_order; ix++){
                        Eyp[ip] += sx[ix]*sy0[iy]*sz[iz]*
                            ey_arr(lo.x+j+ix, lo.y+k0+iy, lo.z+l+iz);
                    }
                }
            }
            // Gather field on particle Ezp[i] from field on grid ez_arr
            for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                for (int iy=0; iy<=depos_order; iy++){
                    for (int ix=0; ix<=depos_order; ix++){
                        Ezp[ip] += sx[ix]*sy[iy]*sz0[iz]*
                            ez_arr(lo.x+j+ix, lo.y+k+iy, lo.z+l0+iz);
                    }
                }
            }
            // Gather field on particle Bzp[i] from field on grid bz_arr
            for (int iz=0; iz<=depos_order; iz++){
                for (int iy=0; iy<=depos_order-lower_in_v; iy++){
                    for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                        Bzp[ip] += sx0[ix]*sy0[iy]*sz[iz]*
                            bz_arr(lo.x+j0+ix, lo.y+k0+iy, lo.z+l+iz);
                    }
                }
            }
            // Gather field on particle Byp[i] from field on grid by_arr
            for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                for (int iy=0; iy<=depos_order; iy++){
                    for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                        Byp[ip] += sx0[ix]*sy[iy]*sz0[iz]*
                            by_arr(lo.x+j0+ix, lo.y+k+iy, lo.z+l0+iz);
                    }
                }
            }
            // Gather field on particle Bxp[i] from field on grid bx_arr
            for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                for (int iy=0; iy<=depos_order-lower_in_v; iy++){
                    for (int ix=0; ix<=depos_order; ix++){
                        Bxp[ip] += sx[ix]*sy0[iy]*sz0[iz]*
                            bx_arr(lo.x+j+ix, lo.y+k0+iy, lo.z+l0+iz);
                    }
                }
            }
#endif
        }
        );
}

#endif // FIELDGATHER_H_
