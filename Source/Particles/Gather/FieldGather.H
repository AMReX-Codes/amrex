#ifndef FIELDGATHER_H_
#define FIELDGATHER_H_

#include "ShapeFactors.H"

using namespace amrex;

/* \brief Field gather for particles handled by thread thread_num
 * /param xp, yp, zp   : Pointer to arrays of particle positions.
 * \param Exp, Eyp, Ezp: Pointer to array of electric field on particles.
 * \param Bxp, Byp, Bzp: Pointer to array of magnetic field on particles.
 * \param ex_arr ey_arr: Array4 of current density, either full array or tile.
 * \param ez_arr bx_arr: Array4 of current density, either full array or tile.
 * \param by_arr bz_arr: Array4 of current density, either full array or tile.
 * \param np_to_gather : Number of particles for which current is deposited.
 * \param dx           : 3D cell size
 * \param xyzmin       : Physical lower bounds of domain.
 * \param lo           : Index lower bounds of domain.
 * \param stagger_shift: 0 if nodal, 0.5 if staggered.
 */
template <int depos_order, int lower_in_v>
void doGatherShapeN(const Real * const xp, 
                    const Real * const yp, 
                    const Real * const zp,
                    Real * const Exp, Real * const Eyp, Real * const Ezp,
                    Real * const Bxp, Real * const Byp, Real * const Bzp,
                    const amrex::Array4<const amrex::Real>& ex_arr, 
                    const amrex::Array4<const amrex::Real>& ey_arr, 
                    const amrex::Array4<const amrex::Real>& ez_arr,
                    const amrex::Array4<const amrex::Real>& bx_arr, 
                    const amrex::Array4<const amrex::Real>& by_arr, 
                    const amrex::Array4<const amrex::Real>& bz_arr,
                    const long np_to_gather, 
                    const std::array<amrex::Real,3>& dx,
                    const std::array<Real, 3> xyzmin,
                    const Dim3 lo,
                    const amrex::Real stagger_shift)
{
    const Real dxi = 1.0/dx[0];
    const Real dzi = 1.0/dx[2];
#if (AMREX_SPACEDIM == 3)
    const Real dyi = 1.0/dx[1];
#endif

    const Real xmin = xyzmin[0];
    const Real ymin = xyzmin[1];
    const Real zmin = xyzmin[2];

    // Loop over particles and gather fields from {e,b}{x,y,z}_arr to {E,B}{xyz}p.
    ParallelFor( np_to_gather,
                 [=] AMREX_GPU_DEVICE (long ip) {
                     // --- Compute shape factors
                     // x direction
                     // Get particle position after 1/2 push back in position
                     const Real xmid = (xp[ip]-xmin)*dxi;
                     // Compute shape factors for node-centered quantities
                     Real AMREX_RESTRICT sx [depos_order + 1];
                     // j: leftmost grid point (node-centered) that the particle touches
                     const int j  = compute_shape_factor<depos_order>(sx, xmid);
                     // Compute shape factors for cell-centered quantities
                     Real AMREX_RESTRICT sx0[depos_order + 1 - lower_in_v];
                     // j0: leftmost grid point (cell-centered) that the particle touches
                     const int j0 = compute_shape_factor<depos_order - lower_in_v>(sx0, xmid-stagger_shift);
                     
#if (AMREX_SPACEDIM == 3)
                     // y direction
                     const Real ymid= (yp[ip]-ymin)*dyi;
                     Real AMREX_RESTRICT sy [depos_order + 1];
                     const int k  = compute_shape_factor<depos_order>(sy,  ymid);
                     Real AMREX_RESTRICT sy0[depos_order + 1 - lower_in_v];
                     const int k0 = compute_shape_factor<depos_order-lower_in_v>(sy0, ymid-stagger_shift);
#endif
                     // z direction
                     const Real zmid= (zp[ip]-zmin)*dzi;
                     Real AMREX_RESTRICT sz [depos_order + 1];
                     const int l  = compute_shape_factor<depos_order>(sz,  zmid);
                     Real AMREX_RESTRICT sz0[depos_order + 1 - lower_in_v];
                     const int l0 = compute_shape_factor<depos_order - lower_in_v>(sz0, zmid-stagger_shift);

                     // Deposit current into jx_arr, jy_arr and jz_arr
                     Exp[ip] = 0;
                     Eyp[ip] = 0;
                     Ezp[ip] = 0;
                     Bxp[ip] = 0;
                     Byp[ip] = 0;
                     Bzp[ip] = 0;
#if (AMREX_SPACEDIM == 2)
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int ix=0; ix<=depos_order; ix++){
                             Eyp[ip] += sx[ix]*sz[iz]*
                                 ey_arr(lo.x+j+ix, lo.y+l+iz, 0);
                         }
                     }
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                             Exp[ip] += sx0[ix]*sz[iz]*
                                 ex_arr(lo.x+j0+ix, lo.y+l +iz, 0);
                             Bzp[ip] += sx0[ix]*sz[iz]*
                                 bz_arr(lo.x+j0+ix, lo.y+l +iz, 0);
                         }
                     }
                     for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                         for (int ix=0; ix<=depos_order; ix++){
                             Ezp[ip] += sx[ix]*sz0[iz]*
                                 ez_arr(lo.x+j+ix, lo.y+l0 +iz, 0);
                             Bxp[ip] += sx[ix]*sz0[iz]*
                                 bx_arr(lo.x+j+ix, lo.y+l0 +iz, 0);
                         }
                     }
                     for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                         for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                             Byp[ip] += sx0[ix]*sz0[iz]*
                                 by_arr(lo.x+j0+ix, lo.y+l0+iz, 0);
                         }
                     }
#else // (AMREX_SPACEDIM == 3)
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int iy=0; iy<=depos_order; iy++){
                             for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                                 Exp[ip] += sx0[ix]*sy[iy]*sz[iz]*
                                     ex_arr(lo.x+j0+ix, lo.y+k+iy, lo.z+l+iz);
                             }
                         }
                     }
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int iy=0; iy<=depos_order-lower_in_v; iy++){
                             for (int ix=0; ix<=depos_order; ix++){
                                 Eyp[ip] += sx[ix]*sy0[iy]*sz[iz]*
                                     ey_arr(lo.x+j+ix, lo.y+k0+iy, lo.z+l+iz);
                             }
                         }
                     }
                     for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                         for (int iy=0; iy<=depos_order; iy++){
                             for (int ix=0; ix<=depos_order; ix++){
                                 Ezp[ip] += sx[ix]*sy[iy]*sz0[iz]*
                                     ez_arr(lo.x+j+ix, lo.y+k+iy, lo.z+l0+iz);
                             }
                         }
                     }

                     for (int iz=0; iz<=depos_order; iz++){
                         for (int iy=0; iy<=depos_order-lower_in_v; iy++){
                             for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                                 Bzp[ip] += sx0[ix]*sy0[iy]*sz[iz]*
                                     bz_arr(lo.x+j0+ix, lo.y+k0+iy, lo.z+l+iz);
                             }
                         }
                     }
                     for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                         for (int iy=0; iy<=depos_order; iy++){
                             for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                                 Byp[ip] += sx0[ix]*sy[iy]*sz0[iz]*
                                     by_arr(lo.x+j0+ix, lo.y+k+iy, lo.z+l0+iz);
                             }
                         }
                     }
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
