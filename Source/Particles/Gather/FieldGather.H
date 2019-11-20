#ifndef FIELDGATHER_H_
#define FIELDGATHER_H_

#include "ShapeFactors.H"
#include <WarpX_Complex.H>

/* \brief Field gather for particles handled by thread thread_num
 * \param xp, yp, zp   : Pointer to arrays of particle positions.
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
 * \param n_rz_azimuthal_modes: Number of azimuthal modes when using RZ geometry
 */
template <int depos_order, int lower_in_v>
void doGatherShapeN(const amrex::ParticleReal * const xp,
                    const amrex::ParticleReal * const yp,
                    const amrex::ParticleReal * const zp,
                    amrex::ParticleReal * const Exp, amrex::ParticleReal * const Eyp,
                    amrex::ParticleReal * const Ezp, amrex::ParticleReal * const Bxp,
                    amrex::ParticleReal * const Byp, amrex::ParticleReal * const Bzp,
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
                    const amrex::Real stagger_shift,
                    const long n_rz_azimuthal_modes)
{
    const amrex::Real dxi = 1.0/dx[0];
    const amrex::Real dzi = 1.0/dx[2];
#if (AMREX_SPACEDIM == 3)
    const amrex::Real dyi = 1.0/dx[1];
#endif

    const amrex::Real xmin = xyzmin[0];
#if (AMREX_SPACEDIM == 3)
    const amrex::Real ymin = xyzmin[1];
#endif
    const amrex::Real zmin = xyzmin[2];

    // Loop over particles and gather fields from
    // {e,b}{x,y,z}_arr to {E,B}{xyz}p.
    amrex::ParallelFor(
        np_to_gather,
        [=] AMREX_GPU_DEVICE (long ip) {
            // --- Compute shape factors
            // x direction
            // Get particle position
#ifdef WARPX_DIM_RZ
            const amrex::Real rp = std::sqrt(xp[ip]*xp[ip] + yp[ip]*yp[ip]);
            const amrex::Real x = (rp - xmin)*dxi;
#else
            const amrex::Real x = (xp[ip]-xmin)*dxi;
#endif
            // Compute shape factors for node-centered quantities
            amrex::Real sx [depos_order + 1];
            // j: leftmost grid point (node-centered) that particle touches
            const int j  = compute_shape_factor<depos_order>(sx, x);
            // Compute shape factors for cell-centered quantities
            amrex::Real sx0[depos_order + 1 - lower_in_v];
            // j0: leftmost grid point (cell-centered) that particle touches
            const int j0 = compute_shape_factor<depos_order - lower_in_v>(
                sx0, x-stagger_shift);
#if (AMREX_SPACEDIM == 3)
            // y direction
            const amrex::Real y = (yp[ip]-ymin)*dyi;
            amrex::Real sy [depos_order + 1];
            const int k  = compute_shape_factor<depos_order>(sy, y);
            amrex::Real sy0[depos_order + 1 - lower_in_v];
            const int k0 = compute_shape_factor<depos_order-lower_in_v>(
                sy0, y-stagger_shift);
#endif
            // z direction
            const amrex::Real z = (zp[ip]-zmin)*dzi;
            amrex::Real sz [depos_order + 1];
            const int l  = compute_shape_factor<depos_order>(sz, z);
            amrex::Real sz0[depos_order + 1 - lower_in_v];
            const int l0 = compute_shape_factor<depos_order - lower_in_v>(
                sz0, z-stagger_shift);

            // Each field is gathered in a separate block of
            // AMREX_SPACEDIM nested loops because the deposition
            // order can differ for each component of each field
            // when lower_in_v is set to 1
#if (AMREX_SPACEDIM == 2)
            // Gather field on particle Eyp[i] from field on grid ey_arr
            for (int iz=0; iz<=depos_order; iz++){
                for (int ix=0; ix<=depos_order; ix++){
                    Eyp[ip] += sx[ix]*sz[iz]*
                        ey_arr(lo.x+j+ix, lo.y+l+iz, 0, 0);
                }
            }
            // Gather field on particle Exp[i] from field on grid ex_arr
            // Gather field on particle Bzp[i] from field on grid bz_arr
            for (int iz=0; iz<=depos_order; iz++){
                for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                    Exp[ip] += sx0[ix]*sz[iz]*
                        ex_arr(lo.x+j0+ix, lo.y+l +iz, 0, 0);
                    Bzp[ip] += sx0[ix]*sz[iz]*
                        bz_arr(lo.x+j0+ix, lo.y+l +iz, 0, 0);
                }
            }
            // Gather field on particle Ezp[i] from field on grid ez_arr
            // Gather field on particle Bxp[i] from field on grid bx_arr
            for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                for (int ix=0; ix<=depos_order; ix++){
                    Ezp[ip] += sx[ix]*sz0[iz]*
                        ez_arr(lo.x+j+ix, lo.y+l0 +iz, 0, 0);
                    Bxp[ip] += sx[ix]*sz0[iz]*
                        bx_arr(lo.x+j+ix, lo.y+l0 +iz, 0, 0);
                }
            }
            // Gather field on particle Byp[i] from field on grid by_arr
            for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                    Byp[ip] += sx0[ix]*sz0[iz]*
                        by_arr(lo.x+j0+ix, lo.y+l0+iz, 0, 0);
                }
            }

#ifdef WARPX_DIM_RZ

            amrex::Real costheta;
            amrex::Real sintheta;
            if (rp > 0.) {
                costheta = xp[ip]/rp;
                sintheta = yp[ip]/rp;
            } else {
                costheta = 1.;
                sintheta = 0.;
            }
            const Complex xy0 = Complex{costheta, -sintheta};
            Complex xy = xy0;

            for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {

                // Gather field on particle Eyp[i] from field on grid ey_arr
                for (int iz=0; iz<=depos_order; iz++){
                    for (int ix=0; ix<=depos_order; ix++){
                        const amrex::Real dEy = (+ ey_arr(lo.x+j+ix, lo.y+l+iz, 0, 2*imode-1)*xy.real()
                                                 - ey_arr(lo.x+j+ix, lo.y+l+iz, 0, 2*imode)*xy.imag());
                        Eyp[ip] += sx[ix]*sz[iz]*dEy;
                    }
                }
                // Gather field on particle Exp[i] from field on grid ex_arr
                // Gather field on particle Bzp[i] from field on grid bz_arr
                for (int iz=0; iz<=depos_order; iz++){
                    for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                        const amrex::Real dEx = (+ ex_arr(lo.x+j0+ix, lo.y+l +iz, 0, 2*imode-1)*xy.real()
                                                 - ex_arr(lo.x+j0+ix, lo.y+l +iz, 0, 2*imode)*xy.imag());
                        Exp[ip] += sx0[ix]*sz[iz]*dEx;
                        const amrex::Real dBz = (+ bz_arr(lo.x+j0+ix, lo.y+l +iz, 0, 2*imode-1)*xy.real()
                                                 - bz_arr(lo.x+j0+ix, lo.y+l +iz, 0, 2*imode)*xy.imag());
                        Bzp[ip] += sx0[ix]*sz[iz]*dBz;
                    }
                }
                // Gather field on particle Ezp[i] from field on grid ez_arr
                // Gather field on particle Bxp[i] from field on grid bx_arr
                for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                    for (int ix=0; ix<=depos_order; ix++){
                        const amrex::Real dEz = (+ ez_arr(lo.x+j+ix, lo.y+l0 +iz, 0, 2*imode-1)*xy.real()
                                                 - ez_arr(lo.x+j+ix, lo.y+l0 +iz, 0, 2*imode)*xy.imag());
                        Ezp[ip] += sx[ix]*sz0[iz]*dEz;
                        const amrex::Real dBx = (+ bx_arr(lo.x+j+ix, lo.y+l0 +iz, 0, 2*imode-1)*xy.real()
                                                 - bx_arr(lo.x+j+ix, lo.y+l0 +iz, 0, 2*imode)*xy.imag());
                        Bxp[ip] += sx[ix]*sz0[iz]*dBx;
                    }
                }
                // Gather field on particle Byp[i] from field on grid by_arr
                for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                    for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                        const amrex::Real dBy = (+ by_arr(lo.x+j0+ix, lo.y+l0+iz, 0, 2*imode-1)*xy.real()
                                                 - by_arr(lo.x+j0+ix, lo.y+l0+iz, 0, 2*imode)*xy.imag());
                        Byp[ip] += sx0[ix]*sz0[iz]*dBy;
                    }
                }
                xy = xy*xy0;
            }

            // Convert Exp and Eyp (which are actually Er and Etheta) to Ex and Ey
            const amrex::Real Exp_save = Exp[ip];
            Exp[ip] = costheta*Exp[ip] - sintheta*Eyp[ip];
            Eyp[ip] = costheta*Eyp[ip] + sintheta*Exp_save;
            const amrex::Real Bxp_save = Bxp[ip];
            Bxp[ip] = costheta*Bxp[ip] - sintheta*Byp[ip];
            Byp[ip] = costheta*Byp[ip] + sintheta*Bxp_save;
#endif

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
