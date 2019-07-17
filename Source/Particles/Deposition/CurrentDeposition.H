#ifndef CURRENTDEPOSITION_H_
#define CURRENTDEPOSITION_H_

#include "ShapeFactors.H"

using namespace amrex;

/* \brief Current Deposition for thread thread_num
 * /param xp, yp, zp   : Pointer to arrays of particle positions.
 * \param wp           : Pointer to array of particle weights.
 * \param uxp uyp uzp  : Pointer to arrays of particle momentum.
 * \param jx_arr       : Array4 of current density, either full array or tile.
 * \param jy_arr       : Array4 of current density, either full array or tile.
 * \param jz_arr       : Array4 of current density, either full array or tile.
 * \param np_to_depose : Number of particles for which current is deposited.
 * \param dt           : Time step for particle level
 * \param dx           : 3D cell size
 * \param xyzmin       : Physical lower bounds of domain.
 * \param lo           : Index lower bounds of domain.
 * \param stagger_shift: 0 if nodal, 0.5 if staggered.
 * /param q            : species charge.
 */
template <int depos_order>
void doDepositionShapeN(const Real * const xp, const Real * const yp, const Real * const zp,
                        const Real * const wp, const Real * const uxp,
                        const Real * const uyp, const Real * const uzp,
                        const amrex::Array4<amrex::Real>& jx_arr, 
                        const amrex::Array4<amrex::Real>& jy_arr, 
                        const amrex::Array4<amrex::Real>& jz_arr,
                        const long np_to_depose, 
                        const amrex::Real dt, const std::array<amrex::Real,3>& dx,
                        const std::array<Real, 3> xyzmin,
                        const Dim3 lo,
                        const amrex::Real stagger_shift, 
                        const amrex::Real q)
{
    const Real dxi = 1.0/dx[0];
    const Real dzi = 1.0/dx[2];
    const Real dts2dx = 0.5*dt*dxi;
    const Real dts2dz = 0.5*dt*dzi;
#if (AMREX_SPACEDIM == 2)
    const Real invvol = dxi*dzi;
#else // (AMREX_SPACEDIM == 3)
    const Real dyi = 1.0/dx[1];
    const Real dts2dy = 0.5*dt*dyi;
    const Real invvol = dxi*dyi*dzi;
#endif

    const Real xmin = xyzmin[0];
    const Real ymin = xyzmin[1];
    const Real zmin = xyzmin[2];
    const Real clightsq = 1.0/PhysConst::c/PhysConst::c;

    // Loop over particles and deposit into jx_arr, jy_arr and jz_arr
    ParallelFor( np_to_depose,
                 [=] AMREX_GPU_DEVICE (long ip) {
                     // --- Get particle quantities
                     const Real gaminv = 1.0/std::sqrt(1.0 + uxp[ip]*uxp[ip]*clightsq
                                                       + uyp[ip]*uyp[ip]*clightsq
                                                       + uzp[ip]*uzp[ip]*clightsq);
                     const Real wq  = q*wp[ip];
                     const Real vx  = uxp[ip]*gaminv;
                     const Real vy  = uyp[ip]*gaminv;
                     const Real vz  = uzp[ip]*gaminv;
                     // wqx, wqy wqz are particle current in each direction 
                     const Real wqx = wq*invvol*vx;
                     const Real wqy = wq*invvol*vy;
                     const Real wqz = wq*invvol*vz;

                     // --- Compute shape factors
                     // x direction
                     // Get particle position after 1/2 push back in position
                     const Real xmid = (xp[ip]-xmin)*dxi-dts2dx*vx;
                     // Compute shape factors for node-centered quantities
                     Real AMREX_RESTRICT sx [depos_order + 1];
                     // j: leftmost grid point (node-centered) that the particle touches
                     const int j  = compute_shape_factor<depos_order>(sx,  xmid);
                     // Compute shape factors for cell-centered quantities
                     Real AMREX_RESTRICT sx0[depos_order + 1];
                     // j0: leftmost grid point (cell-centered) that the particle touches
                     const int j0 = compute_shape_factor<depos_order>(sx0, xmid-stagger_shift);
                     
#if (AMREX_SPACEDIM == 3)
                     // y direction
                     const Real ymid= (yp[ip]-ymin)*dyi-dts2dy*vy;
                     Real AMREX_RESTRICT sy [depos_order + 1];
                     const int k  = compute_shape_factor<depos_order>(sy,  ymid);
                     Real AMREX_RESTRICT sy0[depos_order + 1];
                     const int k0 = compute_shape_factor<depos_order>(sy0, ymid-stagger_shift);
#endif
                     // z direction
                     const Real zmid= (zp[ip]-zmin)*dzi-dts2dz*vz;
                     Real AMREX_RESTRICT sz [depos_order + 1];
                     const int l  = compute_shape_factor<depos_order>(sz,  zmid);
                     Real AMREX_RESTRICT sz0[depos_order + 1];
                     const int l0 = compute_shape_factor<depos_order>(sz0, zmid-stagger_shift);

                     // Deposit current into jx_arr, jy_arr and jz_arr
#if (AMREX_SPACEDIM == 2)
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int ix=0; ix<=depos_order; ix++){
                             amrex::Gpu::Atomic::Add(
                                 &jx_arr(lo.x+j0+ix, lo.y+l +iz, 0), 
                                 sx0[ix]*sz [iz]*wqx);
                             amrex::Gpu::Atomic::Add(
                                 &jy_arr(lo.x+j +ix, lo.y+l +iz, 0), 
                                 sx [ix]*sz [iz]*wqy);
                             amrex::Gpu::Atomic::Add(
                                 &jz_arr(lo.x+j +ix, lo.y+l0+iz, 0), 
                                 sx [ix]*sz0[iz]*wqz);
                         }
                     }
#else // (AMREX_SPACEDIM == 3)
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int iy=0; iy<=depos_order; iy++){
                             for (int ix=0; ix<=depos_order; ix++){
                                 amrex::Gpu::Atomic::Add(
                                     &jx_arr(lo.x+j0+ix, lo.y+k +iy, lo.z+l +iz),
                                     sx0[ix]*sy [iy]*sz [iz]*wqx);
                                 amrex::Gpu::Atomic::Add(
                                     &jy_arr(lo.x+j +ix, lo.y+k0+iy, lo.z+l +iz), 
                                     sx [ix]*sy0[iy]*sz [iz]*wqy);
                                 amrex::Gpu::Atomic::Add(
                                     &jz_arr(lo.x+j +ix, lo.y+k +iy, lo.z+l0+iz),
                                     sx [ix]*sy [iy]*sz0[iz]*wqz);
                             }
                         }
                     }
#endif
                 }
        );
}

#endif // CURRENTDEPOSITION_H_
