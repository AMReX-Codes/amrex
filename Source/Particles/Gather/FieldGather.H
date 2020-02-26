/* Copyright 2019 Axel Huebl, David Grote, Maxence Thevenet
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef FIELDGATHER_H_
#define FIELDGATHER_H_

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/ShapeFactors.H"
#include "Utils/WarpX_Complex.H"


/**
 * \brief Field gather for particles handled by thread thread_num
 * /param GetPosition : A functor for returning the particle position.
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
void doGatherShapeN(const GetParticlePosition& GetPosition,
                    amrex::ParticleReal * const Exp, amrex::ParticleReal * const Eyp,
                    amrex::ParticleReal * const Ezp, amrex::ParticleReal * const Bxp,
                    amrex::ParticleReal * const Byp, amrex::ParticleReal * const Bzp,
                    amrex::FArrayBox const * const exfab,
                    amrex::FArrayBox const * const eyfab,
                    amrex::FArrayBox const * const ezfab,
                    amrex::FArrayBox const * const bxfab,
                    amrex::FArrayBox const * const byfab,
                    amrex::FArrayBox const * const bzfab,
                    const long np_to_gather,
                    const std::array<amrex::Real, 3>& dx,
                    const std::array<amrex::Real, 3> xyzmin,
                    const amrex::Dim3 lo,
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

    amrex::Array4<const amrex::Real> const& ex_arr = exfab->array();
    amrex::Array4<const amrex::Real> const& ey_arr = eyfab->array();
    amrex::Array4<const amrex::Real> const& ez_arr = ezfab->array();
    amrex::Array4<const amrex::Real> const& bx_arr = bxfab->array();
    amrex::Array4<const amrex::Real> const& by_arr = byfab->array();
    amrex::Array4<const amrex::Real> const& bz_arr = bzfab->array();

    amrex::IntVect const ex_type = exfab->box().type();
    amrex::IntVect const ey_type = eyfab->box().type();
    amrex::IntVect const ez_type = ezfab->box().type();
    amrex::IntVect const bx_type = bxfab->box().type();
    amrex::IntVect const by_type = byfab->box().type();
    amrex::IntVect const bz_type = bzfab->box().type();

    constexpr int zdir = (AMREX_SPACEDIM - 1);
    constexpr int NODE = amrex::IndexType::NODE;
    constexpr int CELL = amrex::IndexType::CELL;

    // Loop over particles and gather fields from
    // {e,b}{x,y,z}_arr to {E,B}{xyz}p.
    amrex::ParallelFor(
        np_to_gather,
        [=] AMREX_GPU_DEVICE (long ip) {

            amrex::ParticleReal xp, yp, zp;
            GetPosition(ip, xp, yp, zp);

            // --- Compute shape factors
            // x direction
            // Get particle position
#ifdef WARPX_DIM_RZ
            const amrex::Real rp = std::sqrt(xp*xp + yp*yp);
            const amrex::Real x = (rp - xmin)*dxi;
#else
            const amrex::Real x = (xp-xmin)*dxi;
#endif

            // j_[eb][xyz] leftmost grid point in x that the particle touches for the centering of each current
            // sx_[eb][xyz] shape factor along x for the centering of each current
            // There are only two possible centerings, node or cell centered, so at most only two shape factor
            // arrays will be needed.
            amrex::Real sx_node[depos_order + 1];
            amrex::Real sx_cell[depos_order + 1];
            amrex::Real sx_node_v[depos_order + 1 - lower_in_v];
            amrex::Real sx_cell_v[depos_order + 1 - lower_in_v];
            int j_node;
            int j_cell;
            int j_node_v;
            int j_cell_v;
            if ((ey_type[0] == NODE) || (ez_type[0] == NODE) || (bx_type[0] == NODE)) {
                j_node = compute_shape_factor<depos_order>(sx_node, x);
            }
            if ((ey_type[0] == CELL) || (ez_type[0] == CELL) || (bx_type[0] == CELL)) {
                j_cell = compute_shape_factor<depos_order>(sx_cell, x - 0.5);
            }
            if ((ex_type[0] == NODE) || (by_type[0] == NODE) || (bz_type[0] == NODE)) {
                j_node_v = compute_shape_factor<depos_order-lower_in_v>(sx_node_v, x);
            }
            if ((ex_type[0] == CELL) || (by_type[0] == CELL) || (bz_type[0] == CELL)) {
                j_cell_v = compute_shape_factor<depos_order-lower_in_v>(sx_cell_v, x - 0.5);
            }
            const amrex::Real (&sx_ex)[depos_order + 1 - lower_in_v] = ((ex_type[0] == NODE) ? sx_node_v : sx_cell_v);
            const amrex::Real (&sx_ey)[depos_order + 1             ] = ((ey_type[0] == NODE) ? sx_node   : sx_cell  );
            const amrex::Real (&sx_ez)[depos_order + 1             ] = ((ez_type[0] == NODE) ? sx_node   : sx_cell  );
            const amrex::Real (&sx_bx)[depos_order + 1             ] = ((bx_type[0] == NODE) ? sx_node   : sx_cell  );
            const amrex::Real (&sx_by)[depos_order + 1 - lower_in_v] = ((by_type[0] == NODE) ? sx_node_v : sx_cell_v);
            const amrex::Real (&sx_bz)[depos_order + 1 - lower_in_v] = ((bz_type[0] == NODE) ? sx_node_v : sx_cell_v);
            int const j_ex = ((ex_type[0] == NODE) ? j_node_v : j_cell_v);
            int const j_ey = ((ey_type[0] == NODE) ? j_node   : j_cell  );
            int const j_ez = ((ez_type[0] == NODE) ? j_node   : j_cell  );
            int const j_bx = ((bx_type[0] == NODE) ? j_node   : j_cell  );
            int const j_by = ((by_type[0] == NODE) ? j_node_v : j_cell_v);
            int const j_bz = ((bz_type[0] == NODE) ? j_node_v : j_cell_v);

#if (AMREX_SPACEDIM == 3)
            // y direction
            const amrex::Real y = (yp-ymin)*dyi;
            amrex::Real sy_node[depos_order + 1];
            amrex::Real sy_cell[depos_order + 1];
            amrex::Real sy_node_v[depos_order + 1 - lower_in_v];
            amrex::Real sy_cell_v[depos_order + 1 - lower_in_v];
            int k_node;
            int k_cell;
            int k_node_v;
            int k_cell_v;
            if ((ex_type[1] == NODE) || (ez_type[1] == NODE) || (by_type[1] == NODE)) {
                k_node = compute_shape_factor<depos_order>(sy_node, y);
            }
            if ((ex_type[1] == CELL) || (ez_type[1] == CELL) || (by_type[1] == CELL)) {
                k_cell = compute_shape_factor<depos_order>(sy_cell, y - 0.5);
            }
            if ((ey_type[1] == NODE) || (bx_type[1] == NODE) || (bz_type[1] == NODE)) {
                k_node_v = compute_shape_factor<depos_order-lower_in_v>(sy_node_v, y);
            }
            if ((ey_type[1] == CELL) || (bx_type[1] == CELL) || (bz_type[1] == CELL)) {
                k_cell_v = compute_shape_factor<depos_order-lower_in_v>(sy_cell_v, y - 0.5);
            }
            const amrex::Real (&sy_ex)[depos_order + 1             ] = ((ex_type[1] == NODE) ? sy_node   : sy_cell  );
            const amrex::Real (&sy_ey)[depos_order + 1 - lower_in_v] = ((ey_type[1] == NODE) ? sy_node_v : sy_cell_v);
            const amrex::Real (&sy_ez)[depos_order + 1             ] = ((ez_type[1] == NODE) ? sy_node   : sy_cell  );
            const amrex::Real (&sy_bx)[depos_order + 1 - lower_in_v] = ((bx_type[1] == NODE) ? sy_node_v : sy_cell_v);
            const amrex::Real (&sy_by)[depos_order + 1             ] = ((by_type[1] == NODE) ? sy_node   : sy_cell  );
            const amrex::Real (&sy_bz)[depos_order + 1 - lower_in_v] = ((bz_type[1] == NODE) ? sy_node_v : sy_cell_v);
            int const k_ex = ((ex_type[1] == NODE) ? k_node   : k_cell  );
            int const k_ey = ((ey_type[1] == NODE) ? k_node_v : k_cell_v);
            int const k_ez = ((ez_type[1] == NODE) ? k_node   : k_cell  );
            int const k_bx = ((bx_type[1] == NODE) ? k_node_v : k_cell_v);
            int const k_by = ((by_type[1] == NODE) ? k_node   : k_cell  );
            int const k_bz = ((bz_type[1] == NODE) ? k_node_v : k_cell_v);

#endif
            // z direction
            const amrex::Real z = (zp-zmin)*dzi;
            amrex::Real sz_node[depos_order + 1];
            amrex::Real sz_cell[depos_order + 1];
            amrex::Real sz_node_v[depos_order + 1 - lower_in_v];
            amrex::Real sz_cell_v[depos_order + 1 - lower_in_v];
            int l_node;
            int l_cell;
            int l_node_v;
            int l_cell_v;
            if ((ex_type[zdir] == NODE) || (ey_type[zdir] == NODE) || (bz_type[zdir] == NODE)) {
                l_node = compute_shape_factor<depos_order>(sz_node, z);
            }
            if ((ex_type[zdir] == CELL) || (ey_type[zdir] == CELL) || (bz_type[zdir] == CELL)) {
                l_cell = compute_shape_factor<depos_order>(sz_cell, z - 0.5);
            }
            if ((ez_type[zdir] == NODE) || (bx_type[zdir] == NODE) || (by_type[zdir] == NODE)) {
                l_node_v = compute_shape_factor<depos_order-lower_in_v>(sz_node_v, z);
            }
            if ((ez_type[zdir] == CELL) || (bx_type[zdir] == CELL) || (by_type[zdir] == CELL)) {
                l_cell_v = compute_shape_factor<depos_order-lower_in_v>(sz_cell_v, z - 0.5);
            }
            const amrex::Real (&sz_ex)[depos_order + 1             ] = ((ex_type[zdir] == NODE) ? sz_node   : sz_cell  );
            const amrex::Real (&sz_ey)[depos_order + 1             ] = ((ey_type[zdir] == NODE) ? sz_node   : sz_cell  );
            const amrex::Real (&sz_ez)[depos_order + 1 - lower_in_v] = ((ez_type[zdir] == NODE) ? sz_node_v : sz_cell_v);
            const amrex::Real (&sz_bx)[depos_order + 1 - lower_in_v] = ((bx_type[zdir] == NODE) ? sz_node_v : sz_cell_v);
            const amrex::Real (&sz_by)[depos_order + 1 - lower_in_v] = ((by_type[zdir] == NODE) ? sz_node_v : sz_cell_v);
            const amrex::Real (&sz_bz)[depos_order + 1             ] = ((bz_type[zdir] == NODE) ? sz_node   : sz_cell  );
            int const l_ex = ((ex_type[zdir] == NODE) ? l_node   : l_cell  );
            int const l_ey = ((ey_type[zdir] == NODE) ? l_node   : l_cell  );
            int const l_ez = ((ez_type[zdir] == NODE) ? l_node_v : l_cell_v);
            int const l_bx = ((bx_type[zdir] == NODE) ? l_node_v : l_cell_v);
            int const l_by = ((by_type[zdir] == NODE) ? l_node_v : l_cell_v);
            int const l_bz = ((bz_type[zdir] == NODE) ? l_node   : l_cell  );


            // Each field is gathered in a separate block of
            // AMREX_SPACEDIM nested loops because the deposition
            // order can differ for each component of each field
            // when lower_in_v is set to 1
#if (AMREX_SPACEDIM == 2)
            // Gather field on particle Eyp[i] from field on grid ey_arr
            for (int iz=0; iz<=depos_order; iz++){
                for (int ix=0; ix<=depos_order; ix++){
                    Eyp[ip] += sx_ey[ix]*sz_ey[iz]*
                        ey_arr(lo.x+j_ey+ix, lo.y+l_ey+iz, 0, 0);
                }
            }
            // Gather field on particle Exp[i] from field on grid ex_arr
            // Gather field on particle Bzp[i] from field on grid bz_arr
            for (int iz=0; iz<=depos_order; iz++){
                for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                    Exp[ip] += sx_ex[ix]*sz_ex[iz]*
                        ex_arr(lo.x+j_ex+ix, lo.y+l_ex+iz, 0, 0);
                    Bzp[ip] += sx_bz[ix]*sz_bz[iz]*
                        bz_arr(lo.x+j_bz+ix, lo.y+l_bz+iz, 0, 0);
                }
            }
            // Gather field on particle Ezp[i] from field on grid ez_arr
            // Gather field on particle Bxp[i] from field on grid bx_arr
            for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                for (int ix=0; ix<=depos_order; ix++){
                    Ezp[ip] += sx_ez[ix]*sz_ez[iz]*
                        ez_arr(lo.x+j_ez+ix, lo.y+l_ez+iz, 0, 0);
                    Bxp[ip] += sx_bx[ix]*sz_bx[iz]*
                        bx_arr(lo.x+j_bx+ix, lo.y+l_bx+iz, 0, 0);
                }
            }
            // Gather field on particle Byp[i] from field on grid by_arr
            for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                    Byp[ip] += sx_by[ix]*sz_by[iz]*
                        by_arr(lo.x+j_by+ix, lo.y+l_by+iz, 0, 0);
                }
            }

#ifdef WARPX_DIM_RZ

            amrex::Real costheta;
            amrex::Real sintheta;
            if (rp > 0.) {
                costheta = xp/rp;
                sintheta = yp/rp;
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
                        const amrex::Real dEy = (+ ey_arr(lo.x+j_ey+ix, lo.y+l_ey+iz, 0, 2*imode-1)*xy.real()
                                                 - ey_arr(lo.x+j_ey+ix, lo.y+l_ey+iz, 0, 2*imode)*xy.imag());
                        Eyp[ip] += sx_ey[ix]*sz_ey[iz]*dEy;
                    }
                }
                // Gather field on particle Exp[i] from field on grid ex_arr
                // Gather field on particle Bzp[i] from field on grid bz_arr
                for (int iz=0; iz<=depos_order; iz++){
                    for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                        const amrex::Real dEx = (+ ex_arr(lo.x+j_ex+ix, lo.y+l_ex+iz, 0, 2*imode-1)*xy.real()
                                                 - ex_arr(lo.x+j_ex+ix, lo.y+l_ex+iz, 0, 2*imode)*xy.imag());
                        Exp[ip] += sx_ex[ix]*sz_ex[iz]*dEx;
                        const amrex::Real dBz = (+ bz_arr(lo.x+j_bz+ix, lo.y+l_bz+iz, 0, 2*imode-1)*xy.real()
                                                 - bz_arr(lo.x+j_bz+ix, lo.y+l_bz+iz, 0, 2*imode)*xy.imag());
                        Bzp[ip] += sx_bz[ix]*sz_bz[iz]*dBz;
                    }
                }
                // Gather field on particle Ezp[i] from field on grid ez_arr
                // Gather field on particle Bxp[i] from field on grid bx_arr
                for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                    for (int ix=0; ix<=depos_order; ix++){
                        const amrex::Real dEz = (+ ez_arr(lo.x+j_ez+ix, lo.y+l_ez+iz, 0, 2*imode-1)*xy.real()
                                                 - ez_arr(lo.x+j_ez+ix, lo.y+l_ez+iz, 0, 2*imode)*xy.imag());
                        Ezp[ip] += sx_ez[ix]*sz_ez[iz]*dEz;
                        const amrex::Real dBx = (+ bx_arr(lo.x+j_bx+ix, lo.y+l_bx+iz, 0, 2*imode-1)*xy.real()
                                                 - bx_arr(lo.x+j_bx+ix, lo.y+l_bx+iz, 0, 2*imode)*xy.imag());
                        Bxp[ip] += sx_bx[ix]*sz_bx[iz]*dBx;
                    }
                }
                // Gather field on particle Byp[i] from field on grid by_arr
                for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                    for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                        const amrex::Real dBy = (+ by_arr(lo.x+j_by+ix, lo.y+l_by+iz, 0, 2*imode-1)*xy.real()
                                                 - by_arr(lo.x+j_by+ix, lo.y+l_by+iz, 0, 2*imode)*xy.imag());
                        Byp[ip] += sx_by[ix]*sz_by[iz]*dBy;
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
                        Exp[ip] += sx_ex[ix]*sy_ex[iy]*sz_ex[iz]*
                            ex_arr(lo.x+j_ex+ix, lo.y+k_ex+iy, lo.z+l_ex+iz);
                    }
                }
            }
            // Gather field on particle Eyp[i] from field on grid ey_arr
            for (int iz=0; iz<=depos_order; iz++){
                for (int iy=0; iy<=depos_order-lower_in_v; iy++){
                    for (int ix=0; ix<=depos_order; ix++){
                        Eyp[ip] += sx_ey[ix]*sy_ey[iy]*sz_ey[iz]*
                            ey_arr(lo.x+j_ey+ix, lo.y+k_ey+iy, lo.z+l_ey+iz);
                    }
                }
            }
            // Gather field on particle Ezp[i] from field on grid ez_arr
            for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                for (int iy=0; iy<=depos_order; iy++){
                    for (int ix=0; ix<=depos_order; ix++){
                        Ezp[ip] += sx_ez[ix]*sy_ez[iy]*sz_ez[iz]*
                            ez_arr(lo.x+j_ez+ix, lo.y+k_ez+iy, lo.z+l_ez+iz);
                    }
                }
            }
            // Gather field on particle Bzp[i] from field on grid bz_arr
            for (int iz=0; iz<=depos_order; iz++){
                for (int iy=0; iy<=depos_order-lower_in_v; iy++){
                    for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                        Bzp[ip] += sx_bz[ix]*sy_bz[iy]*sz_bz[iz]*
                            bz_arr(lo.x+j_bz+ix, lo.y+k_bz+iy, lo.z+l_bz+iz);
                    }
                }
            }
            // Gather field on particle Byp[i] from field on grid by_arr
            for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                for (int iy=0; iy<=depos_order; iy++){
                    for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                        Byp[ip] += sx_by[ix]*sy_by[iy]*sz_by[iz]*
                            by_arr(lo.x+j_by+ix, lo.y+k_by+iy, lo.z+l_by+iz);
                    }
                }
            }
            // Gather field on particle Bxp[i] from field on grid bx_arr
            for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                for (int iy=0; iy<=depos_order-lower_in_v; iy++){
                    for (int ix=0; ix<=depos_order; ix++){
                        Bxp[ip] += sx_bx[ix]*sy_bx[iy]*sz_bx[iz]*
                            bx_arr(lo.x+j_bx+ix, lo.y+k_bx+iy, lo.z+l_bx+iz);
                    }
                }
            }
#endif
        }
        );
}

#endif // FIELDGATHER_H_
