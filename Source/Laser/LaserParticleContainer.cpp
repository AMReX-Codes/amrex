/* Copyright 2019-2020 Andrew Myers, Axel Huebl, David Grote
 * Luca Fedeli, Maxence Thevenet, Remi Lehe
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpX_Complex.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/Pusher/GetAndSetPosition.H"

#include <limits>
#include <cmath>
#include <algorithm>
#include <numeric>


using namespace amrex;
using namespace WarpXLaserProfiles;

namespace
{
    Vector<Real> CrossProduct (const Vector<Real>& a, const Vector<Real>& b)
    {
        return { a[1]*b[2]-a[2]*b[1],  a[2]*b[0]-a[0]*b[2],  a[0]*b[1]-a[1]*b[0] };
    }
}

LaserParticleContainer::LaserParticleContainer (AmrCore* amr_core, int ispecies, const std::string& name)
    : WarpXParticleContainer(amr_core, ispecies),
      laser_name(name)
{
    charge = 1.0;
    mass = std::numeric_limits<Real>::max();
    do_back_transformed_diagnostics = 0;

    ParmParse pp(laser_name);

    // Parse the type of laser profile and set the corresponding flag `profile`
    std::string laser_type_s;
    pp.get("profile", laser_type_s);
    std::transform(laser_type_s.begin(), laser_type_s.end(), laser_type_s.begin(), ::tolower);

    // Parse the properties of the antenna
    pp.getarr("position", position);
    pp.getarr("direction", nvec);
    pp.getarr("polarization", p_X);

    pp.query("pusher_algo", pusher_algo);
    pp.get("wavelength", wavelength);
    pp.get("e_max", e_max);
    pp.query("do_continuous_injection", do_continuous_injection);
    pp.query("min_particles_per_mode", min_particles_per_mode);

    //Select laser profile
    if(laser_profiles_dictionary.count(laser_type_s) == 0){
        amrex::Abort(std::string("Unknown laser type: ").append(laser_type_s));
    }
    m_up_laser_profile = laser_profiles_dictionary.at(laser_type_s)();
    //__________

    // Plane normal
    Real s = 1.0_rt / std::sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1] + nvec[2]*nvec[2]);
    nvec = { nvec[0]*s, nvec[1]*s, nvec[2]*s };

    if (WarpX::gamma_boost > 1.) {
        // Check that the laser direction is equal to the boost direction
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(  nvec[0]*WarpX::boost_direction[0]
                                         + nvec[1]*WarpX::boost_direction[1]
                                         + nvec[2]*WarpX::boost_direction[2] - 1. < 1.e-12,
                                           "The Lorentz boost should be in the same direction as the laser propagation");
        // Get the position of the plane, along the boost direction, in the lab frame
        // and convert the position of the antenna to the boosted frame
        Z0_lab = nvec[0]*position[0] + nvec[1]*position[1] + nvec[2]*position[2];
        Real Z0_boost = Z0_lab/WarpX::gamma_boost;
        position[0] += (Z0_boost-Z0_lab)*nvec[0];
        position[1] += (Z0_boost-Z0_lab)*nvec[1];
        position[2] += (Z0_boost-Z0_lab)*nvec[2];
    }

    // The first polarization vector
    s = 1.0_rt / std::sqrt(p_X[0]*p_X[0] + p_X[1]*p_X[1] + p_X[2]*p_X[2]);
    p_X = { p_X[0]*s, p_X[1]*s, p_X[2]*s };

    Real const dp = std::inner_product(nvec.begin(), nvec.end(), p_X.begin(), 0.0);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(dp) < 1.0e-14,
        "Laser plane vector is not perpendicular to the main polarization vector");

    p_Y = CrossProduct(nvec, p_X);   // The second polarization vector

#if (defined WARPX_DIM_3D) || (defined WARPX_DIM_RZ)
    u_X = p_X;
    u_Y = p_Y;
#else
    u_X = CrossProduct({0., 1., 0.}, nvec);
    u_Y = {0., 1., 0.};
#endif

    laser_injection_box= Geom(0).ProbDomain();
    {
        Vector<Real> lo, hi;
        if (pp.queryarr("prob_lo", lo)) {
            laser_injection_box.setLo(lo);
        }
        if (pp.queryarr("prob_hi", hi)) {
            laser_injection_box.setHi(hi);
        }
    }

    if (do_continuous_injection){
        // If laser antenna initially outside of the box, store its theoretical
        // position in z_antenna_th
        updated_position = position;

        // Sanity checks
        int dir = WarpX::moving_window_dir;
        std::vector<Real> windir(3, 0.0);
#if (AMREX_SPACEDIM==2)
        windir[2*dir] = 1.0;
#else
        windir[dir] = 1.0;
#endif
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            (nvec[0]-windir[0]) + (nvec[1]-windir[1]) + (nvec[2]-windir[2])
            < 1.e-12, "do_continous_injection for laser particle only works" +
            " if moving window direction and laser propagation direction are the same");
        if ( WarpX::gamma_boost>1 ){
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                (WarpX::boost_direction[0]-0)*(WarpX::boost_direction[0]-0) +
                (WarpX::boost_direction[1]-0)*(WarpX::boost_direction[1]-0) +
                (WarpX::boost_direction[2]-1)*(WarpX::boost_direction[2]-1) < 1.e-12,
                "do_continous_injection for laser particle only works if " +
                "warpx.boost_direction = z. TODO: all directions.");
        }
    }

    //Init laser profile

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(e_max > 0.,
        "Laser amplitude (e_max) must be positive.");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(wavelength > 0.,
        "Laser wavelength must be positive.");

    CommonLaserParameters common_params;
    common_params.wavelength = wavelength;
    common_params.e_max = e_max;
    common_params.p_X = p_X;
    common_params.nvec = nvec;
    m_up_laser_profile->init(pp, ParmParse{"my_constants"}, common_params);
}

/* \brief Check if laser particles enter the box, and inject if necessary.
 * \param injection_box: a RealBox where particles should be injected.
 */
void
LaserParticleContainer::ContinuousInjection (const RealBox& injection_box)
{
    // Input parameter injection_box contains small box where injection
    // should occur.
    // So far, LaserParticleContainer::laser_injection_box contains the
    // outdated full problem domain at t=0.

    // Convert updated_position to Real* to use RealBox::contains().
#if (AMREX_SPACEDIM == 3)
    const Real* p_pos = updated_position.dataPtr();
#else
    const Real p_pos[2] = {updated_position[0], updated_position[2]};
#endif
    if ( injection_box.contains(p_pos) ){
        // Update laser_injection_box with current value
        laser_injection_box = injection_box;
        // Inject laser particles. LaserParticleContainer::InitData
        // is called only once, when the antenna enters the simulation
        // domain.
        InitData();
    }
}

/* \brief update position of the antenna if running in boosted frame.
 * \param dt time step (level 0).
 * The up-to-date antenna position is stored in updated_position.
 */
void
LaserParticleContainer::UpdateContinuousInjectionPosition (Real dt)
{
    int dir = WarpX::moving_window_dir;
    if (do_continuous_injection and (WarpX::gamma_boost > 1)){
        // In boosted-frame simulations, the antenna has moved since the last
        // call to this function, and injection position needs to be updated
#if ( AMREX_SPACEDIM == 3 )
        updated_position[dir] -= WarpX::beta_boost *
            WarpX::boost_direction[dir] * PhysConst::c * dt;
#elif ( AMREX_SPACEDIM == 2 )
        // In 2D, dir=0 corresponds to x and dir=1 corresponds to z
        // This needs to be converted in order to index `boost_direction`
        // which has 3 components, for both 2D and 3D simulations.
        updated_position[2*dir] -= WarpX::beta_boost *
            WarpX::boost_direction[2*dir] * PhysConst::c * dt;
#endif
    }
}

void
LaserParticleContainer::InitData ()
{
    // Call InitData on max level to inject one laser particle per
    // finest cell.
    InitData(maxLevel());
}

void
LaserParticleContainer::InitData (int lev)
{
    // spacing of laser particles in the laser plane.
    // has to be done after geometry is set up.
    Real S_X, S_Y;
    ComputeSpacing(lev, S_X, S_Y);
    ComputeWeightMobility(S_X, S_Y);

    // LaserParticleContainer::position contains the initial position of the
    // laser antenna. In the boosted frame, the antenna is moving.
    // Update its position with updated_position.
    if (do_continuous_injection){
        position = updated_position;
    }

    auto Transform = [&](int const i, int const j) -> Vector<Real>{
#if (AMREX_SPACEDIM == 3)
        return { position[0] + (S_X*(Real(i)+0.5_rt))*u_X[0] + (S_Y*(Real(j)+0.5_rt))*u_Y[0],
                 position[1] + (S_X*(Real(i)+0.5_rt))*u_X[1] + (S_Y*(Real(j)+0.5_rt))*u_Y[1],
                 position[2] + (S_X*(Real(i)+0.5_rt))*u_X[2] + (S_Y*(Real(j)+0.5_rt))*u_Y[2] };
#else
#   if (defined WARPX_DIM_RZ)
        return { position[0] + (S_X*(Real(i)+0.5)),
                 0.0,
                 position[2]};
#   else
        return { position[0] + (S_X*(Real(i)+0.5))*u_X[0],
                 0.0,
                 position[2] + (S_X*(Real(i)+0.5))*u_X[2] };
#   endif
#endif
    };

    // Given the "lab" frame coordinates, return the real coordinates in the laser plane coordinates
    auto InverseTransform = [&](const Vector<Real>& pos) -> Vector<Real>{
#if (AMREX_SPACEDIM == 3)
        return {u_X[0]*(pos[0]-position[0])+u_X[1]*(pos[1]-position[1])+u_X[2]*(pos[2]-position[2]),
                u_Y[0]*(pos[0]-position[0])+u_Y[1]*(pos[1]-position[1])+u_Y[2]*(pos[2]-position[2])};
#else
#   if (defined WARPX_DIM_RZ)
        return {pos[0]-position[0], 0.0};
#   else
        return {u_X[0]*(pos[0]-position[0])+u_X[2]*(pos[2]-position[2]), 0.0};
#   endif
#endif
    };

    Vector<int> plane_lo(2, std::numeric_limits<int>::max());
    Vector<int> plane_hi(2, std::numeric_limits<int>::min());
    {
        auto compute_min_max = [&](Real x, Real y, Real z){
            const Vector<Real>& pos_plane = InverseTransform({x, y, z});
            auto i = static_cast<int>(pos_plane[0]/S_X);
            auto j = static_cast<int>(pos_plane[1]/S_Y);
            plane_lo[0] = std::min(plane_lo[0], i);
            plane_lo[1] = std::min(plane_lo[1], j);
            plane_hi[0] = std::max(plane_hi[0], i);
            plane_hi[1] = std::max(plane_hi[1], j);
        };

        const Real* prob_lo = laser_injection_box.lo();
        const Real* prob_hi = laser_injection_box.hi();
#if (AMREX_SPACEDIM == 3)
        compute_min_max(prob_lo[0], prob_lo[1], prob_lo[2]);
        compute_min_max(prob_hi[0], prob_lo[1], prob_lo[2]);
        compute_min_max(prob_lo[0], prob_hi[1], prob_lo[2]);
        compute_min_max(prob_hi[0], prob_hi[1], prob_lo[2]);
        compute_min_max(prob_lo[0], prob_lo[1], prob_hi[2]);
        compute_min_max(prob_hi[0], prob_lo[1], prob_hi[2]);
        compute_min_max(prob_lo[0], prob_hi[1], prob_hi[2]);
        compute_min_max(prob_hi[0], prob_hi[1], prob_hi[2]);
#else
        compute_min_max(prob_lo[0], 0.0, prob_lo[1]);
        compute_min_max(prob_hi[0], 0.0, prob_lo[1]);
        compute_min_max(prob_lo[0], 0.0, prob_hi[1]);
        compute_min_max(prob_hi[0], 0.0, prob_hi[1]);
#endif
    }

    const int nprocs = ParallelDescriptor::NProcs();
    const int myproc = ParallelDescriptor::MyProc();

#if (AMREX_SPACEDIM == 3)
    const Box plane_box {IntVect(plane_lo[0],plane_lo[1],0),
                         IntVect(plane_hi[0],plane_hi[1],0)};
    BoxArray plane_ba {plane_box};
    {
        IntVect chunk(plane_box.size());
        const int min_size = 8;
        while (plane_ba.size() < nprocs && chunk[0] > min_size && chunk[1] > min_size)
        {
            for (int j = 1; j >= 0 ; j--)
            {
                chunk[j] /= 2;

                if (plane_ba.size() < nprocs) {
                    plane_ba.maxSize(chunk);
                }
            }
        }
    }
#else
    BoxArray plane_ba { Box {IntVect(plane_lo[0],0), IntVect(plane_hi[0],0)} };
#endif

    RealVector particle_x, particle_y, particle_z, particle_w;

    const DistributionMapping plane_dm {plane_ba, nprocs};
    const Vector<int>& procmap = plane_dm.ProcessorMap();
    for (int i = 0, n = plane_ba.size(); i < n; ++i)
    {
        if (procmap[i] == myproc)
        {
            const Box& bx = plane_ba[i];
            for (IntVect cell = bx.smallEnd(); cell <= bx.bigEnd(); bx.next(cell))
            {
                const Vector<Real>& pos = Transform(cell[0], cell[1]);
#if (AMREX_SPACEDIM == 3)
                const Real* x = pos.dataPtr();
#else
                const Real x[2] = {pos[0], pos[2]};
#endif
                if (laser_injection_box.contains(x))
                {
#ifndef WARPX_DIM_RZ
                    for (int k = 0; k<2; ++k) {
                        particle_x.push_back(pos[0]);
                        particle_y.push_back(pos[1]);
                        particle_z.push_back(pos[2]);
                    }
                    particle_w.push_back( weight);
                    particle_w.push_back(-weight);
#else
                    // Particles are laid out in radial spokes
                    const int n_spokes = (WarpX::n_rz_azimuthal_modes - 1)*min_particles_per_mode;
                    for (int spoke = 0 ; spoke < n_spokes ; spoke++) {
                        const Real phase = 2.*MathConst::pi*spoke/n_spokes;
                        for (int k = 0; k<2; ++k) {
                            particle_x.push_back(pos[0]*std::cos(phase));
                            particle_y.push_back(pos[0]*std::sin(phase));
                            particle_z.push_back(pos[2]);
                        }
                        const Real r_weight = weight*2.*MathConst::pi*pos[0]/n_spokes;
                        particle_w.push_back( r_weight);
                        particle_w.push_back(-r_weight);
                    }
#endif
                }
            }
        }
    }
    const int np = particle_z.size();
    RealVector particle_ux(np, 0.0);
    RealVector particle_uy(np, 0.0);
    RealVector particle_uz(np, 0.0);

    if (Verbose()) amrex::Print() << "Adding laser particles\n";
    // Add particles on level 0. They will be redistributed afterwards
    AddNParticles(0,
                  np, particle_x.dataPtr(), particle_y.dataPtr(), particle_z.dataPtr(),
                  particle_ux.dataPtr(), particle_uy.dataPtr(), particle_uz.dataPtr(),
                  1, particle_w.dataPtr(), 1);
}

void
LaserParticleContainer::Evolve (int lev,
                                const MultiFab&, const MultiFab&, const MultiFab&,
                                const MultiFab&, const MultiFab&, const MultiFab&,
                                MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                MultiFab* cjx, MultiFab* cjy, MultiFab* cjz,
                                MultiFab* rho, MultiFab* crho,
                                const MultiFab*, const MultiFab*, const MultiFab*,
                                const MultiFab*, const MultiFab*, const MultiFab*,
                                Real t, Real dt, DtType /*a_dt_type*/)
{
    WARPX_PROFILE("Laser::Evolve()");
    WARPX_PROFILE_VAR_NS("Laser::Evolve::Copy", blp_copy);
    WARPX_PROFILE_VAR_NS("Laser::ParticlePush", blp_pp);
    WARPX_PROFILE_VAR_NS("Laser::CurrentDepo", blp_cd);
    WARPX_PROFILE_VAR_NS("Laser::Evolve::Accumulate", blp_accumulate);

    Real t_lab = t;
    if (WarpX::gamma_boost > 1) {
        // Convert time from the boosted to the lab-frame
        // (in order to later calculate the amplitude of the field,
        // at the position of the antenna, in the lab-frame)
        t_lab = 1./WarpX::gamma_boost*t + WarpX::beta_boost*Z0_lab/PhysConst::c;
    }

    // Update laser profile
    m_up_laser_profile->update(t);

    BL_ASSERT(OnSameGrids(lev,jx));

    MultiFab* cost = WarpX::getCosts(lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
        int const thread_num = omp_get_thread_num();
#else
        int const thread_num = 0;
#endif

        Gpu::ManagedDeviceVector<Real> plane_Xp, plane_Yp, amplitude_E;

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            Real wt = amrex::second();

            auto& attribs = pti.GetAttribs();

            auto&  wp = attribs[PIdx::w ];
            auto& uxp = attribs[PIdx::ux];
            auto& uyp = attribs[PIdx::uy];
            auto& uzp = attribs[PIdx::uz];

            const long np  = pti.numParticles();
            // For now, laser particles do not take the current buffers into account
            const long np_current = np;

            plane_Xp.resize(np);
            plane_Yp.resize(np);
            amplitude_E.resize(np);

            if (rho) {
                int* AMREX_RESTRICT ion_lev = nullptr;
                DepositCharge(pti, wp, ion_lev, rho, 0, 0,
                              np_current, thread_num, lev, lev);
                if (crho) {
                    DepositCharge(pti, wp, ion_lev, crho, 0, np_current,
                                  np-np_current, thread_num, lev, lev-1);
                }
            }

            //
            // Particle Push
            //
            WARPX_PROFILE_VAR_START(blp_pp);
            // Find the coordinates of the particles in the emission plane
            calculate_laser_plane_coordinates(pti, np,
                                              plane_Xp.dataPtr(),
                                              plane_Yp.dataPtr());

            // Calculate the laser amplitude to be emitted,
            // at the position of the emission plane
            m_up_laser_profile->fill_amplitude(
                np, plane_Xp.dataPtr(), plane_Yp.dataPtr(),
                t_lab, amplitude_E.dataPtr());

            // Calculate the corresponding momentum and position for the particles
            update_laser_particle(pti, np, uxp.dataPtr(), uyp.dataPtr(),
                                  uzp.dataPtr(), wp.dataPtr(),
                                  amplitude_E.dataPtr(), dt);
            WARPX_PROFILE_VAR_STOP(blp_pp);

            //
            // Current Deposition
            //
            // Deposit inside domains
            {
                int* ion_lev = nullptr;
                DepositCurrent(pti, wp, uxp, uyp, uzp, ion_lev, &jx, &jy, &jz,
                               0, np_current, thread_num,
                               lev, lev, dt);

                bool has_buffer = cjx;
                if (has_buffer){
                    // Deposit in buffers
                    DepositCurrent(pti, wp, uxp, uyp, uzp, ion_lev, cjx, cjy, cjz,
                                   np_current, np-np_current, thread_num,
                                   lev, lev-1, dt);
                }
            }

            if (rho) {
                int* AMREX_RESTRICT ion_lev = nullptr;
                DepositCharge(pti, wp, ion_lev, rho, 1, 0,
                              np_current, thread_num, lev, lev);
                if (crho) {
                    DepositCharge(pti, wp, ion_lev, crho, 1, np_current,
                                  np-np_current, thread_num, lev, lev-1);
                }
            }

            if (cost) {
                const Box& tbx = pti.tilebox();
                wt = (amrex::second() - wt) / tbx.d_numPts();
                Array4<Real> const& costarr = cost->array(pti);
                amrex::ParallelFor(tbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    costarr(i,j,k) += wt;
                });
            }
        }
    }
}

void
LaserParticleContainer::PostRestart ()
{
    Real Sx, Sy;
    const int lev = finestLevel();
    ComputeSpacing(lev, Sx, Sy);
    ComputeWeightMobility(Sx, Sy);
}

void
LaserParticleContainer::ComputeSpacing (int lev, Real& Sx, Real& Sy) const
{
    const std::array<Real,3>& dx = WarpX::CellSize(lev);

    const Real eps = dx[0]*1.e-50;
#if (AMREX_SPACEDIM == 3)
    Sx = std::min(std::min(dx[0]/(std::abs(u_X[0])+eps),
                           dx[1]/(std::abs(u_X[1])+eps)),
                           dx[2]/(std::abs(u_X[2])+eps));
    Sy = std::min(std::min(dx[0]/(std::abs(u_Y[0])+eps),
                           dx[1]/(std::abs(u_Y[1])+eps)),
                           dx[2]/(std::abs(u_Y[2])+eps));
#else
#   if (defined WARPX_DIM_RZ)
    Sx = dx[0];
#   else
    Sx = std::min(dx[0]/(std::abs(u_X[0])+eps),
                  dx[2]/(std::abs(u_X[2])+eps));
#   endif
    Sy = 1.0;
#endif
}

void
LaserParticleContainer::ComputeWeightMobility (Real Sx, Real Sy)
{
    constexpr Real eps = 0.01;
    constexpr Real fac = 1.0_rt / (2.0_rt * MathConst::pi * PhysConst::mu0 * PhysConst::c * PhysConst::c * eps);
    weight = fac * wavelength * Sx * Sy / std::min(Sx,Sy) * e_max;

    // The mobility is the constant of proportionality between the field to
    // be emitted, and the corresponding velocity that the particles need to have.
    mobility = (Sx * Sy)/(weight * PhysConst::mu0 * PhysConst::c * PhysConst::c);
    // When running in the boosted-frame, the input parameters (and in particular
    // the amplitude of the field) are given in the lab-frame.
    // Therefore, the mobility needs to be modified by a factor WarpX::gamma_boost.
    mobility = mobility/WarpX::gamma_boost;
}

void
LaserParticleContainer::PushP (int lev, Real dt,
                               const MultiFab&, const MultiFab&, const MultiFab&,
                               const MultiFab&, const MultiFab&, const MultiFab&)
{
    // I don't think we need to do anything.
}

/* \brief compute particles position in laser plane coordinate.
 *
 * \param np: number of laser particles
 * \param thread_num: thread number
 * \param pplane_Xp, pplane_Yp: pointers to arrays of particle positions
 * in laser plane coordinate.
 */
void
LaserParticleContainer::calculate_laser_plane_coordinates (const WarpXParIter& pti, const int np,
                                                           Real * AMREX_RESTRICT const pplane_Xp,
                                                           Real * AMREX_RESTRICT const pplane_Yp)
{
    const auto GetPosition = GetParticlePosition(pti);

    Real tmp_u_X_0 = u_X[0];
    Real tmp_u_X_2 = u_X[2];
    Real tmp_position_0 = position[0];
    Real tmp_position_2 = position[2];
#if (defined WARPX_DIM_3D) || (defined WARPX_DIM_RZ)
    Real tmp_u_X_1 = u_X[1];
    Real tmp_u_Y_0 = u_Y[0];
    Real tmp_u_Y_1 = u_Y[1];
    Real tmp_u_Y_2 = u_Y[2];
    Real tmp_position_1 = position[1];
#endif

    amrex::ParallelFor(
        np,
        [=] AMREX_GPU_DEVICE (int i) {
            ParticleReal x, y, z;
            GetPosition(i, x, y, z);
#if (defined WARPX_DIM_3D) || (defined WARPX_DIM_RZ)
            pplane_Xp[i] =
                tmp_u_X_0 * (x - tmp_position_0) +
                tmp_u_X_1 * (y - tmp_position_1) +
                tmp_u_X_2 * (z - tmp_position_2);
            pplane_Yp[i] =
                tmp_u_Y_0 * (x - tmp_position_0) +
                tmp_u_Y_1 * (y - tmp_position_1) +
                tmp_u_Y_2 * (z - tmp_position_2);
#elif (AMREX_SPACEDIM == 2)
            pplane_Xp[i] =
                tmp_u_X_0 * (x - tmp_position_0) +
                tmp_u_X_2 * (z - tmp_position_2);
            pplane_Yp[i] = 0.;
#endif
        }
        );
}

/* \brief push laser particles, in simulation coordinates.
 *
 * \param pti: Particle iterator
 * \param np: number of laser particles
 * \param puxp, puyp, puzp: pointers to arrays of particle momenta.
 * \param pwp: pointer to array of particle weights.
 * \param amplitude: Electric field amplitude at the position of each particle.
 * \param dt: time step.
 */
void
LaserParticleContainer::update_laser_particle (WarpXParIter& pti,
                                               const int np,
                                               ParticleReal * AMREX_RESTRICT const puxp,
                                               ParticleReal * AMREX_RESTRICT const puyp,
                                               ParticleReal * AMREX_RESTRICT const puzp,
                                               ParticleReal const * AMREX_RESTRICT const pwp,
                                               Real const * AMREX_RESTRICT const amplitude,
                                               const Real dt)
{
    const auto GetPosition = GetParticlePosition(pti);
    auto       SetPosition = SetParticlePosition(pti);

    Real tmp_p_X_0 = p_X[0];
    Real tmp_p_X_1 = p_X[1];
    Real tmp_p_X_2 = p_X[2];
    Real tmp_nvec_0 = nvec[0];
    Real tmp_nvec_1 = nvec[1];
    Real tmp_nvec_2 = nvec[2];

    // Copy member variables to tmp copies for GPU runs.
    Real tmp_mobility = mobility;
    Real gamma_boost = WarpX::gamma_boost;
    Real beta_boost = WarpX::beta_boost;
    amrex::ParallelFor(
        np,
        [=] AMREX_GPU_DEVICE (int i) {
            // Calculate the velocity according to the amplitude of E
            const Real sign_charge = (pwp[i]>0) ? 1 : -1;
            const Real v_over_c = sign_charge * tmp_mobility * amplitude[i];
            // The velocity is along the laser polarization p_X
            Real vx = PhysConst::c * v_over_c * tmp_p_X_0;
            Real vy = PhysConst::c * v_over_c * tmp_p_X_1;
            Real vz = PhysConst::c * v_over_c * tmp_p_X_2;
            // When running in the boosted-frame, their is additional velocity along nvec
            if (gamma_boost > 1.){
                vx -= PhysConst::c * beta_boost * tmp_nvec_0;
                vy -= PhysConst::c * beta_boost * tmp_nvec_1;
                vz -= PhysConst::c * beta_boost * tmp_nvec_2;
            }
            // Get the corresponding momenta
            const Real gamma = gamma_boost/std::sqrt(1. - v_over_c*v_over_c);
            puxp[i] = gamma * vx;
            puyp[i] = gamma * vy;
            puzp[i] = gamma * vz;

            // Push the the particle positions
            ParticleReal x, y, z;
            GetPosition(i, x, y, z);
            x += vx * dt;
#if (defined WARPX_DIM_3D) || (defined WARPX_DIM_RZ)
            y += vy * dt;
#endif
            z += vz * dt;
            SetPosition(i, x, y, z);
        }
        );
}
