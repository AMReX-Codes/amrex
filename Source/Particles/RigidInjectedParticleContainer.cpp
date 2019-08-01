#include <limits>
#include <sstream>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <RigidInjectedParticleContainer.H>
#include <WarpX_f.H>
#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpXAlgorithmSelection.H>
#include <UpdateMomentumBoris.H>
#include <UpdateMomentumVay.H>

using namespace amrex;

RigidInjectedParticleContainer::RigidInjectedParticleContainer (AmrCore* amr_core, int ispecies,
                                                                const std::string& name)
    : PhysicalParticleContainer(amr_core, ispecies, name)
{

    ParmParse pp(species_name);

    pp.get("zinject_plane", zinject_plane);
    pp.query("projected", projected);
    pp.query("focused", focused);
    pp.query("rigid_advance", rigid_advance);

}

void RigidInjectedParticleContainer::InitData()
{
    done_injecting.resize(finestLevel()+1, 0);
    zinject_plane_levels.resize(finestLevel()+1, zinject_plane/WarpX::gamma_boost);

    AddParticles(0); // Note - add on level 0

    // Particles added by AddParticles should already be in the boosted frame
    RemapParticles();

    Redistribute();  // We then redistribute
}

void
RigidInjectedParticleContainer::RemapParticles()
{

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(projected, "ERROR: projected = false is not supported with this particle loading");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!focused, "ERROR: focused = true is not supported with this particle loading");

    // For rigid_advance == false, nothing needs to be done

    if (rigid_advance) {

        // The particle z positions are adjusted to account for the difference between
        // advancing with vzbar and wih vz[i] before injection

        // For now, start with the assumption that this will only happen
        // at the start of the simulation.
        const Real t_lab = 0.;

        const Real uz_boost = WarpX::gamma_boost*WarpX::beta_boost*PhysConst::c;
        const Real csq = PhysConst::c*PhysConst::c;

        vzbeam_ave_boosted = meanParticleVelocity(false)[2];

        for (int lev = 0; lev <= finestLevel(); lev++) {

#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                // Get the average beam velocity in the boosted frame.
                // Note that the particles are already in the boosted frame.
                // This value is saved to advance the particles not injected yet

                Cuda::ManagedDeviceVector<Real> xp, yp, zp;

                for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
                {
                    auto& attribs = pti.GetAttribs();
                    auto& uxp = attribs[PIdx::ux];
                    auto& uyp = attribs[PIdx::uy];
                    auto& uzp = attribs[PIdx::uz];

                    // Copy data from particle container to temp arrays
                    pti.GetPosition(xp, yp, zp);

                    // Loop over particles
                    const long np = pti.numParticles();
                    for (int i=0 ; i < np ; i++) {

                        const Real gammapr = std::sqrt(1. + (uxp[i]*uxp[i] + uyp[i]*uyp[i] + uzp[i]*uzp[i])/csq);
                        const Real vzpr = uzp[i]/gammapr;

                        // Back out the value of z_lab
                        const Real z_lab = (zp[i] + uz_boost*t_lab + WarpX::gamma_boost*t_lab*vzpr)/(WarpX::gamma_boost + uz_boost*vzpr/csq);

                        // Time of the particle in the boosted frame given its position in the lab frame at t=0.
                        const Real tpr = WarpX::gamma_boost*t_lab - uz_boost*z_lab/csq;

                        // Adjust the position, taking away its motion from its own velocity and adding
                        // the motion from the average velocity
                        zp[i] = zp[i] + tpr*vzpr - tpr*vzbeam_ave_boosted;

                    }

                    // Copy the data back to the particle container
                    pti.SetPosition(xp, yp, zp);

                }
            }
        }
    }
}


void
RigidInjectedParticleContainer::BoostandRemapParticles()
{

    // Boost the particles into the boosted frame and map the particles
    // to the t=0 in the boosted frame. If using rigid_advance, the z position
    // is adjusted using vzbar, otherwise using vz[i]

    if (rigid_advance) {
        // Get the average beam velocity in the boosted frame
        // This value is saved to advance the particles not injected yet
        const Real vzbeam_ave_lab = meanParticleVelocity(false)[2];
        vzbeam_ave_boosted = (vzbeam_ave_lab - WarpX::beta_boost*PhysConst::c)/(1. - vzbeam_ave_lab*WarpX::beta_boost/PhysConst::c);
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Cuda::ManagedDeviceVector<Real> xp, yp, zp;

        for (WarpXParIter pti(*this, 0); pti.isValid(); ++pti)
        {

            auto& attribs = pti.GetAttribs();
            auto& uxp = attribs[PIdx::ux];
            auto& uyp = attribs[PIdx::uy];
            auto& uzp = attribs[PIdx::uz];

            // Copy data from particle container to temp arrays
            pti.GetPosition(xp, yp, zp);

            // Loop over particles
            const long np = pti.numParticles();
            for (int i=0 ; i < np ; i++) {

                const Real gamma_lab = std::sqrt(1. + (uxp[i]*uxp[i] + uyp[i]*uyp[i] + uzp[i]*uzp[i])/(PhysConst::c*PhysConst::c));

                const Real vx_lab = uxp[i]/gamma_lab;
                const Real vy_lab = uyp[i]/gamma_lab;
                const Real vz_lab = uzp[i]/gamma_lab;

                // t0_lab is the time in the lab frame that the particles reaches z=0
                // The location and time (z=0, t=0) is a synchronization point between the
                // lab and boosted frames.
                const Real t0_lab = -zp[i]/vz_lab;

                if (!projected) {
                    xp[i] += t0_lab*vx_lab;
                    yp[i] += t0_lab*vy_lab;
                }
                if (focused) {
                    // Correct for focusing effect from shift from z=0 to zinject
                    const Real tfocus = -zinject_plane*WarpX::gamma_boost/vz_lab;
                    xp[i] -= tfocus*vx_lab;
                    yp[i] -= tfocus*vy_lab;
                }

                // Time of the particle in the boosted frame given its position in the lab frame at t=0.
                const Real tpr = -WarpX::gamma_boost*WarpX::beta_boost*zp[i]/PhysConst::c;

                // Position of the particle in the boosted frame given its position in the lab frame at t=0.
                const Real zpr = WarpX::gamma_boost*zp[i];

                // Momentum of the particle in the boosted frame (assuming that it is fixed).
                uzp[i] = WarpX::gamma_boost*(uzp[i] - WarpX::beta_boost*PhysConst::c*gamma_lab);

                // Put the particle at the location in the boosted frame at boost frame t=0,
                if (rigid_advance) {
                    // with the particle moving at the average velocity
                    zp[i] = zpr - vzbeam_ave_boosted*tpr;
                }
                else {
                    // with the particle moving with its own velocity
                    const Real gammapr = std::sqrt(1. + (uxp[i]*uxp[i] + uyp[i]*uyp[i] + uzp[i]*uzp[i])/(PhysConst::c*PhysConst::c));
                    const Real vzpr = uzp[i]/gammapr;
                    zp[i] = zpr - vzpr*tpr;
                }

            }

            // Copy the data back to the particle container
            pti.SetPosition(xp, yp, zp);

        }
    }
}

void
RigidInjectedParticleContainer::PushPX(WarpXParIter& pti,
                                       Cuda::ManagedDeviceVector<Real>& xp,
                                       Cuda::ManagedDeviceVector<Real>& yp,
                                       Cuda::ManagedDeviceVector<Real>& zp,
                                       Cuda::ManagedDeviceVector<Real>& giv,
                                       Real dt)
{

    // This wraps the momentum and position advance so that inheritors can modify the call.
    auto& attribs = pti.GetAttribs();
    auto& uxp = attribs[PIdx::ux];
    auto& uyp = attribs[PIdx::uy];
    auto& uzp = attribs[PIdx::uz];

    // Save the position and momenta, making copies
    Cuda::ManagedDeviceVector<Real> xp_save, yp_save, zp_save;
    RealVector uxp_save, uyp_save, uzp_save;

    Real* const AMREX_RESTRICT x = xp.dataPtr();
    Real* const AMREX_RESTRICT y = yp.dataPtr();
    Real* const AMREX_RESTRICT z = zp.dataPtr();
    Real* const AMREX_RESTRICT gi = giv.dataPtr();
    Real* const AMREX_RESTRICT ux = uxp.dataPtr();
    Real* const AMREX_RESTRICT uy = uyp.dataPtr();
    Real* const AMREX_RESTRICT uz = uzp.dataPtr();
    Real* const AMREX_RESTRICT Exp = attribs[PIdx::Ex].dataPtr();
    Real* const AMREX_RESTRICT Eyp = attribs[PIdx::Ey].dataPtr();
    Real* const AMREX_RESTRICT Ezp = attribs[PIdx::Ez].dataPtr();
    Real* const AMREX_RESTRICT Bxp = attribs[PIdx::Bx].dataPtr();
    Real* const AMREX_RESTRICT Byp = attribs[PIdx::By].dataPtr();
    Real* const AMREX_RESTRICT Bzp = attribs[PIdx::Bz].dataPtr();

    if (!done_injecting_lev) {
        if (!(WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags)) {
            // If the old values are not already saved, create copies here.
            xp_save = xp;
            yp_save = yp;
            zp_save = zp;
            uxp_save = uxp;
            uyp_save = uyp;
            uzp_save = uzp;
        }

        // Scale the fields of particles about to cross the injection plane.
        // This only approximates what should be happening. The particles
        // should by advanced a fraction of a time step instead.
        // Scaling the fields is much easier and may be good enough.
        const Real v_boost = WarpX::beta_boost*PhysConst::c;
        const Real z_plane_previous = zinject_plane_lev_previous;
        const Real vz_ave_boosted = vzbeam_ave_boosted;
        amrex::ParallelFor( pti.numParticles(),
            [=] AMREX_GPU_DEVICE (long i) {
            const Real dtscale = dt - (z_plane_previous - z[i])/(vz_ave_boosted + v_boost);
            if (0. < dtscale && dtscale < dt) {
                Exp[i] *= dtscale;
                Eyp[i] *= dtscale;
                Ezp[i] *= dtscale;
                Bxp[i] *= dtscale;
                Byp[i] *= dtscale;
                Bzp[i] *= dtscale;
            }
        }
        );
    }

    PhysicalParticleContainer::PushPX(pti, xp, yp, zp, giv, dt);

    if (!done_injecting_lev) {

        Real* AMREX_RESTRICT x_save;
        Real* AMREX_RESTRICT y_save;
        Real* AMREX_RESTRICT z_save;
        Real* AMREX_RESTRICT ux_save;
        Real* AMREX_RESTRICT uy_save;
        Real* AMREX_RESTRICT uz_save;
        if (!(WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags)) {
            x_save = xp_save.dataPtr();
            y_save = yp_save.dataPtr();
            z_save = zp_save.dataPtr();
            ux_save = uxp_save.dataPtr();
            uy_save = uyp_save.dataPtr();
            uz_save = uzp_save.dataPtr();
        } else {
            x_save = pti.GetAttribs(particle_comps["xold"]).dataPtr();
            y_save = pti.GetAttribs(particle_comps["yold"]).dataPtr();
            z_save = pti.GetAttribs(particle_comps["zold"]).dataPtr();
            ux_save = pti.GetAttribs(particle_comps["uxold"]).dataPtr();
            uy_save = pti.GetAttribs(particle_comps["uyold"]).dataPtr();
            uz_save = pti.GetAttribs(particle_comps["uzold"]).dataPtr();
        }

        // Undo the push for particles not injected yet.
        // The zp are advanced a fixed amount.
        const Real z_plane_lev = zinject_plane_lev;
        const Real vz_ave_boosted = vzbeam_ave_boosted;
        const bool rigid = rigid_advance;
        const Real inv_csq = 1./(PhysConst::c*PhysConst::c);
        amrex::ParallelFor( pti.numParticles(),
            [=] AMREX_GPU_DEVICE (long i) {
            if (z[i] <= z_plane_lev) {
                ux[i] = ux_save[i];
                uy[i] = uy_save[i];
                uz[i] = uz_save[i];
                gi[i] = 1./std::sqrt(1. + (ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i])*inv_csq);
                x[i] = x_save[i];
                y[i] = y_save[i];
                if (rigid) {
                    z[i] = z_save[i] + dt*vz_ave_boosted;
                }
                else {
                    z[i] = z_save[i] + dt*uz[i]*gi[i];
                }
            }
        }
        );
    }
}

void
RigidInjectedParticleContainer::Evolve (int lev,
                                        const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                        const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                                        MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                        MultiFab* cjx, MultiFab* cjy, MultiFab* cjz,
                                        MultiFab* rho, MultiFab* crho,
                                        const MultiFab* cEx, const MultiFab* cEy, const MultiFab* cEz,
                                        const MultiFab* cBx, const MultiFab* cBy, const MultiFab* cBz,
                                        Real t, Real dt)
{

    // Update location of injection plane in the boosted frame
    zinject_plane_lev_previous = zinject_plane_levels[lev];
    zinject_plane_levels[lev] -= dt*WarpX::beta_boost*PhysConst::c;
    zinject_plane_lev = zinject_plane_levels[lev];

    // Set the done injecting flag whan the inject plane moves out of the
    // simulation domain.
    // It is much easier to do this check, rather than checking if all of the
    // particles have crossed the inject plane.
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();
    done_injecting[lev] = (zinject_plane_levels[lev] < plo[2] || zinject_plane_levels[lev] > phi[2]);
    done_injecting_lev = done_injecting[lev];

    PhysicalParticleContainer::Evolve (lev,
                                       Ex, Ey, Ez,
                                       Bx, By, Bz,
                                       jx, jy, jz,
                                       cjx, cjy, cjz,
                                       rho, crho,
                                       cEx, cEy, cEz,
                                       cBx, cBy, cBz,
                                       t, dt);
}

void
RigidInjectedParticleContainer::PushP (int lev, Real dt,
                                       const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                       const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    BL_PROFILE("RigidInjectedParticleContainer::PushP");

    if (do_not_push) return;

    const std::array<Real,3>& dx = WarpX::CellSize(lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
        int thread_num = omp_get_thread_num();
#else
        int thread_num = 0;
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const Box& box = pti.validbox();

            auto& attribs = pti.GetAttribs();

            auto& uxp = attribs[PIdx::ux];
            auto& uyp = attribs[PIdx::uy];
            auto& uzp = attribs[PIdx::uz];
            auto& Exp = attribs[PIdx::Ex];
            auto& Eyp = attribs[PIdx::Ey];
            auto& Ezp = attribs[PIdx::Ez];
            auto& Bxp = attribs[PIdx::Bx];
            auto& Byp = attribs[PIdx::By];
            auto& Bzp = attribs[PIdx::Bz];

            const long np = pti.numParticles();

            // Data on the grid
            const FArrayBox& exfab = Ex[pti];
            const FArrayBox& eyfab = Ey[pti];
            const FArrayBox& ezfab = Ez[pti];
            const FArrayBox& bxfab = Bx[pti];
            const FArrayBox& byfab = By[pti];
            const FArrayBox& bzfab = Bz[pti];

            Exp.assign(np,0.0);
            Eyp.assign(np,0.0);
            Ezp.assign(np,0.0);
            Bxp.assign(np,WarpX::B_external[0]);
            Byp.assign(np,WarpX::B_external[1]);
            Bzp.assign(np,WarpX::B_external[2]);

            m_giv[thread_num].resize(np);

            //
            // copy data from particle container to temp arrays
            //
            pti.GetPosition(m_xp[thread_num], m_yp[thread_num], m_zp[thread_num]);

            const std::array<Real,3>& xyzmin_grid = WarpX::LowerCorner(box, lev);
            const int* ixyzmin_grid = box.loVect();

            const int ll4symtry          = false;
            long lvect_fieldgathe = 64;

            warpx_geteb_energy_conserving(
                &np,
                m_xp[thread_num].dataPtr(),
                m_yp[thread_num].dataPtr(),
                m_zp[thread_num].dataPtr(),
                Exp.dataPtr(),Eyp.dataPtr(),Ezp.dataPtr(),
                Bxp.dataPtr(),Byp.dataPtr(),Bzp.dataPtr(),
                ixyzmin_grid,
                &xyzmin_grid[0], &xyzmin_grid[1], &xyzmin_grid[2],
                &dx[0], &dx[1], &dx[2],
                &WarpX::nox, &WarpX::noy, &WarpX::noz,
                BL_TO_FORTRAN_ANYD(exfab),
                BL_TO_FORTRAN_ANYD(eyfab),
                BL_TO_FORTRAN_ANYD(ezfab),
                BL_TO_FORTRAN_ANYD(bxfab),
                BL_TO_FORTRAN_ANYD(byfab),
                BL_TO_FORTRAN_ANYD(bzfab),
                &ll4symtry, &WarpX::l_lower_order_in_v, &WarpX::do_nodal,
                &lvect_fieldgathe, &WarpX::field_gathering_algo);

            // Save the position and momenta, making copies
            auto uxp_save = uxp;
            auto uyp_save = uyp;
            auto uzp_save = uzp;

            // This wraps the momentum advance so that inheritors can modify the call.
            // Extract pointers to the different particle quantities
            const Real* const AMREX_RESTRICT zp = m_zp[thread_num].dataPtr();
            Real* const AMREX_RESTRICT gi = m_giv[thread_num].dataPtr();
            Real* const AMREX_RESTRICT uxpp = uxp.dataPtr();
            Real* const AMREX_RESTRICT uypp = uyp.dataPtr();
            Real* const AMREX_RESTRICT uzpp = uzp.dataPtr();
            const Real* const AMREX_RESTRICT Expp = Exp.dataPtr();
            const Real* const AMREX_RESTRICT Eypp = Eyp.dataPtr();
            const Real* const AMREX_RESTRICT Ezpp = Ezp.dataPtr();
            const Real* const AMREX_RESTRICT Bxpp = Bxp.dataPtr();
            const Real* const AMREX_RESTRICT Bypp = Byp.dataPtr();
            const Real* const AMREX_RESTRICT Bzpp = Bzp.dataPtr();

            // Loop over the particles and update their momentum
            const Real q = this->charge;
            const Real m = this->mass;
            if (WarpX::particle_pusher_algo == ParticlePusherAlgo::Boris){
                amrex::ParallelFor( pti.numParticles(),
                    [=] AMREX_GPU_DEVICE (long i) {
                        UpdateMomentumBoris( uxpp[i], uypp[i], uzpp[i], gi[i],
                              Expp[i], Eypp[i], Ezpp[i], Bxpp[i], Bypp[i], Bzpp[i], q, m, dt);
                    }
                );
            } else if (WarpX::particle_pusher_algo == ParticlePusherAlgo::Vay) {
                amrex::ParallelFor( pti.numParticles(),
                    [=] AMREX_GPU_DEVICE (long i) {
                        UpdateMomentumVay( uxpp[i], uypp[i], uzpp[i], gi[i],
                              Expp[i], Eypp[i], Ezpp[i], Bxpp[i], Bypp[i], Bzpp[i], q, m, dt);
                    }
                );
            } else {
              amrex::Abort("Unknown particle pusher");
            };

            // Undo the push for particles not injected yet.
            // It is assumed that PushP will only be called on the first and last steps
            // and that no particles will cross zinject_plane.
            const Real* const AMREX_RESTRICT ux_save = uxp_save.dataPtr();
            const Real* const AMREX_RESTRICT uy_save = uyp_save.dataPtr();
            const Real* const AMREX_RESTRICT uz_save = uzp_save.dataPtr();
            const Real zz = zinject_plane_levels[lev];
            amrex::ParallelFor( pti.numParticles(),
                [=] AMREX_GPU_DEVICE (long i) {
                if (zp[i] <= zz) {
                    uxpp[i] = ux_save[i];
                    uypp[i] = uy_save[i];
                    uzpp[i] = uz_save[i];
                }
            }
            );

        }
    }
}
