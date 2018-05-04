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

using namespace amrex;

RigidInjectedParticleContainer::RigidInjectedParticleContainer (AmrCore* amr_core, int ispecies,
                                                                const std::string& name)
    : PhysicalParticleContainer(amr_core, ispecies, name)
{

    ParmParse pp(species_name);

    pp.get("zinject_plane", zinject_plane);
    pp.query("projected", projected);
    pp.query("focused", focused);

}

void RigidInjectedParticleContainer::InitData()
{
    AddParticles(0); // Note - add on level 0
    BoostandRemapParticles();
    Redistribute();  // We then redistribute
}

void
RigidInjectedParticleContainer::BoostandRemapParticles()
{

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<Real> xp, yp, zp, giv;

        // Get the average beam velocity in the boosted frame
        // This value is saved to advance the particles not injected yet
        const Real vzbeam_ave_lab = meanParticleVelocity(false)[2];
        vzbeam_ave_boosted = (vzbeam_ave_lab - WarpX::beta_boost*PhysConst::c)/(1. - vzbeam_ave_lab*WarpX::beta_boost/PhysConst::c);

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

                // t_lab is the time in the lab frame that the particles reaches z=0
                // The location and time (z=0, t=0) is a synchronization point between the
                // lab and boosted frames.
                const Real t_lab = -zp[i]/vz_lab;

                if (!projected) {
                    xp[i] += t_lab*vx_lab;
                    yp[i] += t_lab*vy_lab;
                }
                if (focused) {
                    // Correct for focusing effect from shift from z=0 to zinject
                    const Real tfocus = -zinject_plane*WarpX::gamma_boost/vz_lab;
                    xp[i] -= tfocus*vx_lab;
                    yp[i] -= tfocus*vy_lab;
                }

                // Time of the particle in the boosted frame given its position in the lab frame at t=0.
                const Real tpr = -WarpX::gamma_boost*WarpX::beta_boost*zp[i]/PhysConst::c;

                // Position of the particle in the boosted from given its position in the lab frame at t=0.
                const Real zpr = WarpX::gamma_boost*zp[i];

                // Momentum of the particle in the boosted frame (assuming that it is fixed).
                uzp[i] = WarpX::gamma_boost*(uzp[i] - WarpX::beta_boost*PhysConst::c*gamma_lab);

                // Put the particle at the location in the boosted frame at boost frame t=0,
                // with the particle moving at the average velocity
                zp[i] = zpr - vzbeam_ave_boosted*tpr;

            }

            // Copy the data back to the particle container
            pti.SetPosition(xp, yp, zp);

        }
    }
}

void
RigidInjectedParticleContainer::PushPX(WarpXParIter& pti,
	                               Vector<Real>& xp, Vector<Real>& yp, Vector<Real>& zp,
                                       Vector<Real>& giv,
                                       Real dt)
{

    // This wraps the call to warpx_particle_pusher so that inheritors can modify the call.
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
    const long np  = pti.numParticles();

    // Save the position and momenta, making copies
    Vector<Real> xp_save, yp_save, zp_save, uxp_save, uyp_save, uzp_save;

    if (!done_injecting) {
        xp_save = xp;
        yp_save = yp;
        zp_save = zp;
        uxp_save = uxp;
        uyp_save = uyp;
        uzp_save = uzp;
        // Scale the fields of particles about to cross the injection plane.
        // This only approximates what should be happening. The particles
        // should by advanced a fraction of a time step instead.
        // Scaling the fields is much easier and may be good enough.
        for (int i=0 ; i < zp.size() ; i++) {
            const Real dtscale = dt - (zinject_plane_previous - zp[i])/(vzbeam_ave_boosted + WarpX::beta_boost*PhysConst::c);
            if (0. < dtscale && dtscale < dt) {
                Exp[i] *= dtscale;
                Eyp[i] *= dtscale;
                Ezp[i] *= dtscale;
                Bxp[i] *= dtscale;
                Byp[i] *= dtscale;
                Bzp[i] *= dtscale;
            }
        }
    }

    warpx_particle_pusher(&np, xp.data(), yp.data(), zp.data(),
                          uxp.data(), uyp.data(), uzp.data(), giv.data(),
                          Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
                          Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
                          &this->charge, &this->mass, &dt,
                          &WarpX::particle_pusher_algo);

    if (!done_injecting) {
        // Undo the push for particles not injected yet.
        // The zp are advanced a fixed amount.
        bool done = true;
        for (int i=0 ; i < zp.size() ; i++) {
            if (zp[i] <= zinject_plane) {
                xp[i] = xp_save[i];
                yp[i] = yp_save[i];
                zp[i] = zp_save[i] + dt*vzbeam_ave_boosted;
                uxp[i] = uxp_save[i];
                uyp[i] = uyp_save[i];
                uzp[i] = uzp_save[i];
                giv[i] = 1./std::sqrt(1. + (uxp[i]*uxp[i] + uyp[i]*uyp[i] + uzp[i]*uzp[i])/(PhysConst::c*PhysConst::c));
                done = false;
            }
        }

#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        done_injecting_temp[tid] = done;
    }

}

void
RigidInjectedParticleContainer::Evolve (int lev,
                                        const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                        const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                                        MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                        MultiFab* rho, MultiFab* rho2,
                                        Real t, Real dt)
{

    // Update location of injection plane in the boosted frame
    zinject_plane_previous = zinject_plane;
    zinject_plane -= dt*WarpX::beta_boost*PhysConst::c;

    // Setup check of whether more particles need to be injected
#ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
#else
    const int nthreads = 1;
#endif
    done_injecting_temp.assign(nthreads, 1); // We do not use bool because vector<bool> is special.

    PhysicalParticleContainer::Evolve (lev,
				       Ex, Ey, Ez,
				       Bx, By, Bz,
				       jx, jy, jz,
                                       rho, rho2,
                                       t, dt);

    // Check if all done_injecting_temp are still true.
    done_injecting = std::all_of(done_injecting_temp.begin(), done_injecting_temp.end(),
                                 [](int i) -> bool { return i; });
}

void
RigidInjectedParticleContainer::PushP (int lev, Real dt,
                                       const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                       const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    if (do_not_push) return;

    const std::array<Real,3>& dx = WarpX::CellSize(lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<Real> xp, yp, zp, giv;

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

            giv.resize(np);

            //
            // copy data from particle container to temp arrays
            //
            pti.GetPosition(xp, yp, zp);

            const std::array<Real,3>& xyzmin_grid = WarpX::LowerCorner(box, lev);
            const int* ixyzmin_grid = box.loVect();

            const int ll4symtry          = false;
            const int l_lower_order_in_v = true;
            long lvect_fieldgathe = 64;
            warpx_geteb_energy_conserving(
                &np, xp.data(), yp.data(), zp.data(),
                Exp.data(),Eyp.data(),Ezp.data(),
                Bxp.data(),Byp.data(),Bzp.data(),
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
                &ll4symtry, &l_lower_order_in_v,
                &lvect_fieldgathe, &WarpX::field_gathering_algo);

            // Save the position and momenta, making copies
            auto uxp_save = uxp;
            auto uyp_save = uyp;
            auto uzp_save = uzp;

            warpx_particle_pusher_momenta(&np, xp.data(), yp.data(), zp.data(),
                                          uxp.data(), uyp.data(), uzp.data(), giv.data(),
                                          Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
                                          Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
                                          &this->charge, &this->mass, &dt,
                                          &WarpX::particle_pusher_algo);

            // Undo the push for particles not injected yet.
            // It is assumed that PushP will only be called on the first and last steps
            // and that no particles will cross zinject_plane.
            for (int i=0 ; i < zp.size() ; i++) {
                if (zp[i] <= zinject_plane) {
                    uxp[i] = uxp_save[i];
                    uyp[i] = uyp_save[i];
                    uzp[i] = uzp_save[i];
                }
            }

        }
    }
}
