#include <limits>
#include <sstream>

#include <ParticleContainer.H>
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
RigidInjectedParticleContainer::Evolve (int lev,
                                        const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                        const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                                        MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                        MultiFab* rho, MultiFab* rho2,
                                        Real t, Real dt)
{
    BL_PROFILE("PPC::Evolve()");
    BL_PROFILE_VAR_NS("PPC::Evolve::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PICSAR::FieldGather", blp_pxr_fg);
    BL_PROFILE_VAR_NS("PICSAR::ParticlePush", blp_pxr_pp);
    BL_PROFILE_VAR_NS("PICSAR::CurrentDeposition", blp_pxr_cd);
    BL_PROFILE_VAR_NS("PPC::Evolve::Accumulate", blp_accumulate);

    const std::array<Real,3>& dx = WarpX::CellSize(lev);

    // WarpX assumes the same number of guard cells for Ex, Ey, Ez, Bx, By, Bz
    long ngE = Ex.nGrow();
    // WarpX assumes the same number of guard cells for Jx, Jy, Jz
    long ngJ = jx.nGrow();
    long ngJDeposit   = (WarpX::use_filter) ? ngJ +1   : ngJ;

    BL_ASSERT(OnSameGrids(lev,Ex));

    MultiFab* cost = WarpX::getCosts(lev);

    // Update location of injection plane in the boosted frame
    const Real zinject_plane_previous = zinject_plane;
    zinject_plane -= dt*WarpX::beta_boost*PhysConst::c;

    // Setup check of whether more particles need to be injected
    // Does something need to be done with OPENMP since this would
    // be a reduction?
    bool done_injecting_temp = true;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<Real> xp, yp, zp, giv;
        FArrayBox local_rho, local_jx, local_jy, local_jz;
        FArrayBox filtered_rho, filtered_jx, filtered_jy, filtered_jz;

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            Real wt = ParallelDescriptor::second();

            const Box& box = pti.validbox();

            auto& attribs = pti.GetAttribs();

            auto&  wp = attribs[PIdx::w];
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
            FArrayBox&       jxfab = jx[pti];
            FArrayBox&       jyfab = jy[pti];
            FArrayBox&       jzfab = jz[pti];

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
            BL_PROFILE_VAR_START(blp_copy);
            pti.GetPosition(xp, yp, zp);
            BL_PROFILE_VAR_STOP(blp_copy);

            const std::array<Real,3>& xyzmin_tile = WarpX::LowerCorner(pti.tilebox(), lev);
            const std::array<Real,3>& xyzmin_grid = WarpX::LowerCorner(box, lev);

            long lvect = 8;

            auto depositCharge = [&] (MultiFab* rhomf)
            {
                long ngRho = rhomf->nGrow();
                long ngRhoDeposit = (WarpX::use_filter) ? ngRho +1 : ngRho;

                Real* data_ptr;
                const int *rholen;
                FArrayBox& rhofab = (*rhomf)[pti];
                Box tile_box = convert(pti.tilebox(), IntVect::TheUnitVector());
                Box grown_box;
                const std::array<Real, 3>& xyzmin = xyzmin_tile;
                tile_box.grow(ngRho);
                if (WarpX::use_filter) {
                    grown_box = tile_box;
                    grown_box.grow(1);
                    local_rho.resize(grown_box);
                } else {
                    local_rho.resize(tile_box);
                }
                local_rho = 0.0;
                data_ptr = local_rho.dataPtr();
                rholen = local_rho.length();

#if (BL_SPACEDIM == 3)
                const long nx = rholen[0]-1-2*ngRhoDeposit;
                const long ny = rholen[1]-1-2*ngRhoDeposit;
                const long nz = rholen[2]-1-2*ngRhoDeposit;
#else
                const long nx = rholen[0]-1-2*ngRhoDeposit;
                const long ny = 0;
                const long nz = rholen[1]-1-2*ngRhoDeposit;
#endif
                const int l_2drz = (const int)Geometry::IsRZ();
                warpx_charge_deposition(data_ptr, &np,
                                        xp.data(), yp.data(), zp.data(), wp.data(),
                                        &this->charge,
                                        &xyzmin[0], &xyzmin[1], &xyzmin[2],
                                        &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
                                        &ngRhoDeposit, &ngRhoDeposit, &ngRhoDeposit,
                                        &WarpX::nox,&WarpX::noy,&WarpX::noz,
                                        &lvect, &WarpX::charge_deposition_algo, &l_2drz);

                const int ncomp = 1;
                if (WarpX::use_filter) {

                    filtered_rho.resize(tile_box);
                    filtered_rho = 0;

                    WRPX_FILTER(local_rho.dataPtr(),
                                local_rho.loVect(),
                                local_rho.hiVect(),
                                filtered_rho.dataPtr(),
                                filtered_rho.loVect(),
                                filtered_rho.hiVect(),
                                ncomp);

                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(filtered_rho),
                                                BL_TO_FORTRAN_3D(rhofab), ncomp);


                } else {
                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_rho),
                                                BL_TO_FORTRAN_3D(rhofab), ncomp);
                }
            };

            if (rho) depositCharge(rho);

            if (! do_not_push)
            {
                //
                // Field Gather of Aux Data (i.e., the full solution)
                //
                const int ll4symtry          = false;
                const int l_lower_order_in_v = true;
                long lvect_fieldgathe = 64;
                const int l_2drz = (const int)Geometry::IsRZ();
                BL_PROFILE_VAR_START(blp_pxr_fg);
                warpx_geteb_energy_conserving(
                    &np, xp.data(), yp.data(), zp.data(),
                    Exp.data(),Eyp.data(),Ezp.data(),
                    Bxp.data(),Byp.data(),Bzp.data(),
                    &xyzmin_grid[0], &xyzmin_grid[1], &xyzmin_grid[2],
                    &dx[0], &dx[1], &dx[2],
                    &WarpX::nox, &WarpX::noy, &WarpX::noz,
                    exfab.dataPtr(), &ngE, exfab.length(),
                    eyfab.dataPtr(), &ngE, eyfab.length(),
                    ezfab.dataPtr(), &ngE, ezfab.length(),
                    bxfab.dataPtr(), &ngE, bxfab.length(),
                    byfab.dataPtr(), &ngE, byfab.length(),
                    bzfab.dataPtr(), &ngE, bzfab.length(),
                    &ll4symtry, &l_2drz, &l_lower_order_in_v,
                    &lvect_fieldgathe, &WarpX::field_gathering_algo);
                BL_PROFILE_VAR_STOP(blp_pxr_fg);

                // Save the position and momenta, making copies
                auto xp_save = xp;
                auto yp_save = yp;
                auto zp_save = zp;
                auto uxp_save = uxp;
                auto uyp_save = uyp;
                auto uzp_save = uzp;

                if (!done_injecting) {
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

                //
                // Particle Push
                //
                BL_PROFILE_VAR_START(blp_pxr_pp);
                warpx_particle_pusher(&np, xp.data(), yp.data(), zp.data(),
                                      uxp.data(), uyp.data(), uzp.data(), giv.data(),
                                      Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
                                      Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
                                      &this->charge, &this->mass, &dt,
                                      &WarpX::particle_pusher_algo);
                BL_PROFILE_VAR_STOP(blp_pxr_pp);

                if (!done_injecting) {
                    // Undo the push for particles not injected yet.
                    // The zp are advanced a fixed amount.
                    for (int i=0 ; i < zp.size() ; i++) {
                        if (zp[i] <= zinject_plane) {
                            xp[i] = xp_save[i];
                            yp[i] = yp_save[i];
                            zp[i] = zp_save[i] + dt*vzbeam_ave_boosted;
                            uxp[i] = uxp_save[i];
                            uyp[i] = uyp_save[i];
                            uzp[i] = uzp_save[i];
                            giv[i] = 1./std::sqrt(1. + (uxp[i]*uxp[i] + uyp[i]*uyp[i] + uzp[i]*uzp[i])/(PhysConst::c*PhysConst::c));
                            done_injecting_temp = false;
                        }
                    }
                }

                //
                // Current Deposition onto fine patch
                //

                BL_PROFILE_VAR_START(blp_pxr_cd);
                Real *jx_ptr, *jy_ptr, *jz_ptr;
                const int  *jxntot, *jyntot, *jzntot;
                Box tbx = convert(pti.tilebox(), WarpX::jx_nodal_flag);
                Box tby = convert(pti.tilebox(), WarpX::jy_nodal_flag);
                Box tbz = convert(pti.tilebox(), WarpX::jz_nodal_flag);
                Box gtbx, gtby, gtbz;

                const std::array<Real, 3>& xyzmin = xyzmin_tile;

                tbx.grow(ngJ);
                tby.grow(ngJ);
                tbz.grow(ngJ);

                if (WarpX::use_filter) {

                    gtbx = tbx;
                    gtbx.grow(1);

                    gtby = tby;
                    gtby.grow(1);

                    gtbz = tbz;
                    gtbz.grow(1);

                    local_jx.resize(gtbx);
                    local_jy.resize(gtby);
                    local_jz.resize(gtbz);
                } else {
                    local_jx.resize(tbx);
                    local_jy.resize(tby);
                    local_jz.resize(tbz);
                }

                local_jx = 0.0;
                local_jy = 0.0;
                local_jz = 0.0;

                jx_ptr = local_jx.dataPtr();
                jy_ptr = local_jy.dataPtr();
                jz_ptr = local_jz.dataPtr();

                jxntot = local_jx.length();
                jyntot = local_jy.length();
                jzntot = local_jz.length();

                warpx_current_deposition(
                    jx_ptr, &ngJDeposit, jxntot,
                    jy_ptr, &ngJDeposit, jyntot,
                    jz_ptr, &ngJDeposit, jzntot,
                    &np, xp.data(), yp.data(), zp.data(),
                    uxp.data(), uyp.data(), uzp.data(),
                    giv.data(), wp.data(), &this->charge,
                    &xyzmin[0], &xyzmin[1], &xyzmin[2],
                    &dt, &dx[0], &dx[1], &dx[2],
                    &WarpX::nox,&WarpX::noy,&WarpX::noz,
                    &lvect,&WarpX::current_deposition_algo,&l_2drz);


                BL_PROFILE_VAR_STOP(blp_pxr_cd);

                BL_PROFILE_VAR_START(blp_accumulate);

                const int ncomp = 1;
                if (WarpX::use_filter) {

                    filtered_jx.resize(tbx);
                    filtered_jx = 0.0;

                    WRPX_FILTER(local_jx.dataPtr(),
                                local_jx.loVect(),
                                local_jx.hiVect(),
                                filtered_jx.dataPtr(),
                                filtered_jx.loVect(),
                                filtered_jx.hiVect(),
                                ncomp);

                    filtered_jy.resize(tby);
                    filtered_jy = 0.0;

                    WRPX_FILTER(local_jy.dataPtr(),
                                local_jy.loVect(),
                                local_jy.hiVect(),
                                filtered_jy.dataPtr(),
                                filtered_jy.loVect(),
                                filtered_jy.hiVect(),
                                ncomp);

                    filtered_jz.resize(tbz);
                    filtered_jz = 0.0;

                    WRPX_FILTER(local_jz.dataPtr(),
                                local_jz.loVect(),
                                local_jz.hiVect(),
                                filtered_jz.dataPtr(),
                                filtered_jz.loVect(),
                                filtered_jz.hiVect(),
                                ncomp);

                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(filtered_jx),
                                                BL_TO_FORTRAN_3D(jxfab), ncomp);

                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(filtered_jy),
                                                BL_TO_FORTRAN_3D(jyfab), ncomp);

                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(filtered_jz),
                                                BL_TO_FORTRAN_3D(jzfab), ncomp);

                } else {

                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_jx),
                                                BL_TO_FORTRAN_3D(jxfab), ncomp);

                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_jy),
                                                BL_TO_FORTRAN_3D(jyfab), ncomp);

                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_jz),
                                                BL_TO_FORTRAN_3D(jzfab), ncomp);
                }
                BL_PROFILE_VAR_STOP(blp_accumulate);

                //
                // copy particle data back
                //
                BL_PROFILE_VAR_START(blp_copy);
                pti.SetPosition(xp, yp, zp);
                BL_PROFILE_VAR_STOP(blp_copy);
            }

            if (rho2) depositCharge(rho2);

            if (cost) {
                const Box& tbx = pti.tilebox();
                wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
                (*cost)[pti].plus(wt, tbx);
            }
        }
    }

    done_injecting = done_injecting_temp;

}

void
RigidInjectedParticleContainer::PushP (int lev, Real dt,
                                       const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                       const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    if (do_not_push) return;

    const std::array<Real,3>& dx = WarpX::CellSize(lev);

    // WarpX assumes the same number of guard cells for Ex, Ey, Ez, Bx, By, Bz
    long ngE = Ex.nGrow();

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

            const int ll4symtry          = false;
            const int l_lower_order_in_v = true;
            long lvect_fieldgathe = 64;
            warpx_geteb_energy_conserving(
                &np, xp.data(), yp.data(), zp.data(),
                Exp.data(),Eyp.data(),Ezp.data(),
                Bxp.data(),Byp.data(),Bzp.data(),
                &xyzmin_grid[0], &xyzmin_grid[1], &xyzmin_grid[2],
                &dx[0], &dx[1], &dx[2],
                &WarpX::nox, &WarpX::noy, &WarpX::noz,
                exfab.dataPtr(), &ngE, exfab.length(),
                eyfab.dataPtr(), &ngE, eyfab.length(),
                ezfab.dataPtr(), &ngE, ezfab.length(),
                bxfab.dataPtr(), &ngE, bxfab.length(),
                byfab.dataPtr(), &ngE, byfab.length(),
                bzfab.dataPtr(), &ngE, bzfab.length(),
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
