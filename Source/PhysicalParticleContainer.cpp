#include <limits>
#include <sstream>

#include <ParticleContainer.H>
#include <WarpX_f.H>
#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpXWrappers.h>


using namespace amrex;

PhysicalParticleContainer::PhysicalParticleContainer (AmrCore* amr_core, int ispecies,
                                                      const std::string& name)
    : WarpXParticleContainer(amr_core, ispecies),
      species_name(name)
{
    plasma_injector.reset(new PlasmaInjector(species_id, species_name));
    charge = plasma_injector->getCharge();
    mass = plasma_injector->getMass();
}

void PhysicalParticleContainer::InitData()
{
    AddParticles(0); // Note - add on level 0
    if (maxLevel() > 0) {
        Redistribute();  // We then redistribute
    }
}

void
PhysicalParticleContainer::AddGaussianBeam(Real x_m, Real y_m, Real z_m,
                                           Real x_rms, Real y_rms, Real z_rms,
                                           Real q_tot, long npart) {

    const Geometry& geom     = m_gdb->Geom(0);
    RealBox containing_bx = geom.ProbDomain();

    std::mt19937_64 mt(0451);
    std::normal_distribution<double> distx(x_m, x_rms);
    std::normal_distribution<double> disty(y_m, y_rms);
    std::normal_distribution<double> distz(z_m, z_rms);

    std::array<Real,PIdx::nattribs> attribs;
    attribs.fill(0.0);

    if (ParallelDescriptor::IOProcessor()) {
       std::array<Real, 3> u;
       Real weight;
       for (long i = 0; i < npart; ++i) {
#if ( BL_SPACEDIM == 3 )
            weight = q_tot/npart/charge;
            Real x = distx(mt);
            Real y = disty(mt);
            Real z = distz(mt);
#elif ( BL_SPACEDIM == 2 )
            weight = q_tot/npart/charge/y_rms;
            Real x = distx(mt);
            Real y = 0.;
            Real z = distz(mt);
#endif
        if (plasma_injector->insideBounds(x, y, z)) {
	    plasma_injector->getMomentum(u, x, y, z);
            attribs[PIdx::ux] = u[0];
            attribs[PIdx::uy] = u[1];
            attribs[PIdx::uz] = u[2];
            attribs[PIdx::w ] = weight;
            AddOneParticle(0, 0, 0, x, y, z, attribs);
            }
        }
    }
    Redistribute();
}

void
PhysicalParticleContainer::AddParticles (int lev)
{
    BL_PROFILE("PhysicalParticleContainer::AddParticles()");

    if (plasma_injector->add_single_particle) {
        AddNParticles(lev, 1,
                      &(plasma_injector->single_particle_pos[0]),
                      &(plasma_injector->single_particle_pos[1]),
                      &(plasma_injector->single_particle_pos[2]),
                      &(plasma_injector->single_particle_vel[0]),
                      &(plasma_injector->single_particle_vel[1]),
                      &(plasma_injector->single_particle_vel[2]),
                      1, &(plasma_injector->single_particle_weight), 0);
        return;
    }

    if (plasma_injector->gaussian_beam) {
        AddGaussianBeam(plasma_injector->x_m,
                        plasma_injector->y_m,
                        plasma_injector->z_m,
                        plasma_injector->x_rms,
                        plasma_injector->y_rms,
                        plasma_injector->z_rms,
                        plasma_injector->q_tot,
                        plasma_injector->npart);

        return;
    }

    if ( plasma_injector->doInjection() ) {
        AddPlasma( lev );
    }
}

/**
 * Create new macroparticles for this species, with a fixed
 * number of particles per cell (in the cells of `part_realbox`).
 * The new particles are only created inside the intersection of `part_realbox`
 * with the local grid for the current proc.
 * @param lev the index of the refinement level
 * @param part_realbox the box in which new particles should be created
 * (this box should correspond to an integer number of cells in each direction,
 * but its boundaries need not be aligned with the actual cells of the simulation)
 */
void
PhysicalParticleContainer::AddPlasma(int lev, RealBox part_realbox )
{

    // If no part_realbox is provided, initialize particles in the whole domain
    const Geometry& geom = Geom(lev);
    if (!part_realbox.ok()) part_realbox = geom.ProbDomain();

    int num_ppc = plasma_injector->num_particles_per_cell;

    const Real* dx = geom.CellSize();

    Real scale_fac;
#if BL_SPACEDIM==3
    scale_fac = dx[0]*dx[1]*dx[2]/num_ppc;
#elif BL_SPACEDIM==2
    scale_fac = dx[0]*dx[1]/num_ppc;
#endif

#ifdef _OPENMP
    // First touch all tiles in the map in serial
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        GetParticles(lev)[std::make_pair(grid_id, tile_id)];
    }
#endif

#ifdef _OPENMP
#pragma omp parallel if (not WarpX::serialize_ics)
#endif
    {
        std::array<Real,PIdx::nattribs> attribs;
        attribs.fill(0.0);

        // Loop through the tiles
        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {

            const Box& tile_box  = mfi.tilebox();
            const RealBox tile_realbox = WarpX::getRealBox(tile_box, lev);

            // Find the cells of part_box that overlap with tile_realbox
            // If there is no overlap, just go to the next tile in the loop
            RealBox overlap_realbox;
            Box overlap_box;
            Real ncells_adjust;
            bool no_overlap = 0;

            for (int dir=0; dir<BL_SPACEDIM; dir++) {
                if ( tile_realbox.lo(dir) < part_realbox.hi(dir) ) {
                    ncells_adjust = std::floor( (tile_realbox.lo(dir) - part_realbox.lo(dir))/dx[dir] );
                    overlap_realbox.setLo( dir, part_realbox.lo(dir) + std::max(ncells_adjust, 0.) * dx[dir]);
                } else {
                    no_overlap = 1; break;
                }
                if ( tile_realbox.hi(dir) > part_realbox.lo(dir) ) {
                    ncells_adjust = std::floor( (part_realbox.hi(dir) - tile_realbox.hi(dir))/dx[dir] );
                    overlap_realbox.setHi( dir, part_realbox.hi(dir) - std::max(ncells_adjust, 0.) * dx[dir]);
                } else {
                    no_overlap = 1; break;
                }
                // Count the number of cells in this direction in overlap_realbox
                overlap_box.setSmall( dir, 0 );
                overlap_box.setBig( dir,
				    int( round((overlap_realbox.hi(dir)-overlap_realbox.lo(dir))/dx[dir] )) - 1);
            }
            if (no_overlap == 1) continue; // Go to the next tile

            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();

            // Loop through the cells of overlap_box and inject
            // the corresponding particles
            const auto& overlap_corner = overlap_realbox.lo();
            for (IntVect iv = overlap_box.smallEnd();
                 iv <= overlap_box.bigEnd(); overlap_box.next(iv)) {
                for (int i_part=0; i_part<num_ppc;i_part++) {
                    std::array<Real, 3> r;
                    plasma_injector->getPositionUnitBox(r, i_part);
#if ( BL_SPACEDIM == 3 )
                    Real x = overlap_corner[0] + (iv[0] + r[0])*dx[0];
                    Real y = overlap_corner[1] + (iv[1] + r[1])*dx[1];
                    Real z = overlap_corner[2] + (iv[2] + r[2])*dx[2];
#elif ( BL_SPACEDIM == 2 )
                    Real x = overlap_corner[0] + (iv[0] + r[0])*dx[0];
                    Real y = 0;
                    Real z = overlap_corner[1] + (iv[1] + r[1])*dx[1];
#endif
                    // If the new particle is not inside the tile box,
                    // go to the next generated particle.
#if ( BL_SPACEDIM == 3 )
                    if(!tile_realbox.contains( RealVect{x, y, z} )) continue;
#elif ( BL_SPACEDIM == 2 )
                    if(!tile_realbox.contains( RealVect{x, z} )) continue;
#endif

                    if (plasma_injector->insideBounds(x, y, z)) {
                        Real dens;
                        std::array<Real, 3> u;
                        if (WarpX::gamma_boost == 1.){
                            // Lab-frame simulation
                            plasma_injector->getMomentum(u, x, y, z);
                            dens = plasma_injector->getDensity(x, y, z);
                        } else {
                            // Boosted-frame simulation: call `getMomentum`
                            // and `getDensity` with lab-frame Parameters
                            // (Assumes that the plasma has a low velocity,
                            // and that the boost is along z)
                            Real t = WarpX::GetInstance().gett_new(lev);
                            Real v_boost = WarpX::beta_boost*PhysConst::c;
                            Real z_lab = WarpX::gamma_boost*( z + v_boost*t );
                            plasma_injector->getMomentum(u, x, y, z_lab);
                            dens = plasma_injector->getDensity(x, y, z_lab);
                            // Perform Lorentz transform
                            // (Assumes that the plasma has a low velocity)
                            u[2] = WarpX::gamma_boost * ( u[2] - v_boost );
                            dens = WarpX::gamma_boost * dens;
                        }
                        attribs[PIdx::w ] = dens * scale_fac;
                        attribs[PIdx::ux] = u[0];
                        attribs[PIdx::uy] = u[1];
                        attribs[PIdx::uz] = u[2];
                        AddOneParticle(lev, grid_id, tile_id, x, y, z, attribs);
                    }
                }
            }
        }
    }
}

void
PhysicalParticleContainer::
FieldGatherES (const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
               const amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks)
{

    const int num_levels = E.size();
    const int ng = E[0][0]->nGrow();

    if (num_levels == 1) {
        const int lev = 0;
        const auto& gm = m_gdb->Geom(lev);
        const auto& ba = m_gdb->ParticleBoxArray(lev);

        BoxArray nba = ba;
        nba.surroundingNodes();

        const Real* dx  = gm.CellSize();
        const Real* plo = gm.ProbLo();

        BL_ASSERT(OnSameGrids(lev, *E[lev][0]));

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti) {
            const Box& box = nba[pti];

            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long np  = pti.numParticles();

            auto& attribs = pti.GetAttribs();
            auto& Exp = attribs[PIdx::Ex];
            auto& Eyp = attribs[PIdx::Ey];
#if BL_SPACEDIM == 3
            auto& Ezp = attribs[PIdx::Ez];
#endif
            Exp.assign(np,0.0);
            Eyp.assign(np,0.0);
#if BL_SPACEDIM == 3
            Ezp.assign(np,0.0);
#endif

            const FArrayBox& exfab = (*E[lev][0])[pti];
            const FArrayBox& eyfab = (*E[lev][1])[pti];
#if BL_SPACEDIM == 3
            const FArrayBox& ezfab = (*E[lev][2])[pti];
#endif

            WRPX_INTERPOLATE_CIC(particles.data(), nstride, np,
                                 Exp.data(), Eyp.data(),
#if BL_SPACEDIM == 3
                                 Ezp.data(),
#endif
                                 exfab.dataPtr(), eyfab.dataPtr(),
#if BL_SPACEDIM == 3
                                 ezfab.dataPtr(),
#endif
                                 box.loVect(), box.hiVect(), plo, dx, &ng);
        }

        return;
    }

    const BoxArray& fine_BA = E[1][0]->boxArray();
    const DistributionMapping& fine_dm = E[1][0]->DistributionMap();
    BoxArray coarsened_fine_BA = fine_BA;
    coarsened_fine_BA.coarsen(IntVect(D_DECL(2,2,2)));

    MultiFab coarse_Ex(coarsened_fine_BA, fine_dm, 1, 1);
    MultiFab coarse_Ey(coarsened_fine_BA, fine_dm, 1, 1);
#if BL_SPACEDIM == 3
    MultiFab coarse_Ez(coarsened_fine_BA, fine_dm, 1, 1);
#endif

    coarse_Ex.copy(*E[0][0], 0, 0, 1, 1, 1);
    coarse_Ey.copy(*E[0][1], 0, 0, 1, 1, 1);
#if BL_SPACEDIM == 3
    coarse_Ez.copy(*E[0][2], 0, 0, 1, 1, 1);
#endif

    for (int lev = 0; lev < num_levels; ++lev) {
        const auto& gm = m_gdb->Geom(lev);
        const auto& ba = m_gdb->ParticleBoxArray(lev);

        BoxArray nba = ba;
        nba.surroundingNodes();

        const Real* dx  = gm.CellSize();
        const Real* plo = gm.ProbLo();

        BL_ASSERT(OnSameGrids(lev, *E[lev][0]));

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti) {
            const Box& box = nba[pti];

            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long np  = pti.numParticles();

            auto& attribs = pti.GetAttribs();
            auto& Exp = attribs[PIdx::Ex];
            auto& Eyp = attribs[PIdx::Ey];
#if BL_SPACEDIM == 3
            auto& Ezp = attribs[PIdx::Ez];
#endif
            Exp.assign(np,0.0);
            Eyp.assign(np,0.0);
#if BL_SPACEDIM == 3
            Ezp.assign(np,0.0);
#endif

            const FArrayBox& exfab = (*E[lev][0])[pti];
            const FArrayBox& eyfab = (*E[lev][1])[pti];
#if BL_SPACEDIM == 3
            const FArrayBox& ezfab = (*E[lev][2])[pti];
#endif

            if (lev == 0) {
                WRPX_INTERPOLATE_CIC(particles.data(), nstride, np,
                                     Exp.data(), Eyp.data(),
#if BL_SPACEDIM == 3
                Ezp.data(),
#endif
                                exfab.dataPtr(), eyfab.dataPtr(),
#if BL_SPACEDIM == 3
                                ezfab.dataPtr(),
#endif
                                box.loVect(), box.hiVect(), plo, dx, &ng);
            } else {

                const FArrayBox& exfab_coarse = coarse_Ex[pti];
                const FArrayBox& eyfab_coarse = coarse_Ey[pti];
#if BL_SPACEDIM == 3
                const FArrayBox& ezfab_coarse = coarse_Ez[pti];
#endif
                const Box& coarse_box = coarsened_fine_BA[pti];
                const Real* coarse_dx = Geom(0).CellSize();

                WRPX_INTERPOLATE_CIC_TWO_LEVELS(particles.data(), nstride, np,
                                                Exp.data(), Eyp.data(),
#if BL_SPACEDIM == 3
                                                Ezp.data(),
#endif
                                                exfab.dataPtr(), eyfab.dataPtr(),
#if BL_SPACEDIM == 3
                                                ezfab.dataPtr(),
#endif
                                                box.loVect(), box.hiVect(), dx,
                                                exfab_coarse.dataPtr(), eyfab_coarse.dataPtr(),
#if BL_SPACEDIM == 3
                                                ezfab_coarse.dataPtr(),
#endif
                                                (*masks[1])[pti].dataPtr(),
                                                coarse_box.loVect(), coarse_box.hiVect(), coarse_dx,
                                                plo, &ng, &lev);
            }
        }
    }
}

void
PhysicalParticleContainer::FieldGather (int lev,
                                        const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                        const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    const std::array<Real,3>& dx = WarpX::CellSize(lev);

    // WarpX assumes the same number of guard cells for Ex, Ey, Ez, Bx, By, Bz
    long ng = Ex.nGrow();

    BL_ASSERT(OnSameGrids(lev,Ex));

    MultiFab* cost = WarpX::getCosts(lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	Vector<Real> xp, yp, zp;

	for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
	{
            Real wt = ParallelDescriptor::second();

	    const Box& box = pti.validbox();

            auto& attribs = pti.GetAttribs();

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
	    Bxp.assign(np,0.0);
	    Byp.assign(np,0.0);
	    Bzp.assign(np,0.0);

	    //
	    // copy data from particle container to temp arrays
	    //
            pti.GetPosition(xp, yp, zp);

            const std::array<Real,3>& xyzmin = WarpX::LowerCorner(box, lev);

	    //
	    // Field Gather
	    //
	    const int ll4symtry          = false;
	    const int l_lower_order_in_v = warpx_l_lower_order_in_v();
            long lvect_fieldgathe = 64;
	    warpx_geteb_energy_conserving(
	       &np, xp.data(), yp.data(), zp.data(),
	       Exp.data(),Eyp.data(),Ezp.data(),
	       Bxp.data(),Byp.data(),Bzp.data(),
	       &xyzmin[0], &xyzmin[1], &xyzmin[2],
	       &dx[0], &dx[1], &dx[2],
	       &WarpX::nox, &WarpX::noy, &WarpX::noz,
	       exfab.dataPtr(), &ng, exfab.length(),
	       eyfab.dataPtr(), &ng, eyfab.length(),
	       ezfab.dataPtr(), &ng, ezfab.length(),
               bxfab.dataPtr(), &ng, bxfab.length(),
	       byfab.dataPtr(), &ng, byfab.length(),
	       bzfab.dataPtr(), &ng, bzfab.length(),
	       &ll4symtry, &l_lower_order_in_v,
	       &lvect_fieldgathe, &WarpX::field_gathering_algo);

            if (cost) {
                const Box& tbx = pti.tilebox();
                wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
                (*cost)[pti].plus(wt, tbx);
            }
        }
    }
}

void
PhysicalParticleContainer::EvolveES (const Vector<std::array<std::unique_ptr<MultiFab>, 3> >& E,
                                           Vector<std::unique_ptr<MultiFab> >& rho,
                                     Real t, Real dt)
{
    BL_PROFILE("PPC::EvolveES()");

    int num_levels = rho.size();
    for (int lev = 0; lev < num_levels; ++lev) {
        BL_ASSERT(OnSameGrids(lev, *rho[lev]));
        const auto& gm = m_gdb->Geom(lev);
        const RealBox& prob_domain = gm.ProbDomain();
	for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti) {
            // Particle structs
            auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long np  = pti.numParticles();

            // Particle attributes
            auto& attribs = pti.GetAttribs();
            auto& uxp = attribs[PIdx::ux];
            auto& uyp = attribs[PIdx::uy];

#if BL_SPACEDIM == 3
            auto& uzp = attribs[PIdx::uz];
#endif

            auto& Exp = attribs[PIdx::Ex];
            auto& Eyp = attribs[PIdx::Ey];

#if BL_SPACEDIM == 3
            auto& Ezp = attribs[PIdx::Ez];
#endif
            //
            // Particle Push
            //
            WRPX_PUSH_LEAPFROG(particles.data(), nstride, np,
                               uxp.data(), uyp.data(),
#if BL_SPACEDIM == 3
                               uzp.data(),
#endif
                               Exp.data(), Eyp.data(),
#if BL_SPACEDIM == 3
                               Ezp.data(),
#endif
                               &this->charge, &this->mass, &dt,
                               prob_domain.lo(), prob_domain.hi());
        }
    }
}

void
PhysicalParticleContainer::Evolve (int lev,
				   const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
				   const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
				   MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                   MultiFab* rho, Real t, Real dt)
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

    long ngRho = (rho) ? rho->nGrow() : 0;

    long ngRhoDeposit = (WarpX::use_filter) ? ngRho +1 : ngRho;
    long ngJDeposit   = (WarpX::use_filter) ? ngJ +1   : ngJ;

    BL_ASSERT(OnSameGrids(lev,Ex));

    MultiFab* cost = WarpX::getCosts(lev);

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

            if (rho)
            {
                Real* data_ptr;
                const int *rholen;
                FArrayBox& rhofab = (*rho)[pti];
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
            	warpx_charge_deposition(data_ptr, &np,
                                        xp.data(), yp.data(), zp.data(), wp.data(),
                                        &this->charge,
                                        &xyzmin[0], &xyzmin[1], &xyzmin[2],
                                        &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
                                        &ngRhoDeposit, &ngRhoDeposit, &ngRhoDeposit,
                                        &WarpX::nox,&WarpX::noy,&WarpX::noz,
                                        &lvect, &WarpX::charge_deposition_algo);

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
            }

            if (! do_not_push)
            {
                //
                // Field Gather of Aux Data (i.e., the full solution)
                //
                const int ll4symtry          = false;
                const int l_lower_order_in_v = warpx_l_lower_order_in_v();
                long lvect_fieldgathe = 64;
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
                    &ll4symtry, &l_lower_order_in_v,
                    &lvect_fieldgathe, &WarpX::field_gathering_algo);
                BL_PROFILE_VAR_STOP(blp_pxr_fg);

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
                    &lvect,&WarpX::current_deposition_algo);

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

            if (cost) {
                const Box& tbx = pti.tilebox();
                wt = (ParallelDescriptor::second() - wt) / tbx.d_numPts();
                (*cost)[pti].plus(wt, tbx);
            }
	}
    }
}

void
PhysicalParticleContainer::PushP (int lev, Real dt,
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
            
            warpx_particle_pusher_momenta(&np, xp.data(), yp.data(), zp.data(),
                                          uxp.data(), uyp.data(), uzp.data(), giv.data(),
                                          Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
                                          Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
                                          &this->charge, &this->mass, &dt,
                                          &WarpX::particle_pusher_algo);
        }
    }
}
