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

    ParmParse pp(species_name);

    pp.query("boost_adjust_transverse_positions", boost_adjust_transverse_positions);
    pp.query("do_backward_propagation", do_backward_propagation);
}

void PhysicalParticleContainer::InitData()
{
    AddParticles(0); // Note - add on level 0
    if (maxLevel() > 0) {
        Redistribute();  // We then redistribute
    }
}

void PhysicalParticleContainer::MapParticletoBoostedFrame(Real& x, Real& y, Real& z, std::array<Real, 3>& u)
{
    // Map the particles from the lab frame to the boosted frame.
    // This boosts the particle to the lab frame and calculates
    // the particle time in the boosted frame. It then maps
    // the position to the time in the boosted frame.

    // For now, start with the assumption that this will only happen
    // at the start of the simulation.
    const Real t_lab = 0.;

    const Real uz_boost = WarpX::gamma_boost*WarpX::beta_boost*PhysConst::c;

    // tpr is the particle's time in the boosted frame
    Real tpr = WarpX::gamma_boost*t_lab - uz_boost*z/(PhysConst::c*PhysConst::c);

    // The particle's transformed location in the boosted frame
    Real xpr = x;
    Real ypr = y;
    Real zpr = WarpX::gamma_boost*z - uz_boost*t_lab;

    // transform u and gamma to the boosted frame
    Real gamma_lab = std::sqrt(1. + (u[0]*u[0] + u[1]*u[1] + u[2]*u[2])/(PhysConst::c*PhysConst::c));
    // u[0] = u[0];
    // u[1] = u[1];
    u[2] = WarpX::gamma_boost*u[2] - uz_boost*gamma_lab;
    Real gammapr = std::sqrt(1. + (u[0]*u[0] + u[1]*u[1] + u[2]*u[2])/(PhysConst::c*PhysConst::c));

    Real vxpr = u[0]/gammapr;
    Real vypr = u[1]/gammapr;
    Real vzpr = u[2]/gammapr;

    if (do_backward_propagation){
        u[2] = -u[2];
    }

    // Move the particles to where they will be at t = 0 in the boosted frame
    if (boost_adjust_transverse_positions) {
      x = xpr - tpr*vxpr;
      y = ypr - tpr*vypr;
    }

    z = zpr - tpr*vzpr;

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
#if ( AMREX_SPACEDIM == 3 )
            weight = q_tot/npart/charge;
            Real x = distx(mt);
            Real y = disty(mt);
            Real z = distz(mt);
#elif ( AMREX_SPACEDIM == 2 )
            weight = q_tot/npart/charge/y_rms;
            Real x = distx(mt);
            Real y = 0.;
            Real z = distz(mt);
#endif
        if (plasma_injector->insideBounds(x, y, z)) {
	    plasma_injector->getMomentum(u, x, y, z);
            if (WarpX::gamma_boost > 1.) {
                MapParticletoBoostedFrame(x, y, z, u);
            }
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
PhysicalParticleContainer::AddPlasma(int lev, RealBox part_realbox)
{
    BL_PROFILE("PhysicalParticleContainer::AddPlasma");

    // If no part_realbox is provided, initialize particles in the whole domain
    const Geometry& geom = Geom(lev);
    if (!part_realbox.ok()) part_realbox = geom.ProbDomain();

    int num_ppc = plasma_injector->num_particles_per_cell;

    const Real* dx = geom.CellSize();

    Real scale_fac;
#if AMREX_SPACEDIM==3
    scale_fac = dx[0]*dx[1]*dx[2]/num_ppc;
#elif AMREX_SPACEDIM==2
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

    MultiFab* cost = WarpX::getCosts(lev);

    if ( (not m_refined_injection_mask) and WarpX::do_moving_window)
    {
        Box mask_box = geom.Domain();
        mask_box.setSmall(WarpX::moving_window_dir, 0);
        mask_box.setBig(WarpX::moving_window_dir, 0);
        m_refined_injection_mask.reset( new IArrayBox(mask_box));
        m_refined_injection_mask->setVal(-1);
    }

    MFItInfo info;
    if (do_tiling) {
        info.EnableTiling(tile_size);
    }
    info.SetDynamic(true);

#ifdef _OPENMP
#pragma omp parallel if (not WarpX::serialize_ics)
#endif
    {
        std::array<Real,PIdx::nattribs> attribs;
        attribs.fill(0.0);

        // Loop through the tiles
        for (MFIter mfi = MakeMFIter(lev, info); mfi.isValid(); ++mfi) {

            Real wt = amrex::second();

            const Box& tile_box = mfi.tilebox();
            const RealBox tile_realbox = WarpX::getRealBox(tile_box, lev);

            // Find the cells of part_box that overlap with tile_realbox
            // If there is no overlap, just go to the next tile in the loop
            RealBox overlap_realbox;
            Box overlap_box;
            Real ncells_adjust;
            bool no_overlap = 0;

            for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
                if ( tile_realbox.lo(dir) <= part_realbox.hi(dir) ) {
                    ncells_adjust = std::floor( (tile_realbox.lo(dir) - part_realbox.lo(dir))/dx[dir] );
                    overlap_realbox.setLo( dir, part_realbox.lo(dir) + std::max(ncells_adjust, 0.) * dx[dir]);
                } else {
                    no_overlap = 1; break;
                }
                if ( tile_realbox.hi(dir) >= part_realbox.lo(dir) ) {
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
            if (no_overlap == 1) {
                continue; // Go to the next tile
            }

            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();

            // Loop through the cells of overlap_box and inject
            // the corresponding particles
            const auto& overlap_corner = overlap_realbox.lo();
            for (IntVect iv = overlap_box.smallEnd(); iv <= overlap_box.bigEnd(); overlap_box.next(iv))
            {
                int fac;
                if (injected) {
#if ( AMREX_SPACEDIM == 3 )
                    Real x = overlap_corner[0] + (iv[0] + 0.5)*dx[0];
                    Real y = overlap_corner[1] + (iv[1] + 0.5)*dx[1];
                    Real z = overlap_corner[2] + (iv[2] + 0.5)*dx[2];
#elif ( AMREX_SPACEDIM == 2 )
                    Real x = overlap_corner[0] + (iv[0] + 0.5)*dx[0];
                    Real y = 0;
                    Real z = overlap_corner[1] + (iv[1] + 0.5)*dx[1];
#endif
                    fac = GetRefineFac(x, y, z);
                } else {
                    fac = 1.0;
                }

                int ref_num_ppc = num_ppc * AMREX_D_TERM(fac, *fac, *fac);
                for (int i_part=0; i_part<ref_num_ppc;i_part++) {
                    std::array<Real, 3> r;
                    plasma_injector->getPositionUnitBox(r, i_part, fac);
#if ( AMREX_SPACEDIM == 3 )
                    Real x = overlap_corner[0] + (iv[0] + r[0])*dx[0];
                    Real y = overlap_corner[1] + (iv[1] + r[1])*dx[1];
                    Real z = overlap_corner[2] + (iv[2] + r[2])*dx[2];
#elif ( AMREX_SPACEDIM == 2 )
                    Real x = overlap_corner[0] + (iv[0] + r[0])*dx[0];
                    Real y = 0;
                    Real z = overlap_corner[1] + (iv[1] + r[1])*dx[1];
#endif
                    // If the new particle is not inside the tile box,
                    // go to the next generated particle.
#if ( AMREX_SPACEDIM == 3 )
                    if(!tile_realbox.contains( RealVect{x, y, z} )) continue;
#elif ( AMREX_SPACEDIM == 2 )
                    if(!tile_realbox.contains( RealVect{x, z} )) continue;
#endif

                    Real dens;
                    std::array<Real, 3> u;
                    if (WarpX::gamma_boost == 1.){
                      // Lab-frame simulation
                      // If the particle is not within the species's
                      // xmin, xmax, ymin, ymax, zmin, zmax, go to
                      // the next generated particle.
                      if (!plasma_injector->insideBounds(x, y, z)) continue;
                      plasma_injector->getMomentum(u, x, y, z);
                      dens = plasma_injector->getDensity(x, y, z);
                    } else {
                      // Boosted-frame simulation
                      Real c = PhysConst::c;
                      Real gamma_boost = WarpX::gamma_boost;
                      Real beta_boost = WarpX::beta_boost;
                      // Since the user provides the density distribution
                      // at t_lab=0 and in the lab-frame coordinates,
                      // we need to find the lab-frame position of this
                      // particle at t_lab=0, from its boosted-frame coordinates
                      // Assuming ballistic motion, this is given by:
                      // z0_lab = gamma*( z_boost*(1-beta*betaz_lab) - ct_boost*(betaz_lab-beta) )
                      // where betaz_lab is the speed of the particle in the lab frame
                      //
                      // In order for this equation to be solvable, betaz_lab
                      // is explicitly assumed to have no dependency on z0_lab
                      plasma_injector->getMomentum(u, x, y, 0.); // No z0_lab dependency
                      // At this point u is the lab-frame momentum
                      // => Apply the above formula for z0_lab
                      Real gamma_lab = std::sqrt( 1 + (u[0]*u[0] + u[1]*u[1] + u[2]*u[2])/(c*c) );
                      Real betaz_lab = u[2]/gamma_lab/c;
                      Real t = WarpX::GetInstance().gett_new(lev);
                      Real z0_lab = gamma_boost * ( z*(1-beta_boost*betaz_lab) - c*t*(betaz_lab-beta_boost) );
                      // If the particle is not within the lab-frame zmin, zmax, etc.
                      // go to the next generated particle.
                      if (!plasma_injector->insideBounds(x, y, z0_lab)) continue;
                      // call `getDensity` with lab-frame parameters
                      dens = plasma_injector->getDensity(x, y, z0_lab);
                      // At this point u and dens are the lab-frame quantities
                      // => Perform Lorentz transform
                      dens = gamma_boost * dens * ( 1 - beta_boost*betaz_lab );
                      u[2] = gamma_boost * ( u[2] -beta_boost*c*gamma_lab );
                    }
                    attribs[PIdx::w ] = dens * scale_fac / (AMREX_D_TERM(fac, *fac, *fac));
                    attribs[PIdx::ux] = u[0];
                    attribs[PIdx::uy] = u[1];
                    attribs[PIdx::uz] = u[2];

#ifdef WARPX_STORE_OLD_PARTICLE_ATTRIBS
                    attribs[PIdx::xold] = x;
                    attribs[PIdx::yold] = y;
                    attribs[PIdx::zold] = z;

                    attribs[PIdx::uxold] = u[0];
                    attribs[PIdx::uyold] = u[1];
                    attribs[PIdx::uzold] = u[2];
#endif

                    AddOneParticle(lev, grid_id, tile_id, x, y, z, attribs);
                }
            }

            if (cost) {
                wt = (amrex::second() - wt) / tile_box.d_numPts();
                (*cost)[mfi].plus(wt, tile_box);
            }
        }
    }
}

#ifdef WARPX_DO_ELECTROSTATIC
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
#if AMREX_SPACEDIM == 3
            auto& Ezp = attribs[PIdx::Ez];
#endif
            Exp.assign(np,0.0);
            Eyp.assign(np,0.0);
#if AMREX_SPACEDIM == 3
            Ezp.assign(np,0.0);
#endif

            const FArrayBox& exfab = (*E[lev][0])[pti];
            const FArrayBox& eyfab = (*E[lev][1])[pti];
#if AMREX_SPACEDIM == 3
            const FArrayBox& ezfab = (*E[lev][2])[pti];
#endif

            WRPX_INTERPOLATE_CIC(particles.data(), nstride, np,
                                 Exp.data(), Eyp.data(),
#if AMREX_SPACEDIM == 3
                                 Ezp.data(),
#endif
                                 exfab.dataPtr(), eyfab.dataPtr(),
#if AMREX_SPACEDIM == 3
                                 ezfab.dataPtr(),
#endif
                                 box.loVect(), box.hiVect(), plo, dx, &ng);
        }

        return;
    }

    const BoxArray& fine_BA = E[1][0]->boxArray();
    const DistributionMapping& fine_dm = E[1][0]->DistributionMap();
    BoxArray coarsened_fine_BA = fine_BA;
    coarsened_fine_BA.coarsen(IntVect(AMREX_D_DECL(2,2,2)));

    MultiFab coarse_Ex(coarsened_fine_BA, fine_dm, 1, 1);
    MultiFab coarse_Ey(coarsened_fine_BA, fine_dm, 1, 1);
#if AMREX_SPACEDIM == 3
    MultiFab coarse_Ez(coarsened_fine_BA, fine_dm, 1, 1);
#endif

    coarse_Ex.copy(*E[0][0], 0, 0, 1, 1, 1);
    coarse_Ey.copy(*E[0][1], 0, 0, 1, 1, 1);
#if AMREX_SPACEDIM == 3
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
#if AMREX_SPACEDIM == 3
            auto& Ezp = attribs[PIdx::Ez];
#endif
            Exp.assign(np,0.0);
            Eyp.assign(np,0.0);
#if AMREX_SPACEDIM == 3
            Ezp.assign(np,0.0);
#endif

            const FArrayBox& exfab = (*E[lev][0])[pti];
            const FArrayBox& eyfab = (*E[lev][1])[pti];
#if AMREX_SPACEDIM == 3
            const FArrayBox& ezfab = (*E[lev][2])[pti];
#endif

            if (lev == 0) {
                WRPX_INTERPOLATE_CIC(particles.data(), nstride, np,
                                     Exp.data(), Eyp.data(),
#if AMREX_SPACEDIM == 3
                Ezp.data(),
#endif
                                exfab.dataPtr(), eyfab.dataPtr(),
#if AMREX_SPACEDIM == 3
                                ezfab.dataPtr(),
#endif
                                box.loVect(), box.hiVect(), plo, dx, &ng);
            } else {

                const FArrayBox& exfab_coarse = coarse_Ex[pti];
                const FArrayBox& eyfab_coarse = coarse_Ey[pti];
#if AMREX_SPACEDIM == 3
                const FArrayBox& ezfab_coarse = coarse_Ez[pti];
#endif
                const Box& coarse_box = coarsened_fine_BA[pti];
                const Real* coarse_dx = Geom(0).CellSize();

                WRPX_INTERPOLATE_CIC_TWO_LEVELS(particles.data(), nstride, np,
                                                Exp.data(), Eyp.data(),
#if AMREX_SPACEDIM == 3
                                                Ezp.data(),
#endif
                                                exfab.dataPtr(), eyfab.dataPtr(),
#if AMREX_SPACEDIM == 3
                                                ezfab.dataPtr(),
#endif
                                                box.loVect(), box.hiVect(), dx,
                                                exfab_coarse.dataPtr(), eyfab_coarse.dataPtr(),
#if AMREX_SPACEDIM == 3
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

#if AMREX_SPACEDIM == 3
            auto& uzp = attribs[PIdx::uz];
#endif

            auto& Exp = attribs[PIdx::Ex];
            auto& Eyp = attribs[PIdx::Ey];

#if AMREX_SPACEDIM == 3
            auto& Ezp = attribs[PIdx::Ez];
#endif
            //
            // Particle Push
            //
            WRPX_PUSH_LEAPFROG(particles.data(), nstride, np,
                               uxp.data(), uyp.data(),
#if AMREX_SPACEDIM == 3
                               uzp.data(),
#endif
                               Exp.data(), Eyp.data(),
#if AMREX_SPACEDIM == 3
                               Ezp.data(),
#endif
                               &this->charge, &this->mass, &dt,
                               prob_domain.lo(), prob_domain.hi());
        }
    }
}
#endif // WARPX_DO_ELECTROSTATIC

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
            Real wt = amrex::second();

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
            const int* ixyzmin = box.loVect();

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
               ixyzmin,
               &xyzmin[0], &xyzmin[1], &xyzmin[2],
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

            if (cost) {
                const Box& tbx = pti.tilebox();
                wt = (amrex::second() - wt) / tbx.d_numPts();
                (*cost)[pti].plus(wt, tbx);
            }
        }
    }
}

void
PhysicalParticleContainer::Evolve (int lev,
				   const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
				   const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
				   MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                   MultiFab* cjx, MultiFab* cjy, MultiFab* cjz,
                                   MultiFab* rho,
                                   const MultiFab* cEx, const MultiFab* cEy, const MultiFab* cEz,
                                   const MultiFab* cBx, const MultiFab* cBy, const MultiFab* cBz,
                                   Real t, Real dt)
{
    BL_PROFILE("PPC::Evolve()");
    BL_PROFILE_VAR_NS("PPC::Evolve::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PICSAR::FieldGather", blp_pxr_fg);
    BL_PROFILE_VAR_NS("PICSAR::ParticlePush", blp_pxr_pp);
    BL_PROFILE_VAR_NS("PICSAR::CurrentDeposition", blp_pxr_cd);
    BL_PROFILE_VAR_NS("PPC::Evolve::Accumulate", blp_accumulate);
    BL_PROFILE_VAR_NS("PPC::Evolve::partition", blp_partition);

    const std::array<Real,3>& dx = WarpX::CellSize(lev);
    const std::array<Real,3>& cdx = WarpX::CellSize(std::max(lev-1,0));

    const auto& mypc = WarpX::GetInstance().GetPartContainer();
    const int nstencilz_fdtd_nci_corr = mypc.nstencilz_fdtd_nci_corr;

    // WarpX assumes the same number of guard cells for Jx, Jy, Jz
    long ngJ = jx.nGrow();

    BL_ASSERT(OnSameGrids(lev,Ex));

    MultiFab* cost = WarpX::getCosts(lev);

    const iMultiFab* current_masks = WarpX::CurrentBufferMasks(lev);
    const iMultiFab* gather_masks = WarpX::GatherBufferMasks(lev);

    bool has_buffer = cEx || cjx;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	Vector<Real> xp, yp, zp, giv;
        FArrayBox local_rho, local_jx, local_jy, local_jz;
        FArrayBox filtered_Ex, filtered_Ey, filtered_Ez;
        FArrayBox filtered_Bx, filtered_By, filtered_Bz;
        std::vector<bool> inexflag;
        Vector<long> pid;
        Vector<Real> tmp;
        Vector<ParticleType> particle_tmp;

	for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
	{
            Real wt = amrex::second();

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
            FArrayBox const* exfab = &(Ex[pti]);
            FArrayBox const* eyfab = &(Ey[pti]);
            FArrayBox const* ezfab = &(Ez[pti]);
            FArrayBox const* bxfab = &(Bx[pti]);
            FArrayBox const* byfab = &(By[pti]);
            FArrayBox const* bzfab = &(Bz[pti]);

            if (warpx_use_fdtd_nci_corr())
            {
#if (AMREX_SPACEDIM == 2)
                const Box& tbox = amrex::grow(pti.tilebox(),{static_cast<int>(WarpX::nox),
                            static_cast<int>(WarpX::noz)});
#else
                const Box& tbox = amrex::grow(pti.tilebox(),{static_cast<int>(WarpX::nox),
                            static_cast<int>(WarpX::noy),
                            static_cast<int>(WarpX::noz)});
#endif

                // both 2d and 3d
                filtered_Ex.resize(amrex::convert(tbox,WarpX::Ex_nodal_flag));
                WRPX_PXR_GODFREY_FILTER(BL_TO_FORTRAN_BOX(filtered_Ex),
                                        BL_TO_FORTRAN_ANYD(filtered_Ex),
                                        BL_TO_FORTRAN_ANYD(Ex[pti]),
                                        mypc.fdtd_nci_stencilz_ex.data(),
                                        &nstencilz_fdtd_nci_corr);
                exfab = &filtered_Ex;

                filtered_Ez.resize(amrex::convert(tbox,WarpX::Ez_nodal_flag));
                WRPX_PXR_GODFREY_FILTER(BL_TO_FORTRAN_BOX(filtered_Ez),
                                        BL_TO_FORTRAN_ANYD(filtered_Ez),
                                        BL_TO_FORTRAN_ANYD(Ez[pti]),
                                        mypc.fdtd_nci_stencilz_by.data(),
                                        &nstencilz_fdtd_nci_corr);
                ezfab = &filtered_Ez;

                filtered_By.resize(amrex::convert(tbox,WarpX::By_nodal_flag));
                WRPX_PXR_GODFREY_FILTER(BL_TO_FORTRAN_BOX(filtered_By),
                                        BL_TO_FORTRAN_ANYD(filtered_By),
                                        BL_TO_FORTRAN_ANYD(By[pti]),
                                        mypc.fdtd_nci_stencilz_by.data(),
                                        &nstencilz_fdtd_nci_corr);
                byfab = &filtered_By;

#if (AMREX_SPACEDIM == 3)
                filtered_Ey.resize(amrex::convert(tbox,WarpX::Ey_nodal_flag));
                WRPX_PXR_GODFREY_FILTER(BL_TO_FORTRAN_BOX(filtered_Ey),
                                        BL_TO_FORTRAN_ANYD(filtered_Ey),
                                        BL_TO_FORTRAN_ANYD(Ey[pti]),
                                        mypc.fdtd_nci_stencilz_ex.data(),
                                        &nstencilz_fdtd_nci_corr);
                eyfab = &filtered_Ey;

                filtered_Bx.resize(amrex::convert(tbox,WarpX::Bx_nodal_flag));
                WRPX_PXR_GODFREY_FILTER(BL_TO_FORTRAN_BOX(filtered_Bx),
                                        BL_TO_FORTRAN_ANYD(filtered_Bx),
                                        BL_TO_FORTRAN_ANYD(Bx[pti]),
                                        mypc.fdtd_nci_stencilz_by.data(),
                                        &nstencilz_fdtd_nci_corr);
                bxfab = &filtered_Bx;

                filtered_Bz.resize(amrex::convert(tbox,WarpX::Bz_nodal_flag));
                WRPX_PXR_GODFREY_FILTER(BL_TO_FORTRAN_BOX(filtered_Bz),
                                        BL_TO_FORTRAN_ANYD(filtered_Bz),
                                        BL_TO_FORTRAN_ANYD(Bz[pti]),
                                        mypc.fdtd_nci_stencilz_ex.data(),
                                        &nstencilz_fdtd_nci_corr);
                bzfab = &filtered_Bz;
#endif
            }

	    FArrayBox& jxfab = jx[pti];
	    FArrayBox& jyfab = jy[pti];
	    FArrayBox& jzfab = jz[pti];

	    Exp.assign(np,0.0);
	    Eyp.assign(np,0.0);
	    Ezp.assign(np,0.0);
	    Bxp.assign(np,WarpX::B_external[0]);
	    Byp.assign(np,WarpX::B_external[1]);
	    Bzp.assign(np,WarpX::B_external[2]);

	    giv.resize(np);

            long nfine_current = np;
            long nfine_gather = np;
            if (has_buffer && !do_not_push)
            {
                BL_PROFILE_VAR_START(blp_partition);
                inexflag.resize(np);
                auto& aos = pti.GetArrayOfStructs();
                // We need to partition the large buffer first
                iMultiFab const* bmasks = (WarpX::n_field_gather_buffer >= WarpX::n_current_deposition_buffer) ?
                    gather_masks : current_masks;
                int i = 0;
                const auto& msk = (*bmasks)[pti];
                for (const auto& p : aos) {
                    const IntVect& iv = Index(p, lev);
                    inexflag[i++] = msk(iv);
                }

                pid.resize(np);
                std::iota(pid.begin(), pid.end(), 0L);

                auto sep = std::stable_partition(pid.begin(), pid.end(),
                                                 [&inexflag](long id) { return inexflag[id]; });

                if (WarpX::n_current_deposition_buffer == WarpX::n_field_gather_buffer) {
                    nfine_current = nfine_gather = std::distance(pid.begin(), sep);
                } else if (sep != pid.end()) {
                    int n_buf;
                    if (bmasks == gather_masks) {
                        nfine_gather = std::distance(pid.begin(), sep);
                        bmasks = current_masks;
                        n_buf = WarpX::n_current_deposition_buffer;
                    } else {
                        nfine_current = std::distance(pid.begin(), sep);
                        bmasks = gather_masks;
                        n_buf = WarpX::n_field_gather_buffer;
                    }
                    if (n_buf > 0)
                    {
                        const auto& msk2 = (*bmasks)[pti];
                        for (auto it = sep; it != pid.end(); ++it) {
                            const long id = *it;
                            const IntVect& iv = Index(aos[id], lev);
                            inexflag[id] = msk2(iv);
                        }

                        auto sep2 = std::stable_partition(sep, pid.end(),
                                                          [&inexflag](long id) {return inexflag[id];});
                        if (bmasks == gather_masks) {
                            nfine_gather = std::distance(pid.begin(), sep2);
                        } else {
                            nfine_current = std::distance(pid.begin(), sep2);
                        }
                    }
                }

                if (deposit_on_main_grid && lev > 0) {
                    nfine_current = 0;
                }

                if (nfine_current != np || nfine_gather != np)
                {
                    particle_tmp.resize(np);
                    for (long ip = 0; ip < np; ++ip) {
                        particle_tmp[ip] = aos[pid[ip]];
                    }
                    std::swap(aos(), particle_tmp);

                    tmp.resize(np);
                    for (long ip = 0; ip < np; ++ip) {
                        tmp[ip] = wp[pid[ip]];
                    }
                    std::swap(wp, tmp);

                    for (long ip = 0; ip < np; ++ip) {
                        tmp[ip] = uxp[pid[ip]];
                    }
                    std::swap(uxp, tmp);

                    for (long ip = 0; ip < np; ++ip) {
                        tmp[ip] = uyp[pid[ip]];
                    }
                    std::swap(uyp, tmp);

                    for (long ip = 0; ip < np; ++ip) {
                        tmp[ip] = uzp[pid[ip]];
                    }
                    std::swap(uzp, tmp);
                }
                BL_PROFILE_VAR_STOP(blp_partition);
            }

	    //
	    // copy data from particle container to temp arrays
	    //
	    BL_PROFILE_VAR_START(blp_copy);
            pti.GetPosition(xp, yp, zp);
	    BL_PROFILE_VAR_STOP(blp_copy);

            const std::array<Real,3>& xyzmin_tile = WarpX::LowerCorner(pti.tilebox(), lev);
            const std::array<Real,3>& xyzmin_grid = WarpX::LowerCorner(box, lev);
            const int* ixyzmin_grid = box.loVect();

	    long lvect = 8;

            auto depositCharge = [&] (MultiFab* rhomf, int icomp)
            {
                long ngRho = rhomf->nGrow();

                Real* data_ptr;
                const int *rholen;
                FArrayBox& rhofab = (*rhomf)[pti];
                Box tile_box = convert(pti.tilebox(), IntVect::TheUnitVector());
                Box grown_box;
                const std::array<Real, 3>& xyzmin = xyzmin_tile;
                tile_box.grow(ngRho);
                local_rho.resize(tile_box);
                local_rho = 0.0;
                data_ptr = local_rho.dataPtr();
                rholen = local_rho.length();

#if (AMREX_SPACEDIM == 3)
                const long nx = rholen[0]-1-2*ngRho;
                const long ny = rholen[1]-1-2*ngRho;
                const long nz = rholen[2]-1-2*ngRho;
#else
                const long nx = rholen[0]-1-2*ngRho;
                const long ny = 0;
                const long nz = rholen[1]-1-2*ngRho;
#endif
            	warpx_charge_deposition(data_ptr, &np,
                                        xp.data(), yp.data(), zp.data(), wp.data(),
                                        &this->charge,
                                        &xyzmin[0], &xyzmin[1], &xyzmin[2],
                                        &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
                                        &ngRho, &ngRho, &ngRho,
                                        &WarpX::nox,&WarpX::noy,&WarpX::noz,
                                        &lvect, &WarpX::charge_deposition_algo);

                const int ncomp = 1;
                amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_rho),
                                            BL_TO_FORTRAN_N_3D(rhofab,icomp), ncomp);
            };

            if (rho) depositCharge(rho,0);

            if (! do_not_push)
            {
                //
                // Field Gather of Aux Data (i.e., the full solution)
                //
                const int ll4symtry          = false;
                const int l_lower_order_in_v = warpx_l_lower_order_in_v();
                long lvect_fieldgathe = 64;

                const long np_gather = (cEx) ? nfine_gather : np;

                BL_PROFILE_VAR_START(blp_pxr_fg);

                warpx_geteb_energy_conserving(
                    &np_gather, xp.data(), yp.data(), zp.data(),
                    Exp.data(),Eyp.data(),Ezp.data(),
                    Bxp.data(),Byp.data(),Bzp.data(),
                    ixyzmin_grid,
                    &xyzmin_grid[0], &xyzmin_grid[1], &xyzmin_grid[2],
                    &dx[0], &dx[1], &dx[2],
                    &WarpX::nox, &WarpX::noy, &WarpX::noz,
                    BL_TO_FORTRAN_ANYD(*exfab),
                    BL_TO_FORTRAN_ANYD(*eyfab),
                    BL_TO_FORTRAN_ANYD(*ezfab),
                    BL_TO_FORTRAN_ANYD(*bxfab),
                    BL_TO_FORTRAN_ANYD(*byfab),
                    BL_TO_FORTRAN_ANYD(*bzfab),
                    &ll4symtry, &l_lower_order_in_v,
                    &lvect_fieldgathe, &WarpX::field_gathering_algo);

                if (np_gather < np)
                {
                    const IntVect& ref_ratio = WarpX::RefRatio(lev-1);
                    const Box& cbox = amrex::coarsen(box,ref_ratio);
                    const std::array<Real,3>& cxyzmin_grid = WarpX::LowerCorner(cbox, lev-1);
                    const int* cixyzmin_grid = cbox.loVect();

                    const FArrayBox& cexfab = (*cEx)[pti];
                    const FArrayBox& ceyfab = (*cEy)[pti];
                    const FArrayBox& cezfab = (*cEz)[pti];
                    const FArrayBox& cbxfab = (*cBx)[pti];
                    const FArrayBox& cbyfab = (*cBy)[pti];
                    const FArrayBox& cbzfab = (*cBz)[pti];

                    long ncrse = np - nfine_gather;
                    warpx_geteb_energy_conserving(
                        &ncrse, xp.data()+nfine_gather, yp.data()+nfine_gather, zp.data()+nfine_gather,
                        Exp.data()+nfine_gather, Eyp.data()+nfine_gather, Ezp.data()+nfine_gather,
                        Bxp.data()+nfine_gather, Byp.data()+nfine_gather, Bzp.data()+nfine_gather,
                        cixyzmin_grid,
                        &cxyzmin_grid[0], &cxyzmin_grid[1], &cxyzmin_grid[2],
                        &cdx[0], &cdx[1], &cdx[2],
                        &WarpX::nox, &WarpX::noy, &WarpX::noz,
                        BL_TO_FORTRAN_ANYD(cexfab),
                        BL_TO_FORTRAN_ANYD(ceyfab),
                        BL_TO_FORTRAN_ANYD(cezfab),
                        BL_TO_FORTRAN_ANYD(cbxfab),
                        BL_TO_FORTRAN_ANYD(cbyfab),
                        BL_TO_FORTRAN_ANYD(cbzfab),
                        &ll4symtry, &l_lower_order_in_v,
                        &lvect_fieldgathe, &WarpX::field_gathering_algo);
                }

                BL_PROFILE_VAR_STOP(blp_pxr_fg);

                //
                // Particle Push
                //
                BL_PROFILE_VAR_START(blp_pxr_pp);
                PushPX(pti, xp, yp, zp, giv, dt);
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

                const long np_current = (cjx) ? nfine_current : np;

                const std::array<Real, 3>& xyzmin = xyzmin_tile;

                if (np_current > 0)
                {
                    tbx.grow(ngJ);
                    tby.grow(ngJ);
                    tbz.grow(ngJ);
                    
                    local_jx.resize(tbx);
                    local_jy.resize(tby);
                    local_jz.resize(tbz);

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
                        jx_ptr, &ngJ, jxntot,
                        jy_ptr, &ngJ, jyntot,
                        jz_ptr, &ngJ, jzntot,
                        &np_current, xp.data(), yp.data(), zp.data(),
                        uxp.data(), uyp.data(), uzp.data(),
                        giv.data(), wp.data(), &this->charge,
                        &xyzmin[0], &xyzmin[1], &xyzmin[2],
                        &dt, &dx[0], &dx[1], &dx[2],
                        &WarpX::nox,&WarpX::noy,&WarpX::noz,
                        &lvect,&WarpX::current_deposition_algo);

                    BL_PROFILE_VAR_STOP(blp_pxr_cd);

                    BL_PROFILE_VAR_START(blp_accumulate);
                    const int ncomp = 1;
                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_jx),
                                                BL_TO_FORTRAN_3D(jxfab), ncomp);
                    
                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_jy),
                                                BL_TO_FORTRAN_3D(jyfab), ncomp);

                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_jz),
                                                BL_TO_FORTRAN_3D(jzfab), ncomp);
                    BL_PROFILE_VAR_STOP(blp_accumulate);
                }

                if (np_current < np)
                {
                    const IntVect& ref_ratio = WarpX::RefRatio(lev-1);
                    const Box& ctilebox = amrex::coarsen(pti.tilebox(),ref_ratio);
                    const std::array<Real,3>& cxyzmin_tile = WarpX::LowerCorner(ctilebox, lev-1);

                    tbx = amrex::convert(ctilebox, WarpX::jx_nodal_flag);
                    tby = amrex::convert(ctilebox, WarpX::jy_nodal_flag);
                    tbz = amrex::convert(ctilebox, WarpX::jz_nodal_flag);
                    tbx.grow(ngJ);
                    tby.grow(ngJ);
                    tbz.grow(ngJ);

                    local_jx.resize(tbx);
                    local_jy.resize(tby);
                    local_jz.resize(tbz);

                    local_jx = 0.0;
                    local_jy = 0.0;
                    local_jz = 0.0;

                    jx_ptr = local_jx.dataPtr();
                    jy_ptr = local_jy.dataPtr();
                    jz_ptr = local_jz.dataPtr();

                    jxntot = local_jx.length();
                    jyntot = local_jy.length();
                    jzntot = local_jz.length();

                    long ncrse = np - nfine_current;
                    warpx_current_deposition(
                        jx_ptr, &ngJ, jxntot,
                        jy_ptr, &ngJ, jyntot,
                        jz_ptr, &ngJ, jzntot,
                        &ncrse, xp.data()+nfine_current, yp.data()+nfine_current, zp.data()+nfine_current,
                        uxp.data()+nfine_current, uyp.data()+nfine_current, uzp.data()+nfine_current,
                        giv.data()+nfine_current, wp.data()+nfine_current, &this->charge,
                        &cxyzmin_tile[0], &cxyzmin_tile[1], &cxyzmin_tile[2],
                        &dt, &cdx[0], &cdx[1], &cdx[2],
                        &WarpX::nox,&WarpX::noy,&WarpX::noz,
                        &lvect,&WarpX::current_deposition_algo);

                    FArrayBox& cjxfab = (*cjx)[pti];
                    FArrayBox& cjyfab = (*cjy)[pti];
                    FArrayBox& cjzfab = (*cjz)[pti];

                    const int ncomp = 1;
                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_jx),
                                                BL_TO_FORTRAN_3D(cjxfab), ncomp);
                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_jy),
                                                BL_TO_FORTRAN_3D(cjyfab), ncomp);
                    amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_jz),
                                                BL_TO_FORTRAN_3D(cjzfab), ncomp);
                }

                //
                // copy particle data back
                //
                BL_PROFILE_VAR_START(blp_copy);
                pti.SetPosition(xp, yp, zp);
                BL_PROFILE_VAR_STOP(blp_copy);
            }

            if (rho) depositCharge(rho,1);

            if (cost) {
                const Box& tbx = pti.tilebox();
                wt = (amrex::second() - wt) / tbx.d_numPts();
                (*cost)[pti].plus(wt, tbx);
            }
        }
    }
}

void
PhysicalParticleContainer::PushPX(WarpXParIter& pti,
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

#ifdef WARPX_STORE_OLD_PARTICLE_ATTRIBS
    auto& xpold  = attribs[PIdx::xold];
    auto& ypold  = attribs[PIdx::yold];
    auto& zpold  = attribs[PIdx::zold];
    auto& uxpold = attribs[PIdx::uxold];
    auto& uypold = attribs[PIdx::uyold];
    auto& uzpold = attribs[PIdx::uzold];

    warpx_copy_attribs(&np, xp.data(), yp.data(), zp.data(),
                       uxp.data(), uyp.data(), uzp.data(),
                       xpold.data(), ypold.data(), zpold.data(),
                       uxpold.data(), uypold.data(), uzpold.data());

#endif

    warpx_particle_pusher(&np, xp.data(), yp.data(), zp.data(),
                          uxp.data(), uyp.data(), uzp.data(), giv.data(),
                          Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
                          Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
                          &this->charge, &this->mass, &dt,
                          &WarpX::particle_pusher_algo);

}

void
PhysicalParticleContainer::PushP (int lev, Real dt,
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

            warpx_particle_pusher_momenta(&np, xp.data(), yp.data(), zp.data(),
                                          uxp.data(), uyp.data(), uzp.data(), giv.data(),
                                          Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
                                          Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
                                          &this->charge, &this->mass, &dt,
                                          &WarpX::particle_pusher_algo);
        }
    }
}

void PhysicalParticleContainer::GetParticleSlice(const int direction, const Real z_old,
                                                 const Real z_new, const Real t_boost,
                                                 const Real t_lab, const Real dt,
                                                 DiagnosticParticles& diagnostic_particles)
{
    BL_PROFILE("PhysicalParticleContainer::GetParticleSlice");

#ifdef WARPX_STORE_OLD_PARTICLE_ATTRIBS
    // Assume that the boost in the positive z direction.
#if (AMREX_SPACEDIM == 2)
    AMREX_ALWAYS_ASSERT(direction == 1);
#else
    AMREX_ALWAYS_ASSERT(direction == 2);
#endif

    // Note the the slice should always move in the negative boost direction.
    AMREX_ALWAYS_ASSERT(z_new < z_old);

    const int nlevs = std::max(0, finestLevel()+1);

    // we figure out a box for coarse-grained rejection. If the RealBox corresponding to a
    // given tile doesn't intersect with this, there is no need to check any particles.
    const Real* base_dx = Geom(0).CellSize();
    const Real z_min = z_new - base_dx[direction];
    const Real z_max = z_old + base_dx[direction];

    RealBox slice_box = Geom(0).ProbDomain();
    slice_box.setLo(direction, z_min);
    slice_box.setHi(direction, z_max);

    for (int lev = 0; lev < nlevs; ++lev) {

        const Real* dx  = Geom(lev).CellSize();
        const Real* plo = Geom(lev).ProbLo();

        // first we touch each map entry in serial
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
            diagnostic_particles[index];
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            Vector<Real> xp_new, yp_new, zp_new;

            for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
            {
                const Box& box = pti.validbox();

                auto index = std::make_pair(pti.index(), pti.LocalTileIndex());

                const RealBox tile_real_box(box, dx, plo);

                if ( !slice_box.intersects(tile_real_box) ) continue;

                pti.GetPosition(xp_new, yp_new, zp_new);

                auto& attribs = pti.GetAttribs();

                auto& wp = attribs[PIdx::w ];

                auto& uxp_new = attribs[PIdx::ux   ];
                auto& uyp_new = attribs[PIdx::uy   ];
                auto& uzp_new = attribs[PIdx::uz   ];

                auto&  xp_old = attribs[PIdx::xold ];
                auto&  yp_old = attribs[PIdx::yold ];
                auto&  zp_old = attribs[PIdx::zold ];
                auto& uxp_old = attribs[PIdx::uxold];
                auto& uyp_old = attribs[PIdx::uyold];
                auto& uzp_old = attribs[PIdx::uzold];

                const long np = pti.numParticles();

                Real uzfrm = -WarpX::gamma_boost*WarpX::beta_boost*PhysConst::c;
                Real inv_c2 = 1.0/PhysConst::c/PhysConst::c;

                for (long i = 0; i < np; ++i) {

                    // if the particle did not cross the plane of z_boost in the last
                    // timestep, skip it.
                    if ( not (((zp_new[i] >= z_new) && (zp_old[i] <= z_old)) ||
                              ((zp_new[i] <= z_new) && (zp_old[i] >= z_old))) ) continue;

                    // Lorentz transform particles to lab frame
                    Real gamma_new_p = std::sqrt(1.0 + inv_c2*(uxp_new[i]*uxp_new[i] + uyp_new[i]*uyp_new[i] + uzp_new[i]*uzp_new[i]));
                    Real t_new_p = WarpX::gamma_boost*t_boost - uzfrm*zp_new[i]*inv_c2;
                    Real z_new_p = WarpX::gamma_boost*(zp_new[i] + WarpX::beta_boost*PhysConst::c*t_boost);
                    Real uz_new_p = WarpX::gamma_boost*uzp_new[i] - gamma_new_p*uzfrm;

                    Real gamma_old_p = std::sqrt(1.0 + inv_c2*(uxp_old[i]*uxp_old[i] + uyp_old[i]*uyp_old[i] + uzp_old[i]*uzp_old[i]));
                    Real t_old_p = WarpX::gamma_boost*(t_boost - dt) - uzfrm*zp_old[i]*inv_c2;
                    Real z_old_p = WarpX::gamma_boost*(zp_old[i] + WarpX::beta_boost*PhysConst::c*(t_boost-dt));
                    Real uz_old_p = WarpX::gamma_boost*uzp_old[i] - gamma_old_p*uzfrm;

                    // interpolate in time to t_lab
                    Real weight_old = (t_new_p - t_lab) / (t_new_p - t_old_p);
                    Real weight_new = (t_lab - t_old_p) / (t_new_p - t_old_p);

                    Real xp = xp_old[i]*weight_old + xp_new[i]*weight_new;
                    Real yp = yp_old[i]*weight_old + yp_new[i]*weight_new;
                    Real zp = z_old_p  *weight_old + z_new_p  *weight_new;

                    Real uxp = uxp_old[i]*weight_old + uxp_new[i]*weight_new;
                    Real uyp = uyp_old[i]*weight_old + uyp_new[i]*weight_new;
                    Real uzp = uz_old_p  *weight_old + uz_new_p  *weight_new;

                    diagnostic_particles[index].GetRealData(DiagIdx::w).push_back(wp[i]);

                    diagnostic_particles[index].GetRealData(DiagIdx::x).push_back(xp);
                    diagnostic_particles[index].GetRealData(DiagIdx::y).push_back(yp);
                    diagnostic_particles[index].GetRealData(DiagIdx::z).push_back(zp);

                    diagnostic_particles[index].GetRealData(DiagIdx::ux).push_back(uxp);
                    diagnostic_particles[index].GetRealData(DiagIdx::uy).push_back(uyp);
                    diagnostic_particles[index].GetRealData(DiagIdx::uz).push_back(uzp);
                }
            }
        }
    }
#else
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( false ,
"ERROR: WarpX must be compiled with STORE_OLD_PARTICLE_ATTRIBS=TRUE to use the back-transformed diagnostics");
#endif
}

int PhysicalParticleContainer::GetRefineFac(const Real x, const Real y, const Real z)
{
    if (finestLevel() == 0) return 1;
    if (not WarpX::refine_plasma) return 1;

    IntVect iv;
    const Geometry& geom = Geom(0);

    std::array<Real, 3> offset;

#if ( AMREX_SPACEDIM == 3)
    offset[0] = geom.ProbLo(0);
    offset[1] = geom.ProbLo(1);
    offset[2] = geom.ProbLo(2);
#elif ( AMREX_SPACEDIM == 2 )
    offset[0] = geom.ProbLo(0);
    offset[1] = 0.0;
    offset[2] = geom.ProbLo(1);
#endif

    AMREX_D_TERM(iv[0]=static_cast<int>(floor((x-offset[0])*geom.InvCellSize(0)));,
                 iv[1]=static_cast<int>(floor((y-offset[1])*geom.InvCellSize(1)));,
                 iv[2]=static_cast<int>(floor((z-offset[2])*geom.InvCellSize(2))););

    iv += geom.Domain().smallEnd();

    const int dir = WarpX::moving_window_dir;

    IntVect iv2 = iv;
    iv2[dir] = 0;

    if ( (*m_refined_injection_mask)(iv2) != -1) return (*m_refined_injection_mask)(iv2);

    int ref_fac = 1;
    for (int lev = 0; lev < finestLevel(); ++lev)
    {
        const IntVect rr = m_gdb->refRatio(lev);
        const BoxArray& fine_ba = this->ParticleBoxArray(lev+1);
        const int num_boxes = fine_ba.size();
        Vector<Box> stretched_boxes;
        const int safety_factor = 4;
        for (int i = 0; i < num_boxes; ++i)
        {
            Box bx = fine_ba[i];
            bx.coarsen(ref_fac*rr[dir]);
            bx.setSmall(dir, std::numeric_limits<int>::min()/safety_factor);
            bx.setBig(dir, std::numeric_limits<int>::max()/safety_factor);
            stretched_boxes.push_back(bx);
        }

        BoxArray stretched_ba(stretched_boxes.data(), stretched_boxes.size());

        const int num_ghost = 0;
        if ( stretched_ba.intersects(Box(iv, iv), num_ghost) )
        {
            ref_fac *= rr[dir];
        }
        else
        {
            break;
        }
    }

    (*m_refined_injection_mask)(iv2) = ref_fac;

    return ref_fac;
}
