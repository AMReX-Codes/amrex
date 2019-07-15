#include <limits>
#include <sstream>

#include <MultiParticleContainer.H>
#include <WarpX_f.H>
#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpXWrappers.h>


using namespace amrex;

long PhysicalParticleContainer::
NumParticlesToAdd(const Box& overlap_box, const RealBox& overlap_realbox,
		  const RealBox& tile_realbox, const RealBox& particle_real_box)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    int num_ppc = plasma_injector->num_particles_per_cell;
    const Real* dx = geom.CellSize();

    long np = 0;
    const auto& overlap_corner = overlap_realbox.lo();
    for (IntVect iv = overlap_box.smallEnd(); iv <= overlap_box.bigEnd(); overlap_box.next(iv))
    {
        int fac;
	if (do_continuous_injection) {
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
	    ++np;
	}
    }
    
    return np;
}

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
    pp.query("do_splitting", do_splitting);
    pp.query("split_type", split_type);
    pp.query("do_continuous_injection", do_continuous_injection);
    // Whether to plot back-transformed (lab-frame) diagnostics 
    // for this species.
    pp.query("do_boosted_frame_diags", do_boosted_frame_diags);

    pp.query("plot_species", plot_species);
    int do_user_plot_vars;
    do_user_plot_vars = pp.queryarr("plot_vars", plot_vars);
    if (not do_user_plot_vars){
        // By default, all particle variables are dumped to plotfiles,
        // including {x,y,z,ux,uy,uz}old variables when running in a 
        // boosted frame
        if (WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags){
            plot_flags.resize(PIdx::nattribs + 6, 1);
        } else {
            plot_flags.resize(PIdx::nattribs, 1);
        }
    } else {
        // Set plot_flag to 0 for all attribs
        if (WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags){
            plot_flags.resize(PIdx::nattribs + 6, 0);
        } else {
            plot_flags.resize(PIdx::nattribs, 0);
        }
        // If not none, set plot_flags values to 1 for elements in plot_vars.
        if (plot_vars[0] != "none"){
            for (const auto& var : plot_vars){
                // Return error if var not in PIdx. 
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE( 
                    ParticleStringNames::to_index.count(var), 
                    "plot_vars argument not in ParticleStringNames");
                plot_flags[ParticleStringNames::to_index.at(var)] = 1;
            }
        }
    }
}

PhysicalParticleContainer::PhysicalParticleContainer (AmrCore* amr_core)
    : WarpXParticleContainer(amr_core, 0)
{
    plasma_injector.reset(new PlasmaInjector());
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
                                           Real q_tot, long npart, 
                                           int do_symmetrize) {

    const Geometry& geom     = m_gdb->Geom(0);
    RealBox containing_bx = geom.ProbDomain();

    std::mt19937_64 mt(0451);
    std::normal_distribution<double> distx(x_m, x_rms);
    std::normal_distribution<double> disty(y_m, y_rms);
    std::normal_distribution<double> distz(z_m, z_rms);

    if (ParallelDescriptor::IOProcessor()) {
        std::array<Real, 3> u;
        Real weight;
        // If do_symmetrize, create 4x fewer particles, and 
        // Replicate each particle 4 times (x,y) (-x,y) (x,-y) (-x,-y)
        if (do_symmetrize){
            npart /= 4;
        }
        for (long i = 0; i < npart; ++i) {
#if ( AMREX_SPACEDIM == 3 | WARPX_RZ)
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
                if (do_symmetrize){
                    std::array<Real, 3> u_tmp;
                    Real x_tmp, y_tmp;
                    // Add four particles to the beam:
                    // (x,ux,y,uy) (-x,-ux,y,uy) (x,ux,-y,-uy) (-x,-ux,-y,-uy)
                    for (int ix=0; ix<2; ix++){
                        for (int iy=0; iy<2; iy++){
                            u_tmp = u;
                            x_tmp     = x*std::pow(-1,ix);
                            u_tmp[0] *= std::pow(-1,ix);
                            y_tmp     = y*std::pow(-1,iy);
                            u_tmp[1] *= std::pow(-1,iy);
                            CheckAndAddParticle(x_tmp, y_tmp, z, 
                                                u_tmp, weight/4);
                        }
                    }
                } else {
                    CheckAndAddParticle(x, y, z, u, weight);
                }
            }
        }
    }
    Redistribute();
}

void
PhysicalParticleContainer::CheckAndAddParticle(Real x, Real y, Real z,
                                               std::array<Real, 3> u,
                                               Real weight)
{
    std::array<Real,PIdx::nattribs> attribs;
    attribs.fill(0.0);

    // update attribs with input arguments
    if (WarpX::gamma_boost > 1.) {
        MapParticletoBoostedFrame(x, y, z, u);
    }
    attribs[PIdx::ux] = u[0];
    attribs[PIdx::uy] = u[1];
    attribs[PIdx::uz] = u[2];
    attribs[PIdx::w ] = weight;

    if (WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags)
    {
        // need to create old values
        auto& particle_tile = DefineAndReturnParticleTile(0, 0, 0);
        particle_tile.push_back_real(particle_comps["xold"], x);
        particle_tile.push_back_real(particle_comps["yold"], y);
        particle_tile.push_back_real(particle_comps["zold"], z);
                
        particle_tile.push_back_real(particle_comps["uxold"], u[0]);
        particle_tile.push_back_real(particle_comps["uyold"], u[1]);
        particle_tile.push_back_real(particle_comps["uzold"], u[2]);
    }
    // add particle
    AddOneParticle(0, 0, 0, x, y, z, attribs);
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
                        plasma_injector->npart,
                        plasma_injector->do_symmetrize);


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
PhysicalParticleContainer::AddPlasma (int lev, RealBox part_realbox)
{
#ifdef AMREX_USE_GPU
  AddPlasmaGPU(lev, part_realbox);
#else
  AddPlasmaCPU(lev, part_realbox);
#endif
}

void
PhysicalParticleContainer::AddPlasmaCPU (int lev, RealBox part_realbox)
{
    BL_PROFILE("PhysicalParticleContainer::AddPlasmaCPU");

    // If no part_realbox is provided, initialize particles in the whole domain
    const Geometry& geom = Geom(lev);
    if (!part_realbox.ok()) part_realbox = geom.ProbDomain();

    int num_ppc = plasma_injector->num_particles_per_cell;
#ifdef WARPX_RZ
    Real rmax = std::min(plasma_injector->xmax, part_realbox.hi(0));
#endif

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
                if (do_continuous_injection) {
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

                    // Save the x and y values to use in the insideBounds checks.
                    // This is needed with WARPX_RZ since x and y are modified.
                    Real xb = x;
                    Real yb = y;

#ifdef WARPX_RZ
                    // Replace the x and y, choosing the angle randomly.
                    // These x and y are used to get the momentum and density
                    Real theta = 2.*MathConst::pi*amrex::Random();
                    y = x*std::sin(theta);
                    x = x*std::cos(theta);
#endif

                    Real dens;
                    std::array<Real, 3> u;
                    if (WarpX::gamma_boost == 1.){
                      // Lab-frame simulation
                      // If the particle is not within the species's
                      // xmin, xmax, ymin, ymax, zmin, zmax, go to
                      // the next generated particle.
                      if (!plasma_injector->insideBounds(xb, yb, z)) continue;
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
                      if (!plasma_injector->insideBounds(xb, yb, z0_lab)) continue;
                      // call `getDensity` with lab-frame parameters
                      dens = plasma_injector->getDensity(x, y, z0_lab);
                      // At this point u and dens are the lab-frame quantities
                      // => Perform Lorentz transform
                      dens = gamma_boost * dens * ( 1 - beta_boost*betaz_lab );
                      u[2] = gamma_boost * ( u[2] -beta_boost*c*gamma_lab );
                    }
                    Real weight = dens * scale_fac / (AMREX_D_TERM(fac, *fac, *fac));
#ifdef WARPX_RZ
                    if (plasma_injector->radially_weighted) {
                      weight *= 2*MathConst::pi*xb;
                    } else {
                      // This is not correct since it might shift the particle
                      // out of the local grid
                      x = std::sqrt(xb*rmax);
                      weight *= dx[0];
                    }
#endif
                    attribs[PIdx::w ] = weight;
                    attribs[PIdx::ux] = u[0];
                    attribs[PIdx::uy] = u[1];
                    attribs[PIdx::uz] = u[2];
                    
                    if (WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags)
                    {
                        auto& particle_tile = DefineAndReturnParticleTile(lev, grid_id, tile_id);
                        particle_tile.push_back_real(particle_comps["xold"], x);
                        particle_tile.push_back_real(particle_comps["yold"], y);
                        particle_tile.push_back_real(particle_comps["zold"], z);

                        particle_tile.push_back_real(particle_comps["uxold"], u[0]);
                        particle_tile.push_back_real(particle_comps["uyold"], u[1]);
                        particle_tile.push_back_real(particle_comps["uzold"], u[2]);
                    }

		    AddOneParticle(lev, grid_id, tile_id, x, y, z, attribs);
                }
            }

            if (cost) {
	        wt = (amrex::second() - wt) / tile_box.d_numPts();
                Array4<Real> const& costarr = cost->array(mfi);
                amrex::ParallelFor(tile_box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    costarr(i,j,k) += wt;
                });
            }
        }
    }
}

#ifdef AMREX_USE_GPU
void
PhysicalParticleContainer::AddPlasmaGPU (int lev, RealBox part_realbox)
{
    BL_PROFILE("PhysicalParticleContainer::AddPlasmaGPU");

    // If no part_realbox is provided, initialize particles in the whole domain
    const Geometry& geom = Geom(lev);
    if (!part_realbox.ok()) part_realbox = geom.ProbDomain();

    int num_ppc = plasma_injector->num_particles_per_cell;
#ifdef WARPX_RZ
    Real rmax = std::min(plasma_injector->xmax, part_realbox.hi(0));
#endif

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

	    Cuda::HostVector<ParticleType> host_particles;
	    std::array<Cuda::HostVector<Real>, PIdx::nattribs> host_attribs;
	    
            // Loop through the cells of overlap_box and inject
            // the corresponding particles
            const auto& overlap_corner = overlap_realbox.lo();
            for (IntVect iv = overlap_box.smallEnd(); iv <= overlap_box.bigEnd(); overlap_box.next(iv))
            {
                int fac;
                if (do_continuous_injection) {
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

                    // Save the x and y values to use in the insideBounds checks.
                    // This is needed with WARPX_RZ since x and y are modified.
                    Real xb = x;
                    Real yb = y;

#ifdef WARPX_RZ
                    // Replace the x and y, choosing the angle randomly.
                    // These x and y are used to get the momentum and density
                    Real theta = 2.*MathConst::pi*amrex::Random();
                    x = xb*std::cos(theta);
                    y = xb*std::sin(theta);
#endif

                    Real dens;
                    std::array<Real, 3> u;
                    if (WarpX::gamma_boost == 1.){
                      // Lab-frame simulation
                      // If the particle is not within the species's
                      // xmin, xmax, ymin, ymax, zmin, zmax, go to
                      // the next generated particle.
                      if (!plasma_injector->insideBounds(xb, yb, z)) continue;
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
                      if (!plasma_injector->insideBounds(xb, yb, z0_lab)) continue;
                      // call `getDensity` with lab-frame parameters
                      dens = plasma_injector->getDensity(x, y, z0_lab);
                      // At this point u and dens are the lab-frame quantities
                      // => Perform Lorentz transform
                      dens = gamma_boost * dens * ( 1 - beta_boost*betaz_lab );
                      u[2] = gamma_boost * ( u[2] -beta_boost*c*gamma_lab );
                    }
                    Real weight = dens * scale_fac / (AMREX_D_TERM(fac, *fac, *fac));
#ifdef WARPX_RZ
                    if (plasma_injector->radially_weighted) {
                      weight *= 2*MathConst::pi*xb;
                    } else {
                      // This is not correct since it might shift the particle
                      // out of the local grid
                      x = std::sqrt(xb*rmax);
                      weight *= dx[0];
                    }
#endif
                    attribs[PIdx::w ] = weight;
                    attribs[PIdx::ux] = u[0];
                    attribs[PIdx::uy] = u[1];
                    attribs[PIdx::uz] = u[2];

                    // note - this will be slow on the GPU, need to revisit
                    if (WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags)
                    {
                        auto& particle_tile = DefineAndReturnParticleTile(lev, grid_id, tile_id);
                        particle_tile.push_back_real(particle_comps["xold"], x);
                        particle_tile.push_back_real(particle_comps["yold"], y);
                        particle_tile.push_back_real(particle_comps["zold"], z);

                        particle_tile.push_back_real(particle_comps["uxold"], u[0]);
                        particle_tile.push_back_real(particle_comps["uyold"], u[1]);
                        particle_tile.push_back_real(particle_comps["uzold"], u[2]);
                    }

		    ParticleType p;
		    p.id()  = ParticleType::NextID();
		    p.cpu() = ParallelDescriptor::MyProc();
#if (AMREX_SPACEDIM == 3)
		    p.pos(0) = x;
		    p.pos(1) = y;
		    p.pos(2) = z;
#elif (AMREX_SPACEDIM == 2)
#ifdef WARPX_RZ
                    attribs[PIdx::theta] = theta;
#endif
		    p.pos(0) = xb;
		    p.pos(1) = z;
#endif

		    host_particles.push_back(p);
		    for (int kk = 0; kk < PIdx::nattribs; ++kk)
		      host_attribs[kk].push_back(attribs[kk]);
                }
            }

	    auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
            auto old_size = particle_tile.GetArrayOfStructs().size();
            auto new_size = old_size + host_particles.size();
	    particle_tile.resize(new_size);

            Cuda::thrust_copy(host_particles.begin(),
                              host_particles.end(),
                              particle_tile.GetArrayOfStructs().begin() + old_size);

	    for (int kk = 0; kk < PIdx::nattribs; ++kk) {
                Cuda::thrust_copy(host_attribs[kk].begin(),
                                  host_attribs[kk].end(),
                                  particle_tile.GetStructOfArrays().GetRealData(kk).begin() + old_size);
	    }
	    			 
            if (cost) {
	        wt = (amrex::second() - wt) / tile_box.d_numPts();
                Array4<Real> const& costarr = cost->array(mfi);
                amrex::ParallelFor(tile_box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    costarr(i,j,k) += wt;
                });
            }
        }		
    }
}
#endif

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

            WRPX_INTERPOLATE_CIC(particles.dataPtr(), nstride, np,
                                 Exp.dataPtr(), Eyp.dataPtr(),
#if AMREX_SPACEDIM == 3
                                 Ezp.dataPtr(),
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
                WRPX_INTERPOLATE_CIC(particles.dataPtr(), nstride, np,
                                     Exp.dataPtr(), Eyp.dataPtr(),
#if AMREX_SPACEDIM == 3
                Ezp.dataPtr(),
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

                WRPX_INTERPOLATE_CIC_TWO_LEVELS(particles.dataPtr(), nstride, np,
                                                Exp.dataPtr(), Eyp.dataPtr(),
#if AMREX_SPACEDIM == 3
                                                Ezp.dataPtr(),
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
            WRPX_PUSH_LEAPFROG(particles.dataPtr(), nstride, np,
                               uxp.dataPtr(), uyp.dataPtr(),
#if AMREX_SPACEDIM == 3
                               uzp.dataPtr(),
#endif
                               Exp.dataPtr(), Eyp.dataPtr(),
#if AMREX_SPACEDIM == 3
                               Ezp.dataPtr(),
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

    BL_ASSERT(OnSameGrids(lev,Ex));

    MultiFab* cost = WarpX::getCosts(lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Cuda::ManagedDeviceVector<Real> xp, yp, zp;

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
            long lvect_fieldgathe = 64;
	    warpx_geteb_energy_conserving(
	       &np,
               xp.dataPtr(),
               yp.dataPtr(),
               zp.dataPtr(),
	       Exp.dataPtr(),Eyp.dataPtr(),Ezp.dataPtr(),
	       Bxp.dataPtr(),Byp.dataPtr(),Bzp.dataPtr(),
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
	       &ll4symtry, &WarpX::l_lower_order_in_v, &WarpX::do_nodal,
	       &lvect_fieldgathe, &WarpX::field_gathering_algo);

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
PhysicalParticleContainer::Evolve (int lev,
				   const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
				   const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
				   MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                   MultiFab* cjx, MultiFab* cjy, MultiFab* cjz,
                                   MultiFab* rho, MultiFab* crho,
                                   const MultiFab* cEx, const MultiFab* cEy, const MultiFab* cEz,
                                   const MultiFab* cBx, const MultiFab* cBy, const MultiFab* cBz,
                                   Real t, Real dt)
{
    BL_PROFILE("PPC::Evolve()");
    BL_PROFILE_VAR_NS("PPC::Evolve::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PICSAR::FieldGather", blp_pxr_fg);
    BL_PROFILE_VAR_NS("PICSAR::ParticlePush", blp_pxr_pp);
    BL_PROFILE_VAR_NS("PPC::Evolve::partition", blp_partition);
    
    const std::array<Real,3>& dx = WarpX::CellSize(lev);
    const std::array<Real,3>& cdx = WarpX::CellSize(std::max(lev-1,0));

    // Get instances of NCI Godfrey filters 
    const auto& nci_godfrey_filter_exeybz = WarpX::GetInstance().nci_godfrey_filter_exeybz;
    const auto& nci_godfrey_filter_bxbyez = WarpX::GetInstance().nci_godfrey_filter_bxbyez;

    BL_ASSERT(OnSameGrids(lev,jx));

    MultiFab* cost = WarpX::getCosts(lev);

    const iMultiFab* current_masks = WarpX::CurrentBufferMasks(lev);
    const iMultiFab* gather_masks = WarpX::GatherBufferMasks(lev);

    bool has_buffer = cEx || cjx;

#ifdef _OPENMP
#pragma omp parallel 
#endif
    {
#ifdef _OPENMP
        int thread_num = omp_get_thread_num();
#else
        int thread_num = 0;
#endif

        FArrayBox filtered_Ex, filtered_Ey, filtered_Ez;
        FArrayBox filtered_Bx, filtered_By, filtered_Bz;
        std::vector<bool> inexflag;
        Vector<long> pid;
        RealVector tmp;
        ParticleVector particle_tmp;

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

            Elixir exeli, eyeli, ezeli, bxeli, byeli, bzeli;
            if (WarpX::use_fdtd_nci_corr)
            {
#if (AMREX_SPACEDIM == 2)
                const Box& tbox = amrex::grow(pti.tilebox(),{static_cast<int>(WarpX::nox),
                            static_cast<int>(WarpX::noz)});
#else
                const Box& tbox = amrex::grow(pti.tilebox(),{static_cast<int>(WarpX::nox),
                            static_cast<int>(WarpX::noy),
                            static_cast<int>(WarpX::noz)});
#endif

                // Filter Ex (Both 2D and 3D)
                filtered_Ex.resize(amrex::convert(tbox,WarpX::Ex_nodal_flag));
                // Safeguard for GPU
                exeli = filtered_Ex.elixir();
                // Apply filter on Ex, result stored in filtered_Ex
                nci_godfrey_filter_exeybz[lev]->ApplyStencil(filtered_Ex, Ex[pti], filtered_Ex.box());
                // Update exfab reference
                exfab = &filtered_Ex;

                // Filter Ez
                filtered_Ez.resize(amrex::convert(tbox,WarpX::Ez_nodal_flag));
                ezeli = filtered_Ez.elixir();
                nci_godfrey_filter_bxbyez[lev]->ApplyStencil(filtered_Ez, Ez[pti], filtered_Ez.box());
                ezfab = &filtered_Ez;

                // Filter By
                filtered_By.resize(amrex::convert(tbox,WarpX::By_nodal_flag));
                byeli = filtered_By.elixir();
                nci_godfrey_filter_bxbyez[lev]->ApplyStencil(filtered_By, By[pti], filtered_By.box());
                byfab = &filtered_By;
#if (AMREX_SPACEDIM == 3)
                // Filter Ey
                filtered_Ey.resize(amrex::convert(tbox,WarpX::Ey_nodal_flag));
                eyeli = filtered_Ey.elixir();
                nci_godfrey_filter_exeybz[lev]->ApplyStencil(filtered_Ey, Ey[pti], filtered_Ey.box());
                eyfab = &filtered_Ey;

                // Filter Bx
                filtered_Bx.resize(amrex::convert(tbox,WarpX::Bx_nodal_flag));
                bxeli = filtered_Bx.elixir();
                nci_godfrey_filter_bxbyez[lev]->ApplyStencil(filtered_Bx, Bx[pti], filtered_Bx.box());
                bxfab = &filtered_Bx;

                // Filter Bz
                filtered_Bz.resize(amrex::convert(tbox,WarpX::Bz_nodal_flag));
                bzeli = filtered_Bz.elixir();
                nci_godfrey_filter_exeybz[lev]->ApplyStencil(filtered_Bz, Bz[pti], filtered_Bz.box());
                bzfab = &filtered_Bz;
#endif
            }

	    Exp.assign(np,0.0);
	    Eyp.assign(np,0.0);
	    Ezp.assign(np,0.0);
	    Bxp.assign(np,WarpX::B_external[0]);
	    Byp.assign(np,WarpX::B_external[1]);
	    Bzp.assign(np,WarpX::B_external[2]);

	    m_giv[thread_num].resize(np);

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

            const long np_current = (cjx) ? nfine_current : np;
            
	    //
	    // copy data from particle container to temp arrays
	    //
	    BL_PROFILE_VAR_START(blp_copy);
            pti.GetPosition(m_xp[thread_num], m_yp[thread_num], m_zp[thread_num]);
	    BL_PROFILE_VAR_STOP(blp_copy);

            if (rho) DepositCharge(pti, wp, rho, crho, 0, np_current, np, thread_num, lev);
            
            if (! do_not_push)
            {
                //
                // Field Gather of Aux Data (i.e., the full solution)
                //
                const int ll4symtry          = false;
                long lvect_fieldgathe = 64;

                const std::array<Real,3>& xyzmin_grid = WarpX::LowerCorner(box, lev);
                const int* ixyzmin_grid = box.loVect();
                
                const long np_gather = (cEx) ? nfine_gather : np;

                BL_PROFILE_VAR_START(blp_pxr_fg);

                warpx_geteb_energy_conserving(
                    &np_gather,
                    m_xp[thread_num].dataPtr(),
                    m_yp[thread_num].dataPtr(),
                    m_zp[thread_num].dataPtr(),
                    Exp.dataPtr(),Eyp.dataPtr(),Ezp.dataPtr(),
                    Bxp.dataPtr(),Byp.dataPtr(),Bzp.dataPtr(),
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
                    &ll4symtry, &WarpX::l_lower_order_in_v, &WarpX::do_nodal,
                    &lvect_fieldgathe, &WarpX::field_gathering_algo);

                if (np_gather < np)
                {
                    const IntVect& ref_ratio = WarpX::RefRatio(lev-1);
                    const Box& cbox = amrex::coarsen(box,ref_ratio);
                    const std::array<Real,3>& cxyzmin_grid = WarpX::LowerCorner(cbox, lev-1);
                    const int* cixyzmin_grid = cbox.loVect();

                    const FArrayBox* cexfab = &(*cEx)[pti];
                    const FArrayBox* ceyfab = &(*cEy)[pti];
                    const FArrayBox* cezfab = &(*cEz)[pti];
                    const FArrayBox* cbxfab = &(*cBx)[pti];
                    const FArrayBox* cbyfab = &(*cBy)[pti];
                    const FArrayBox* cbzfab = &(*cBz)[pti];

                    if (WarpX::use_fdtd_nci_corr)
                    {
#if (AMREX_SPACEDIM == 2)
                        const Box& tbox = amrex::grow(cbox,{static_cast<int>(WarpX::nox),
                                    static_cast<int>(WarpX::noz)});
#else
                        const Box& tbox = amrex::grow(cbox,{static_cast<int>(WarpX::nox),
                                    static_cast<int>(WarpX::noy),
                                    static_cast<int>(WarpX::noz)});
#endif

                        // Filter Ex (both 2D and 3D)
                        filtered_Ex.resize(amrex::convert(tbox,WarpX::Ex_nodal_flag));
                        // Safeguard for GPU
                        exeli = filtered_Ex.elixir();
                        // Apply filter on Ex, result stored in filtered_Ex
                        nci_godfrey_filter_exeybz[lev-1]->ApplyStencil(filtered_Ex, (*cEx)[pti], filtered_Ex.box());
                        // Update exfab reference
                        cexfab = &filtered_Ex;

                        // Filter Ez
                        filtered_Ez.resize(amrex::convert(tbox,WarpX::Ez_nodal_flag));
                        ezeli = filtered_Ez.elixir();
                        nci_godfrey_filter_bxbyez[lev-1]->ApplyStencil(filtered_Ez, (*cEz)[pti], filtered_Ez.box());
                        cezfab = &filtered_Ez;

                        // Filter By
                        filtered_By.resize(amrex::convert(tbox,WarpX::By_nodal_flag));
                        byeli = filtered_By.elixir();
                        nci_godfrey_filter_bxbyez[lev-1]->ApplyStencil(filtered_By, (*cBy)[pti], filtered_By.box());
                        cbyfab = &filtered_By;
#if (AMREX_SPACEDIM == 3)
                        // Filter Ey
                        filtered_Ey.resize(amrex::convert(tbox,WarpX::Ey_nodal_flag));
                        eyeli = filtered_Ey.elixir();
                        nci_godfrey_filter_exeybz[lev-1]->ApplyStencil(filtered_Ey, (*cEy)[pti], filtered_Ey.box());
                        ceyfab = &filtered_Ey;
                        
                        // Filter Bx
                        filtered_Bx.resize(amrex::convert(tbox,WarpX::Bx_nodal_flag));
                        bxeli = filtered_Bx.elixir();
                        nci_godfrey_filter_bxbyez[lev-1]->ApplyStencil(filtered_Bx, (*cBx)[pti], filtered_Bx.box());
                        cbxfab = &filtered_Bx;
                        
                        // Filter Bz
                        filtered_Bz.resize(amrex::convert(tbox,WarpX::Bz_nodal_flag));
                        bzeli = filtered_Bz.elixir();
                        nci_godfrey_filter_exeybz[lev-1]->ApplyStencil(filtered_Bz, (*cBz)[pti], filtered_Bz.box());
                        cbzfab = &filtered_Bz;
#endif
                    }
                    
                    long ncrse = np - nfine_gather;
                    warpx_geteb_energy_conserving(
                        &ncrse,
                        m_xp[thread_num].dataPtr()+nfine_gather,
                        m_yp[thread_num].dataPtr()+nfine_gather,
                        m_zp[thread_num].dataPtr()+nfine_gather,
                        Exp.dataPtr()+nfine_gather, Eyp.dataPtr()+nfine_gather, Ezp.dataPtr()+nfine_gather,
                        Bxp.dataPtr()+nfine_gather, Byp.dataPtr()+nfine_gather, Bzp.dataPtr()+nfine_gather,
                        cixyzmin_grid,
                        &cxyzmin_grid[0], &cxyzmin_grid[1], &cxyzmin_grid[2],
                        &cdx[0], &cdx[1], &cdx[2],
                        &WarpX::nox, &WarpX::noy, &WarpX::noz,
                        BL_TO_FORTRAN_ANYD(*cexfab),
                        BL_TO_FORTRAN_ANYD(*ceyfab),
                        BL_TO_FORTRAN_ANYD(*cezfab),
                        BL_TO_FORTRAN_ANYD(*cbxfab),
                        BL_TO_FORTRAN_ANYD(*cbyfab),
                        BL_TO_FORTRAN_ANYD(*cbzfab),
                        &ll4symtry, &WarpX::l_lower_order_in_v, &WarpX::do_nodal,
                        &lvect_fieldgathe, &WarpX::field_gathering_algo);
                }

                BL_PROFILE_VAR_STOP(blp_pxr_fg);

                //
                // Particle Push
                //
                BL_PROFILE_VAR_START(blp_pxr_pp);
                PushPX(pti, m_xp[thread_num], m_yp[thread_num], m_zp[thread_num], 
                       m_giv[thread_num], dt);
                BL_PROFILE_VAR_STOP(blp_pxr_pp);

                //
                // Current Deposition
                //
                if (WarpX::use_picsar_deposition) {
                    // Deposit inside domains
                    DepositCurrentFortran(pti, wp, uxp, uyp, uzp, &jx, &jy, &jz,
                                          0, np_current, thread_num,
                                          lev, lev, dt);
                    if (has_buffer){
                        // Deposit in buffers
                        DepositCurrentFortran(pti, wp, uxp, uyp, uzp, cjx, cjy, cjz,
                                              np_current, np-np_current, thread_num,
                                              lev, lev-1, dt);
                    }
                } else {
                    // Deposit inside domains
                    DepositCurrent(pti, wp, uxp, uyp, uzp, &jx, &jy, &jz,
                                   0, np_current, thread_num,
                                   lev, lev, dt);
                    if (has_buffer){
                        // Deposit in buffers
                        DepositCurrent(pti, wp, uxp, uyp, uzp, cjx, cjy, cjz,
                                       np_current, np-np_current, thread_num,
                                       lev, lev-1, dt);
                    }
                }

                //
                // copy particle data back
                //
                BL_PROFILE_VAR_START(blp_copy);
                pti.SetPosition(m_xp[thread_num], m_yp[thread_num], m_zp[thread_num]);
                BL_PROFILE_VAR_STOP(blp_copy);
            }
            
            if (rho) DepositCharge(pti, wp, rho, crho, 1, np_current, np, thread_num, lev);

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
    // Split particles
    if (do_splitting){ SplitParticles(lev); }
}

// Loop over all particles in the particle container and
// split particles tagged with p.id()=DoSplitParticleID
void
PhysicalParticleContainer::SplitParticles(int lev)
{
    auto& mypc = WarpX::GetInstance().GetPartContainer();
    auto& pctmp_split = mypc.GetPCtmp();
    Cuda::ManagedDeviceVector<Real> xp, yp, zp;
    RealVector psplit_x, psplit_y, psplit_z, psplit_w;
    RealVector psplit_ux, psplit_uy, psplit_uz;
    long np_split_to_add = 0;
    long np_split;
    if(split_type==0)
    {
        np_split = pow(2, AMREX_SPACEDIM);
    } else {
        np_split = 2*AMREX_SPACEDIM;
    }

    // Loop over particle interator
    for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        pti.GetPosition(xp, yp, zp);
        const std::array<Real,3>& dx = WarpX::CellSize(lev);
	// particle Array Of Structs data
        auto& particles = pti.GetArrayOfStructs();
	// particle Struct Of Arrays data
        auto& attribs = pti.GetAttribs();
        auto& wp  = attribs[PIdx::w ];
        auto& uxp = attribs[PIdx::ux];
        auto& uyp = attribs[PIdx::uy];
        auto& uzp = attribs[PIdx::uz];
        const long np = pti.numParticles();
        for(int i=0; i<np; i++){
            auto& p = particles[i];
            if (p.id() == DoSplitParticleID){
		// If particle is tagged, split it and put the 
		// split particles in local arrays psplit_x etc.
                np_split_to_add += np_split;
#if (AMREX_SPACEDIM==2)
                if (split_type==0){
		    // Split particle in two along each axis
		    // 4 particles in 2d
                    for (int ishift = -1; ishift < 2; ishift +=2 ){
                        for (int kshift = -1; kshift < 2; kshift +=2 ){
                            // Add one particle with offset in x and z
                            psplit_x.push_back( xp[i] + ishift*dx[0]/2 );
                            psplit_y.push_back( yp[i] );
                            psplit_z.push_back( zp[i] + kshift*dx[2]/2 );
                            psplit_ux.push_back( uxp[i] );
                            psplit_uy.push_back( uyp[i] );
                            psplit_uz.push_back( uzp[i] );
                            psplit_w.push_back( wp[i]/np_split );
                        }
                    }
                } else {
		    // Split particle in two along each diagonal
		    // 4 particles in 2d
                    for (int ishift = -1; ishift < 2; ishift +=2 ){
                        // Add one particle with offset in x
                        psplit_x.push_back( xp[i] + ishift*dx[0]/2 );
                        psplit_y.push_back( yp[i] );
                        psplit_z.push_back( zp[i] );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                        // Add one particle with offset in z
                        psplit_x.push_back( xp[i] );
                        psplit_y.push_back( yp[i] );
                        psplit_z.push_back( zp[i] + ishift*dx[2]/2 );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                    }
                }
#elif (AMREX_SPACEDIM==3)
		if (split_type==0){
		    // Split particle in two along each axis
		    // 6 particles in 2d
		    for (int ishift = -1; ishift < 2; ishift +=2 ){
			for (int jshift = -1; jshift < 2; jshift +=2 ){
			    for (int kshift = -1; kshift < 2; kshift +=2 ){
				// Add one particle with offset in x, y and z
				psplit_x.push_back( xp[i] + ishift*dx[0]/2 );
				psplit_y.push_back( yp[i] + jshift*dx[1]/2 );
				psplit_z.push_back( zp[i] + jshift*dx[2]/2 );
				psplit_ux.push_back( uxp[i] );
				psplit_uy.push_back( uyp[i] );
				psplit_uz.push_back( uzp[i] );
				psplit_w.push_back( wp[i]/np_split );
			    }
			}
		    }
		} else {
		    // Split particle in two along each diagonal
		    // 8 particles in 3d
                    for (int ishift = -1; ishift < 2; ishift +=2 ){
                        // Add one particle with offset in x
                        psplit_x.push_back( xp[i] + ishift*dx[0]/2 );
                        psplit_y.push_back( yp[i] );
                        psplit_z.push_back( zp[i] );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                        // Add one particle with offset in y
                        psplit_x.push_back( xp[i] );
                        psplit_y.push_back( yp[i] + ishift*dx[1]/2 );
                        psplit_z.push_back( zp[i] );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                        // Add one particle with offset in z
                        psplit_x.push_back( xp[i] );
                        psplit_y.push_back( yp[i] );
                        psplit_z.push_back( zp[i] + ishift*dx[2]/2 );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                    }
		}
#endif
		// invalidate the particle
                p.m_idata.id = -p.m_idata.id;
            }
        }
    }
	// Add local arrays psplit_x etc. to the temporary
	// particle container pctmp_split. Split particles
	// are tagged with p.id()=NoSplitParticleID so that 
	// they are not re-split when entering a higher level
	// AddNParticles calls Redistribute, so that particles
	// in pctmp_split are in the proper grids and tiles
	pctmp_split.AddNParticles(lev, 
                                  np_split_to_add,
                                  psplit_x.dataPtr(),
                                  psplit_y.dataPtr(),
                                  psplit_z.dataPtr(),
                                  psplit_ux.dataPtr(),
                                  psplit_uy.dataPtr(),
                                  psplit_uz.dataPtr(),
                                  1,
                                  psplit_w.dataPtr(),
                                  1, NoSplitParticleID);
	// Copy particles from tmp to current particle container
    addParticles(pctmp_split,1);
	// Clear tmp container
    pctmp_split.clearParticles();
}

void
PhysicalParticleContainer::PushPX(WarpXParIter& pti,
	                          Cuda::ManagedDeviceVector<Real>& xp,
                                  Cuda::ManagedDeviceVector<Real>& yp,
                                  Cuda::ManagedDeviceVector<Real>& zp,
                                  Cuda::ManagedDeviceVector<Real>& giv,
                                  Real dt)
{

    if (WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags)
    {
        copy_attribs(pti, xp.dataPtr(), yp.dataPtr(), zp.dataPtr());
    }

    // The following attributes should be included in CPP version of warpx_particle_pusher
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
    
    warpx_particle_pusher(&np,
                          xp.dataPtr(),
                          yp.dataPtr(),
                          zp.dataPtr(),
                          uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(),
                          giv.dataPtr(),
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
    BL_PROFILE("PhysicalParticleContainer::PushP");

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

            warpx_particle_pusher_momenta(&np,
                                          m_xp[thread_num].dataPtr(),
                                          m_yp[thread_num].dataPtr(),
                                          m_zp[thread_num].dataPtr(),
                                          uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(),
                                          m_giv[thread_num].dataPtr(),
                                          Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
                                          Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
                                          &this->charge, &this->mass, &dt,
                                          &WarpX::particle_pusher_algo);
        }
    }
}

void PhysicalParticleContainer::copy_attribs(WarpXParIter& pti,const Real* xp,
                        const Real* yp, const Real* zp)
{

    auto& attribs = pti.GetAttribs();
    
    Real* AMREX_RESTRICT uxp = attribs[PIdx::ux].dataPtr();
    Real* AMREX_RESTRICT uyp = attribs[PIdx::uy].dataPtr();
    Real* AMREX_RESTRICT uzp = attribs[PIdx::uz].dataPtr();
    
    Real* AMREX_RESTRICT xpold = pti.GetAttribs(particle_comps["xold"]).dataPtr();
    Real* AMREX_RESTRICT ypold = pti.GetAttribs(particle_comps["yold"]).dataPtr();
    Real* AMREX_RESTRICT zpold = pti.GetAttribs(particle_comps["zold"]).dataPtr();
    Real* AMREX_RESTRICT uxpold = pti.GetAttribs(particle_comps["uxold"]).dataPtr();
    Real* AMREX_RESTRICT uypold = pti.GetAttribs(particle_comps["uyold"]).dataPtr();
    Real* AMREX_RESTRICT uzpold = pti.GetAttribs(particle_comps["uzold"]).dataPtr();
    
    const long np = pti.numParticles();
    
    ParallelFor( np,
        [=] AMREX_GPU_DEVICE (long i) {
            xpold[i]=xp[i];
            ypold[i]=yp[i];
            zpold[i]=zp[i];
            
            uxpold[i]=uxp[i];
            uypold[i]=uyp[i];
            uzpold[i]=uzp[i];
        }
    );
}

void PhysicalParticleContainer::GetParticleSlice(const int direction, const Real z_old,
                                                 const Real z_new, const Real t_boost,
                                                 const Real t_lab, const Real dt,
                                                 DiagnosticParticles& diagnostic_particles)
{
    BL_PROFILE("PhysicalParticleContainer::GetParticleSlice");

    // Assume that the boost in the positive z direction.
#if (AMREX_SPACEDIM == 2)
    AMREX_ALWAYS_ASSERT(direction == 1);
#else
    AMREX_ALWAYS_ASSERT(direction == 2);
#endif

    // Note the the slice should always move in the negative boost direction.
    AMREX_ALWAYS_ASSERT(z_new < z_old);

    AMREX_ALWAYS_ASSERT(do_boosted_frame_diags == 1);

    const int nlevs = std::max(0, finestLevel()+1);

    // we figure out a box for coarse-grained rejection. If the RealBox corresponding to a
    // given tile doesn't intersect with this, there is no need to check any particles.
    const Real* base_dx = Geom(0).CellSize();
    const Real z_min = z_new - base_dx[direction];
    const Real z_max = z_old + base_dx[direction];

    RealBox slice_box = Geom(0).ProbDomain();
    slice_box.setLo(direction, z_min);
    slice_box.setHi(direction, z_max);

    diagnostic_particles.resize(finestLevel()+1);
    
    for (int lev = 0; lev < nlevs; ++lev) {

        const Real* dx  = Geom(lev).CellSize();
        const Real* plo = Geom(lev).ProbLo();

        // first we touch each map entry in serial
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
            diagnostic_particles[lev][index];
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            RealVector xp_new, yp_new, zp_new;

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

                auto&  xp_old = pti.GetAttribs(particle_comps["xold"]);
                auto&  yp_old = pti.GetAttribs(particle_comps["yold"]);
                auto&  zp_old = pti.GetAttribs(particle_comps["zold"]);
                auto& uxp_old = pti.GetAttribs(particle_comps["uxold"]);
                auto& uyp_old = pti.GetAttribs(particle_comps["uyold"]);
                auto& uzp_old = pti.GetAttribs(particle_comps["uzold"]);

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

                    diagnostic_particles[lev][index].GetRealData(DiagIdx::w).push_back(wp[i]);
                    
                    diagnostic_particles[lev][index].GetRealData(DiagIdx::x).push_back(xp);
                    diagnostic_particles[lev][index].GetRealData(DiagIdx::y).push_back(yp);
                    diagnostic_particles[lev][index].GetRealData(DiagIdx::z).push_back(zp);

                    diagnostic_particles[lev][index].GetRealData(DiagIdx::ux).push_back(uxp);
                    diagnostic_particles[lev][index].GetRealData(DiagIdx::uy).push_back(uyp);
                    diagnostic_particles[lev][index].GetRealData(DiagIdx::uz).push_back(uzp);
                }
            }
        }
    }
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

        BoxArray stretched_ba(stretched_boxes.dataPtr(), stretched_boxes.size());

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

/* \brief Inject particles during the simulation
 * \param injection_box: domain where particles should be injected.
 */
void
PhysicalParticleContainer::ContinuousInjection(const RealBox& injection_box)
{
    // Inject plasma on level 0. Paticles will be redistributed.
    const int lev=0;
    AddPlasma(lev, injection_box);
}
