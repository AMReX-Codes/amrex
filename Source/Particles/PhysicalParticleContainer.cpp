#include <limits>
#include <sstream>

#include <MultiParticleContainer.H>
#include <WarpX_f.H>
#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpXWrappers.h>
#include <IonizationEnergiesTable.H>
#include <FieldGather.H>

#include <WarpXAlgorithmSelection.H>

// Import low-level single-particle kernels
#include <UpdatePosition.H>
#include <UpdateMomentumBoris.H>
#include <UpdateMomentumVay.H>

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
    pp.query("do_splitting", do_splitting);
    pp.query("split_type", split_type);
    pp.query("do_continuous_injection", do_continuous_injection);
    // Whether to plot back-transformed (lab-frame) diagnostics 
    // for this species.
    pp.query("do_boosted_frame_diags", do_boosted_frame_diags);

    pp.query("do_field_ionization", do_field_ionization);
#ifdef AMREX_USE_GPU
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        do_field_ionization == 0,
        "Field ionization does not work on GPU so far, because the current "
        "version of Redistribute in AMReX does not work with runtime parameters");
#endif

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
    // Init ionization module here instead of in the PhysicalParticleContainer
    // constructor because dt is required
    if (do_field_ionization) {InitIonizationModule();}
    AddParticles(0); // Note - add on level 0
    Redistribute();  // We then redistribute
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
        // If do_symmetrize, create 4x fewer particles, and 
        // Replicate each particle 4 times (x,y) (-x,y) (x,-y) (-x,-y)
        if (do_symmetrize){
            npart /= 4;
        }
        for (long i = 0; i < npart; ++i) {
#if (defined WARPX_DIM_3D) || (WARPX_DIM_RZ)
            Real weight = q_tot/npart/charge;
            Real x = distx(mt);
            Real y = disty(mt);
            Real z = distz(mt);
#elif (defined WARPX_DIM_XZ)
            Real weight = q_tot/npart/charge/y_rms;
            Real x = distx(mt);
            Real y = 0.;
            Real z = distz(mt);
#endif
            if (plasma_injector->insideBounds(x, y, z)) {
                XDim3 u = plasma_injector->getMomentum(x, y, z);
                u.x *= PhysConst::c;
                u.y *= PhysConst::c;
                u.z *= PhysConst::c;
                if (do_symmetrize){
                    // Add four particles to the beam:
                    CheckAndAddParticle( x, y, z, { u.x, u.y, u.z}, weight/4. );
                    CheckAndAddParticle( x,-y, z, { u.x,-u.y, u.z}, weight/4. );
                    CheckAndAddParticle(-x, y, z, {-u.x, u.y, u.z}, weight/4. );
                    CheckAndAddParticle(-x,-y, z, {-u.x,-u.y, u.z}, weight/4. );
                } else {
                    CheckAndAddParticle(x, y, z, {u.x,u.y,u.z}, weight);
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

    if ( (NumRuntimeRealComps()>0) || (NumRuntimeIntComps()>0) )
    {
        // need to create old values
        auto& particle_tile = DefineAndReturnParticleTile(0, 0, 0);
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
    BL_PROFILE("PhysicalParticleContainer::AddPlasma");

    // If no part_realbox is provided, initialize particles in the whole domain
    const Geometry& geom = Geom(lev);
    if (!part_realbox.ok()) part_realbox = geom.ProbDomain();

    int num_ppc = plasma_injector->num_particles_per_cell;
#ifdef WARPX_DIM_RZ
    Real rmax = std::min(plasma_injector->xmax, part_realbox.hi(0));
#endif

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    Real scale_fac;
#if AMREX_SPACEDIM==3
    scale_fac = dx[0]*dx[1]*dx[2]/num_ppc;
#elif AMREX_SPACEDIM==2
    scale_fac = dx[0]*dx[1]/num_ppc;
#endif

#ifdef _OPENMP
    // First touch all tiles in the map in serial
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
        GetParticles(lev)[index];
        tmp_particle_data.resize(finestLevel()+1);
        // Create map entry if not there
        tmp_particle_data[lev][index];
        if ( (NumRuntimeRealComps()>0) || (NumRuntimeIntComps()>0) ) {
            DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
        }
    }
#endif

    MultiFab* cost = WarpX::getCosts(lev);

    const int nlevs = numLevels();
    static bool refine_injection = false;
    static Box fine_injection_box;
    static int rrfac = 1;
    // This does not work if the mesh is dynamic.  But in that case, we should
    // not use refined injected either.  We also assume there is only one fine level.
    if (WarpX::do_moving_window and WarpX::refine_plasma
        and do_continuous_injection and nlevs == 2)
    {
        refine_injection = true;
        fine_injection_box = ParticleBoxArray(1).minimalBox();
        fine_injection_box.setSmall(WarpX::moving_window_dir, std::numeric_limits<int>::lowest());
        fine_injection_box.setBig(WarpX::moving_window_dir, std::numeric_limits<int>::max());
        rrfac = m_gdb->refRatio(0)[0];
        fine_injection_box.coarsen(rrfac);
    }

    InjectorPosition* inj_pos = plasma_injector->getInjectorPosition();
    InjectorDensity*  inj_rho = plasma_injector->getInjectorDensity();
    InjectorMomentum* inj_mom = plasma_injector->getInjectorMomentum();
    Real gamma_boost = WarpX::gamma_boost;
    Real beta_boost = WarpX::beta_boost;
    Real t = WarpX::GetInstance().gett_new(lev);
    Real density_min = plasma_injector->density_min;
    Real density_max = plasma_injector->density_max;

#ifdef WARPX_DIM_RZ
    const long nmodes = WarpX::n_rz_azimuthal_modes;
    bool radially_weighted = plasma_injector->radially_weighted;
#endif

    MFItInfo info;
    if (do_tiling && Gpu::notInLaunchRegion()) {
        info.EnableTiling(tile_size);
    }
#ifdef _OPENMP
    info.SetDynamic(true);
#pragma omp parallel if (not WarpX::serialize_ics)
#endif
    for (MFIter mfi = MakeMFIter(lev, info); mfi.isValid(); ++mfi)
    {
        Real wt = amrex::second();

        const Box& tile_box = mfi.tilebox();
        const RealBox tile_realbox = WarpX::getRealBox(tile_box, lev);

        // Find the cells of part_box that overlap with tile_realbox
        // If there is no overlap, just go to the next tile in the loop
        RealBox overlap_realbox;
        Box overlap_box;
        IntVect shifted;
        bool no_overlap = false;

        for (int dir=0; dir<AMREX_SPACEDIM; dir++) {
            if ( tile_realbox.lo(dir) <= part_realbox.hi(dir) ) {
                Real ncells_adjust = std::floor( (tile_realbox.lo(dir) - part_realbox.lo(dir))/dx[dir] );
                overlap_realbox.setLo( dir, part_realbox.lo(dir) + std::max(ncells_adjust, 0.) * dx[dir]);
            } else {
                no_overlap = true; break;
            }
            if ( tile_realbox.hi(dir) >= part_realbox.lo(dir) ) {
                Real ncells_adjust = std::floor( (part_realbox.hi(dir) - tile_realbox.hi(dir))/dx[dir] );
                overlap_realbox.setHi( dir, part_realbox.hi(dir) - std::max(ncells_adjust, 0.) * dx[dir]);
            } else {
                no_overlap = true; break;
            }
            // Count the number of cells in this direction in overlap_realbox
            overlap_box.setSmall( dir, 0 );
            overlap_box.setBig( dir,
                int( std::round((overlap_realbox.hi(dir)-overlap_realbox.lo(dir))
                                /dx[dir] )) - 1);
            shifted[dir] = std::round((overlap_realbox.lo(dir)-problo[dir])/dx[dir]);
            // shifted is exact in non-moving-window direction.  That's all we care.
        }
        if (no_overlap == 1) {
            continue; // Go to the next tile
        }

        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();

        // Max number of new particles, if particles are created in the whole
        // overlap_box. All of them are created, and invalid ones are then 
        // discaded
        int max_new_particles = overlap_box.numPts() * num_ppc;

        // If refine injection, build pointer dp_cellid that holds pointer to 
        // array of refined cell IDs.
        Vector<int> cellid_v;
        if (refine_injection and lev == 0)
        {
            // then how many new particles will be injected is not that simple
            // We have to shift fine_injection_box because overlap_box has been shifted.
            Box fine_overlap_box = overlap_box & amrex::shift(fine_injection_box,shifted);
            if (fine_overlap_box.ok()) {
                max_new_particles += fine_overlap_box.numPts() * num_ppc
                    * (AMREX_D_TERM(rrfac,*rrfac,*rrfac)-1);
                for (int icell = 0, ncells = overlap_box.numPts(); icell < ncells; ++icell) {
                    IntVect iv = overlap_box.atOffset(icell);
                    int r = (fine_overlap_box.contains(iv)) ? AMREX_D_TERM(rrfac,*rrfac,*rrfac) : 1;
                    for (int ipart = 0; ipart < r; ++ipart) {
                        cellid_v.push_back(icell);
                        cellid_v.push_back(ipart);
                    }
                }
            }
        }
        int const* hp_cellid = (cellid_v.empty()) ? nullptr : cellid_v.data();
        amrex::AsyncArray<int> cellid_aa(hp_cellid, cellid_v.size());
        int const* dp_cellid = cellid_aa.data();

        // Update NextID to include particles created in this function
        int pid;
#pragma omp critical (add_plasma_nextid)
        {
            pid = ParticleType::NextID();
            ParticleType::NextID(pid+max_new_particles);
        }
        const int cpuid = ParallelDescriptor::MyProc();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        if ( (NumRuntimeRealComps()>0) || (NumRuntimeIntComps()>0) ) {
            DefineAndReturnParticleTile(lev, grid_id, tile_id);
        }

        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + max_new_particles;
        particle_tile.resize(new_size);

        ParticleType* pp = particle_tile.GetArrayOfStructs()().data() + old_size;
        auto& soa = particle_tile.GetStructOfArrays();
        GpuArray<Real*,PIdx::nattribs> pa;
        for (int ia = 0; ia < PIdx::nattribs; ++ia) {
            pa[ia] = soa.GetRealData(ia).data() + old_size;
        }

        int* pi;
        if (do_field_ionization) {
            pi = soa.GetIntData(particle_icomps["ionization_level"]).data() + old_size;
        }
        
        const GpuArray<Real,AMREX_SPACEDIM> overlap_corner
            {AMREX_D_DECL(overlap_realbox.lo(0),
                          overlap_realbox.lo(1),
                          overlap_realbox.lo(2))};

        std::size_t shared_mem_bytes = plasma_injector->sharedMemoryNeeded();
        int lrrfac = rrfac;

	bool loc_do_field_ionization = do_field_ionization;
	int loc_ionization_initial_level = ionization_initial_level;

        // Loop over all new particles and inject them (creates too many 
        // particles, in particular does not consider xmin, xmax etc.).
        // The invalid ones are given negative ID and are deleted during the 
        // next redistribute.
        amrex::For(max_new_particles, [=] AMREX_GPU_DEVICE (int ip) noexcept
        {
            ParticleType& p = pp[ip];
            p.id() = pid+ip;
            p.cpu() = cpuid;

            int cellid, i_part;
            Real fac;
            if (dp_cellid == nullptr) {
                cellid = ip/num_ppc;
                i_part = ip - cellid*num_ppc;
                fac = 1.0;
            } else {
                cellid = dp_cellid[2*ip];
                i_part = dp_cellid[2*ip+1];
                fac = lrrfac;
            }

            IntVect iv = overlap_box.atOffset(cellid);

            const XDim3 r = inj_pos->getPositionUnitBox(i_part, fac);
#if (AMREX_SPACEDIM == 3)
            Real x = overlap_corner[0] + (iv[0]+r.x)*dx[0];
            Real y = overlap_corner[1] + (iv[1]+r.y)*dx[1];
            Real z = overlap_corner[2] + (iv[2]+r.z)*dx[2];
#else
            Real x = overlap_corner[0] + (iv[0]+r.x)*dx[0];
            Real y = 0.0;
#if   defined WARPX_DIM_XZ
            Real z = overlap_corner[1] + (iv[1]+r.y)*dx[1];
#elif defined WARPX_DIM_RZ
            // Note that for RZ, r.y will be theta
            Real z = overlap_corner[1] + (iv[1]+r.z)*dx[1];
#endif
#endif

#if (AMREX_SPACEDIM == 3)
            if (!tile_realbox.contains(XDim3{x,y,z})) {
                p.id() = -1;
                return;
            }
#else
            if (!tile_realbox.contains(XDim3{x,z,0.0})) {
                p.id() = -1;
                return;
            }
#endif

            // Save the x and y values to use in the insideBounds checks.
            // This is needed with WARPX_DIM_RZ since x and y are modified.
            Real xb = x;
            Real yb = y;

#ifdef WARPX_DIM_RZ
            // Replace the x and y, setting an angle theta.
            // These x and y are used to get the momentum and density
            Real theta;
            if (nmodes == 1) {
                // With only 1 mode, the angle doesn't matter so
                // choose it randomly.
                theta = 2.*MathConst::pi*amrex::Random();
            } else {
                theta = 2.*MathConst::pi*r.y;
            }
            x = xb*std::cos(theta);
            y = xb*std::sin(theta);
#endif

            Real dens;
            XDim3 u;
            if (gamma_boost == 1.) {
                // Lab-frame simulation
                // If the particle is not within the species's
                // xmin, xmax, ymin, ymax, zmin, zmax, go to
                // the next generated particle.
                if (!inj_pos->insideBounds(xb, yb, z)) {
                    p.id() = -1;
                    return;
                }
                u = inj_mom->getMomentum(x, y, z);
                dens = inj_rho->getDensity(x, y, z);
                // Remove particle if density below threshold
                if ( dens < density_min ){
                    p.id() = -1;
                    return;
                }
                // Cut density if above threshold
                dens = amrex::min(dens, density_max);
            } else {
                // Boosted-frame simulation
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
                u = inj_mom->getMomentum(x, y, 0.); // No z0_lab dependency
                // At this point u is the lab-frame momentum
                // => Apply the above formula for z0_lab
                Real gamma_lab = std::sqrt( 1.+(u.x*u.x+u.y*u.y+u.z*u.z) );
                Real betaz_lab = u.z/(gamma_lab);
                Real z0_lab = gamma_boost * ( z*(1-beta_boost*betaz_lab)
                                              - PhysConst::c*t*(betaz_lab-beta_boost) );
                // If the particle is not within the lab-frame zmin, zmax, etc.
                // go to the next generated particle.
                if (!inj_pos->insideBounds(xb, yb, z0_lab)) {
                    p.id() = -1;
                    return;
                }
                // call `getDensity` with lab-frame parameters
                dens = inj_rho->getDensity(x, y, z0_lab);
                // Remove particle if density below threshold
                if ( dens < density_min ){
                    p.id() = -1;
                    return;
                }
                // Cut density if above threshold
                dens = amrex::min(dens, density_max);
                // At this point u and dens are the lab-frame quantities
                // => Perform Lorentz transform
                dens = gamma_boost * dens * ( 1.0 - beta_boost*betaz_lab );
                u.z = gamma_boost * ( u.z -beta_boost*gamma_lab );
            }

            if (loc_do_field_ionization) {
                pi[ip] = loc_ionization_initial_level;
            }

            u.x *= PhysConst::c;
            u.y *= PhysConst::c;
            u.z *= PhysConst::c;

            // Real weight = dens * scale_fac / (AMREX_D_TERM(fac, *fac, *fac));
            Real weight = dens * scale_fac;
#ifdef WARPX_DIM_RZ
            if (radially_weighted) {
                weight *= 2.*MathConst::pi*xb;
            } else {
                // This is not correct since it might shift the particle
                // out of the local grid
                x = std::sqrt(xb*rmax);
                weight *= dx[0];
            }
#endif
            pa[PIdx::w ][ip] = weight;
            pa[PIdx::ux][ip] = u.x;
            pa[PIdx::uy][ip] = u.y;
            pa[PIdx::uz][ip] = u.z;

#if (AMREX_SPACEDIM == 3)
            p.pos(0) = x;
            p.pos(1) = y;
            p.pos(2) = z;
#elif (AMREX_SPACEDIM == 2)
#ifdef WARPX_DIM_RZ
            pa[PIdx::theta][ip] = theta;
#endif
            p.pos(0) = xb;
            p.pos(1) = z;
#endif
        }, shared_mem_bytes);
    			 
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

    // The function that calls this is responsible for redistributing particles.
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
#ifdef _OPENMP
        int thread_num = omp_get_thread_num();
#else
        int thread_num = 0;
#endif
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
            pti.GetPosition(m_xp[thread_num], m_yp[thread_num], m_zp[thread_num]);

            //
            // Field Gather
            //
            int e_is_nodal = Ex.is_nodal() and Ey.is_nodal() and Ez.is_nodal();
            FieldGather(pti, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                        &exfab, &eyfab, &ezfab, &bxfab, &byfab, &bzfab, 
                        Ex.nGrow(), e_is_nodal,
                        0, np, thread_num, lev, lev);

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
    BL_PROFILE_VAR_NS("PPC::FieldGather", blp_fg);
    BL_PROFILE_VAR_NS("PPC::ParticlePush", blp_ppc_pp);
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

    if (WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags)
    {
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const auto np = pti.numParticles();
            const auto lev = pti.GetLevel();
            const auto index = pti.GetPairIndex();
            tmp_particle_data.resize(finestLevel()+1);
            for (int i = 0; i < TmpIdx::nattribs; ++i)
                tmp_particle_data[lev][index][i].resize(np);
        }
    }
    
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

            long nfine_current = np; //! number of particles depositing to fine grid
            long nfine_gather = np;  //! number of particles gathering from fine grid
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

                // only deposit / gather to coarsest grid
                if (m_deposit_on_main_grid && lev > 0) {
                    nfine_current = 0;
                }
                if (m_gather_from_main_grid && lev > 0) {
                    nfine_gather = 0;
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

            if (rho) {
                // Deposit charge before particle push, in component 0 of MultiFab rho.
                int* AMREX_RESTRICT ion_lev;
                if (do_field_ionization){
                    ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
                } else {
                    ion_lev = nullptr;
                }
                DepositCharge(pti, wp, ion_lev, rho, 0, 0, 
                              np_current, thread_num, lev, lev);
                if (has_buffer){
                    DepositCharge(pti, wp, ion_lev, crho, 0, np_current,
                                  np-np_current, thread_num, lev, lev-1);
                }
            }
            
            if (! do_not_push)
            {
                const long np_gather = (cEx) ? nfine_gather : np;

                int e_is_nodal = Ex.is_nodal() and Ey.is_nodal() and Ez.is_nodal();

                //
                // Field Gather of Aux Data (i.e., the full solution)
                //
                BL_PROFILE_VAR_START(blp_fg);
                FieldGather(pti, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                            exfab, eyfab, ezfab, bxfab, byfab, bzfab, 
                            Ex.nGrow(), e_is_nodal,
                            0, np_gather, thread_num, lev, lev);

                if (np_gather < np)
                {
                    const IntVect& ref_ratio = WarpX::RefRatio(lev-1);
                    const Box& cbox = amrex::coarsen(box,ref_ratio);

                    // Data on the grid
                    FArrayBox const* cexfab = &(*cEx)[pti];
                    FArrayBox const* ceyfab = &(*cEy)[pti];
                    FArrayBox const* cezfab = &(*cEz)[pti];
                    FArrayBox const* cbxfab = &(*cBx)[pti];
                    FArrayBox const* cbyfab = &(*cBy)[pti];
                    FArrayBox const* cbzfab = &(*cBz)[pti];
                    
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
                    
                    // Field gather for particles in gather buffers
                    e_is_nodal = cEx->is_nodal() and cEy->is_nodal() and cEz->is_nodal();
                    FieldGather(pti, Exp, Eyp, Ezp, Bxp, Byp, Bzp, 
                                cexfab, ceyfab, cezfab,
                                cbxfab, cbyfab, cbzfab,
                                cEx->nGrow(), e_is_nodal, 
                                nfine_gather, np-nfine_gather, 
                                thread_num, lev, lev-1);
                }

                BL_PROFILE_VAR_STOP(blp_fg);

                //
                // Particle Push
                //
                BL_PROFILE_VAR_START(blp_ppc_pp);
                PushPX(pti, m_xp[thread_num], m_yp[thread_num], m_zp[thread_num], 
                       m_giv[thread_num], dt);
                BL_PROFILE_VAR_STOP(blp_ppc_pp);

                //
                // Current Deposition
                //

                int* AMREX_RESTRICT ion_lev;
                if (do_field_ionization){
                    ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
                } else {
                    ion_lev = nullptr;
                }
                
                // Deposit inside domains
                DepositCurrent(pti, wp, uxp, uyp, uzp, ion_lev, &jx, &jy, &jz,
                               0, np_current, thread_num,
                               lev, lev, dt);
                if (has_buffer){
                    // Deposit in buffers
                    DepositCurrent(pti, wp, uxp, uyp, uzp, ion_lev, cjx, cjy, cjz,
                                   np_current, np-np_current, thread_num,
                                   lev, lev-1, dt);
                }


                //
                // copy particle data back
                //
                BL_PROFILE_VAR_START(blp_copy);
                pti.SetPosition(m_xp[thread_num], m_yp[thread_num], m_zp[thread_num]);
                BL_PROFILE_VAR_STOP(blp_copy);
            }
            
            if (rho) {
                // Deposit charge after particle push, in component 1 of MultiFab rho.
                int* AMREX_RESTRICT ion_lev;
                if (do_field_ionization){
                    ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
                } else {
                    ion_lev = nullptr;
                }
                DepositCharge(pti, wp, ion_lev, rho, 1, 0,
                              np_current, thread_num, lev, lev);
                if (has_buffer){
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

    // This wraps the momentum and position advance so that inheritors can modify the call.
    auto& attribs = pti.GetAttribs();
    // Extract pointers to the different particle quantities
    Real* const AMREX_RESTRICT x = xp.dataPtr();
    Real* const AMREX_RESTRICT y = yp.dataPtr();
    Real* const AMREX_RESTRICT z = zp.dataPtr();
    Real* const AMREX_RESTRICT gi = giv.dataPtr();
    Real* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    Real* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    Real* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
    const Real* const AMREX_RESTRICT Ex = attribs[PIdx::Ex].dataPtr();
    const Real* const AMREX_RESTRICT Ey = attribs[PIdx::Ey].dataPtr();
    const Real* const AMREX_RESTRICT Ez = attribs[PIdx::Ez].dataPtr();
    const Real* const AMREX_RESTRICT Bx = attribs[PIdx::Bx].dataPtr();
    const Real* const AMREX_RESTRICT By = attribs[PIdx::By].dataPtr();
    const Real* const AMREX_RESTRICT Bz = attribs[PIdx::Bz].dataPtr();

    if (WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags)
    {
        copy_attribs(pti, x, y, z);
    }

    int* AMREX_RESTRICT ion_lev = nullptr;
    if (do_field_ionization){
        ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
    }
    
    // Loop over the particles and update their momentum
    const Real q = this->charge;
    const Real m = this-> mass;
    if (WarpX::particle_pusher_algo == ParticlePusherAlgo::Boris){
        amrex::ParallelFor( 
            pti.numParticles(),
            [=] AMREX_GPU_DEVICE (long i) {
                Real qp = q;
                if (ion_lev){ qp *= ion_lev[i]; }
                UpdateMomentumBoris( ux[i], uy[i], uz[i], gi[i],
                                     Ex[i], Ey[i], Ez[i], Bx[i],
                                     By[i], Bz[i], qp, m, dt);
                UpdatePosition( x[i], y[i], z[i],
                      ux[i], uy[i], uz[i], dt );
            }
        );
    } else if (WarpX::particle_pusher_algo == ParticlePusherAlgo::Vay) {
        amrex::ParallelFor(
            pti.numParticles(),
            [=] AMREX_GPU_DEVICE (long i) {
                Real qp = q;
                if (ion_lev){ qp *= ion_lev[i]; }
                UpdateMomentumVay( ux[i], uy[i], uz[i], gi[i],
                                   Ex[i], Ey[i], Ez[i], Bx[i],
                                   By[i], Bz[i], qp, m, dt);
                UpdatePosition( x[i], y[i], z[i],
                                ux[i], uy[i], uz[i], dt );
            }
        );
    } else {
      amrex::Abort("Unknown particle pusher");
    };
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

            int e_is_nodal = Ex.is_nodal() and Ey.is_nodal() and Ez.is_nodal();
            FieldGather(pti, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                        &exfab, &eyfab, &ezfab, &bxfab, &byfab, &bzfab, 
                        Ex.nGrow(), e_is_nodal,
                        0, np, thread_num, lev, lev);

            // This wraps the momentum advance so that inheritors can modify the call.
            // Extract pointers to the different particle quantities
            Real* const AMREX_RESTRICT gi = m_giv[thread_num].dataPtr();
            Real* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
            Real* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
            Real* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
            const Real* const AMREX_RESTRICT Expp = Exp.dataPtr();
            const Real* const AMREX_RESTRICT Eypp = Eyp.dataPtr();
            const Real* const AMREX_RESTRICT Ezpp = Ezp.dataPtr();
            const Real* const AMREX_RESTRICT Bxpp = Bxp.dataPtr();
            const Real* const AMREX_RESTRICT Bypp = Byp.dataPtr();
            const Real* const AMREX_RESTRICT Bzpp = Bzp.dataPtr();

            // Loop over the particles and update their momentum
            const Real q = this->charge;
            const Real m = this-> mass;
            if (WarpX::particle_pusher_algo == ParticlePusherAlgo::Boris){
                amrex::ParallelFor( pti.numParticles(),
                    [=] AMREX_GPU_DEVICE (long i) {
                        UpdateMomentumBoris( ux[i], uy[i], uz[i], gi[i],
                              Expp[i], Eypp[i], Ezpp[i], Bxpp[i], Bypp[i], Bzpp[i], q, m, dt);
                    }
                );
            } else if (WarpX::particle_pusher_algo == ParticlePusherAlgo::Vay) {
                amrex::ParallelFor( pti.numParticles(),
                    [=] AMREX_GPU_DEVICE (long i) {
                        UpdateMomentumVay( ux[i], uy[i], uz[i], gi[i],
                              Expp[i], Eypp[i], Ezpp[i], Bxpp[i], Bypp[i], Bzpp[i], q, m, dt);
                    }
                );
            } else {
              amrex::Abort("Unknown particle pusher");
            };
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
    
    const auto np = pti.numParticles();
    const auto lev = pti.GetLevel();
    const auto index = pti.GetPairIndex();
    Real* AMREX_RESTRICT xpold  = tmp_particle_data[lev][index][TmpIdx::xold ].dataPtr();
    Real* AMREX_RESTRICT ypold  = tmp_particle_data[lev][index][TmpIdx::yold ].dataPtr();
    Real* AMREX_RESTRICT zpold  = tmp_particle_data[lev][index][TmpIdx::zold ].dataPtr();
    Real* AMREX_RESTRICT uxpold = tmp_particle_data[lev][index][TmpIdx::uxold].dataPtr();
    Real* AMREX_RESTRICT uypold = tmp_particle_data[lev][index][TmpIdx::uyold].dataPtr();
    Real* AMREX_RESTRICT uzpold = tmp_particle_data[lev][index][TmpIdx::uzold].dataPtr();

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

                auto&  xp_old = tmp_particle_data[lev][index][TmpIdx::xold];
                auto&  yp_old = tmp_particle_data[lev][index][TmpIdx::yold];
                auto&  zp_old = tmp_particle_data[lev][index][TmpIdx::zold];
                auto& uxp_old = tmp_particle_data[lev][index][TmpIdx::uxold];
                auto& uyp_old = tmp_particle_data[lev][index][TmpIdx::uyold];
                auto& uzp_old = tmp_particle_data[lev][index][TmpIdx::uzold];

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

/* \brief Gather fields from FArrayBox exfab, eyfab, ezfab, bxfab, byfab, 
 * bzfab into arrays of fields on particles Exp, Eyp, Ezp, Bxp, Byp, Bzp.
 * \param Exp-Bzp: fields on particles.
 * \param exfab-bzfab: FAB of electric and magnetic fields for particles in pti
 * \param ngE: number of guard cells for E
 * \param e_is_nodal: 0 if E is staggered, 1 if E is nodal
 * \param offset: index of first particle for which fields are gathered
 * \param np_to_gather: number of particles onto which fields are gathered
 * \param thread_num: if using OpenMP, thread number
 * \param lev: level on which particles are located
 * \param gather_lev: level from which particles gather fields (lev-1) for 
          particles in buffers.
 */
void
PhysicalParticleContainer::FieldGather (WarpXParIter& pti,
                                        RealVector& Exp,
                                        RealVector& Eyp,
                                        RealVector& Ezp,
                                        RealVector& Bxp,
                                        RealVector& Byp,
                                        RealVector& Bzp,
                                        FArrayBox const * exfab,
                                        FArrayBox const * eyfab,
                                        FArrayBox const * ezfab,
                                        FArrayBox const * bxfab,
                                        FArrayBox const * byfab,
                                        FArrayBox const * bzfab,
                                        const int ngE, const int e_is_nodal,
                                        const long offset,
                                        const long np_to_gather,
                                        int thread_num,
                                        int lev,
                                        int gather_lev)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE((gather_lev==(lev-1)) ||
                                     (gather_lev==(lev  )),
                                     "Gather buffers only work for lev-1");
    
    // If no particles, do not do anything
    if (np_to_gather == 0) return;
    // Get cell size on gather_lev
    const std::array<Real,3>& dx = WarpX::CellSize(std::max(gather_lev,0));
    // Set staggering shift depending on e_is_nodal
    const Real stagger_shift = e_is_nodal ? 0.0 : 0.5;
    
    // Get box from which field is gathered.
    // If not gathering from the finest level, the box is coarsened.
    Box box;
    if (lev == gather_lev) {
        box = pti.tilebox();
    } else {
        const IntVect& ref_ratio = WarpX::RefRatio(gather_lev);
        box = amrex::coarsen(pti.tilebox(),ref_ratio);
    }
    
    // Add guard cells to the box.
    box.grow(ngE);
    
    const Array4<const Real>& ex_arr = exfab->array();
    const Array4<const Real>& ey_arr = eyfab->array();
    const Array4<const Real>& ez_arr = ezfab->array();
    const Array4<const Real>& bx_arr = bxfab->array();
    const Array4<const Real>& by_arr = byfab->array();
    const Array4<const Real>& bz_arr = bzfab->array();
    
    const Real * const AMREX_RESTRICT xp = m_xp[thread_num].dataPtr() + offset;
    const Real * const AMREX_RESTRICT zp = m_zp[thread_num].dataPtr() + offset;
    const Real * const AMREX_RESTRICT yp = m_yp[thread_num].dataPtr() + offset;
    
    // Lower corner of tile box physical domain
    const std::array<Real, 3>& xyzmin = WarpX::LowerCorner(box, gather_lev);
    
    const Dim3 lo = lbound(box);
    
    // Depending on l_lower_in_v and WarpX::nox, call
    // different versions of template function doGatherShapeN
    if (WarpX::l_lower_order_in_v){
        if        (WarpX::nox == 1){
            doGatherShapeN<1,1>(xp, yp, zp,
                                Exp.dataPtr() + offset, Eyp.dataPtr() + offset,
                                Ezp.dataPtr() + offset, Bxp.dataPtr() + offset,
                                Byp.dataPtr() + offset, Bzp.dataPtr() + offset,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                np_to_gather, dx,
                                xyzmin, lo, stagger_shift, WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 2){
            doGatherShapeN<2,1>(xp, yp, zp,
                                Exp.dataPtr() + offset, Eyp.dataPtr() + offset,
                                Ezp.dataPtr() + offset, Bxp.dataPtr() + offset,
                                Byp.dataPtr() + offset, Bzp.dataPtr() + offset,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                np_to_gather, dx,
                                xyzmin, lo, stagger_shift, WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 3){
            doGatherShapeN<3,1>(xp, yp, zp,
                                Exp.dataPtr() + offset, Eyp.dataPtr() + offset,
                                Ezp.dataPtr() + offset, Bxp.dataPtr() + offset,
                                Byp.dataPtr() + offset, Bzp.dataPtr() + offset,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                np_to_gather, dx,
                                xyzmin, lo, stagger_shift, WarpX::n_rz_azimuthal_modes);
        }
    } else {
        if        (WarpX::nox == 1){
            doGatherShapeN<1,0>(xp, yp, zp,
                                Exp.dataPtr() + offset, Eyp.dataPtr() + offset,
                                Ezp.dataPtr() + offset, Bxp.dataPtr() + offset,
                                Byp.dataPtr() + offset, Bzp.dataPtr() + offset,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                np_to_gather, dx,
                                xyzmin, lo, stagger_shift, WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 2){
            doGatherShapeN<2,0>(xp, yp, zp,
                                Exp.dataPtr() + offset, Eyp.dataPtr() + offset,
                                Ezp.dataPtr() + offset, Bxp.dataPtr() + offset,
                                Byp.dataPtr() + offset, Bzp.dataPtr() + offset,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                np_to_gather, dx,
                                xyzmin, lo, stagger_shift, WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 3){
            doGatherShapeN<3,0>(xp, yp, zp,
                                Exp.dataPtr() + offset, Eyp.dataPtr() + offset,
                                Ezp.dataPtr() + offset, Bxp.dataPtr() + offset,
                                Byp.dataPtr() + offset, Bzp.dataPtr() + offset,
                                ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                                np_to_gather, dx,
                                xyzmin, lo, stagger_shift, WarpX::n_rz_azimuthal_modes);
        }
    }
}


void PhysicalParticleContainer::InitIonizationModule ()
{
    if (!do_field_ionization) return;
    ParmParse pp(species_name);
    if (charge != PhysConst::q_e){
        amrex::Warning(
            "charge != q_e for ionizable species: overriding user value and setting charge = q_e.");
        charge = PhysConst::q_e;
    }
    pp.query("ionization_initial_level", ionization_initial_level);
    pp.get("ionization_product_species", ionization_product_name);
    pp.get("physical_element", physical_element);
    // Add runtime integer component for ionization level
    AddIntComp("ionization_level");
    // Get atomic number and ionization energies from file
    int ion_element_id = ion_map_ids[physical_element];
    ion_atomic_number = ion_atomic_numbers[ion_element_id];
    ionization_energies.resize(ion_atomic_number);
    int offset = ion_energy_offsets[ion_element_id];
    for(int i=0; i<ion_atomic_number; i++){
        ionization_energies[i] = table_ionization_energies[i+offset];
    }
    // Compute ADK prefactors (See Chen, JCP 236 (2013), equation (2))
    // For now, we assume l=0 and m=0.
    // The approximate expressions are used, 
    // without Gamma function
    Real wa = std::pow(PhysConst::alpha,3) * PhysConst::c / PhysConst::r_e;
    Real Ea = PhysConst::m_e * PhysConst::c*PhysConst::c /PhysConst::q_e * 
        std::pow(PhysConst::alpha,4)/PhysConst::r_e;
    Real UH = table_ionization_energies[0];
    Real l_eff = std::sqrt(UH/ionization_energies[0]) - 1.;

    const Real dt = WarpX::GetInstance().getdt(0);

    adk_power.resize(ion_atomic_number);
    adk_prefactor.resize(ion_atomic_number);
    adk_exp_prefactor.resize(ion_atomic_number);
    for (int i=0; i<ion_atomic_number; ++i){
        Real n_eff = (i+1) * std::sqrt(UH/ionization_energies[i]);
        Real C2 = std::pow(2,2*n_eff)/(n_eff*tgamma(n_eff+l_eff+1)*tgamma(n_eff-l_eff));
        adk_power[i] = -(2*n_eff - 1);
        Real Uion = ionization_energies[i];
        adk_prefactor[i] = dt * wa * C2 * ( Uion/(2*UH) ) 
            * std::pow(2*std::pow((Uion/UH),3./2)*Ea,2*n_eff - 1);
        adk_exp_prefactor[i] = -2./3 * std::pow( Uion/UH,3./2) * Ea;
    }
}

/* \brief create mask of ionized particles (1 if ionized, 0 otherwise)
 * 
 * \param mfi: tile or grid
 * \param lev: MR level
 * \param ionization_mask: Array with as many elements as particles in mfi.
 * This function initialized the array, and set each element to 1 or 0 
 * depending on whether the particle is ionized or not.
 */
void
PhysicalParticleContainer::buildIonizationMask (const amrex::MFIter& mfi, const int lev,
                                                amrex::Gpu::ManagedDeviceVector<int>& ionization_mask)
{
    BL_PROFILE("PPC::buildIonizationMask");
    // Get pointers to ionization data from pc_source
    const Real * const AMREX_RESTRICT p_ionization_energies = ionization_energies.dataPtr();
    const Real * const AMREX_RESTRICT p_adk_prefactor = adk_prefactor.dataPtr();
    const Real * const AMREX_RESTRICT p_adk_exp_prefactor = adk_exp_prefactor.dataPtr();
    const Real * const AMREX_RESTRICT p_adk_power = adk_power.dataPtr();

    // Current tile info
    const int grid_id = mfi.index();
    const int tile_id = mfi.LocalTileIndex();

    // Get GPU-friendly arrays of particle data
    auto& ptile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
    // Only need attribs (i.e., SoA data)
    auto& soa = ptile.GetStructOfArrays();
    const int np = ptile.GetArrayOfStructs().size();

    // If no particle, nothing to do.
    if (np == 0) return;
    // Otherwise, resize ionization_mask, and get poiters to attribs arrays.
    ionization_mask.resize(np);
    int * const AMREX_RESTRICT p_ionization_mask = ionization_mask.data();
    const Real * const AMREX_RESTRICT ux = soa.GetRealData(PIdx::ux).data();
    const Real * const AMREX_RESTRICT uy = soa.GetRealData(PIdx::uy).data();
    const Real * const AMREX_RESTRICT uz = soa.GetRealData(PIdx::uz).data();
    const Real * const AMREX_RESTRICT ex = soa.GetRealData(PIdx::Ex).data();
    const Real * const AMREX_RESTRICT ey = soa.GetRealData(PIdx::Ey).data();
    const Real * const AMREX_RESTRICT ez = soa.GetRealData(PIdx::Ez).data();
    const Real * const AMREX_RESTRICT bx = soa.GetRealData(PIdx::Bx).data();
    const Real * const AMREX_RESTRICT by = soa.GetRealData(PIdx::By).data();
    const Real * const AMREX_RESTRICT bz = soa.GetRealData(PIdx::Bz).data();
    int* ion_lev = soa.GetIntData(particle_icomps["ionization_level"]).data();

    Real c = PhysConst::c;
    Real c2_inv = 1./c/c;
    int atomic_number = ion_atomic_number;

    // Loop over all particles in grid/tile. If ionized, set mask to 1
    // and increment ionization level.
    ParallelFor( 
        np,
        [=] AMREX_GPU_DEVICE (long i) {
            // Get index of ionization_level
            p_ionization_mask[i] = 0;
            if ( ion_lev[i]<atomic_number ){
                Real random_draw = amrex::Random();
                // Compute electric field amplitude in the particle's frame of
                // reference (particularly important when in boosted frame).
                Real ga = std::sqrt(1. + (ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i]) * c2_inv);
                Real E = std::sqrt(
                    - ( ux[i]*ex[i] + uy[i]*ey[i] + uz[i]*ez[i] ) * ( ux[i]*ex[i] + uy[i]*ey[i] + uz[i]*ez[i] ) * c2_inv
                    + ( ga   *ex[i] + uy[i]*bz[i] - uz[i]*by[i] ) * ( ga   *ex[i] + uy[i]*bz[i] - uz[i]*by[i] )
                    + ( ga   *ey[i] + uz[i]*bx[i] - ux[i]*bz[i] ) * ( ga   *ey[i] + uz[i]*bx[i] - ux[i]*bz[i] )
                    + ( ga   *ez[i] + ux[i]*by[i] - uy[i]*bx[i] ) * ( ga   *ez[i] + ux[i]*by[i] - uy[i]*bx[i] )
                    );
                // Compute probability of ionization p
                Real w_dtau = 1./ ga * p_adk_prefactor[ion_lev[i]] *
                    std::pow(E,p_adk_power[ion_lev[i]]) *
                    std::exp( p_adk_exp_prefactor[ion_lev[i]]/E );
                Real p = 1. - std::exp( - w_dtau );

                if (random_draw < p){
                    // increment particle's ionization level
                    ion_lev[i] += 1;
                    // update mask
                    p_ionization_mask[i] = 1;
                }
            }
        }
    );
}
