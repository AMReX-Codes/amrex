/* Copyright 2019-2020 Andrew Myers, Aurore Blelly, Axel Huebl
 * David Grote, Glenn Richardson, Jean-Luc Vay
 * Ligia Diana Amorim, Luca Fedeli, Maxence Thevenet
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 * Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PhysicalParticleContainer.H"

#include "MultiParticleContainer.H"
#include "FortranInterface/WarpX_f.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"
#include "Python/WarpXWrappers.h"
#include "Utils/IonizationEnergiesTable.H"
#include "Particles/Gather/FieldGather.H"
#include "Particles/Pusher/GetAndSetPosition.H"

#include "Utils/WarpXAlgorithmSelection.H"

// Import low-level single-particle kernels
#include "Particles/Pusher/UpdatePosition.H"
#include "Particles/Pusher/UpdateMomentumBoris.H"
#include "Particles/Pusher/UpdateMomentumVay.H"
#include "Particles/Pusher/UpdateMomentumBorisWithRadiationReaction.H"
#include "Particles/Pusher/UpdateMomentumHigueraCary.H"

#include <limits>
#include <sstream>
#include <string>


using namespace amrex;

PhysicalParticleContainer::PhysicalParticleContainer (AmrCore* amr_core, int ispecies,
                                                      const std::string& name)
    : WarpXParticleContainer(amr_core, ispecies),
      species_name(name)
{
    plasma_injector.reset(new PlasmaInjector(species_id, species_name));
    physical_species = plasma_injector->getPhysicalSpecies();
    charge = plasma_injector->getCharge();
    mass = plasma_injector->getMass();

    ParmParse pp(species_name);

    pp.query("boost_adjust_transverse_positions", boost_adjust_transverse_positions);
    pp.query("do_backward_propagation", do_backward_propagation);

    // Initialize splitting
    pp.query("do_splitting", do_splitting);
    pp.query("split_type", split_type);
    pp.query("do_not_deposit", do_not_deposit);
    pp.query("do_not_gather", do_not_gather);
    pp.query("do_not_push", do_not_push);

    pp.query("do_continuous_injection", do_continuous_injection);
    pp.query("initialize_self_fields", initialize_self_fields);
    pp.query("self_fields_required_precision", self_fields_required_precision);
    // Whether to plot back-transformed (lab-frame) diagnostics
    // for this species.
    pp.query("do_back_transformed_diagnostics", do_back_transformed_diagnostics);

    pp.query("do_field_ionization", do_field_ionization);

    //check if Radiation Reaction is enabled and do consistency checks
    pp.query("do_classical_radiation_reaction", do_classical_radiation_reaction);
    //if the species is not a lepton, do_classical_radiation_reaction
    //should be false
    WarpXUtilMsg::AlwaysAssert(
        !(do_classical_radiation_reaction &&
        !(AmIA<PhysicalSpecies::electron>() ||
        AmIA<PhysicalSpecies::positron>() )),
        "ERROR: can't enable classical radiation reaction for non lepton species '"
        + species_name + "'."
    );

    //Only Boris pusher is compatible with radiation reaction
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !(do_classical_radiation_reaction &&
        WarpX::particle_pusher_algo != ParticlePusherAlgo::Boris),
        "Radiation reaction can be enabled only if Boris pusher is used");
    //_____________________________

#ifdef WARPX_QED
    pp.query("do_qed", m_do_qed);
    if(m_do_qed){
        //If do_qed is enabled, find out if Quantum Synchrotron process is enabled
        pp.query("do_qed_quantum_sync", m_do_qed_quantum_sync);
        if (m_do_qed_quantum_sync)
            AddRealComp("optical_depth_QSR");
        pp.query("do_qed_breit_wheeler", m_do_qed_breit_wheeler);
        if (m_do_qed_breit_wheeler)
            AddRealComp("optical_depth_BW");
    }

    if(m_do_qed_quantum_sync){
        pp.get("qed_quantum_sync_phot_product_species",
            m_qed_quantum_sync_phot_product_name);
    }


#endif

    //variable to set plot_flags size
    int plot_flag_size = PIdx::nattribs;
    if(WarpX::do_back_transformed_diagnostics && do_back_transformed_diagnostics)
        plot_flag_size += 6;

#ifdef WARPX_QED
    if(m_do_qed){
        // plot_flag will have an entry for the optical depth
        plot_flag_size++;
    }
#endif
    //_______________________________

    pp.query("plot_species", plot_species);
    int do_user_plot_vars;
    do_user_plot_vars = pp.queryarr("plot_vars", plot_vars);
    if (not do_user_plot_vars){
        // By default, all particle variables are dumped to plotfiles,
        // including {x,y,z,ux,uy,uz}old variables when running in a
        // boosted frame
        plot_flags.resize(plot_flag_size, 1);
    } else {
        // Set plot_flag to 0 for all attribs
        plot_flags.resize(plot_flag_size, 0);

        // If not none, set plot_flags values to 1 for elements in plot_vars.
        if (plot_vars[0] != "none"){
            for (const auto& var : plot_vars){
                // Return error if var not in PIdx.
                WarpXUtilMsg::AlwaysAssert(
                    ParticleStringNames::to_index.count(var),
                    "ERROR: plot_vars argument '" + var +
                    "' not in ParticleStringNames"
                );
                plot_flags[ParticleStringNames::to_index.at(var)] = 1;
            }
        }

#ifdef WARPX_DIM_RZ
        // Always write out theta, whether or not it's requested,
        // to be consistent with always writing out r and z.
        plot_flags[ParticleStringNames::to_index.at("theta")] = 1;
#endif

    }

    // Parse galilean velocity
    ParmParse ppsatd("psatd");
    ppsatd.query("v_galilean", v_galilean);
    // Scale the velocity by the speed of light
    for (int i=0; i<3; i++) v_galilean[i] *= PhysConst::c;

    #ifdef WARPX_QED
        if(m_do_qed){
            //Optical depths is always plotted if QED is on
            plot_flags[plot_flag_size-1] = 1;
        }
    #endif
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

    // Allocate temporary vectors on the CPU
    Gpu::HostVector<ParticleReal> particle_x;
    Gpu::HostVector<ParticleReal> particle_y;
    Gpu::HostVector<ParticleReal> particle_z;
    Gpu::HostVector<ParticleReal> particle_ux;
    Gpu::HostVector<ParticleReal> particle_uy;
    Gpu::HostVector<ParticleReal> particle_uz;
    Gpu::HostVector<ParticleReal> particle_w;
    int np = 0;

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
                    CheckAndAddParticle(x, y, z, { u.x, u.y, u.z}, weight/4.,
                                        particle_x,  particle_y,  particle_z,
                                        particle_ux, particle_uy, particle_uz,
                                        particle_w);
                    CheckAndAddParticle(x, -y, z, { u.x, -u.y, u.z}, weight/4.,
                                        particle_x,  particle_y,  particle_z,
                                        particle_ux, particle_uy, particle_uz,
                                        particle_w);
                    CheckAndAddParticle(-x, y, z, { -u.x, u.y, u.z}, weight/4.,
                                        particle_x,  particle_y,  particle_z,
                                        particle_ux, particle_uy, particle_uz,
                                        particle_w);
                    CheckAndAddParticle(-x, -y, z, { -u.x, -u.y, u.z}, weight/4.,
                                        particle_x,  particle_y,  particle_z,
                                        particle_ux, particle_uy, particle_uz,
                                        particle_w);
                } else {
                    CheckAndAddParticle(x, y, z, { u.x, u.y, u.z}, weight,
                                        particle_x,  particle_y,  particle_z,
                                        particle_ux, particle_uy, particle_uz,
                                        particle_w);
                }
            }
        }
    }
    // Add the temporary CPU vectors to the particle structure
    np = particle_z.size();
    AddNParticles(0,np,
                  particle_x.dataPtr(),  particle_y.dataPtr(),  particle_z.dataPtr(),
                  particle_ux.dataPtr(), particle_uy.dataPtr(), particle_uz.dataPtr(),
                  1, particle_w.dataPtr(),1);
}

void
PhysicalParticleContainer::CheckAndAddParticle(Real x, Real y, Real z,
                                               std::array<Real, 3> u,
                                               Real weight,
                                               Gpu::HostVector<ParticleReal>& particle_x,
                                               Gpu::HostVector<ParticleReal>& particle_y,
                                               Gpu::HostVector<ParticleReal>& particle_z,
                                               Gpu::HostVector<ParticleReal>& particle_ux,
                                               Gpu::HostVector<ParticleReal>& particle_uy,
                                               Gpu::HostVector<ParticleReal>& particle_uz,
                                               Gpu::HostVector<ParticleReal>& particle_w)
{
    if (WarpX::gamma_boost > 1.) {
        MapParticletoBoostedFrame(x, y, z, u);
    }
    particle_x.push_back(x);
    particle_y.push_back(y);
    particle_z.push_back(z);
    particle_ux.push_back(u[0]);
    particle_uy.push_back(u[1]);
    particle_uz.push_back(u[2]);
    particle_w.push_back(weight);
}

void
PhysicalParticleContainer::AddParticles (int lev)
{
    WARPX_PROFILE("PhysicalParticleContainer::AddParticles()");

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
    WARPX_PROFILE("PhysicalParticleContainer::AddPlasma");

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

    defineAllParticleTiles();

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
                overlap_realbox.setLo( dir, part_realbox.lo(dir) + std::max(ncells_adjust, 0._rt) * dx[dir]);
            } else {
                no_overlap = true; break;
            }
            if ( tile_realbox.hi(dir) >= part_realbox.lo(dir) ) {
                Real ncells_adjust = std::floor( (part_realbox.hi(dir) - tile_realbox.hi(dir))/dx[dir] );
                overlap_realbox.setHi( dir, part_realbox.hi(dir) - std::max(ncells_adjust, 0._rt) * dx[dir]);
            } else {
                no_overlap = true; break;
            }
            // Count the number of cells in this direction in overlap_realbox
            overlap_box.setSmall( dir, 0 );
            overlap_box.setBig( dir,
                int( std::round((overlap_realbox.hi(dir)-overlap_realbox.lo(dir))
                                /dx[dir] )) - 1);
            shifted[dir] =
                static_cast<int>(std::round((overlap_realbox.lo(dir)-problo[dir])/dx[dir]));
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
#ifdef _OPENMP
#pragma omp critical (add_plasma_nextid)
#endif
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
        GpuArray<ParticleReal*,PIdx::nattribs> pa;
        for (int ia = 0; ia < PIdx::nattribs; ++ia) {
            pa[ia] = soa.GetRealData(ia).data() + old_size;
        }

        int* pi;
        if (do_field_ionization) {
            pi = soa.GetIntData(particle_icomps["ionization_level"]).data() + old_size;
        }

#ifdef WARPX_QED
        //Pointer to the optical depth component
        amrex::Real* p_optical_depth_QSR;
        amrex::Real* p_optical_depth_BW;

        // If a QED effect is enabled, the corresponding optical depth
        // has to be initialized
        bool loc_has_quantum_sync = has_quantum_sync();
        bool loc_has_breit_wheeler = has_breit_wheeler();
        if (loc_has_quantum_sync)
            p_optical_depth_QSR = soa.GetRealData(
                particle_comps["optical_depth_QSR"]).data() + old_size;
        if(loc_has_breit_wheeler)
            p_optical_depth_BW = soa.GetRealData(
                particle_comps["optical_depth_BW"]).data() + old_size;

        //If needed, get the appropriate functors from the engines
        QuantumSynchrotronGetOpticalDepth quantum_sync_get_opt;
        BreitWheelerGetOpticalDepth breit_wheeler_get_opt;
        if(loc_has_quantum_sync){
            quantum_sync_get_opt =
                m_shr_p_qs_engine->build_optical_depth_functor();
        }
        if(loc_has_breit_wheeler){
            breit_wheeler_get_opt =
                m_shr_p_bw_engine->build_optical_depth_functor();
        }
#endif

        const GpuArray<Real,AMREX_SPACEDIM> overlap_corner
            {AMREX_D_DECL(overlap_realbox.lo(0),
                          overlap_realbox.lo(1),
                          overlap_realbox.lo(2))};

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

            const XDim3 r =
                inj_pos->getPositionUnitBox(i_part, static_cast<int>(fac));
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

                // include ballistic correction for plasma species with bulk motion
                const XDim3 u_bulk = inj_mom->getBulkMomentum(x, y, z);
                const Real gamma_bulk = std::sqrt(1.+(u_bulk.x*u_bulk.x+u_bulk.y*u_bulk.y+u_bulk.z*u_bulk.z));
                const Real betaz_bulk = u_bulk.z/gamma_bulk;
                const Real z0 = z - PhysConst::c*t*betaz_bulk;

                if (!inj_pos->insideBounds(xb, yb, z0)) {
                    p.id() = -1;
                    return;
                }

                u = inj_mom->getMomentum(x, y, z0);
                dens = inj_rho->getDensity(x, y, z0);
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
                //
                // Note that we use the bulk momentum to perform the ballastic correction
                const XDim3 u_bulk = inj_mom->getBulkMomentum(x, y, 0.); // No z0_lab dependency
                // At this point u is the lab-frame momentum
                // => Apply the above formula for z0_lab
                const Real gamma_lab_bulk = std::sqrt(1.+(u_bulk.x*u_bulk.x+u_bulk.y*u_bulk.y+u_bulk.z*u_bulk.z));
                const Real betaz_lab_bulk = u_bulk.z/(gamma_lab_bulk);
                const Real z0_lab = gamma_boost * ( z*(1-beta_boost*betaz_lab_bulk)
                                              - PhysConst::c*t*(betaz_lab_bulk-beta_boost) );
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

                // get the full momentum, including thermal motion
                u = inj_mom->getMomentum(x, y, 0.);
                const Real gamma_lab = std::sqrt( 1.+(u.x*u.x+u.y*u.y+u.z*u.z) );
                const Real betaz_lab = u.z/(gamma_lab);

                // At this point u and dens are the lab-frame quantities
                // => Perform Lorentz transform
                dens = gamma_boost * dens * ( 1.0 - beta_boost*betaz_lab );
                u.z = gamma_boost * ( u.z -beta_boost*gamma_lab );
            }

            if (loc_do_field_ionization) {
                pi[ip] = loc_ionization_initial_level;
            }

#ifdef WARPX_QED
            if(loc_has_quantum_sync){
                p_optical_depth_QSR[ip] = quantum_sync_get_opt();
            }

            if(loc_has_breit_wheeler){
                p_optical_depth_BW[ip] = breit_wheeler_get_opt();
            }
#endif

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
        });

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

void
PhysicalParticleContainer::AssignExternalFieldOnParticles(WarpXParIter& pti,
                           RealVector& Exp, RealVector& Eyp, RealVector& Ezp,
                           RealVector& Bxp, RealVector& Byp, RealVector& Bzp, int lev)
{
   const long np = pti.numParticles();
    /// get WarpX class object
    auto & warpx = WarpX::GetInstance();
    /// get MultiParticleContainer class object
    auto & mypc = warpx.GetPartContainer();
   if (mypc.m_E_ext_particle_s=="constant" ||
       mypc.m_E_ext_particle_s=="default") {
       Exp.assign(np,mypc.m_E_external_particle[0]);
       Eyp.assign(np,mypc.m_E_external_particle[1]);
       Ezp.assign(np,mypc.m_E_external_particle[2]);
   }
   if (mypc.m_B_ext_particle_s=="constant" ||
       mypc.m_B_ext_particle_s=="default") {
       Bxp.assign(np,mypc.m_B_external_particle[0]);
       Byp.assign(np,mypc.m_B_external_particle[1]);
       Bzp.assign(np,mypc.m_B_external_particle[2]);
   }
   if (mypc.m_E_ext_particle_s=="parse_e_ext_particle_function") {
      const auto GetPosition = GetParticlePosition(pti);
      Real* const AMREX_RESTRICT Exp_data = Exp.dataPtr();
      Real* const AMREX_RESTRICT Eyp_data = Eyp.dataPtr();
      Real* const AMREX_RESTRICT Ezp_data = Ezp.dataPtr();
      ParserWrapper<4> *xfield_partparser = mypc.m_Ex_particle_parser.get();
      ParserWrapper<4> *yfield_partparser = mypc.m_Ey_particle_parser.get();
      ParserWrapper<4> *zfield_partparser = mypc.m_Ez_particle_parser.get();
      Real time = warpx.gett_new(lev);
      amrex::ParallelFor(pti.numParticles(),
                         [=] AMREX_GPU_DEVICE (long i) {
                             ParticleReal x, y, z;
                             GetPosition(i, x, y, z);
                             Exp_data[i] = (*xfield_partparser)(x, y, z, time);
                             Eyp_data[i] = (*yfield_partparser)(x, y, z, time);
                             Ezp_data[i] = (*zfield_partparser)(x, y, z, time);
                         });
   }
   if (mypc.m_B_ext_particle_s=="parse_b_ext_particle_function") {
      const auto GetPosition = GetParticlePosition(pti);
      Real* const AMREX_RESTRICT Bxp_data = Bxp.dataPtr();
      Real* const AMREX_RESTRICT Byp_data = Byp.dataPtr();
      Real* const AMREX_RESTRICT Bzp_data = Bzp.dataPtr();
      ParserWrapper<4> *xfield_partparser = mypc.m_Bx_particle_parser.get();
      ParserWrapper<4> *yfield_partparser = mypc.m_By_particle_parser.get();
      ParserWrapper<4> *zfield_partparser = mypc.m_Bz_particle_parser.get();
      Real time = warpx.gett_new(lev);
      amrex::ParallelFor(pti.numParticles(),
            [=] AMREX_GPU_DEVICE (long i) {
                             ParticleReal x, y, z;
                             GetPosition(i, x, y, z);
                             Bxp_data[i] = (*xfield_partparser)(x, y, z, time);
                             Byp_data[i] = (*yfield_partparser)(x, y, z, time);
                             Bzp_data[i] = (*zfield_partparser)(x, y, z, time);
                         });
   }
}



void
PhysicalParticleContainer::FieldGather (int lev,
                                        const amrex::MultiFab& Ex,
                                        const amrex::MultiFab& Ey,
                                        const amrex::MultiFab& Ez,
                                        const amrex::MultiFab& Bx,
                                        const amrex::MultiFab& By,
                                        const amrex::MultiFab& Bz)
{
    const std::array<Real,3>& dx = WarpX::CellSize(lev);

    BL_ASSERT(OnSameGrids(lev,Ex));

    MultiFab* cost = WarpX::getCosts(lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
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

            //
            // Field Gather
            //
            int e_is_nodal = Ex.is_nodal() and Ey.is_nodal() and Ez.is_nodal();
            FieldGather(pti, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                        &exfab, &eyfab, &ezfab, &bxfab, &byfab, &bzfab,
                        Ex.nGrow(), e_is_nodal,
                        0, np, lev, lev);

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
                                   Real /*t*/, Real dt, DtType a_dt_type)
{
    WARPX_PROFILE("PPC::Evolve()");
    WARPX_PROFILE_VAR_NS("PPC::Evolve::Copy", blp_copy);
    WARPX_PROFILE_VAR_NS("PPC::FieldGather", blp_fg);
    WARPX_PROFILE_VAR_NS("PPC::EvolveOpticalDepth", blp_ppc_qed_ev);
    WARPX_PROFILE_VAR_NS("PPC::ParticlePush", blp_ppc_pp);

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

    if (WarpX::do_back_transformed_diagnostics && do_back_transformed_diagnostics)
    {
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const auto np = pti.numParticles();
            const auto t_lev = pti.GetLevel();
            const auto index = pti.GetPairIndex();
            tmp_particle_data.resize(finestLevel()+1);
            for (int i = 0; i < TmpIdx::nattribs; ++i)
                tmp_particle_data[t_lev][index][i].resize(np);
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
                // Filter arrays Ex[pti], store the result in
                // filtered_Ex and update pointer exfab so that it
                // points to filtered_Ex (and do the same for all
                // components of E and B).
                applyNCIFilter(lev, pti.tilebox(), exeli, eyeli, ezeli, bxeli, byeli, bzeli,
                               filtered_Ex, filtered_Ey, filtered_Ez,
                               filtered_Bx, filtered_By, filtered_Bz,
                               Ex[pti], Ey[pti], Ez[pti], Bx[pti], By[pti], Bz[pti],
                               exfab, eyfab, ezfab, bxfab, byfab, bzfab);
            }

            // Determine which particles deposit/gather in the buffer, and
            // which particles deposit/gather in the fine patch
            long nfine_current = np;
            long nfine_gather = np;
            if (has_buffer && !do_not_push) {
                // - Modify `nfine_current` and `nfine_gather` (in place)
                //    so that they correspond to the number of particles
                //    that deposit/gather in the fine patch respectively.
                // - Reorder the particle arrays,
                //    so that the `nfine_current`/`nfine_gather` first particles
                //    deposit/gather in the fine patch
                //    and (thus) the `np-nfine_current`/`np-nfine_gather` last particles
                //    deposit/gather in the buffer
                PartitionParticlesInBuffers( nfine_current, nfine_gather, np,
                    pti, lev, current_masks, gather_masks, uxp, uyp, uzp, wp );
            }

            const long np_current = (cjx) ? nfine_current : np;

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
                WARPX_PROFILE_VAR_START(blp_fg);
                FieldGather(pti, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                            exfab, eyfab, ezfab, bxfab, byfab, bzfab,
                            Ex.nGrow(), e_is_nodal,
                            0, np_gather, lev, lev);

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
                        // Filter arrays (*cEx)[pti], store the result in
                        // filtered_Ex and update pointer cexfab so that it
                        // points to filtered_Ex (and do the same for all
                        // components of E and B)
                        applyNCIFilter(lev-1, cbox, exeli, eyeli, ezeli, bxeli, byeli, bzeli,
                                       filtered_Ex, filtered_Ey, filtered_Ez,
                                       filtered_Bx, filtered_By, filtered_Bz,
                                       (*cEx)[pti], (*cEy)[pti], (*cEz)[pti],
                                       (*cBx)[pti], (*cBy)[pti], (*cBz)[pti],
                                       cexfab, ceyfab, cezfab, cbxfab, cbyfab, cbzfab);
                    }

                    // Field gather for particles in gather buffers
                    e_is_nodal = cEx->is_nodal() and cEy->is_nodal() and cEz->is_nodal();
                    FieldGather(pti, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                cexfab, ceyfab, cezfab,
                                cbxfab, cbyfab, cbzfab,
                                cEx->nGrow(), e_is_nodal,
                                nfine_gather, np-nfine_gather,
                                lev, lev-1);
                }

                WARPX_PROFILE_VAR_STOP(blp_fg);

#ifdef WARPX_QED
                //
                //Evolve Optical Depth
                //
                WARPX_PROFILE_VAR_START(blp_ppc_qed_ev);
                EvolveOpticalDepth(pti, dt);
                WARPX_PROFILE_VAR_STOP(blp_ppc_qed_ev);
#endif

                //
                // Particle Push
                //
                WARPX_PROFILE_VAR_START(blp_ppc_pp);
                PushPX(pti, dt, a_dt_type);
                WARPX_PROFILE_VAR_STOP(blp_ppc_pp);

                //
                // Current Deposition (only needed for electromagnetic solver)
                //
                if (!WarpX::do_electrostatic) {
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
                } // end of "if !do_electrostatic"
            } // end of "if do_not_push"

            if (rho) {
                // Deposit charge after particle push, in component 1 of MultiFab rho.
                // (Skipped for electrostatic solver, as this may lead to out-of-bounds)
                if (!WarpX::do_electrostatic) {
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
    // Split particles at the end of the timestep.
    // When subcycling is ON, the splitting is done on the last call to
    // PhysicalParticleContainer::Evolve on the finest level, i.e., at the
    // end of the large timestep. Otherwise, the pushes on different levels
    // are not consistent, and the call to Redistribute (inside
    // SplitParticles) may result in split particles to deposit twice on the
    // coarse level.
    if (do_splitting && (a_dt_type == DtType::SecondHalf || a_dt_type == DtType::Full) ){
        SplitParticles(lev);
    }
}

void
PhysicalParticleContainer::applyNCIFilter (
    int lev, const Box& box,
    Elixir& exeli, Elixir& eyeli, Elixir& ezeli,
    Elixir& bxeli, Elixir& byeli, Elixir& bzeli,
    FArrayBox& filtered_Ex, FArrayBox& filtered_Ey, FArrayBox& filtered_Ez,
    FArrayBox& filtered_Bx, FArrayBox& filtered_By, FArrayBox& filtered_Bz,
    const FArrayBox& Ex, const FArrayBox& Ey, const FArrayBox& Ez,
    const FArrayBox& Bx, const FArrayBox& By, const FArrayBox& Bz,
    FArrayBox const * & ex_ptr, FArrayBox const * & ey_ptr,
    FArrayBox const * & ez_ptr, FArrayBox const * & bx_ptr,
    FArrayBox const * & by_ptr, FArrayBox const * & bz_ptr)
{

    // Get instances of NCI Godfrey filters
    const auto& nci_godfrey_filter_exeybz = WarpX::GetInstance().nci_godfrey_filter_exeybz;
    const auto& nci_godfrey_filter_bxbyez = WarpX::GetInstance().nci_godfrey_filter_bxbyez;

#if (AMREX_SPACEDIM == 2)
    const Box& tbox = amrex::grow(box,{static_cast<int>(WarpX::nox),
                static_cast<int>(WarpX::noz)});
#else
    const Box& tbox = amrex::grow(box,{static_cast<int>(WarpX::nox),
                static_cast<int>(WarpX::noy),
                static_cast<int>(WarpX::noz)});
#endif

    // Filter Ex (Both 2D and 3D)
    filtered_Ex.resize(amrex::convert(tbox,Ex.box().ixType()));
    // Safeguard for GPU
    exeli = filtered_Ex.elixir();
    // Apply filter on Ex, result stored in filtered_Ex

    nci_godfrey_filter_exeybz[lev]->ApplyStencil(filtered_Ex, Ex, filtered_Ex.box());
    // Update ex_ptr reference
    ex_ptr = &filtered_Ex;

    // Filter Ez
    filtered_Ez.resize(amrex::convert(tbox,Ez.box().ixType()));
    ezeli = filtered_Ez.elixir();
    nci_godfrey_filter_bxbyez[lev]->ApplyStencil(filtered_Ez, Ez, filtered_Ez.box());
    ez_ptr = &filtered_Ez;

    // Filter By
    filtered_By.resize(amrex::convert(tbox,By.box().ixType()));
    byeli = filtered_By.elixir();
    nci_godfrey_filter_bxbyez[lev]->ApplyStencil(filtered_By, By, filtered_By.box());
    by_ptr = &filtered_By;
#if (AMREX_SPACEDIM == 3)
    // Filter Ey
    filtered_Ey.resize(amrex::convert(tbox,Ey.box().ixType()));
    eyeli = filtered_Ey.elixir();
    nci_godfrey_filter_exeybz[lev]->ApplyStencil(filtered_Ey, Ey, filtered_Ey.box());
    ey_ptr = &filtered_Ey;

    // Filter Bx
    filtered_Bx.resize(amrex::convert(tbox,Bx.box().ixType()));
    bxeli = filtered_Bx.elixir();
    nci_godfrey_filter_bxbyez[lev]->ApplyStencil(filtered_Bx, Bx, filtered_Bx.box());
    bx_ptr = &filtered_Bx;

    // Filter Bz
    filtered_Bz.resize(amrex::convert(tbox,Bz.box().ixType()));
    bzeli = filtered_Bz.elixir();
    nci_godfrey_filter_exeybz[lev]->ApplyStencil(filtered_Bz, Bz, filtered_Bz.box());
    bz_ptr = &filtered_Bz;
#endif
}

// Loop over all particles in the particle container and
// split particles tagged with p.id()=DoSplitParticleID
void
PhysicalParticleContainer::SplitParticles(int lev)
{
    auto& mypc = WarpX::GetInstance().GetPartContainer();
    auto& pctmp_split = mypc.GetPCtmp();
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
        const auto GetPosition = GetParticlePosition(pti);

        const amrex::Vector<int> ppc_nd = plasma_injector->num_particles_per_cell_each_dim;
        const std::array<Real,3>& dx = WarpX::CellSize(lev);
        amrex::Vector<Real> split_offset = {dx[0]/2._rt,
                                            dx[1]/2._rt,
                                            dx[2]/2._rt};
        if (ppc_nd[0] > 0){
            // offset for split particles is computed as a function of cell size
            // and number of particles per cell, so that a uniform distribution
            // before splitting results in a uniform distribution after splitting
            split_offset[0] /= ppc_nd[0];
            split_offset[1] /= ppc_nd[1];
            split_offset[2] /= ppc_nd[2];
        }
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
            ParticleReal xp, yp, zp;
            GetPosition(i, xp, yp, zp);
            auto& p = particles[i];
            if (p.id() == DoSplitParticleID){
                // If particle is tagged, split it and put the
                // split particles in local arrays psplit_x etc.
                np_split_to_add += np_split;
#if (AMREX_SPACEDIM==2)
                if (split_type==0){
                    // Split particle in two along each diagonals
                    // 4 particles in 2d
                    for (int ishift = -1; ishift < 2; ishift +=2 ){
                        for (int kshift = -1; kshift < 2; kshift +=2 ){
                            // Add one particle with offset in x and z
                            psplit_x.push_back( xp + ishift*split_offset[0] );
                            psplit_y.push_back( yp );
                            psplit_z.push_back( zp + kshift*split_offset[2] );
                            psplit_ux.push_back( uxp[i] );
                            psplit_uy.push_back( uyp[i] );
                            psplit_uz.push_back( uzp[i] );
                            psplit_w.push_back( wp[i]/np_split );
                        }
                    }
                } else {
                    // Split particle in two along each axis
                    // 4 particles in 2d
                    for (int ishift = -1; ishift < 2; ishift +=2 ){
                        // Add one particle with offset in x
                        psplit_x.push_back( xp + ishift*split_offset[0] );
                        psplit_y.push_back( yp );
                        psplit_z.push_back( zp );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                        // Add one particle with offset in z
                        psplit_x.push_back( xp );
                        psplit_y.push_back( yp );
                        psplit_z.push_back( zp + ishift*split_offset[2] );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                    }
                }
#elif (AMREX_SPACEDIM==3)
                if (split_type==0){
                    // Split particle in two along each diagonals
                    // 8 particles in 3d
                    for (int ishift = -1; ishift < 2; ishift +=2 ){
                        for (int jshift = -1; jshift < 2; jshift +=2 ){
                            for (int kshift = -1; kshift < 2; kshift +=2 ){
                                // Add one particle with offset in x, y and z
                                psplit_x.push_back( xp + ishift*split_offset[0] );
                                psplit_y.push_back( yp + jshift*split_offset[1] );
                                psplit_z.push_back( zp + kshift*split_offset[2] );
                                psplit_ux.push_back( uxp[i] );
                                psplit_uy.push_back( uyp[i] );
                                psplit_uz.push_back( uzp[i] );
                                psplit_w.push_back( wp[i]/np_split );
                            }
                        }
                    }
                } else {
                    // Split particle in two along each axis
                    // 6 particles in 3d
                    for (int ishift = -1; ishift < 2; ishift +=2 ){
                        // Add one particle with offset in x
                        psplit_x.push_back( xp + ishift*split_offset[0] );
                        psplit_y.push_back( yp );
                        psplit_z.push_back( zp );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                        // Add one particle with offset in y
                        psplit_x.push_back( xp );
                        psplit_y.push_back( yp + ishift*split_offset[1] );
                        psplit_z.push_back( zp );
                        psplit_ux.push_back( uxp[i] );
                        psplit_uy.push_back( uyp[i] );
                        psplit_uz.push_back( uzp[i] );
                        psplit_w.push_back( wp[i]/np_split );
                        // Add one particle with offset in z
                        psplit_x.push_back( xp );
                        psplit_y.push_back( yp );
                        psplit_z.push_back( zp + ishift*split_offset[2] );
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
PhysicalParticleContainer::PushPX (WarpXParIter& pti, Real dt, DtType a_dt_type)
{

    // This wraps the momentum and position advance so that inheritors can modify the call.
    auto& attribs = pti.GetAttribs();
    // Extract pointers to the different particle quantities

    const auto GetPosition = GetParticlePosition(pti);
          auto SetPosition = SetParticlePosition(pti);

    ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Ex = attribs[PIdx::Ex].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Ey = attribs[PIdx::Ey].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Ez = attribs[PIdx::Ez].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Bx = attribs[PIdx::Bx].dataPtr();
    const ParticleReal* const AMREX_RESTRICT By = attribs[PIdx::By].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Bz = attribs[PIdx::Bz].dataPtr();

    if (WarpX::do_back_transformed_diagnostics && do_back_transformed_diagnostics && (a_dt_type!=DtType::SecondHalf))
    {
        copy_attribs(pti);
    }

    int* AMREX_RESTRICT ion_lev = nullptr;
    if (do_field_ionization){
        ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
    }

    // Loop over the particles and update their momentum
    const Real q = this->charge;
    const Real m = this-> mass;

#ifdef WARPX_QED
    if(do_classical_radiation_reaction){
        if(m_do_qed_quantum_sync){
            const auto t_chi_max = m_shr_p_qs_engine->get_ref_ctrl().chi_part_min;
            amrex::ParallelFor(
                pti.numParticles(),
                [=] AMREX_GPU_DEVICE (long i) {
                    auto chi = QedUtils::chi_lepton(m*ux[i], m*uy[i], m*uz[i],
                         Ex[i], Ey[i], Ez[i],
                         Bx[i], By[i], Bz[i]);
                    if(chi < t_chi_max){
                        UpdateMomentumBorisWithRadiationReaction( ux[i], uy[i], uz[i],
                                           Ex[i], Ey[i], Ez[i], Bx[i],
                                           By[i], Bz[i], q, m, dt);
                    }
                    else{
                        UpdateMomentumBoris( ux[i], uy[i], uz[i],
                                           Ex[i], Ey[i], Ez[i], Bx[i],
                                           By[i], Bz[i], q, m, dt);
                    }
                    ParticleReal x, y, z;
                    GetPosition(i, x, y, z);
                    UpdatePosition(x, y, z, ux[i], uy[i], uz[i], dt );
                    SetPosition(i, x, y, z);
                }
            );
        }else{
            amrex::ParallelFor(
                pti.numParticles(),
                [=] AMREX_GPU_DEVICE (long i) {
                    UpdateMomentumBorisWithRadiationReaction( ux[i], uy[i], uz[i],
                                       Ex[i], Ey[i], Ez[i], Bx[i],
                                       By[i], Bz[i], q, m, dt);
                    ParticleReal x, y, z;
                    GetPosition(i, x, y, z);
                    UpdatePosition(x, y, z, ux[i], uy[i], uz[i], dt );
                    SetPosition(i, x, y, z);
                }
            );
        }
#else
    if(do_classical_radiation_reaction){
        amrex::ParallelFor(
            pti.numParticles(),
            [=] AMREX_GPU_DEVICE (long i) {
                Real qp = q;
                if (ion_lev){ qp *= ion_lev[i]; }
                UpdateMomentumBorisWithRadiationReaction( ux[i], uy[i], uz[i],
                                   Ex[i], Ey[i], Ez[i], Bx[i],
                                   By[i], Bz[i], qp, m, dt);
                ParticleReal x, y, z;
                GetPosition(i, x, y, z);
                UpdatePosition(x, y, z, ux[i], uy[i], uz[i], dt );
                SetPosition(i, x, y, z);
            }
        );
#endif
    } else if (WarpX::particle_pusher_algo == ParticlePusherAlgo::Boris){
        amrex::ParallelFor(
            pti.numParticles(),
            [=] AMREX_GPU_DEVICE (long i) {
                Real qp = q;
                if (ion_lev){ qp *= ion_lev[i]; }
                UpdateMomentumBoris( ux[i], uy[i], uz[i],
                                     Ex[i], Ey[i], Ez[i], Bx[i],
                                     By[i], Bz[i], qp, m, dt);
                ParticleReal x, y, z;
                GetPosition(i, x, y, z);
                UpdatePosition(x, y, z, ux[i], uy[i], uz[i], dt );
                SetPosition(i, x, y, z);
            }
        );
    } else if (WarpX::particle_pusher_algo == ParticlePusherAlgo::Vay) {
        amrex::ParallelFor(
            pti.numParticles(),
            [=] AMREX_GPU_DEVICE (long i) {
                Real qp = q;
                if (ion_lev){ qp *= ion_lev[i]; }
                UpdateMomentumVay( ux[i], uy[i], uz[i],
                                   Ex[i], Ey[i], Ez[i], Bx[i],
                                   By[i], Bz[i], qp, m, dt);
                ParticleReal x, y, z;
                GetPosition(i, x, y, z);
                UpdatePosition(x, y, z, ux[i], uy[i], uz[i], dt );
                SetPosition(i, x, y, z);
            }
        );
    } else if (WarpX::particle_pusher_algo == ParticlePusherAlgo::HigueraCary) {
        amrex::ParallelFor(
            pti.numParticles(),
            [=] AMREX_GPU_DEVICE (long i) {
                Real qp = q;
                if (ion_lev){ qp *= ion_lev[i]; }
                UpdateMomentumHigueraCary( ux[i], uy[i], uz[i],
                                   Ex[i], Ey[i], Ez[i], Bx[i],
                                   By[i], Bz[i], qp, m, dt);
                ParticleReal x, y, z;
                GetPosition(i, x, y, z);
                UpdatePosition(x, y, z, ux[i], uy[i], uz[i], dt );
                SetPosition(i, x, y, z);
            }
        );
    } else {
      amrex::Abort("Unknown particle pusher");
    };
}

#ifdef WARPX_QED
void PhysicalParticleContainer::EvolveOpticalDepth(
    WarpXParIter& pti, amrex::Real dt)
{
    if(!has_quantum_sync())
        return;

    QuantumSynchrotronEvolveOpticalDepth evolve_opt =
        m_shr_p_qs_engine->build_evolve_functor();

    auto& attribs = pti.GetAttribs();
    const ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
    const ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
    const ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Ex = attribs[PIdx::Ex].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Ey = attribs[PIdx::Ey].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Ez = attribs[PIdx::Ez].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Bx = attribs[PIdx::Bx].dataPtr();
    const ParticleReal* const AMREX_RESTRICT By = attribs[PIdx::By].dataPtr();
    const ParticleReal* const AMREX_RESTRICT Bz = attribs[PIdx::Bz].dataPtr();

    ParticleReal* const AMREX_RESTRICT p_optical_depth_QSR =
        pti.GetAttribs(particle_comps["optical_depth_QSR"]).dataPtr();

    const ParticleReal m = this->mass;

    amrex::ParallelFor(pti.numParticles(),
            [=] AMREX_GPU_DEVICE (long i) {
                const ParticleReal px = m * ux[i];
                const ParticleReal py = m * uy[i];
                const ParticleReal pz = m * uz[i];

                bool has_event_happened = evolve_opt(
                    px, py, pz,
                    Ex[i], Ey[i], Ez[i],
                    Bx[i], By[i], Bz[i],
                    dt, p_optical_depth_QSR[i]);
            }
    );

}
#endif

void
PhysicalParticleContainer::PushP (int lev, Real dt,
                                  const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                  const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    WARPX_PROFILE("PhysicalParticleContainer::PushP");

    if (do_not_push) return;

    const std::array<Real,3>& dx = WarpX::CellSize(lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
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

            int e_is_nodal = Ex.is_nodal() and Ey.is_nodal() and Ez.is_nodal();
            FieldGather(pti, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                        &exfab, &eyfab, &ezfab, &bxfab, &byfab, &bzfab,
                        Ex.nGrow(), e_is_nodal,
                        0, np, lev, lev);

            // This wraps the momentum advance so that inheritors can modify the call.
            // Extract pointers to the different particle quantities
            ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
            ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
            ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();
            const ParticleReal* const AMREX_RESTRICT Expp = Exp.dataPtr();
            const ParticleReal* const AMREX_RESTRICT Eypp = Eyp.dataPtr();
            const ParticleReal* const AMREX_RESTRICT Ezpp = Ezp.dataPtr();
            const ParticleReal* const AMREX_RESTRICT Bxpp = Bxp.dataPtr();
            const ParticleReal* const AMREX_RESTRICT Bypp = Byp.dataPtr();
            const ParticleReal* const AMREX_RESTRICT Bzpp = Bzp.dataPtr();

            // Loop over the particles and update their momentum
            const Real q = this->charge;
            const Real m = this-> mass;

            int* AMREX_RESTRICT ion_lev = nullptr;
            if (do_field_ionization){
                ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
            }

            //Assumes that all consistency checks have been done at initialization
            if(do_classical_radiation_reaction){
                 amrex::ParallelFor(pti.numParticles(),
                    [=] AMREX_GPU_DEVICE (long i){
                        Real qp = q;
                        if (ion_lev){ qp *= ion_lev[i]; }
                        UpdateMomentumBorisWithRadiationReaction(
                            ux[i], uy[i], uz[i],
                            Expp[i], Eypp[i], Ezpp[i],
                            Bxpp[i], Bypp[i], Bzpp[i],
                            qp, m, dt);
                    }
                );
            } else if (WarpX::particle_pusher_algo == ParticlePusherAlgo::Boris){
                amrex::ParallelFor(pti.numParticles(),
                    [=] AMREX_GPU_DEVICE (long i) {
                        Real qp = q;
                        if (ion_lev){ qp *= ion_lev[i]; }
                        UpdateMomentumBoris(
                            ux[i], uy[i], uz[i],
                            Expp[i], Eypp[i], Ezpp[i],
                            Bxpp[i], Bypp[i], Bzpp[i],
                            qp, m, dt);
                    }
                );
            } else if (WarpX::particle_pusher_algo == ParticlePusherAlgo::Vay){
                amrex::ParallelFor(pti.numParticles(),
                    [=] AMREX_GPU_DEVICE (long i) {
                        Real qp = q;
                        if (ion_lev){ qp *= ion_lev[i]; }
                        UpdateMomentumVay(
                            ux[i], uy[i], uz[i],
                            Expp[i], Eypp[i], Ezpp[i],
                            Bxpp[i], Bypp[i], Bzpp[i],
                            qp, m, dt);
                    }
                );
            } else if (WarpX::particle_pusher_algo == ParticlePusherAlgo::HigueraCary){
                amrex::ParallelFor(pti.numParticles(),
                    [=] AMREX_GPU_DEVICE (long i) {
                        UpdateMomentumHigueraCary( ux[i], uy[i], uz[i],
                            Expp[i], Eypp[i], Ezpp[i],
                            Bxpp[i], Bypp[i], Bzpp[i],
                            q, m, dt);
                    }
                );
            } else {
              amrex::Abort("Unknown particle pusher");
            }

        }
    }
}

void PhysicalParticleContainer::copy_attribs (WarpXParIter& pti)
{
    auto& attribs = pti.GetAttribs();
    ParticleReal* AMREX_RESTRICT uxp = attribs[PIdx::ux].dataPtr();
    ParticleReal* AMREX_RESTRICT uyp = attribs[PIdx::uy].dataPtr();
    ParticleReal* AMREX_RESTRICT uzp = attribs[PIdx::uz].dataPtr();

    const auto np = pti.numParticles();
    const auto lev = pti.GetLevel();
    const auto index = pti.GetPairIndex();
    ParticleReal* AMREX_RESTRICT xpold  = tmp_particle_data[lev][index][TmpIdx::xold ].dataPtr();
    ParticleReal* AMREX_RESTRICT ypold  = tmp_particle_data[lev][index][TmpIdx::yold ].dataPtr();
    ParticleReal* AMREX_RESTRICT zpold  = tmp_particle_data[lev][index][TmpIdx::zold ].dataPtr();
    ParticleReal* AMREX_RESTRICT uxpold = tmp_particle_data[lev][index][TmpIdx::uxold].dataPtr();
    ParticleReal* AMREX_RESTRICT uypold = tmp_particle_data[lev][index][TmpIdx::uyold].dataPtr();
    ParticleReal* AMREX_RESTRICT uzpold = tmp_particle_data[lev][index][TmpIdx::uzold].dataPtr();

    const auto GetPosition = GetParticlePosition(pti);

    ParallelFor( np,
                 [=] AMREX_GPU_DEVICE (long i) {
                     ParticleReal x, y, z;
                     GetPosition(i, x, y, z);
                     xpold[i]=x;
                     ypold[i]=y;
                     zpold[i]=z;

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
    WARPX_PROFILE("PhysicalParticleContainer::GetParticleSlice");

    // Assume that the boost in the positive z direction.
#if (AMREX_SPACEDIM == 2)
    AMREX_ALWAYS_ASSERT(direction == 1);
#else
    AMREX_ALWAYS_ASSERT(direction == 2);
#endif

    // Note the the slice should always move in the negative boost direction.
    AMREX_ALWAYS_ASSERT(z_new < z_old);

    AMREX_ALWAYS_ASSERT(do_back_transformed_diagnostics == 1);

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
            for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
            {
                const Box& box = pti.validbox();
                auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
                const RealBox tile_real_box(box, dx, plo);

                if ( !slice_box.intersects(tile_real_box) ) continue;

                const auto GetPosition = GetParticlePosition(pti);

                auto& attribs = pti.GetAttribs();
                Real* const AMREX_RESTRICT wpnew = attribs[PIdx::w].dataPtr();
                Real* const AMREX_RESTRICT uxpnew = attribs[PIdx::ux].dataPtr();
                Real* const AMREX_RESTRICT uypnew = attribs[PIdx::uy].dataPtr();
                Real* const AMREX_RESTRICT uzpnew = attribs[PIdx::uz].dataPtr();

                Real* const AMREX_RESTRICT
                  xpold = tmp_particle_data[lev][index][TmpIdx::xold].dataPtr();
                Real* const AMREX_RESTRICT
                  ypold = tmp_particle_data[lev][index][TmpIdx::yold].dataPtr();
                Real* const AMREX_RESTRICT
                  zpold = tmp_particle_data[lev][index][TmpIdx::zold].dataPtr();
                Real* const AMREX_RESTRICT
                  uxpold = tmp_particle_data[lev][index][TmpIdx::uxold].dataPtr();
                Real* const AMREX_RESTRICT
                  uypold = tmp_particle_data[lev][index][TmpIdx::uyold].dataPtr();
                Real* const AMREX_RESTRICT
                  uzpold = tmp_particle_data[lev][index][TmpIdx::uzold].dataPtr();

                const long np = pti.numParticles();

                Real uzfrm = -WarpX::gamma_boost*WarpX::beta_boost*PhysConst::c;
                Real inv_c2 = 1.0/PhysConst::c/PhysConst::c;

                // temporary arrays to store copy_flag and copy_index
                // for particles that cross the z-slice
                amrex::Gpu::ManagedDeviceVector<int> FlagForPartCopy(np);
                amrex::Gpu::ManagedDeviceVector<int> IndexForPartCopy(np);

                int* const AMREX_RESTRICT Flag = FlagForPartCopy.dataPtr();
                int* const AMREX_RESTRICT IndexLocation = IndexForPartCopy.dataPtr();

                //Flag particles that need to be copied if they cross the z_slice
                amrex::ParallelFor(np,
                [=] AMREX_GPU_DEVICE(int i)
                {
                    ParticleReal xp, yp, zp;
                    GetPosition(i, xp, yp, zp);
                    Flag[i] = 0;
                    if ( (((zp >= z_new) && (zpold[i] <= z_old)) ||
                          ((zp <= z_new) && (zpold[i] >= z_old))) )
                    {
                        Flag[i] = 1;
                    }
                });

                // exclusive scan to obtain location indices using flag values
                // These location indices are used to copy data from
                // src to dst when the copy-flag is set to 1.
                amrex::Gpu::exclusive_scan(Flag,Flag+np,IndexLocation);

                const int total_partdiag_size = IndexLocation[np-1] + Flag[np-1];

                // allocate array size for diagnostic particle array
                diagnostic_particles[lev][index].resize(total_partdiag_size);

                amrex::Real gammaboost = WarpX::gamma_boost;
                amrex::Real betaboost = WarpX::beta_boost;
                amrex::Real Phys_c = PhysConst::c;

                Real* const AMREX_RESTRICT diag_wp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::w).data();
                Real* const AMREX_RESTRICT diag_xp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::x).data();
                Real* const AMREX_RESTRICT diag_yp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::y).data();
                Real* const AMREX_RESTRICT diag_zp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::z).data();
                Real* const AMREX_RESTRICT diag_uxp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::ux).data();
                Real* const AMREX_RESTRICT diag_uyp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::uy).data();
                Real* const AMREX_RESTRICT diag_uzp =
                diagnostic_particles[lev][index].GetRealData(DiagIdx::uz).data();

                // Copy particle data to diagnostic particle array on the GPU
                //  using flag and index values
                amrex::ParallelFor(np,
                [=] AMREX_GPU_DEVICE(int i)
                {
                    ParticleReal xp_new, yp_new, zp_new;
                    GetPosition(i, xp_new, yp_new, zp_new);
                    if (Flag[i] == 1)
                    {
                         // Lorentz Transform particles to lab-frame
                         const Real gamma_new_p = std::sqrt(1.0 + inv_c2*
                                                  (uxpnew[i]*uxpnew[i]
                                                 + uypnew[i]*uypnew[i]
                                                 + uzpnew[i]*uzpnew[i]));
                         const Real t_new_p = gammaboost*t_boost - uzfrm*zp_new*inv_c2;
                         const Real z_new_p = gammaboost*(zp_new + betaboost*Phys_c*t_boost);
                         const Real uz_new_p = gammaboost*uzpnew[i] - gamma_new_p*uzfrm;

                         const Real gamma_old_p = std::sqrt(1.0 + inv_c2*
                                                  (uxpold[i]*uxpold[i]
                                                 + uypold[i]*uypold[i]
                                                 + uzpold[i]*uzpold[i]));
                         const Real t_old_p = gammaboost*(t_boost - dt)
                                              - uzfrm*zpold[i]*inv_c2;
                         const Real z_old_p = gammaboost*(zpold[i]
                                              + betaboost*Phys_c*(t_boost-dt));
                         const Real uz_old_p = gammaboost*uzpold[i]
                                              - gamma_old_p*uzfrm;

                         // interpolate in time to t_lab
                         const Real weight_old = (t_new_p - t_lab)
                                               / (t_new_p - t_old_p);
                         const Real weight_new = (t_lab - t_old_p)
                                               / (t_new_p - t_old_p);

                         const Real xp = xpold[i]*weight_old + xp_new*weight_new;
                         const Real yp = ypold[i]*weight_old + yp_new*weight_new;
                         const Real zp = z_old_p*weight_old  + z_new_p*weight_new;

                         const Real uxp = uxpold[i]*weight_old
                                        + uxpnew[i]*weight_new;
                         const Real uyp = uypold[i]*weight_old
                                        + uypnew[i]*weight_new;
                         const Real uzp = uz_old_p*weight_old
                                        + uz_new_p  *weight_new;

                         const int loc = IndexLocation[i];
                         diag_wp[loc] = wpnew[i];
                         diag_xp[loc] = xp;
                         diag_yp[loc] = yp;
                         diag_zp[loc] = zp;
                         diag_uxp[loc] = uxp;
                         diag_uyp[loc] = uyp;
                         diag_uzp[loc] = uzp;
                    }
                });
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
                                        amrex::FArrayBox const * exfab,
                                        amrex::FArrayBox const * eyfab,
                                        amrex::FArrayBox const * ezfab,
                                        amrex::FArrayBox const * bxfab,
                                        amrex::FArrayBox const * byfab,
                                        amrex::FArrayBox const * bzfab,
                                        const int ngE, const int /*e_is_nodal*/,
                                        const long offset,
                                        const long np_to_gather,
                                        int lev,
                                        int gather_lev)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE((gather_lev==(lev-1)) ||
                                     (gather_lev==(lev  )),
                                     "Gather buffers only work for lev-1");
    // If no particles, do not do anything
    // If do_not_gather = 1 by user, do not do anything
    if (np_to_gather == 0 || do_not_gather) return;

    // initializing the field value to the externally applied field before
    // gathering fields from the grid to the particles.
    AssignExternalFieldOnParticles(pti, Exp, Eyp, Ezp, Bxp, Byp, Bzp, lev);


    // Get cell size on gather_lev
    const std::array<Real,3>& dx = WarpX::CellSize(std::max(gather_lev,0));

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

    const auto GetPosition = GetParticlePosition(pti, offset);

    // Lower corner of tile box physical domain (take into account Galilean shift)
    Real cur_time = WarpX::GetInstance().gett_new(lev);
    const auto& time_of_last_gal_shift = WarpX::GetInstance().time_of_last_gal_shift;
    Real time_shift = (cur_time - time_of_last_gal_shift);
    amrex::Array<amrex::Real,3> galilean_shift = { v_galilean[0]*time_shift, v_galilean[1]*time_shift, v_galilean[2]*time_shift };
    const std::array<Real, 3>& xyzmin = WarpX::LowerCorner(box, galilean_shift, gather_lev);

    const Dim3 lo = lbound(box);

    // Depending on l_lower_in_v and WarpX::nox, call
    // different versions of template function doGatherShapeN
    if (WarpX::l_lower_order_in_v){
        if        (WarpX::nox == 1){
            doGatherShapeN<1,1>(GetPosition,
                                Exp.dataPtr() + offset, Eyp.dataPtr() + offset,
                                Ezp.dataPtr() + offset, Bxp.dataPtr() + offset,
                                Byp.dataPtr() + offset, Bzp.dataPtr() + offset,
                                exfab, eyfab, ezfab, bxfab, byfab, bzfab,
                                np_to_gather, dx,
                                xyzmin, lo, WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 2){
            doGatherShapeN<2,1>(GetPosition,
                                Exp.dataPtr() + offset, Eyp.dataPtr() + offset,
                                Ezp.dataPtr() + offset, Bxp.dataPtr() + offset,
                                Byp.dataPtr() + offset, Bzp.dataPtr() + offset,
                                exfab, eyfab, ezfab, bxfab, byfab, bzfab,
                                np_to_gather, dx,
                                xyzmin, lo, WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 3){
            doGatherShapeN<3,1>(GetPosition,
                                Exp.dataPtr() + offset, Eyp.dataPtr() + offset,
                                Ezp.dataPtr() + offset, Bxp.dataPtr() + offset,
                                Byp.dataPtr() + offset, Bzp.dataPtr() + offset,
                                exfab, eyfab, ezfab, bxfab, byfab, bzfab,
                                np_to_gather, dx,
                                xyzmin, lo, WarpX::n_rz_azimuthal_modes);
        }
    } else {
        if        (WarpX::nox == 1){
            doGatherShapeN<1,0>(GetPosition,
                                Exp.dataPtr() + offset, Eyp.dataPtr() + offset,
                                Ezp.dataPtr() + offset, Bxp.dataPtr() + offset,
                                Byp.dataPtr() + offset, Bzp.dataPtr() + offset,
                                exfab, eyfab, ezfab, bxfab, byfab, bzfab,
                                np_to_gather, dx,
                                xyzmin, lo, WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 2){
            doGatherShapeN<2,0>(GetPosition,
                                Exp.dataPtr() + offset, Eyp.dataPtr() + offset,
                                Ezp.dataPtr() + offset, Bxp.dataPtr() + offset,
                                Byp.dataPtr() + offset, Bzp.dataPtr() + offset,
                                exfab, eyfab, ezfab, bxfab, byfab, bzfab,
                                np_to_gather, dx,
                                xyzmin, lo, WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 3){
            doGatherShapeN<3,0>(GetPosition,
                                Exp.dataPtr() + offset, Eyp.dataPtr() + offset,
                                Ezp.dataPtr() + offset, Bxp.dataPtr() + offset,
                                Byp.dataPtr() + offset, Bzp.dataPtr() + offset,
                                exfab, eyfab, ezfab, bxfab, byfab, bzfab,
                                np_to_gather, dx,
                                xyzmin, lo, WarpX::n_rz_azimuthal_modes);
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

IonizationFilterFunc
PhysicalParticleContainer::getIonizationFunc ()
{
    WARPX_PROFILE("PPC::getIonizationFunc");

    return IonizationFilterFunc{ionization_energies.dataPtr(),
                                adk_prefactor.dataPtr(),
                                adk_exp_prefactor.dataPtr(),
                                adk_power.dataPtr(),
                                particle_icomps["ionization_level"],
                                ion_atomic_number};
}

#ifdef WARPX_QED


bool PhysicalParticleContainer::has_quantum_sync()
{
    return m_do_qed_quantum_sync;
}

bool PhysicalParticleContainer::has_breit_wheeler()
{
    return m_do_qed_breit_wheeler;
}

void
PhysicalParticleContainer::
set_breit_wheeler_engine_ptr(std::shared_ptr<BreitWheelerEngine> ptr)
{
    m_shr_p_bw_engine = ptr;
}

void
PhysicalParticleContainer::
set_quantum_sync_engine_ptr(std::shared_ptr<QuantumSynchrotronEngine> ptr)
{
    m_shr_p_qs_engine = ptr;
}

PhotonEmissionFilterFunc
PhysicalParticleContainer::getPhotonEmissionFilterFunc ()
{
    WARPX_PROFILE("PPC::getPhotonEmissionFunc");
    return PhotonEmissionFilterFunc{particle_runtime_comps["optical_depth_QSR"]};
}

PairGenerationFilterFunc
PhysicalParticleContainer::getPairGenerationFilterFunc ()
{
    WARPX_PROFILE("PPC::getPairGenerationFunc");
    return PairGenerationFilterFunc{particle_runtime_comps["optical_depth_BW"]};
}

#endif
