#include <limits>
#include <algorithm>
#include <string>

#include <MultiParticleContainer.H>
#include <WarpX_f.H>
#include <WarpX.H>

using namespace amrex;

MultiParticleContainer::MultiParticleContainer (AmrCore* amr_core)
{
    ReadParameters();

    allcontainers.resize(nspecies + nlasers);
    for (int i = 0; i < nspecies; ++i) {
        if (species_types[i] == PCTypes::Physical) {
            allcontainers[i].reset(new PhysicalParticleContainer(amr_core, i, species_names[i]));
        }
        else if (species_types[i] == PCTypes::RigidInjected) {
            allcontainers[i].reset(new RigidInjectedParticleContainer(amr_core, i, species_names[i]));
        }
        allcontainers[i]->deposit_on_main_grid = deposit_on_main_grid[i];
    }
    
    for (int i = nspecies; i < nspecies+nlasers; ++i) {
        allcontainers[i].reset(new LaserParticleContainer(amr_core,i, lasers_names[i-nspecies]));
    }

    pc_tmp.reset(new PhysicalParticleContainer(amr_core));

    // Compute the number of species for which lab-frame data is dumped
    // nspecies_lab_frame_diags, and map their ID to MultiParticleContainer
    // particle IDs in map_species_lab_diags.
    map_species_boosted_frame_diags.resize(nspecies);
    nspecies_boosted_frame_diags = 0;
    for (int i=0; i<nspecies; i++){
        auto& pc = allcontainers[i];
        if (pc->do_boosted_frame_diags){
            map_species_boosted_frame_diags[nspecies_boosted_frame_diags] = i;
            do_boosted_frame_diags = 1;
            nspecies_boosted_frame_diags += 1;
        }
    }
}

void
MultiParticleContainer::ReadParameters ()
{
    static bool initialized = false;
    if (!initialized)
    {
        ParmParse pp("particles");

        pp.query("nspecies", nspecies);
        BL_ASSERT(nspecies >= 0);

        if (nspecies > 0) {
            pp.getarr("species_names", species_names);
            BL_ASSERT(species_names.size() == nspecies);

            deposit_on_main_grid.resize(nspecies, 0);
            std::vector<std::string> tmp;
            pp.queryarr("deposit_on_main_grid", tmp);
            for (auto const& name : tmp) {
                auto it = std::find(species_names.begin(), species_names.end(), name);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != species_names.end(), "ERROR: species in particles.deposit_on_main_grid must be part of particles.species_names");
                int i = std::distance(species_names.begin(), it);
                deposit_on_main_grid[i] = 1;
            }

            species_types.resize(nspecies, PCTypes::Physical);

            std::vector<std::string> rigid_injected_species;
            pp.queryarr("rigid_injected_species", rigid_injected_species);

            if (!rigid_injected_species.empty()) {
                for (auto const& name : rigid_injected_species) {
                    auto it = std::find(species_names.begin(), species_names.end(), name);
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != species_names.end(), "ERROR: species in particles.rigid_injected_species must be part of particles.species_names");
                    int i = std::distance(species_names.begin(), it);
                    species_types[i] = PCTypes::RigidInjected;
                }
            }
        }

        pp.query("use_fdtd_nci_corr", WarpX::use_fdtd_nci_corr);
        pp.query("l_lower_order_in_v", WarpX::l_lower_order_in_v);

        ParmParse ppl("lasers");
        ppl.query("nlasers", nlasers);
        BL_ASSERT(nlasers >= 0);
        if (nlasers > 0) {
            ppl.getarr("names", lasers_names);
            BL_ASSERT(lasers_names.size() == nlasers);
        }

        initialized = true;
    }
}

void
MultiParticleContainer::AllocData ()
{
    for (auto& pc : allcontainers) {
        pc->AllocData();
    }
    pc_tmp->AllocData();
}

void
MultiParticleContainer::InitData ()
{
    for (auto& pc : allcontainers) {
        pc->InitData();
    }
    pc_tmp->InitData();
    // For each species, get the ID of its product species.
    // This is used for ionization and pair creation processes.
    mapSpeciesProduct();
}


#ifdef WARPX_DO_ELECTROSTATIC
void
MultiParticleContainer::FieldGatherES (const Vector<std::array<std::unique_ptr<MultiFab>, 3> >& E,
                                       const amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks)
{
    for (auto& pc : allcontainers) {
        pc->FieldGatherES(E, masks);
    }
}

void
MultiParticleContainer::EvolveES (const Vector<std::array<std::unique_ptr<MultiFab>, 3> >& E,
                                        Vector<std::unique_ptr<MultiFab> >& rho,
                                  Real t, Real dt)
{

    int nlevs = rho.size();
    int ng = rho[0]->nGrow();

    for (unsigned i = 0; i < nlevs; i++) {
        rho[i]->setVal(0.0, ng);
    }

    for (auto& pc : allcontainers) {
        pc->EvolveES(E, rho, t, dt);
    }

    for (unsigned i = 0; i < nlevs; i++) {
        const Geometry& gm = allcontainers[0]->Geom(i);
        rho[i]->SumBoundary(gm.periodicity());
    }
}

void
MultiParticleContainer::PushXES (Real dt)
{
    for (auto& pc : allcontainers) {
        pc->PushXES(dt);
    }
}

void
MultiParticleContainer::
DepositCharge (Vector<std::unique_ptr<MultiFab> >& rho, bool local)
{
    int nlevs = rho.size();
    int ng = rho[0]->nGrow();

    for (unsigned i = 0; i < nlevs; i++) {
        rho[i]->setVal(0.0, ng);
    }

    for (unsigned i = 0, n = allcontainers.size(); i < n; ++i) {
        allcontainers[i]->DepositCharge(rho, true);
    }

    if (!local) {
        for (unsigned i = 0; i < nlevs; i++) {
            const Geometry& gm = allcontainers[0]->Geom(i);
            rho[i]->SumBoundary(gm.periodicity());
        }
    }
}

amrex::Real
MultiParticleContainer::sumParticleCharge (bool local)
{
    amrex::Real total_charge = allcontainers[0]->sumParticleCharge(local);
    for (unsigned i = 1, n = allcontainers.size(); i < n; ++i) {
        total_charge += allcontainers[i]->sumParticleCharge(local);
    }
    return total_charge;
}

#endif // WARPX_DO_ELECTROSTATIC

void
MultiParticleContainer::FieldGather (int lev,
                                     const MultiFab& Ex, const MultiFab& Ey,
                                     const MultiFab& Ez, const MultiFab& Bx,
                                     const MultiFab& By, const MultiFab& Bz)
{
    for (auto& pc : allcontainers) {
        pc->FieldGather(lev, Ex, Ey, Ez, Bx, By, Bz);
    }
}

void
MultiParticleContainer::Evolve (int lev,
                                const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                                MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                MultiFab* cjx,  MultiFab* cjy, MultiFab* cjz,
                                MultiFab* rho, MultiFab* crho,
                                const MultiFab* cEx, const MultiFab* cEy, const MultiFab* cEz,
                                const MultiFab* cBx, const MultiFab* cBy, const MultiFab* cBz,
                                Real t, Real dt)
{
    jx.setVal(0.0);
    jy.setVal(0.0);
    jz.setVal(0.0);
    if (cjx) cjx->setVal(0.0);
    if (cjy) cjy->setVal(0.0);
    if (cjz) cjz->setVal(0.0);
    if (rho) rho->setVal(0.0);
    if (crho) crho->setVal(0.0);
    for (auto& pc : allcontainers) {
	pc->Evolve(lev, Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, cjx, cjy, cjz,
               rho, crho, cEx, cEy, cEz, cBx, cBy, cBz, t, dt);
    }
}

void
MultiParticleContainer::PushX (Real dt)
{
    for (auto& pc : allcontainers) {
        pc->PushX(dt);
    }
}

void
MultiParticleContainer::PushP (int lev, Real dt,
                               const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                               const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    for (auto& pc : allcontainers) {
        pc->PushP(lev, dt, Ex, Ey, Ez, Bx, By, Bz);
    }
}

std::unique_ptr<MultiFab>
MultiParticleContainer::GetChargeDensity (int lev, bool local)
{
    std::unique_ptr<MultiFab> rho = allcontainers[0]->GetChargeDensity(lev, true);
    for (unsigned i = 1, n = allcontainers.size(); i < n; ++i) {
        std::unique_ptr<MultiFab> rhoi = allcontainers[i]->GetChargeDensity(lev, true);
        MultiFab::Add(*rho, *rhoi, 0, 0, 1, rho->nGrow());
    }
    if (!local) {
        const Geometry& gm = allcontainers[0]->Geom(lev);
        rho->SumBoundary(gm.periodicity());
    }
    return rho;
}

void
MultiParticleContainer::SortParticlesByCell ()
{
    for (auto& pc : allcontainers) {
        pc->SortParticlesByCell();
    }
}

void
MultiParticleContainer::Redistribute ()
{
    for (auto& pc : allcontainers) {
        pc->Redistribute();
    }
}

void
MultiParticleContainer::RedistributeLocal (const int num_ghost)
{
    for (auto& pc : allcontainers) {
        pc->Redistribute(0, 0, 0, num_ghost);
    }
}

Vector<long>
MultiParticleContainer::NumberOfParticlesInGrid (int lev) const
{
    const bool only_valid=true, only_local=true;
    Vector<long> r = allcontainers[0]->NumberOfParticlesInGrid(lev,only_valid,only_local);
    for (unsigned i = 1, n = allcontainers.size(); i < n; ++i) {
        const auto& ri = allcontainers[i]->NumberOfParticlesInGrid(lev,only_valid,only_local);
        for (unsigned j=0, m=ri.size(); j<m; ++j) {
            r[j] += ri[j];
        }
    }
    ParallelDescriptor::ReduceLongSum(r.data(),r.size());
    return r;
}

void
MultiParticleContainer::Increment (MultiFab& mf, int lev)
{
    for (auto& pc : allcontainers) {
        pc->Increment(mf,lev);
    }
}

void
MultiParticleContainer::SetParticleBoxArray (int lev, BoxArray& new_ba)
{
    for (auto& pc : allcontainers) {
        pc->SetParticleBoxArray(lev,new_ba);
    }
}

void
MultiParticleContainer::SetParticleDistributionMap (int lev, DistributionMapping& new_dm)
{
    for (auto& pc : allcontainers) {
        pc->SetParticleDistributionMap(lev,new_dm);
    }
}

void
MultiParticleContainer::PostRestart ()
{
    for (auto& pc : allcontainers) {
        pc->PostRestart();
    }
    pc_tmp->PostRestart();
}

void
MultiParticleContainer
::GetLabFrameData(const std::string& snapshot_name,
                  const int i_lab, const int direction,
                  const Real z_old, const Real z_new,
                  const Real t_boost, const Real t_lab, const Real dt,
                  Vector<WarpXParticleContainer::DiagnosticParticleData>& parts) const
{

    BL_PROFILE("MultiParticleContainer::GetLabFrameData");
    
    // Loop over particle species
    for (int i = 0; i < nspecies_boosted_frame_diags; ++i){
        int isp = map_species_boosted_frame_diags[i];
        WarpXParticleContainer* pc = allcontainers[isp].get();
        WarpXParticleContainer::DiagnosticParticles diagnostic_particles;
        pc->GetParticleSlice(direction, z_old, z_new, t_boost, t_lab, dt, diagnostic_particles);
        // Here, diagnostic_particles[lev][index] is a WarpXParticleContainer::DiagnosticParticleData
        // where "lev" is the AMR level and "index" is a [grid index][tile index] pair.
        
        // Loop over AMR levels
        for (int lev = 0; lev <= pc->finestLevel(); ++lev){
            // Loop over [grid index][tile index] pairs 
            // and Fills parts[species number i] with particle data from all grids and 
            // tiles in diagnostic_particles. parts contains particles from all 
            // AMR levels indistinctly.
            for (auto it = diagnostic_particles[lev].begin(); it != diagnostic_particles[lev].end(); ++it){
                // it->first is the [grid index][tile index] key
                // it->second is the corresponding 
                // WarpXParticleContainer::DiagnosticParticleData value
                parts[i].GetRealData(DiagIdx::w).insert(  parts[i].GetRealData(DiagIdx::w  ).end(),
                                                          it->second.GetRealData(DiagIdx::w  ).begin(),
                                                          it->second.GetRealData(DiagIdx::w  ).end());
                
                parts[i].GetRealData(DiagIdx::x).insert(  parts[i].GetRealData(DiagIdx::x  ).end(),
                                                          it->second.GetRealData(DiagIdx::x  ).begin(),
                                                          it->second.GetRealData(DiagIdx::x  ).end());
                
                parts[i].GetRealData(DiagIdx::y).insert(  parts[i].GetRealData(DiagIdx::y  ).end(),
                                                          it->second.GetRealData(DiagIdx::y  ).begin(),
                                                          it->second.GetRealData(DiagIdx::y  ).end());
                
                parts[i].GetRealData(DiagIdx::z).insert(  parts[i].GetRealData(DiagIdx::z  ).end(),
                                                          it->second.GetRealData(DiagIdx::z  ).begin(),
                                                          it->second.GetRealData(DiagIdx::z  ).end());
                
                parts[i].GetRealData(DiagIdx::ux).insert(  parts[i].GetRealData(DiagIdx::ux).end(),
                                                           it->second.GetRealData(DiagIdx::ux).begin(),
                                                           it->second.GetRealData(DiagIdx::ux).end());
                
                parts[i].GetRealData(DiagIdx::uy).insert(  parts[i].GetRealData(DiagIdx::uy).end(),
                                                           it->second.GetRealData(DiagIdx::uy).begin(),
                                                           it->second.GetRealData(DiagIdx::uy).end());
                
                parts[i].GetRealData(DiagIdx::uz).insert(  parts[i].GetRealData(DiagIdx::uz).end(),
                                                           it->second.GetRealData(DiagIdx::uz).begin(),
                                                           it->second.GetRealData(DiagIdx::uz).end());
            }
        }
    }
}

/* \brief Continuous injection for particles initially outside of the domain.
 * \param injection_box: Domain where new particles should be injected.
 * Loop over all WarpXParticleContainer in MultiParticleContainer and 
 * calls virtual function ContinuousInjection.
 */
void
MultiParticleContainer::ContinuousInjection(const RealBox& injection_box) const
{
    for (int i=0; i<nspecies+nlasers; i++){
        auto& pc = allcontainers[i];
        if (pc->do_continuous_injection){
            pc->ContinuousInjection(injection_box);
        }
    }
}

/* \brief Update position of continuous injection parameters.
 * \param dt: simulation time step (level 0)
 * All classes inherited from WarpXParticleContainer do not have 
 * a position to update (PhysicalParticleContainer does not do anything).
 */
void
MultiParticleContainer::UpdateContinuousInjectionPosition(Real dt) const
{
    for (int i=0; i<nspecies+nlasers; i++){
        auto& pc = allcontainers[i];
        if (pc->do_continuous_injection){
            pc->UpdateContinuousInjectionPosition(dt);
        }
    }
}

int
MultiParticleContainer::doContinuousInjection () const
{
    int warpx_do_continuous_injection = 0;
    for (int i=0; i<nspecies+nlasers; i++){
        auto& pc = allcontainers[i];
        if (pc->do_continuous_injection){
            warpx_do_continuous_injection = 1;
        }
    }
    return warpx_do_continuous_injection;
}

/* \brief Get ID of product species of each species.
 * The users specifies the name of the product species, 
 * this routine get its ID.
 */
void
MultiParticleContainer::mapSpeciesProduct ()
{
    for (int i=0; i<nspecies; i++){
        auto& pc = allcontainers[i];
        // If species pc has ionization on, find species with name 
        // pc->ionization_product_name and store its ID into 
        // pc->ionization_product.
        if (pc->do_field_ionization){
            int i_product = getSpeciesID(pc->ionization_product_name);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                i != i_product,
                "ERROR: ionization product cannot be the same species");
            pc->ionization_product = i_product;
        }
    }
}

/* \brief Given a species name, return its ID.
 */
int
MultiParticleContainer::getSpeciesID (std::string product_str)
{
    int i_product;
    bool found = 0;
    // Loop over species
    for (int i=0; i<nspecies; i++){
        // If species name matches, store its ID
        // into i_product
        if (species_names[i] == product_str){
            found = 1;
            i_product = i;
        }
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        found != 0,
        "ERROR: could not find product species ID for ionization. Wrong name?");
    return i_product;
}

namespace
{
    // For particle i in mfi, if is_ionized[i]=1, copy particle
    // particle i from container pc_source into pc_product
    void createIonizedParticles (
        int lev, const MFIter& mfi,
        std::unique_ptr< WarpXParticleContainer>& pc_source,
        std::unique_ptr< WarpXParticleContainer>& pc_product,
        amrex::Gpu::ManagedDeviceVector<int>& is_ionized)
    {
        BL_PROFILE("createIonizedParticles");

        const int * const AMREX_RESTRICT p_is_ionized = is_ionized.dataPtr();

        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();

        // Get source particle data
        auto& ptile_source = pc_source->GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        const int np_source = ptile_source.GetArrayOfStructs().size();
        if (np_source == 0) return;
        // --- source AoS particle data
        WarpXParticleContainer::ParticleType* particles_source = ptile_source.GetArrayOfStructs()().data();
        // --- source SoA particle data
        auto& soa_source = ptile_source.GetStructOfArrays();
        GpuArray<Real*,PIdx::nattribs> attribs_source;
        for (int ia = 0; ia < PIdx::nattribs; ++ia) {
            attribs_source[ia] = soa_source.GetRealData(ia).data();
        }
        // --- source runtime attribs
        GpuArray<Real*,3> runtime_uold_source;
        // Prepare arrays for boosted frame diagnostics.
        runtime_uold_source[0] = soa_source.GetRealData(PIdx::ux).data();
        runtime_uold_source[1] = soa_source.GetRealData(PIdx::uy).data();
        runtime_uold_source[2] = soa_source.GetRealData(PIdx::uz).data();

        // Indices of product particle for each ionized source particle.
        // i_product[i]-1 is the location in product tile of product particle
        // from source particle i.
        amrex::Gpu::ManagedDeviceVector<int> i_product;
        i_product.resize(np_source);
        // 0<i<np_source
        // 0<i_product<np_ionized
        // Strictly speaking, i_product should be an exclusive_scan of
        // is_ionized. However, for indices where is_ionized is 1, the
        // inclusive scan gives the same result with an offset of 1.
        // The advantage of inclusive_scan is that the sum of is_ionized
        // is in the last element, so no other reduction is required to get
        // number of particles.
        // Gpu::inclusive_scan runs on the current GPU stream, and synchronizes
        // with the CPU, so that the next line (executed by the CPU) has the
        // updated values of i_product
        amrex::Gpu::inclusive_scan(is_ionized.begin(), is_ionized.end(), i_product.begin());
        int np_ionized = i_product[np_source-1];
        if (np_ionized == 0) return;
        int* AMREX_RESTRICT p_i_product = i_product.dataPtr();

        // Get product particle data
        auto& ptile_product = pc_product->GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        // old and new (i.e., including ionized particles) number of particles
        // for product species
        const int np_product_old = ptile_product.GetArrayOfStructs().size();
        const int np_product_new = np_product_old + np_ionized;
        // Allocate extra space in product species for ionized particles.
        ptile_product.resize(np_product_new);
        // --- product AoS particle data
        // First element is the first newly-created product particle
        WarpXParticleContainer::ParticleType* particles_product = ptile_product.GetArrayOfStructs()().data() + np_product_old;
        // --- product SoA particle data
        auto& soa_product = ptile_product.GetStructOfArrays();
        GpuArray<Real*,PIdx::nattribs> attribs_product;
        for (int ia = 0; ia < PIdx::nattribs; ++ia) {
            // First element is the first newly-created product particle
            attribs_product[ia] = soa_product.GetRealData(ia).data() + np_product_old;
        }
        // --- product runtime attribs
        GpuArray<Real*,6> runtime_attribs_product;
        bool do_boosted_product = WarpX::do_boosted_frame_diagnostic 
            && pc_product->DoBoostedFrameDiags();
        if (do_boosted_product) {
            std::map<std::string, int> comps_product = pc_product->getParticleComps();
            runtime_attribs_product[0] = soa_product.GetRealData(comps_product[ "xold"]).data() + np_product_old;
            runtime_attribs_product[1] = soa_product.GetRealData(comps_product[ "yold"]).data() + np_product_old;
            runtime_attribs_product[2] = soa_product.GetRealData(comps_product[ "zold"]).data() + np_product_old;
            runtime_attribs_product[3] = soa_product.GetRealData(comps_product["uxold"]).data() + np_product_old;
            runtime_attribs_product[4] = soa_product.GetRealData(comps_product["uyold"]).data() + np_product_old;
            runtime_attribs_product[5] = soa_product.GetRealData(comps_product["uzold"]).data() + np_product_old;
        }

        int pid_product;
#pragma omp critical (doFieldIonization_nextid)
        {
            // ID of first newly-created product particle
            pid_product = pc_product->NextID();
            // Update NextID to include particles created in this function
            pc_product->setNextID(pid_product+np_ionized);
        }
        const int cpuid = ParallelDescriptor::MyProc();

        // Loop over all source particles. If is_ionized, copy particle data
        // to corresponding product particle.
        amrex::For(
            np_source, [=] AMREX_GPU_DEVICE (int is) noexcept
            {
                if(p_is_ionized[is]){
                    // offset of 1 due to inclusive scan
                    int ip = p_i_product[is]-1;
                    // is: index of ionized particle in source species
                    // ip: index of corresponding new particle in product species
                    WarpXParticleContainer::ParticleType& p_product = particles_product[ip];
                    WarpXParticleContainer::ParticleType& p_source  = particles_source[is];
                    // Copy particle from source to product: AoS
                    p_product.id() = pid_product + ip;
                    p_product.cpu() = cpuid;
                    p_product.pos(0) = p_source.pos(0);
                    p_product.pos(1) = p_source.pos(1);
#if (AMREX_SPACEDIM == 3)
                    p_product.pos(2) = p_source.pos(2);
#endif
                    // Copy particle from source to product: SoA
                    for (int ia = 0; ia < PIdx::nattribs; ++ia) {
                        attribs_product[ia][ip] = attribs_source[ia][is];
                    }
                    // Update xold etc. if boosted frame diagnostics required
                    // for product species. Fill runtime attribs with a copy of 
                    // current properties (xold = x etc.).
                    if (do_boosted_product) {
                        runtime_attribs_product[0][ip] = p_source.pos(0);
                        runtime_attribs_product[1][ip] = p_source.pos(1);
                        runtime_attribs_product[2][ip] = p_source.pos(2);
                        runtime_attribs_product[3][ip] = runtime_uold_source[0][ip];
                        runtime_attribs_product[4][ip] = runtime_uold_source[1][ip];
                        runtime_attribs_product[5][ip] = runtime_uold_source[2][ip];
                    }
                }
            }
        );
    }
}

void
MultiParticleContainer::doFieldIonization ()
{
    BL_PROFILE("MPC::doFieldIonization");
    // Loop over all species.
    // Ionized particles in pc_source create particles in pc_product
    for (auto& pc_source : allcontainers){
    
        // Skip if not ionizable
        if (!pc_source->do_field_ionization){ continue; }

        // Get product species
        auto& pc_product = allcontainers[pc_source->ionization_product];

        for (int lev = 0; lev <= pc_source->finestLevel(); ++lev){

            // When using runtime components, AMReX requires to touch all tiles
            // in serial and create particles tiles with runtime components if
            // they do not exist (or if they were defined by default, i.e.,
            // without runtime component).
#ifdef _OPENMP
            // Touch all tiles of source species in serial if runtime attribs 
            for (MFIter mfi = pc_source->MakeMFIter(lev); mfi.isValid(); ++mfi) {
                const int grid_id = mfi.index();
                const int tile_id = mfi.LocalTileIndex();
                pc_source->GetParticles(lev)[std::make_pair(grid_id,tile_id)];
                if ( (pc_source->NumRuntimeRealComps()>0) || (pc_source->NumRuntimeIntComps()>0) ) {
                    pc_source->DefineAndReturnParticleTile(lev, grid_id, tile_id);
                }
            }
#endif
            // Touch all tiles of product species in serial
            for (MFIter mfi = pc_source->MakeMFIter(lev); mfi.isValid(); ++mfi) {
                const int grid_id = mfi.index();
                const int tile_id = mfi.LocalTileIndex();
                pc_product->GetParticles(lev)[std::make_pair(grid_id,tile_id)];
                pc_product->DefineAndReturnParticleTile(lev, grid_id, tile_id);
            }

            // Enable tiling
            MFItInfo info;
            if (pc_source->do_tiling && Gpu::notInLaunchRegion()) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    pc_product->do_tiling, 
                    "For ionization, either all or none of the "
                    "particle species must use tiling.");
                info.EnableTiling(pc_source->tile_size);
            }

#ifdef _OPENMP
            info.SetDynamic(true);
#pragma omp parallel
#endif
            // Loop over all grids (if not tiling) or grids and tiles (if tiling)
            for (MFIter mfi = pc_source->MakeMFIter(lev, info); mfi.isValid(); ++mfi)
            {
                // Ionization mask: one element per source particles.
                // 0 if not ionized, 1 if ionized.
                amrex::Gpu::ManagedDeviceVector<int> is_ionized;
                pc_source->buildIonizationMask(mfi, lev, is_ionized);
                // Create particles in pc_product
                createIonizedParticles(lev, mfi, pc_source, pc_product, is_ionized);
            }
        } // lev
    } // pc_source
}
