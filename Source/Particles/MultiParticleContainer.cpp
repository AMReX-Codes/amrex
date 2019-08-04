#include <limits>
#include <algorithm>
#include <string>
#include <math.h>

#include <AMReX_Particles.H>  
#include <AMReX_AmrCore.H>
#include <AMReX_Gpu.H>
#include <MultiParticleContainer.H>
#include <WarpX_f.H>
#include <WarpX.H>
#include <WarpXConst.H>

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

    // For each species, get the ID of its product species.
    // This is used for ionization and pair creation processes.
    mapSpeciesProduct();
    
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

    if (WarpX::do_boosted_frame_diagnostic && do_boosted_frame_diags)
    {
        for (int i = 0; i < nspecies_boosted_frame_diags; ++i)
        {
            int is = map_species_boosted_frame_diags[i];
            allcontainers[is]->AddRealComp("xold");
            allcontainers[is]->AddRealComp("yold");
            allcontainers[is]->AddRealComp("zold");
            allcontainers[is]->AddRealComp("uxold");
            allcontainers[is]->AddRealComp("uyold");
            allcontainers[is]->AddRealComp("uzold");
        }
        pc_tmp->AddRealComp("xold");
        pc_tmp->AddRealComp("yold");
        pc_tmp->AddRealComp("zold");
        pc_tmp->AddRealComp("uxold");
        pc_tmp->AddRealComp("uyold");
        pc_tmp->AddRealComp("uzold");
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
MultiParticleContainer::doContinuousInjection() const
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
MultiParticleContainer::mapSpeciesProduct()
{
    for (int i=0; i<nspecies; i++){
        auto& pc = allcontainers[i];
        // If species pc has ionization on, find species with name 
        // pc->ionization_product_name and store its ID into 
        // pc->ionization_product.
        if (pc->do_field_ionization){
            int i_product = getSpeciesID(pc->ionization_product_name);
            pc->ionization_product = i_product;
        }
    }
}

/* \brief Given a species name, return its ID.
 */
int
MultiParticleContainer::getSpeciesID(std::string product_str)
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
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(found != 0,
                                     "ERROR: could not find ID for species");
    return i_product;
}

void
MultiParticleContainer::doFieldIonization()
{
    /*
      amrex::Gpu::SharedMemory<unsigned short int> grid_ids;
      amrex::Gpu::SharedMemory<unsigned short int> tile_ids;
      amrex::Gpu::SharedMemory<bool> is_ionized;
    */

    for (auto& pc : allcontainers){
        if (!pc->do_field_ionization){ continue; }
        const Real * const AMREX_RESTRICT p_ionization_energies = pc->ionization_energies.dataPtr();
        const Real * const AMREX_RESTRICT p_adk_prefactor = pc->adk_prefactor.dataPtr();
        const Real * const AMREX_RESTRICT p_adk_exp_prefactor = pc->adk_exp_prefactor.dataPtr();
        const Real * const AMREX_RESTRICT p_adk_power = pc->adk_power.dataPtr();
        for (int lev = 0; lev <= pc->finestLevel(); ++lev){

#ifdef _OPENMP
            // First touch all tiles in the map in serial
            for (MFIter mfi = pc->MakeMFIter(lev); mfi.isValid(); ++mfi) {
                const int grid_id = mfi.index();
                const int tile_id = mfi.LocalTileIndex();
                pc->GetParticles(lev)[std::make_pair(grid_id,tile_id)];
                if ( (pc->NumRuntimeRealComps()>0) || (pc->NumRuntimeIntComps()>0) ) {
                    pc->DefineAndReturnParticleTile(lev, grid_id, tile_id);
                }
            }
#endif
            MFItInfo info;
            if (pc->do_tiling && Gpu::notInLaunchRegion()) {
                info.EnableTiling(pc->tile_size);
            }
#ifdef _OPENMP
            info.SetDynamic(true);
#pragma omp parallel
#endif
            for (MFIter mfi = pc->MakeMFIter(lev, info); mfi.isValid(); ++mfi)
            {
                const int grid_id = mfi.index();
                const int tile_id = mfi.LocalTileIndex();
                auto& particle_tile = pc->GetParticles(lev)[std::make_pair(grid_id,tile_id)];
                auto& soa = particle_tile.GetStructOfArrays();
                const int np = particle_tile.GetArrayOfStructs().size();
                if (np == 0) break;
                amrex::Gpu::ManagedVector<int> is_ionized;
                // const int np_lev = pc->NumberOfParticlesAtLevel(lev);
                is_ionized.resize(np);
                int * const AMREX_RESTRICT p_is_ionized = is_ionized.dataPtr();
                const Real * const AMREX_RESTRICT ux = soa.GetRealData(PIdx::ux).data();
                const Real * const AMREX_RESTRICT uy = soa.GetRealData(PIdx::uy).data();
                const Real * const AMREX_RESTRICT uz = soa.GetRealData(PIdx::uz).data();
                const Real * const AMREX_RESTRICT ex = soa.GetRealData(PIdx::Ex).data();
                const Real * const AMREX_RESTRICT ey = soa.GetRealData(PIdx::Ey).data();
                const Real * const AMREX_RESTRICT ez = soa.GetRealData(PIdx::Ez).data();
                const Real * const AMREX_RESTRICT bx = soa.GetRealData(PIdx::Bx).data();
                const Real * const AMREX_RESTRICT by = soa.GetRealData(PIdx::By).data();
                const Real * const AMREX_RESTRICT bz = soa.GetRealData(PIdx::Bz).data();
                Real * const AMREX_RESTRICT ilev_real = soa.GetRealData(pc->particle_comps["ionization_level"]).data();
                Real c = PhysConst::c;
                Real c2_inv = 1./c/c;

                ParallelFor( 
                    np,
                    [=] AMREX_GPU_DEVICE (long i) {
                        Real random_draw = Random();
                        Real ga = std::sqrt(1. +
                                            (ux[i]*ux[i] + 
                                             uy[i]*uy[i] + 
                                             uz[i]*uz[i]) * c2_inv);
                        Real E = std::sqrt(
                            - ( ux[i]*ex[i] + uy[i]*ey[i] + uz[i]*ez[i] ) * ( ux[i]*ex[i] + uy[i]*ey[i] + uz[i]*ez[i] ) * c2_inv
                            + ( ga   *ex[i] + uy[i]*bz[i] - uz[i]*by[i] ) * ( ga   *ex[i] + uy[i]*bz[i] - uz[i]*by[i] )
                            + ( ga   *ey[i] + uz[i]*bx[i] - ux[i]*bz[i] ) * ( ga   *ey[i] + uz[i]*bx[i] - ux[i]*bz[i] )
                            + ( ga   *ez[i] + ux[i]*by[i] - uy[i]*bx[i] ) * ( ga   *ez[i] + ux[i]*by[i] - uy[i]*bx[i] )
                            );
                        int ilev = (int) round(ilev_real[i]);
                        // int ilev = static_cast<int>(round(ilev_real[i]));
                        Real p;
                        Real w_dtau;
                        Print()<<p_ionization_energies[0]<<" IE \n";
                        if (E<1.e-100*(p_ionization_energies[0])){
                            p = 0.;
                        } else {
                            w_dtau = 1./ ga * p_adk_prefactor[ilev] * 
                                std::pow(E,p_adk_power[ilev]) * 
                                std::exp( p_adk_exp_prefactor[ilev]/E );
                            p = 1. - std::exp( - w_dtau );
                        }
                        p_is_ionized[i] = 0;
                        if (random_draw < p){
                            ilev_real[i] += 1.;
                            p_is_ionized[i] = 1;
                        }
                    }
                ); // ParallelFor

                amrex::Gpu::ManagedVector<int> is_ionized_cumsum_vector;
                is_ionized_cumsum_vector.resize(np);
                int np_ionized = p_is_ionized[0];
                for(int i=1; i<np; ++i){
                    np_ionized += p_is_ionized[i];
                    is_ionized_cumsum_vector[i] = is_ionized_cumsum_vector[i-1]+p_is_ionized[i-1];
                }
                if (np_ionized == 0){
                    break;
                }

                auto& prod_pc = allcontainers[pc->ionization_product];
                auto& prod_particle_tile = prod_pc->GetParticles(lev)[std::make_pair(grid_id,tile_id)];
                const int np_old = prod_particle_tile.GetArrayOfStructs().size();
                const int np_new = np_old + np_ionized;
                prod_particle_tile.resize(np_new);
                auto& prod_soa = prod_particle_tile.GetStructOfArrays();

                GpuArray<Real*,PIdx::nattribs> prod_pa;
                for (int ia = 0; ia < PIdx::nattribs; ++ia) {
                    prod_pa[ia] = prod_soa.GetRealData(ia).data() + np_old;
                }
                WarpXParticleContainer::ParticleType* prod_pp = prod_particle_tile.GetArrayOfStructs()().data() + np_old;
                WarpXParticleContainer::ParticleType* pp = particle_tile.GetArrayOfStructs()().data();

                GpuArray<Real*,PIdx::nattribs> pa;
                for (int ia = 0; ia < PIdx::nattribs; ++ia) {
                    pa[ia] = soa.GetRealData(ia).data();
                }

                const int cpuid = ParallelDescriptor::MyProc();
                int* AMREX_RESTRICT is_ionized_cumsum = is_ionized_cumsum_vector.dataPtr();
                amrex::For(
                    np, [=] AMREX_GPU_DEVICE (int ip) noexcept
                    {
                        if(is_ionized[ip]){
                            int i = is_ionized_cumsum[ip];
                            WarpXParticleContainer::ParticleType& prod_p = prod_pp[i];
                            WarpXParticleContainer::ParticleType& p = pp[ip];
                            prod_p.id() = 9999;
                            prod_p.cpu() = cpuid;
                            prod_p.pos(0) = p.pos(0);
                            prod_p.pos(1) = p.pos(1);
#if (AMREX_SPACEDIM == 3)
                            prod_p.pos(20 = p.pos(2);
#endif
                            for (int ia = 0; ia < PIdx::nattribs; ++ia) {
                                prod_pa[ia][i] = pa[ia][ip];
                            }
                        }
                    }
                );
            } // MFIter
        }
    }
}
