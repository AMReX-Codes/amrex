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
MultiParticleContainer::Evolve (int lev,
                                const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                                MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                MultiFab* cjx,  MultiFab* cjy, MultiFab* cjz, 
                                MultiFab* rho,
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
    for (auto& pc : allcontainers) {
	pc->Evolve(lev, Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, cjx, cjy, cjz,
               rho, cEx, cEy, cEz, cBx, cBy, cBz, t, dt);
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
MultiParticleContainer::FieldGatherFortran (int lev,
                                            const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                            const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    for (auto& pc : allcontainers) {
        pc->FieldGatherFortran(lev, Ex, Ey, Ez, Bx, By, Bz);
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
MultiParticleContainer::NumberOfParticlesInGrid(int lev) const
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
