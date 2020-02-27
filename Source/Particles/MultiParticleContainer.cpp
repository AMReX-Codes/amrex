/* Copyright 2019-2020 Andrew Myers, Ann Almgren, Axel Huebl
 * David Grote, Jean-Luc Vay, Luca Fedeli
 * Mathieu Lobet, Maxence Thevenet, Remi Lehe
 * Revathi Jambunathan, Weiqun Zhang, Yinjian Zhao
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "MultiParticleContainer.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

#include <AMReX_Vector.H>

#include <limits>
#include <algorithm>
#include <string>


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
        else if (species_types[i] == PCTypes::Photon) {
            allcontainers[i].reset(new PhotonParticleContainer(amr_core, i, species_names[i]));
        }
        allcontainers[i]->m_deposit_on_main_grid = m_deposit_on_main_grid[i];
        allcontainers[i]->m_gather_from_main_grid = m_gather_from_main_grid[i];
    }

    for (int i = nspecies; i < nspecies+nlasers; ++i) {
        allcontainers[i].reset(new LaserParticleContainer(amr_core, i, lasers_names[i-nspecies]));
    }

    pc_tmp.reset(new PhysicalParticleContainer(amr_core));

    // Compute the number of species for which lab-frame data is dumped
    // nspecies_lab_frame_diags, and map their ID to MultiParticleContainer
    // particle IDs in map_species_lab_diags.
    map_species_back_transformed_diagnostics.resize(nspecies);
    nspecies_back_transformed_diagnostics = 0;
    for (int i=0; i<nspecies; i++){
        auto& pc = allcontainers[i];
        if (pc->do_back_transformed_diagnostics){
            map_species_back_transformed_diagnostics[nspecies_back_transformed_diagnostics] = i;
            do_back_transformed_diagnostics = 1;
            nspecies_back_transformed_diagnostics += 1;
        }
    }

    // collision
    allcollisions.resize(ncollisions);
    for (int i = 0; i < ncollisions; ++i) {
        allcollisions[i].reset
            (new CollisionType(species_names, collision_names[i]));
    }

}

void
MultiParticleContainer::ReadParameters ()
{
    static bool initialized = false;
    if (!initialized)
    {
        ParmParse pp("particles");

        // allocating and initializing default values of external fields for particles
        m_E_external_particle.resize(3);
        m_B_external_particle.resize(3);
        // initialize E and B fields to 0.0
        for (int idim = 0; idim < 3; ++idim) {
            m_E_external_particle[idim] = 0.0;
            m_B_external_particle[idim] = 0.0;
        }
        // default values of E_external_particle and B_external_particle
        // are used to set the E and B field when "constant" or "parser"
        // is not explicitly used in the input
        pp.query("B_ext_particle_init_style", m_B_ext_particle_s);
        std::transform(m_B_ext_particle_s.begin(),
                       m_B_ext_particle_s.end(),
                       m_B_ext_particle_s.begin(),
                       ::tolower);
        pp.query("E_ext_particle_init_style", m_E_ext_particle_s);
        std::transform(m_E_ext_particle_s.begin(),
                       m_E_ext_particle_s.end(),
                       m_E_ext_particle_s.begin(),
                       ::tolower);
        // if the input string for B_external on particles is "constant"
        // then the values for the external B on particles must
        // be provided in the input file.
        if (m_B_ext_particle_s == "constant")
            pp.getarr("B_external_particle", m_B_external_particle);

        // if the input string for E_external on particles is "constant"
        // then the values for the external E on particles must
        // be provided in the input file.
        if (m_E_ext_particle_s == "constant")
            pp.getarr("E_external_particle", m_E_external_particle);

        // if the input string for B_ext_particle_s is
        // "parse_b_ext_particle_function" then the mathematical expression
        // for the Bx_, By_, Bz_external_particle_function(x,y,z)
        // must be provided in the input file.
        if (m_B_ext_particle_s == "parse_b_ext_particle_function") {
           // store the mathematical expression as string
           std::string str_Bx_ext_particle_function;
           std::string str_By_ext_particle_function;
           std::string str_Bz_ext_particle_function;
           Store_parserString(pp, "Bx_external_particle_function(x,y,z,t)",
                                      str_Bx_ext_particle_function);
           Store_parserString(pp, "By_external_particle_function(x,y,z,t)",
                                      str_By_ext_particle_function);
           Store_parserString(pp, "Bz_external_particle_function(x,y,z,t)",
                                      str_Bz_ext_particle_function);

           // Parser for B_external on the particle
           m_Bx_particle_parser.reset(new ParserWrapper<4>(
                                    makeParser(str_Bx_ext_particle_function,{"x","y","z","t"})));
           m_By_particle_parser.reset(new ParserWrapper<4>(
                                    makeParser(str_By_ext_particle_function,{"x","y","z","t"})));
           m_Bz_particle_parser.reset(new ParserWrapper<4>(
                                    makeParser(str_Bz_ext_particle_function,{"x","y","z","t"})));

        }

        // if the input string for E_ext_particle_s is
        // "parse_e_ext_particle_function" then the mathematical expression
        // for the Ex_, Ey_, Ez_external_particle_function(x,y,z)
        // must be provided in the input file.
        if (m_E_ext_particle_s == "parse_e_ext_particle_function") {
           // store the mathematical expression as string
           std::string str_Ex_ext_particle_function;
           std::string str_Ey_ext_particle_function;
           std::string str_Ez_ext_particle_function;
           Store_parserString(pp, "Ex_external_particle_function(x,y,z,t)",
                                      str_Ex_ext_particle_function);
           Store_parserString(pp, "Ey_external_particle_function(x,y,z,t)",
                                      str_Ey_ext_particle_function);
           Store_parserString(pp, "Ez_external_particle_function(x,y,z,t)",
                                      str_Ez_ext_particle_function);
           // Parser for E_external on the particle
           m_Ex_particle_parser.reset(new ParserWrapper<4>(
                                    makeParser(str_Ex_ext_particle_function,{"x","y","z","t"})));
           m_Ey_particle_parser.reset(new ParserWrapper<4>(
                                    makeParser(str_Ey_ext_particle_function,{"x","y","z","t"})));
           m_Ez_particle_parser.reset(new ParserWrapper<4>(
                                    makeParser(str_Ez_ext_particle_function,{"x","y","z","t"})));

        }




        pp.query("nspecies", nspecies);
        BL_ASSERT(nspecies >= 0);

        if (nspecies > 0) {
            // Get species names
            pp.getarr("species_names", species_names);
            BL_ASSERT(species_names.size() == nspecies);

            // Get species to deposit on main grid
            m_deposit_on_main_grid.resize(nspecies, false);
            std::vector<std::string> tmp;
            pp.queryarr("deposit_on_main_grid", tmp);
            for (auto const& name : tmp) {
                auto it = std::find(species_names.begin(), species_names.end(), name);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != species_names.end(), "ERROR: species in particles.deposit_on_main_grid must be part of particles.species_names");
                int i = std::distance(species_names.begin(), it);
                m_deposit_on_main_grid[i] = true;
            }

            m_gather_from_main_grid.resize(nspecies, false);
            std::vector<std::string> tmp_gather;
            pp.queryarr("gather_from_main_grid", tmp_gather);
            for (auto const& name : tmp_gather) {
                auto it = std::find(species_names.begin(), species_names.end(), name);
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != species_names.end(), "ERROR: species in particles.gather_from_main_grid must be part of particles.species_names");
                int i = std::distance(species_names.begin(), it);
                m_gather_from_main_grid.at(i) = true;
            }

            species_types.resize(nspecies, PCTypes::Physical);

            // Get rigid-injected species
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
            // Get photon species
            std::vector<std::string> photon_species;
            pp.queryarr("photon_species", photon_species);
            if (!photon_species.empty()) {
                for (auto const& name : photon_species) {
                    auto it = std::find(species_names.begin(), species_names.end(), name);
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                        it != species_names.end(),
                        "ERROR: species in particles.rigid_injected_species must be part of particles.species_names");
                    int i = std::distance(species_names.begin(), it);
                    species_types[i] = PCTypes::Photon;
                }
            }

            // collision
            ParmParse pc("collisions");
            pc.query("ncollisions", ncollisions);
            BL_ASSERT(ncollisions >= 0);
            if (ncollisions > 0) {
                pc.getarr("collision_names", collision_names);
                BL_ASSERT(collision_names.size() == ncollisions);
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

#ifdef WARPX_QED
    InitQED();
#endif

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
                                Real t, Real dt, DtType a_dt_type)
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
                   rho, crho, cEx, cEy, cEz, cBx, cBy, cBz, t, dt, a_dt_type);
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
        MultiFab::Add(*rho, *rhoi, 0, 0, rho->nComp(), rho->nGrow());
    }
    if (!local) {
        const Geometry& gm = allcontainers[0]->Geom(lev);
        rho->SumBoundary(gm.periodicity());
    }
    return rho;
}

void
MultiParticleContainer::SortParticlesByBin (amrex::IntVect bin_size)
{
    for (auto& pc : allcontainers) {
        pc->SortParticlesByBin(bin_size);
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
::GetLabFrameData (const std::string& /*snapshot_name*/,
                   const int /*i_lab*/, const int direction,
                   const Real z_old, const Real z_new,
                   const Real t_boost, const Real t_lab, const Real dt,
                   Vector<WarpXParticleContainer::DiagnosticParticleData>& parts) const
{

    WARPX_PROFILE("MultiParticleContainer::GetLabFrameData");

    // Loop over particle species
    for (int i = 0; i < nspecies_back_transformed_diagnostics; ++i){
        int isp = map_species_back_transformed_diagnostics[i];
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
MultiParticleContainer::ContinuousInjection (const RealBox& injection_box) const
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
MultiParticleContainer::UpdateContinuousInjectionPosition (Real dt) const
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

void
MultiParticleContainer::doFieldIonization ()
{
    WARPX_PROFILE("MPC::doFieldIonization");

    // Loop over all species.
    // Ionized particles in pc_source create particles in pc_product
    for (auto& pc_source : allcontainers)
    {
        if (!pc_source->do_field_ionization){ continue; }

        auto& pc_product = allcontainers[pc_source->ionization_product];

        SmartCopyFactory copy_factory(*pc_source, *pc_product);
        auto phys_pc_ptr = static_cast<PhysicalParticleContainer*>(pc_source.get());

        auto Filter    = phys_pc_ptr->getIonizationFunc();
        auto Copy      = copy_factory.getSmartCopy();
        auto Transform = IonizationTransformFunc();

        pc_source ->defineAllParticleTiles();
        pc_product->defineAllParticleTiles();

        for (int lev = 0; lev <= pc_source->finestLevel(); ++lev)
        {
            auto info = getMFItInfo(*pc_source, *pc_product);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi = pc_source->MakeMFIter(lev, info); mfi.isValid(); ++mfi)
            {
                auto& src_tile = pc_source ->ParticlesAt(lev, mfi);
                auto& dst_tile = pc_product->ParticlesAt(lev, mfi);

                auto np_dst = dst_tile.numParticles();
                auto num_added = filterCopyTransformParticles<1>(dst_tile, src_tile, np_dst,
                                                                 Filter, Copy, Transform);

                setNewParticleIDs(dst_tile, np_dst, num_added);
            }
        }
    }
}

void
MultiParticleContainer::doCoulombCollisions ()
{
    WARPX_PROFILE("MPC::doCoulombCollisions");

    for (int i = 0; i < ncollisions; ++i)
    {
        auto& species1 = allcontainers[ allcollisions[i]->m_species1_index ];
        auto& species2 = allcontainers[ allcollisions[i]->m_species2_index ];

        // Enable tiling
        MFItInfo info;
        if (Gpu::notInLaunchRegion()) info.EnableTiling(species1->tile_size);

        // Loop over refinement levels
        for (int lev = 0; lev <= species1->finestLevel(); ++lev){

            // Loop over all grids/tiles at this level
#ifdef _OPENMP
            info.SetDynamic(true);
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi = species1->MakeMFIter(lev, info); mfi.isValid(); ++mfi){

                CollisionType::doCoulombCollisionsWithinTile
                    ( lev, mfi, species1, species2,
                      allcollisions[i]->m_isSameSpecies,
                      allcollisions[i]->m_CoulombLog );

            }
        }
    }
}

MFItInfo MultiParticleContainer::getMFItInfo (const WarpXParticleContainer& pc_src,
                                              const WarpXParticleContainer& pc_dst) const noexcept
{
    MFItInfo info;

    if (pc_src.do_tiling && Gpu::notInLaunchRegion()) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pc_dst.do_tiling,
                                         "For ionization, either all or none of the "
                                         "particle species must use tiling.");
        info.EnableTiling(pc_src.tile_size);
    }

#ifdef _OPENMP
    info.SetDynamic(true);
#endif

    return info;
}

#ifdef WARPX_QED
void MultiParticleContainer::InitQED ()
{
    m_shr_p_qs_engine = std::make_shared<QuantumSynchrotronEngine>();
    m_shr_p_bw_engine = std::make_shared<BreitWheelerEngine>();

    m_nspecies_quantum_sync = 0;
    m_nspecies_breit_wheeler = 0;

    for (auto& pc : allcontainers) {
        if(pc->has_quantum_sync()){
            pc->set_quantum_sync_engine_ptr
                (m_shr_p_qs_engine);
            m_nspecies_quantum_sync++;
        }
        if(pc->has_breit_wheeler()){
            pc->set_breit_wheeler_engine_ptr
                (m_shr_p_bw_engine);
            m_nspecies_breit_wheeler++;
        }
    }

    if(m_nspecies_quantum_sync != 0)
        InitQuantumSync();

    if(m_nspecies_breit_wheeler !=0)
        InitBreitWheeler();

}

void MultiParticleContainer::InitQuantumSync ()
{
    std::string lookup_table_mode;
    ParmParse pp("qed_qs");
    pp.query("lookup_table_mode", lookup_table_mode);
    if(lookup_table_mode.empty()){
        amrex::Abort("Quantum Synchrotron table mode should be provided");
    }

    if(lookup_table_mode == "generate"){
        amrex::Print() << "Quantum Synchrotron table will be generated. \n" ;
#ifndef WARPX_QED_TABLE_GEN
        amrex::Error("Error: Compile with QED_TABLE_GEN=TRUE to enable table generation!\n");
#else
        QuantumSyncGenerateTable();
#endif
    }
    else if(lookup_table_mode == "load"){
        amrex::Print() << "Quantum Synchrotron table will be read from file. \n" ;
        std::string load_table_name;
        pp.query("load_table_from", load_table_name);
        if(load_table_name.empty()){
            amrex::Abort("Quantum Synchrotron table name should be provided");
        }
        Vector<char> table_data;
        ParallelDescriptor::ReadAndBcastFile(load_table_name, table_data);
        ParallelDescriptor::Barrier();
        m_shr_p_qs_engine->init_lookup_tables_from_raw_data(table_data);
    }
    else if(lookup_table_mode == "dummy_builtin"){
        amrex::Print() << "Built-in Quantum Synchrotron dummy table will be used. \n" ;
        m_shr_p_qs_engine->init_dummy_tables();
    }
    else{
        amrex::Abort("Unknown Quantum Synchrotron table mode");
    }

    if(!m_shr_p_qs_engine->are_lookup_tables_initialized()){
        amrex::Abort("Table initialization has failed!");
    }
}

void MultiParticleContainer::InitBreitWheeler ()
{
    std::string lookup_table_mode;
    ParmParse pp("qed_bw");
    pp.query("lookup_table_mode", lookup_table_mode);
    if(lookup_table_mode.empty()){
        amrex::Abort("Breit Wheeler table mode should be provided");
    }

    if(lookup_table_mode == "generate"){
        amrex::Print() << "Breit Wheeler table will be generated. \n" ;
#ifndef WARPX_QED_TABLE_GEN
        amrex::Error("Error: Compile with QED_TABLE_GEN=TRUE to enable table generation!\n");
#else
        BreitWheelerGenerateTable();
#endif
    }
    else if(lookup_table_mode == "load"){
        amrex::Print() << "Breit Wheeler table will be read from file. \n" ;
        std::string load_table_name;
        pp.query("load_table_from", load_table_name);
        if(load_table_name.empty()){
            amrex::Abort("Breit Wheeler table name should be provided");
        }
        Vector<char> table_data;
        ParallelDescriptor::ReadAndBcastFile(load_table_name, table_data);
        ParallelDescriptor::Barrier();
        m_shr_p_bw_engine->init_lookup_tables_from_raw_data(table_data);
    }
    else if(lookup_table_mode == "dummy_builtin"){
        amrex::Print() << "Built-in Breit Wheeler dummy table will be used. \n" ;
        m_shr_p_bw_engine->init_dummy_tables();
    }
    else{
        amrex::Abort("Unknown Breit Wheeler table mode");
    }

    if(!m_shr_p_bw_engine->are_lookup_tables_initialized()){
        amrex::Abort("Table initialization has failed!");
    }
}

void
MultiParticleContainer::QuantumSyncGenerateTable ()
{
    ParmParse pp("qed_qs");
    std::string table_name;
    pp.query("save_table_in", table_name);
    if(table_name.empty())
        amrex::Abort("qed_qs.save_table_in should be provided!");

    if(ParallelDescriptor::IOProcessor()){
        PicsarQuantumSynchrotronCtrl ctrl;
        int t_int;

        // Engine paramenter: chi_part_min is the minium chi parameter to be
        // considered by the engine. If a lepton has chi < chi_part_min,
        // the optical depth is not evolved and photon generation is ignored
        if(!pp.query("chi_min", ctrl.chi_part_min))
            amrex::Abort("qed_qs.chi_min should be provided!");

        //==Table parameters==

        //--- sub-table 1 (1D)
        //These parameters are used to pre-compute a function
        //which appears in the evolution of the optical depth

        //Minimun chi for the table. If a lepton has chi < chi_part_tdndt_min,
        //chi is considered as it were equal to chi_part_tdndt_min
        if(!pp.query("tab_dndt_chi_min", ctrl.chi_part_tdndt_min))
            amrex::Abort("qed_qs.tab_dndt_chi_min should be provided!");

        //Maximum chi for the table. If a lepton has chi > chi_part_tdndt_max,
        //chi is considered as it were equal to chi_part_tdndt_max
        if(!pp.query("tab_dndt_chi_max", ctrl.chi_part_tdndt_max))
            amrex::Abort("qed_qs.tab_dndt_chi_max should be provided!");

        //How many points should be used for chi in the table
        if(!pp.query("tab_dndt_how_many", t_int))
            amrex::Abort("qed_qs.tab_dndt_how_many should be provided!");
        ctrl.chi_part_tdndt_how_many = t_int;
        //------

        //--- sub-table 2 (2D)
        //These parameters are used to pre-compute a function
        //which is used to extract the properties of the generated
        //photons.

        //Minimun chi for the table. If a lepton has chi < chi_part_tem_min,
        //chi is considered as it were equal to chi_part_tem_min
        if(!pp.query("tab_em_chi_min", ctrl.chi_part_tem_min))
            amrex::Abort("qed_qs.tab_em_chi_min should be provided!");

        //Maximum chi for the table. If a lepton has chi > chi_part_tem_max,
        //chi is considered as it were equal to chi_part_tem_max
        if(!pp.query("tab_em_chi_max", ctrl.chi_part_tem_max))
            amrex::Abort("qed_qs.tab_em_chi_max should be provided!");

        //How many points should be used for chi in the table
        if(!pp.query("tab_em_chi_how_many", t_int))
            amrex::Abort("qed_qs.tab_em_chi_how_many should be provided!");
        ctrl.chi_part_tem_how_many = t_int;

        //The other axis of the table is a cumulative probability distribution
        //(corresponding to different energies of the generated particles)
        //This parameter is the number of different points to consider
        if(!pp.query("tab_em_prob_how_many", t_int))
            amrex::Abort("qed_qs.tab_em_prob_how_many should be provided!");
        ctrl.prob_tem_how_many = t_int;
        //====================

        m_shr_p_qs_engine->compute_lookup_tables(ctrl);
        WarpXUtilIO::WriteBinaryDataOnFile(table_name,
            m_shr_p_qs_engine->export_lookup_tables_data());
    }

    ParallelDescriptor::Barrier();
    Vector<char> table_data;
    ParallelDescriptor::ReadAndBcastFile(table_name, table_data);
    ParallelDescriptor::Barrier();

    //No need to initialize from raw data for the processor that
    //has just generated the table
    if(!ParallelDescriptor::IOProcessor()){
        m_shr_p_qs_engine->init_lookup_tables_from_raw_data(table_data);
    }
}

void
MultiParticleContainer::BreitWheelerGenerateTable ()
{
    ParmParse pp("qed_bw");
    std::string table_name;
    pp.query("save_table_in", table_name);
    if(table_name.empty())
        amrex::Abort("qed_bw.save_table_in should be provided!");

    if(ParallelDescriptor::IOProcessor()){
        PicsarBreitWheelerCtrl ctrl;
        int t_int;

        // Engine paramenter: chi_phot_min is the minium chi parameter to be
        // considered by the engine. If a photon has chi < chi_phot_min,
        // the optical depth is not evolved and pair generation is ignored
        if(!pp.query("chi_min", ctrl.chi_phot_min))
            amrex::Abort("qed_bw.chi_min should be provided!");

        //==Table parameters==

        //--- sub-table 1 (1D)
        //These parameters are used to pre-compute a function
        //which appears in the evolution of the optical depth

        //Minimun chi for the table. If a photon has chi < chi_phot_tdndt_min,
        //an analytical approximation is used.
        if(!pp.query("tab_dndt_chi_min", ctrl.chi_phot_tdndt_min))
            amrex::Abort("qed_bw.tab_dndt_chi_min should be provided!");

        //Maximum chi for the table. If a photon has chi > chi_phot_tdndt_min,
        //an analytical approximation is used.
        if(!pp.query("tab_dndt_chi_max", ctrl.chi_phot_tdndt_max))
            amrex::Abort("qed_bw.tab_dndt_chi_max should be provided!");

        //How many points should be used for chi in the table
        if(!pp.query("tab_dndt_how_many", t_int))
            amrex::Abort("qed_bw.tab_dndt_how_many should be provided!");
        ctrl.chi_phot_tdndt_how_many = t_int;
        //------

        //--- sub-table 2 (2D)
        //These parameters are used to pre-compute a function
        //which is used to extract the properties of the generated
        //particles.

        //Minimun chi for the table. If a photon has chi < chi_phot_tpair_min
        //chi is considered as it were equal to chi_phot_tpair_min
        if(!pp.query("tab_pair_chi_min", ctrl.chi_phot_tpair_min))
            amrex::Abort("qed_bw.tab_pair_chi_min should be provided!");

        //Maximum chi for the table. If a photon has chi > chi_phot_tpair_max
        //chi is considered as it were equal to chi_phot_tpair_max
        if(!pp.query("tab_pair_chi_max", ctrl.chi_phot_tpair_max))
            amrex::Abort("qed_bw.tab_pair_chi_max should be provided!");

        //How many points should be used for chi in the table
        if(!pp.query("tab_pair_chi_how_many", t_int))
            amrex::Abort("qed_bw.tab_pair_chi_how_many should be provided!");
        ctrl.chi_phot_tpair_how_many = t_int;

        //The other axis of the table is the fraction of the initial energy
        //'taken away' by the most energetic particle of the pair.
        //This parameter is the number of different fractions to consider
        if(!pp.query("tab_pair_frac_how_many", t_int))
            amrex::Abort("qed_bw.tab_pair_frac_how_many should be provided!");
        ctrl.chi_frac_tpair_how_many = t_int;
        //====================

        m_shr_p_bw_engine->compute_lookup_tables(ctrl);
        WarpXUtilIO::WriteBinaryDataOnFile(table_name,
            m_shr_p_bw_engine->export_lookup_tables_data());
    }

    ParallelDescriptor::Barrier();
    Vector<char> table_data;
    ParallelDescriptor::ReadAndBcastFile(table_name, table_data);
    ParallelDescriptor::Barrier();

    //No need to initialize from raw data for the processor that
    //has just generated the table
    if(!ParallelDescriptor::IOProcessor()){
        m_shr_p_bw_engine->init_lookup_tables_from_raw_data(table_data);
    }
}
#endif
