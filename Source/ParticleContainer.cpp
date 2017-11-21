#include <limits>

#include <ParticleContainer.H>
#include <WarpX_f.H>
#include <WarpX.H>

using namespace amrex;

MultiParticleContainer::MultiParticleContainer (AmrCore* amr_core)
{
    ReadParameters();

    int n = WarpX::use_laser ? nspecies+1 : nspecies;
    allcontainers.resize(n);
    for (int i = 0; i < nspecies; ++i) {
	allcontainers[i].reset(new PhysicalParticleContainer(amr_core, i, species_names[i]));
    }
    if (WarpX::use_laser) {
	allcontainers[n-1].reset(new LaserParticleContainer(amr_core,n-1));
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
}

void
MultiParticleContainer::InitData ()
{
    for (auto& pc : allcontainers) {
	pc->InitData();
    }
}

void
MultiParticleContainer::FieldGatherES (const Vector<std::array<std::unique_ptr<MultiFab>, 3> >& E,
                                       const amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks)
{
    for (auto& pc : allcontainers) {
	pc->FieldGatherES(E, masks);
    }
}

void
MultiParticleContainer::FieldGather (int lev,
                                     const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                     const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    for (auto& pc : allcontainers) {
	pc->FieldGather(lev, Ex, Ey, Ez, Bx, By, Bz);
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
MultiParticleContainer::PushX (Real dt)
{
    for (auto& pc : allcontainers) {
	pc->PushX(dt);
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

amrex::Real
MultiParticleContainer::sumParticleCharge (bool local)
{
    amrex::Real total_charge = allcontainers[0]->sumParticleCharge(local);
    for (unsigned i = 1, n = allcontainers.size(); i < n; ++i) {
        total_charge += allcontainers[i]->sumParticleCharge(local);
    }
    return total_charge;
}

void
MultiParticleContainer::Redistribute ()
{
    for (auto& pc : allcontainers) {
	pc->Redistribute();
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
}
