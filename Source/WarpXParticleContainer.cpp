
#include <limits>

#include <ParticleContainer.H>
#include <WarpXParticleContainer.H>
#include <AMReX_AmrParGDB.H>
#include <WarpX_f.H>
#include <WarpX.H>

using namespace amrex;

#if (BL_SPACEDIM == 2)
void
WarpXParIter::GetPosition (Array<Real>& x, Array<Real>& y, Array<Real>& z) const
{
    amrex::ParIter<0,0,PIdx::nattribs>::GetPosition(x, z);
    y.resize(x.size(), std::numeric_limits<Real>::quiet_NaN());
}

void
WarpXParIter::SetPosition (const Array<Real>& x, const Array<Real>& y, const Array<Real>& z)
{
    amrex::ParIter<0,0,PIdx::nattribs>::SetPosition(x, z);
}
#endif

WarpXParticleContainer::WarpXParticleContainer (AmrCore* amr_core, int ispecies)
    : ParticleContainer<0,0,PIdx::nattribs>(amr_core->GetParGDB())
    , species_id(ispecies)
{
    for (unsigned int i = PIdx::Ex; i < PIdx::nattribs; ++i) {
        communicate_real_comp[i] = false; // Don't need to communicate E and B.
    }
    ReadParameters();
}

void
WarpXParticleContainer::ReadParameters ()
{
    static bool initialized = false;
    if (!initialized)
    {
	ParmParse pp("particles");

        do_tiling = true;  // because the default in amrex is false
	pp.query("do_tiling",  do_tiling);

	initialized = true;
    }
}

void
WarpXParticleContainer::AllocData ()
{
    // have to resize here, not in the constructor because grids have not
    // been built when constructor was called.
    reserveData();
    resizeData();
}

void
WarpXParticleContainer::AddOneParticle (int lev, int grid, int tile,
                                        Real x, Real y, Real z,
                                        const std::array<Real,PIdx::nattribs>& attribs)
{
    auto& particle_tile = GetParticles(lev)[std::make_pair(grid,tile)];
    AddOneParticle(particle_tile, x, y, z, attribs); 
}

void
WarpXParticleContainer::AddOneParticle (ParticleTileType& particle_tile,
                                        Real x, Real y, Real z,
                                        const std::array<Real,PIdx::nattribs>& attribs)
{
    ParticleType p;
    p.id()  = ParticleType::NextID();
    p.cpu() = ParallelDescriptor::MyProc();
#if (BL_SPACEDIM == 3)
    p.pos(0) = x;
    p.pos(1) = y;
    p.pos(2) = z;
#elif (BL_SPACEDIM == 2)
    p.pos(0) = x;
    p.pos(1) = z;
#endif
    
    particle_tile.push_back(p);
    particle_tile.push_back_real(attribs);
}

void
WarpXParticleContainer::AddNParticles (int n, const Real* x, const Real* y, const Real* z,
				       const Real* vx, const Real* vy, const Real* vz,
				       int nattr, const Real* attr, int uniqueparticles)
{
    const int lev = 0;

    BL_ASSERT(nattr == 1);
    const Real* weight = attr;

    auto npart_before = TotalNumberOfParticles();  // xxxxx move this into if (verbose > xxx)

    int ibegin, iend;
    if (uniqueparticles) {
	ibegin = 0;
	iend = n;
    } else {
	int myproc = ParallelDescriptor::MyProc();
	int nprocs = ParallelDescriptor::NProcs();
	int navg = n/nprocs;
	int nleft = n - navg * nprocs;
	if (myproc < nleft) {
	    ibegin = myproc*(navg+1);
	    iend = ibegin + navg+1;
	} else {
	    ibegin = myproc*navg + nleft;
	    iend = ibegin + navg;
	}
    }

    //  Add to grid 0 and tile 0
    // Redistribute() will move them to proper places.
    std::pair<int,int> key {0,0};
    auto& particle_tile = GetParticles(lev)[key];

    for (int i = ibegin; i < iend; ++i)
    {
        ParticleType p;
        p.id()  = ParticleType::NextID();
        p.cpu() = ParallelDescriptor::MyProc();
#if (BL_SPACEDIM == 3)
        p.pos(0) = x[i];
        p.pos(1) = y[i];
        p.pos(2) = z[i];
#elif (BL_SPACEDIM == 2)
        p.pos(0) = x[i];
        p.pos(1) = z[i];
#endif
        particle_tile.push_back(p);
    }

    particle_tile.push_back_real(PIdx::w , weight + ibegin, weight + iend);
    particle_tile.push_back_real(PIdx::ux,     vx + ibegin,     vx + iend);
    particle_tile.push_back_real(PIdx::uy,     vy + ibegin,     vy + iend);
    particle_tile.push_back_real(PIdx::uz,     vz + ibegin,     vz + iend);

    std::size_t np = iend-ibegin;
    for (int comp = PIdx::uz+1; comp < PIdx::nattribs; ++comp)
    {
        particle_tile.push_back_real(comp, np, 0.0);
    }

    Redistribute();

    auto npart_after = TotalNumberOfParticles();  // xxxxx move this into if (verbose > xxx)
    amrex::Print() << "Total number of particles added: " << npart_after - npart_before << "\n";
}

void
WarpXParticleContainer::DepositCharge (Array<std::unique_ptr<MultiFab> >& rho, bool local)
{

    const int lev = 0;

    const auto& gm = m_gdb->Geom(lev);
    const auto& ba = m_gdb->ParticleBoxArray(lev);
    const auto& dm = m_gdb->DistributionMap(lev);
    BoxArray nba = ba;
    nba.surroundingNodes();

    const Real* dx  = gm.CellSize();
    const Real* plo = gm.ProbLo();
    const int ng = 1;

    for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const Box& box = nba[pti];

        auto& wp = pti.GetAttribs(PIdx::w);
        const auto& particles = pti.GetArrayOfStructs();
        int nstride = particles.dataShape().first;           
        const long np  = pti.numParticles();

	// Data on the grid
	FArrayBox& rhofab = (*rho[0])[pti];

        warpx_deposit_cic(particles.data(), nstride, np,
                          wp.data(), &this->charge, 
                          rhofab.dataPtr(), box.loVect(), box.hiVect(), plo, dx);
    }

    if (!local) rho[0]->SumBoundary(gm.periodicity());
}

std::unique_ptr<MultiFab>
WarpXParticleContainer::GetChargeDensity (int lev, bool local)
{
    const auto& gm = m_gdb->Geom(lev);
    const auto& ba = m_gdb->ParticleBoxArray(lev);
    const auto& dm = m_gdb->DistributionMap(lev);
    BoxArray nba = ba;
    nba.surroundingNodes();

    const std::array<Real,3>& dx = WarpX::CellSize(lev);

    const int ng = WarpX::nox;

    auto rho = std::unique_ptr<MultiFab>(new MultiFab(nba,dm,1,ng));
    rho->setVal(0.0);

    Array<Real> xp, yp, zp;

    for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const Box& box = pti.validbox();

        auto& wp = pti.GetAttribs(PIdx::w);
            
        const long np  = pti.numParticles();

	// Data on the grid
	FArrayBox& rhofab = (*rho)[pti];

        pti.GetPosition(xp, yp, zp);

#if (BL_SPACEDIM == 3)
	long nx = box.length(0);
	long ny = box.length(1);
	long nz = box.length(2); 
#elif (BL_SPACEDIM == 2)
	long nx = box.length(0);
	long ny = 0;
	long nz = box.length(1); 
#endif

        const std::array<Real,3>& xyzmin = WarpX::LowerCorner(box, lev);

	long nxg = ng;
	long nyg = ng;
	long nzg = ng;
	long lvect = 8;

	warpx_charge_deposition(rhofab.dataPtr(), 
				&np, xp.data(), yp.data(), zp.data(), wp.data(),
				&this->charge, &xyzmin[0], &xyzmin[1], &xyzmin[2], 
				&dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				&nxg, &nyg, &nzg, &WarpX::nox,&WarpX::noy,&WarpX::noz,
				&lvect, &WarpX::charge_deposition_algo);
				
    }

    if (!local) rho->SumBoundary(gm.periodicity());
    
    return rho;
}

Real WarpXParticleContainer::sumParticleCharge(int lev, bool local) {

    amrex::Real total_charge = 0.0;
    for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& wp = pti.GetAttribs(PIdx::w);
        for (unsigned long i = 0; i < wp.size(); i++) {
            total_charge += wp[i];
        }
    }

    if (!local) ParallelDescriptor::ReduceRealSum(total_charge);
    total_charge *= this->charge;
    return total_charge;
}

void
WarpXParticleContainer::PushXES (Real dt)
{
    const int lev = 0;
    BL_PROFILE("WPC::PushXES()");

    for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        int nstride = particles.dataShape().first;           
        const long np  = pti.numParticles();
        
        auto& attribs = pti.GetAttribs();        
        auto& uxp = attribs[PIdx::ux];
        auto& uyp = attribs[PIdx::uy];
        auto& uzp = attribs[PIdx::uz];
        
        warpx_push_leapfrog_positions(particles.data(), nstride, np,
                                      uxp.data(), uyp.data(), uzp.data(), &dt);
    }
}

void
WarpXParticleContainer::PushX (int lev,
                               Real dt)
{
    BL_PROFILE("WPC::PushX()");
    BL_PROFILE_VAR_NS("WPC::PushX::Copy", blp_copy);
    BL_PROFILE_VAR_NS("WPC:PushX::Push", blp_pxr_pp);

    Array<Real> xp, yp, zp, giv;

    for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
    {

        auto& attribs = pti.GetAttribs();

        auto& uxp = attribs[PIdx::ux];
        auto& uyp = attribs[PIdx::uy];
        auto& uzp = attribs[PIdx::uz];
        
        const long np = pti.numParticles();

        giv.resize(np);

        //
        // copy data from particle container to temp arrays
        //
        BL_PROFILE_VAR_START(blp_copy);
        pti.GetPosition(xp, yp, zp);
        BL_PROFILE_VAR_STOP(blp_copy);

        //
        // Particle Push
        //
        BL_PROFILE_VAR_START(blp_pxr_pp);
        warpx_particle_pusher_positions(&np, xp.data(), yp.data(), zp.data(),
                                        uxp.data(), uyp.data(), uzp.data(), giv.data(), &dt);
        BL_PROFILE_VAR_STOP(blp_pxr_pp);

        //
        // copy particle data back
        //
        BL_PROFILE_VAR_START(blp_copy);
        pti.SetPosition(xp, yp, zp);
        BL_PROFILE_VAR_STOP(blp_copy);
    }
}
