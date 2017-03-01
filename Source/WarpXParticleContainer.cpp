
#include <limits>

#include <ParticleContainer.H>
#include <ParticleIterator.H>
#include <WarpX_f.H>
#include <WarpX.H>

bool    WarpXParticleContainer::do_tiling = 0;
IntVect WarpXParticleContainer::tile_size   { D_DECL(1024000,8,8) };

WarpXParticleContainer::WarpXParticleContainer (AmrCore* amr_core, int ispecies)
    : ParticleContainer<PIdx::nattribs,0,std::vector<Particle<PIdx::nattribs,0> > >
      (amr_core->GetParGDB())
    , species_id(ispecies)
{
    this->SetVerbose(0);

    m_particles.reserve(m_gdb->maxLevel()+1);
    m_particles.resize (m_gdb->finestLevel()+1);

    ReadParameters();
}

void
WarpXParticleContainer::ReadParameters ()
{
    static bool initialized = false;
    if (!initialized)
    {
	ParmParse pp("particles");

	pp.query("do_tiling",  do_tiling);

	Array<int> ts(BL_SPACEDIM);
	if (pp.queryarr("tile_size", ts)) {
	    tile_size = IntVect(ts);
	}

	initialized = true;
    }
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

    int gid = 0;

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

    for (int i = ibegin; i < iend; ++i)
    {
	ParticleType p;
	p.m_id  = ParticleBase::NextID();
	p.m_cpu = ParallelDescriptor::MyProc();
	p.m_lev = lev;
	p.m_grid = gid; 

#if (BL_SPACEDIM == 3)	
	p.m_pos[0] = x[i];
	p.m_pos[1] = y[i];
	p.m_pos[2] = z[i];
#else
	p.m_pos[0] = x[i];
	p.m_pos[1] = z[i];
#endif
	
	p.m_data[PIdx::w] = weight[i];
	
	for (int j = 1; j < PIdx::nattribs; ++j) {
	    p.m_data[j] = 0;
	}
	
	p.m_data[PIdx::ux] = vx[i];
	p.m_data[PIdx::uy] = vy[i];
	p.m_data[PIdx::uz] = vz[i];
	
	if (ParticleBase::Where(p,m_gdb)) // this will update m_lev, m_grid, and m_cell
	{
	    gid = p.m_grid;
	    m_particles[p.m_lev][p.m_grid].push_back(p);
	}
    }

    Redistribute(true);

    auto npart_after = TotalNumberOfParticles();  // xxxxx move this into if (verbose > xxx)
    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Total number of particles added: " << npart_after - npart_before << std::endl;
    }
}

std::unique_ptr<MultiFab>
WarpXParticleContainer::GetChargeDensity (int lev, bool local)
{
    const Geometry& gm = m_gdb->Geom(lev);
    const BoxArray& ba = m_gdb->ParticleBoxArray(lev);
    BoxArray nba = ba;
    nba.surroundingNodes();

#if (BL_SPACEDIM == 3)
    const Real* dx = gm.CellSize();
#elif (BL_SPACEDIM == 2)
    Real dx[3] = { gm.CellSize(0), std::numeric_limits<Real>::quiet_NaN(), gm.CellSize(1) };
#endif

    const int ng = WarpX::nox;

    auto rho = std::unique_ptr<MultiFab>(new MultiFab(nba,1,ng));
    rho->setVal(0.0);

    Array<Real> xp, yp, zp, wp;

    PartIterInfo info {lev, do_tiling, tile_size};
    for (PartIter pti(*this, info); pti.isValid(); ++pti)
    {
	const int  gid = pti.index();
	const Box& vbx = pti.validbox();
	const long np  = pti.numParticles();
	
	// Data on the grid
	FArrayBox& rhofab = (*rho)[gid];

	xp.resize(np);
	yp.resize(np);
	zp.resize(np);
	wp.resize(np);

	pti.foreach([&](int i, ParticleType& p) {
#if (BL_SPACEDIM == 3)
		xp[i] = p.m_pos[0];
		yp[i] = p.m_pos[1];
		zp[i] = p.m_pos[2];
#elif (BL_SPACEDIM == 2)
		xp[i] = p.m_pos[0];
		yp[i] = std::numeric_limits<Real>::quiet_NaN();
		zp[i] = p.m_pos[1];
#endif
		wp[i]  = p.m_data[PIdx::w]; 
	    });

	const Box& box = BoxLib::enclosedCells(ba[gid]);
	BL_ASSERT(box == vbx);
#if (BL_SPACEDIM == 3)
	long nx = box.length(0);
	long ny = box.length(1);
	long nz = box.length(2); 
#elif (BL_SPACEDIM == 2)
	long nx = box.length(0);
	long ny = 0;
	long nz = box.length(1); 
#endif
	RealBox grid_box = RealBox( box, gm.CellSize(), gm.ProbLo() );
#if (BL_SPACEDIM == 3)
	const Real* xyzmin = grid_box.lo();
#elif (BL_SPACEDIM == 2)
	Real xyzmin[3] = { grid_box.lo(0), std::numeric_limits<Real>::quiet_NaN(), grid_box.lo(1) };
#endif

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

void
WarpXParticleContainer::PushX (int lev,
                                  Real dt)
{
    BL_PROFILE("PPC::Evolve()");
    BL_PROFILE_VAR_NS("PPC::Evolve::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PICSAR::ParticlePush", blp_pxr_pp);

    //xxxxx not using m_pardata for now. auto& partleveldata = m_partdata[lev];

    {
	Array<Real> xp, yp, zp, wp, uxp, uyp, uzp, giv;

	PartIterInfo info {lev, do_tiling, tile_size};
	for (PartIter pti(*this, info); pti.isValid(); ++pti)
	{
	    const int  gid = pti.index();
	    const Box& tbx = pti.tilebox();
	    const Box& vbx = pti.validbox();
	    const long np  = pti.numParticles();

	    xp.resize(np);
	    yp.resize(np);
	    zp.resize(np);
	    wp.resize(np);
	    uxp.resize(np);
	    uyp.resize(np);
	    uzp.resize(np);
	    giv.resize(np);

	    //
	    // copy data from particle container to temp arrays
	    //
	    BL_PROFILE_VAR_START(blp_copy);
	    pti.foreach([&](int i, ParticleType& p) {
#if (BL_SPACEDIM == 3)
		    xp[i] = p.m_pos[0];
		    yp[i] = p.m_pos[1];
		    zp[i] = p.m_pos[2];
#elif (BL_SPACEDIM == 2)
		    xp[i] = p.m_pos[0];
		    yp[i] = std::numeric_limits<Real>::quiet_NaN();
		    zp[i] = p.m_pos[1];
#endif
		    wp[i]  = p.m_data[PIdx::w];
		    uxp[i] = p.m_data[PIdx::ux];
		    uyp[i] = p.m_data[PIdx::uy];
		    uzp[i] = p.m_data[PIdx::uz];
		});
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
	    pti.foreach([&](int i, ParticleType& p) {
                    BL_ASSERT(p.m_id > 0);
#if (BL_SPACEDIM == 3)
		    p.m_pos[0] = xp[i];
		    p.m_pos[1] = yp[i];
		    p.m_pos[2] = zp[i];
#elif (BL_SPACEDIM == 2)
		    p.m_pos[0] = xp[i];
		    p.m_pos[1] = zp[i];
#endif
                });
            BL_PROFILE_VAR_STOP(blp_copy);
	}
    }
}
