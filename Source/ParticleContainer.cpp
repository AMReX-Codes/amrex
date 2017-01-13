#include <limits>

#include <ParticleContainer.H>
#include <ParticleIterator.H>
#include <WarpX_f.H>
#include <WarpX.H>

int      MultiSpeciesContainer::nspecies  = 1;
bool    SingleSpeciesContainer::do_tiling = 0;
IntVect SingleSpeciesContainer::tile_size   { D_DECL(1024000,8,8) };

SingleSpeciesContainer::SingleSpeciesContainer (AmrCore* amr_core, int ispecies)
    : ParticleContainer<PIdx::nattribs,0,std::vector<Particle<PIdx::nattribs,0> > >
      (amr_core->GetParGDB())
    , species_id(ispecies)
{
    this->SetVerbose(0);

    m_particles.reserve(m_gdb->maxLevel()+1);
    m_particles.resize (m_gdb->finestLevel()+1);

    m_partdata.reserve(m_gdb->maxLevel()+1);
    m_partdata.resize (m_gdb->finestLevel()+1);
}

void
SingleSpeciesContainer::AllocData ()
{
    m_particles.resize(m_gdb->finestLevel()+1);
    m_partdata.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev <= m_gdb->finestLevel(); ++ lev)
    {
	auto& partleveldata = m_partdata[lev];
	const BoxArray& ba = m_gdb->ParticleBoxArray(lev);
	const DistributionMapping& dm = m_gdb->ParticleDistributionMap(lev);

	MultiFab foo(ba, 1, 0, dm, Fab_noallocate);
	for (MFIter mfi(foo); mfi.isValid(); ++mfi)
	{
	    int i = mfi.index();
	    partleveldata[i] = Array<std::unique_ptr<Array<Real> > > (PIdx::npartdata);
	    for (auto& d : partleveldata[i])
	    {
		d.reset(new Array<Real>());
	    }
	}
    }
}

MultiSpeciesContainer::MultiSpeciesContainer (AmrCore* amr_core)
{
    ReadStaticParameters();

    species.resize(nspecies);
    for (int i = 0; i < nspecies; ++i) {
	species[i].reset(new SingleSpeciesContainer(amr_core, i));
    }
}

void
MultiSpeciesContainer::ReadStaticParameters ()
{
    static bool initialized = false;
    if (!initialized)
    {
	ParmParse pp("particles");

	pp.query("nspecies", nspecies);
	BL_ASSERT(nspecies >= 1);

	pp.query("do_tiling",  SingleSpeciesContainer::do_tiling);

	Array<int> ts(BL_SPACEDIM);
	if (pp.queryarr("tile_size", ts)) {
	    SingleSpeciesContainer::tile_size = IntVect(ts);
	}

	initialized = true;
    }
}

void
MultiSpeciesContainer::AllocData ()
{
    for (auto& spec : species) {
	spec->AllocData();
    }
}

void
MultiSpeciesContainer::InitData ()
{
    for (auto& spec : species) {
	spec->InitData();
    }
}

void
MultiSpeciesContainer::Evolve (int lev,
			     const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
			     const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
			     MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt)
{
    jx.setVal(0.0);
    jy.setVal(0.0);
    jz.setVal(0.0);

    for (auto& spec : species) {
	spec->Evolve(lev, Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, dt);
    }    

    const Geometry& gm = species[0]->m_gdb->Geom(lev);
    jx.SumBoundary(gm.periodicity());
    jy.SumBoundary(gm.periodicity());
    jz.SumBoundary(gm.periodicity());
}

void
SingleSpeciesContainer::Evolve (int lev,
				const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
				const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
				MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt)
{
    BL_PROFILE("SPC::Evolve()");
    BL_PROFILE_VAR_NS("SPC::Evolve::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PICSAR::FieldGather", blp_pxr_fg);
    BL_PROFILE_VAR_NS("PICSAR::ParticlePush", blp_pxr_pp);
    BL_PROFILE_VAR_NS("PICSAR::CurrentDeposition", blp_pxr_cd);

    const Geometry& gm  = m_gdb->Geom(lev);
    const BoxArray& ba  = Ex.boxArray();

#if (BL_SPACEDIM == 3)
    const Real* dx = gm.CellSize();
#elif (BL_SPACEDIM == 2)
    Real dx[3] = { gm.CellSize(0), std::numeric_limits<Real>::quiet_NaN(), gm.CellSize(1) };
#endif

#if (BL_SPACEDIM == 3)
    long ngx_eb = Ex.nGrow();
    long ngy_eb = ngx_eb;
    long ngz_eb = ngx_eb;
    long ngx_j  = jx.nGrow();
    long ngy_j  = ngx_j;
    long ngz_j  = ngx_j;
#elif (BL_SPACEDIM == 2)
    long ngx_eb = Ex.nGrow();
    long ngy_eb = 0;
    long ngz_eb = ngx_eb;
    long ngx_j  = jx.nGrow();;
    long ngy_j  = 0;
    long ngz_j  = ngx_j;
#endif
    
    BL_ASSERT(OnSameGrids(lev,Ex));

    //xxxxx not using m_pardata for now. auto& partleveldata = m_partdata[lev];

    {
	Array<Real> xp, yp, zp, wp, uxp, uyp, uzp, giv;
	Array<Real> Exp, Eyp, Ezp, Bxp, Byp, Bzp;

	PartIterInfo info {lev, do_tiling, tile_size};
	for (PartIter pti(*this, info); pti.isValid(); ++pti)
	{
	    const int  gid = pti.index();
	    const Box& tbx = pti.tilebox();
	    const Box& vbx = pti.validbox();
	    const long np  = pti.numParticles();

	    // Data on the grid
	    const FArrayBox& exfab = Ex[gid];
	    const FArrayBox& eyfab = Ey[gid];
	    const FArrayBox& ezfab = Ez[gid];
	    const FArrayBox& bxfab = Bx[gid];
	    const FArrayBox& byfab = By[gid];
	    const FArrayBox& bzfab = Bz[gid];
	    FArrayBox&       jxfab = jx[gid];
	    FArrayBox&       jyfab = jy[gid];
	    FArrayBox&       jzfab = jz[gid];

	    Exp.assign(np,0.0);
	    Eyp.assign(np,0.0);
	    Ezp.assign(np,0.0);
	    Bxp.assign(np,0.0);
	    Byp.assign(np,0.0);
	    Bzp.assign(np,0.0);

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
	    
	    //
	    // Field Gather
	    //
	    const int ll4symtry          = false;
	    const int l_lower_order_in_v = true;
            long lvect_fieldgathe = 64;
	    BL_PROFILE_VAR_START(blp_pxr_fg);
	    warpx_geteb_energy_conserving(&np, xp.data(), yp.data(), zp.data(),
					  Exp.data(),Eyp.data(),Ezp.data(),
					  Bxp.data(),Byp.data(),Bzp.data(),
					  &xyzmin[0], &xyzmin[1], &xyzmin[2],
					  &dx[0], &dx[1], &dx[2],
					  &nx, &ny, &nz, &ngx_eb, &ngy_eb, &ngz_eb, 
					  &WarpX::nox, &WarpX::noy, &WarpX::noz, 
					  exfab.dataPtr(), eyfab.dataPtr(), ezfab.dataPtr(),
					  bxfab.dataPtr(), byfab.dataPtr(), bzfab.dataPtr(),
					  &ll4symtry, &l_lower_order_in_v,
                                          &lvect_fieldgathe,
		                          &WarpX::field_gathering_algo);
	    BL_PROFILE_VAR_STOP(blp_pxr_fg);

	    //
	    // Particle Push
	    //
	    BL_PROFILE_VAR_START(blp_pxr_pp);
	    warpx_particle_pusher(&np, xp.data(), yp.data(), zp.data(),
				  uxp.data(), uyp.data(), uzp.data(), giv.data(),
				  Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
				  Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
				  &this->charge, &this->mass, &dt, 
                                  &WarpX::particle_pusher_algo);
	    BL_PROFILE_VAR_STOP(blp_pxr_pp);

	    //
	    // Current Deposition
	    // xxxxx this part needs to be thread safe if we have OpenMP over tiles
	    //
	    long lvect = 8;
	    BL_PROFILE_VAR_START(blp_pxr_cd);
	    warpx_current_deposition(jxfab.dataPtr(), jyfab.dataPtr(), jzfab.dataPtr(),
				     &np, xp.data(), yp.data(), zp.data(), 
				     uxp.data(), uyp.data(), uzp.data(), 
				     giv.data(), wp.data(), &this->charge, 
				     &xyzmin[0], &xyzmin[1], &xyzmin[2], 
				     &dt, &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				     &ngx_j, &ngy_j, &ngz_j, 
                                     &WarpX::nox,&WarpX::noy,&WarpX::noz,
				     &lvect,&WarpX::current_deposition_algo);
	    BL_PROFILE_VAR_STOP(blp_pxr_cd);

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
		    p.m_data[PIdx::ux] = uxp[i];
		    p.m_data[PIdx::uy] = uyp[i];
		    p.m_data[PIdx::uz] = uzp[i];
                });
            BL_PROFILE_VAR_STOP(blp_copy);
	}
    }
}

std::unique_ptr<MultiFab>
MultiSpeciesContainer::GetChargeDensity (int lev, bool local)
{
    std::unique_ptr<MultiFab> rho = species[0]->GetChargeDensity(lev, true);
    for (int i = 1; i < nspecies; ++i) {
	std::unique_ptr<MultiFab> rhoi = species[i]->GetChargeDensity(lev, true);
	MultiFab::Add(*rho, *rhoi, 0, 0, 1, rho->nGrow());
    }
    if (!local) {
	const Geometry& gm = species[0]->m_gdb->Geom(lev);
	rho->SumBoundary(gm.periodicity());
    }
}

std::unique_ptr<MultiFab>
SingleSpeciesContainer::GetChargeDensity (int lev, bool local)
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
SingleSpeciesContainer::AddNParticles (int n, const Real* x, const Real* y, const Real* z,
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
	
	p.m_pos[0] = x[i];
	p.m_pos[1] = y[i];
	p.m_pos[2] = z[i];
	
	p.m_data[PIdx::w] = weight[i];
	
	for (int i = 1; i < PIdx::nattribs; i++) {
	    p.m_data[i] = 0;
	}
	
	p.m_data[PIdx::ux] = vx[i];
	p.m_data[PIdx::uy] = vy[i];
	p.m_data[PIdx::uz] = vz[i];
	
	if (!ParticleBase::Where(p,m_gdb)) // this will update m_lev, m_grid, and m_cell
	{
	    BoxLib::Abort("Invalid particle in ParticleContainer::AddNParticles()");
	}

	gid = p.m_grid;

	m_particles[p.m_lev][p.m_grid].push_back(p);
    }

    Redistribute(true);

    auto npart_after = TotalNumberOfParticles();  // xxxxx move this into if (verbose > xxx)
    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Total number of particles injected: " << npart_after - npart_before << std::endl;
    }
}

void
MultiSpeciesContainer::Checkpoint (const std::string& dir, const std::string& name)
{
    for (int i = 0; i < nspecies; ++i) {
	std::string diri = dir + std::to_string(i);
	Checkpoint(diri, name);
    }
}

void
MultiSpeciesContainer::Redistribute (bool where_called)
{
    for (auto& spec : species) {
	spec->Redistribute(where_called,true);
    }
}
