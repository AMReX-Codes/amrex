
#include <ParticleContainer.H>
#include <PICSAR_f.H>
#include <WarpX.H>

MyParticleContainer::MyParticleContainer (AmrCore* amr_core)
    : ParticleContainer<PIdx::nattribs,0,std::vector<Particle<PIdx::nattribs,0> > >
      (amr_core->GetParGDB())
{
    this->SetVerbose(0);

    m_particles.reserve(m_gdb->maxLevel()+1);
    m_particles.resize (m_gdb->finestLevel()+1);

    m_partdata.reserve(m_gdb->maxLevel()+1);
    m_partdata.resize (m_gdb->finestLevel()+1);
}

void
MyParticleContainer::AllocData ()
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
		d = std::unique_ptr<Array<Real> >(new Array<Real>());
	    }
	}
    }
}

void
MyParticleContainer::Evolve (int lev,
			     const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
			     const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
			     MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt)
{
    BL_PROFILE("MyPC::Evolve()");
    BL_PROFILE_VAR_NS("MyPC::Evolve::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PICSAR::FieldGather", blp_pxr_fg);
    BL_PROFILE_VAR_NS("PICSAR::ParticlePush", blp_pxr_pp);
    BL_PROFILE_VAR_NS("PICSAR::CurrentDeposition", blp_pxr_cd);

    const Geometry& gm  = m_gdb->Geom(lev);
    const BoxArray& ba  = Ex.boxArray();
    const Real*     dx  = gm.CellSize();

    const long ng_eb = Ex.nGrow();
    const long ng_j  = jx.nGrow();
    
    BL_ASSERT(OnSameGrids(lev,Ex));

    PMap& pmap = m_particles[lev];
    auto& partleveldata = m_partdata[lev];

    Array<Real> xp, yp, zp, wp, uxp, uyp, uzp, giv;

    jx.setVal(0.0);
    jy.setVal(0.0);
    jz.setVal(0.0);

    for (auto& kv : pmap)
    {
	const int gid = kv.first;
        PBox&     pbx = kv.second;

	const long np = pbx.size();

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

	auto& pdata = partleveldata[gid];
	auto& Exp = pdata[PIdx::Ex];
	auto& Eyp = pdata[PIdx::Ey];
	auto& Ezp = pdata[PIdx::Ez];
	auto& Bxp = pdata[PIdx::Bx];
	auto& Byp = pdata[PIdx::By];
	auto& Bzp = pdata[PIdx::Bz];

	Exp->resize(np);
	Eyp->resize(np);
	Ezp->resize(np);
	Bxp->resize(np);
	Byp->resize(np);
	Bzp->resize(np);

	xp.resize(np);
	yp.resize(np);
	zp.resize(np);
	wp.resize(np);
	uxp.resize(np);
	uyp.resize(np);
	uzp.resize(np);
	giv.resize(np);

	// copy data from particle container to temp arrays
	BL_PROFILE_VAR_START(blp_copy);
        for (auto i = 0; i < np; ++i)
        {
	    const auto& p = pbx[i];
	    BL_ASSERT(p.m_id > 0);
	    xp[i] = p.m_pos[0];
	    yp[i] = p.m_pos[1];
	    zp[i] = p.m_pos[2];
 	    wp[i]  = p.m_data[PIdx::w]; 
 	    uxp[i] = p.m_data[PIdx::ux]; 
 	    uyp[i] = p.m_data[PIdx::uy]; 
 	    uzp[i] = p.m_data[PIdx::uz]; 
	}
	BL_PROFILE_VAR_STOP(blp_copy);
	
	const Box& box = BoxLib::enclosedCells(ba[gid]);
	long nx = box.length(0);
	long ny = box.length(1);
	long nz = box.length(2); 

	RealBox grid_box = RealBox( box, dx, gm.ProbLo() );
	const Real* xyzmin = grid_box.lo();

	{       // Field Gather
	    const long order = 1;
	    const int ll4symtry          = false;
	    const int l_lower_order_in_v = true;

	    BL_PROFILE_VAR_START(blp_pxr_fg);
	    warpx_geteb3d_energy_conserving(&np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(),
					    Exp->dataPtr(),Eyp->dataPtr(),Ezp->dataPtr(),
					    Bxp->dataPtr(),Byp->dataPtr(),Bzp->dataPtr(),
					    &xyzmin[0], &xyzmin[1], &xyzmin[2],
					    &dx[0], &dx[1], &dx[2],
					    &nx, &ny, &nz, &ng_eb, &ng_eb, &ng_eb, 
					    &order, &order, &order, 
					    exfab.dataPtr(), eyfab.dataPtr(), ezfab.dataPtr(),
					    bxfab.dataPtr(), byfab.dataPtr(), bzfab.dataPtr(),
					    &ll4symtry, &l_lower_order_in_v, &WarpX::field_gathering_algo);
	    BL_PROFILE_VAR_STOP(blp_pxr_fg);
	}

	{       // Particle Push
	    
	    BL_PROFILE_VAR_START(blp_pxr_pp);
	    warpx_particle_pusher(&np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(),
				  uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), giv.dataPtr(),
				  Exp->dataPtr(), Eyp->dataPtr(), Ezp->dataPtr(),
				  Bxp->dataPtr(), Byp->dataPtr(), Bzp->dataPtr(),
				  &this->charge, &this->mass, &dt, &WarpX::particle_pusher_algo);
	    BL_PROFILE_VAR_STOP(blp_pxr_pp);
	}
	
	{    // Current Deposition
	    long nox = 1;
	    long noy = 1;
	    long noz = 1;
	    long lvect = 8;

	    BL_PROFILE_VAR_START(blp_pxr_cd);
	    warpx_current_deposition(jxfab.dataPtr(), jyfab.dataPtr(), jzfab.dataPtr(),
				     &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
				     uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), 
				     giv.dataPtr(), wp.dataPtr(), &this->charge, 
				     &xyzmin[0], &xyzmin[1], &xyzmin[2], 
				     &dt, &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				     &ng_j, &ng_j, &ng_j, &nox,&noy,&noz,&lvect,&WarpX::current_deposition_algo);
	    BL_PROFILE_VAR_STOP(blp_pxr_cd);
	}

	// copy particle data back
	BL_PROFILE_VAR_START(blp_copy);
	for (auto i = 0; i < np; ++i)
        {
	    auto& p = pbx[i];
	    BL_ASSERT(p.m_id > 0);
	    p.m_pos[0] = xp[i];
	    p.m_pos[1] = yp[i];
	    p.m_pos[2] = zp[i];
 	    p.m_data[PIdx::ux] = uxp[i];
 	    p.m_data[PIdx::uy] = uyp[i];
 	    p.m_data[PIdx::uz] = uzp[i];
        }
	BL_PROFILE_VAR_STOP(blp_copy);
    }

    jx.SumBoundary(gm.periodicity());
    jy.SumBoundary(gm.periodicity());
    jz.SumBoundary(gm.periodicity());
}

