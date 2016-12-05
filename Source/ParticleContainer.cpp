#include <limits>

#include <ParticleContainer.H>
#include <ParticleIterator.H>
#include <WarpX_f.H>
#include <WarpX.H>

int     MyParticleContainer::do_tiling = 0;
IntVect MyParticleContainer::tile_size   { D_DECL(1024000,8,8) };

MyParticleContainer::MyParticleContainer (AmrCore* amr_core)
    : ParticleContainer<PIdx::nattribs,0,std::vector<Particle<PIdx::nattribs,0> > >
      (amr_core->GetParGDB())
{
    ReadStaticParameters();

    this->SetVerbose(0);

    m_particles.reserve(m_gdb->maxLevel()+1);
    m_particles.resize (m_gdb->finestLevel()+1);

    m_partdata.reserve(m_gdb->maxLevel()+1);
    m_partdata.resize (m_gdb->finestLevel()+1);
}

void
MyParticleContainer::ReadStaticParameters ()
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
		d.reset(new Array<Real>());
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

    jx.setVal(0.0);
    jy.setVal(0.0);
    jz.setVal(0.0);

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

	    Exp.resize(np);
	    Eyp.resize(np);
	    Ezp.resize(np);
	    Bxp.resize(np);
	    Byp.resize(np);
	    Bzp.resize(np);

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

    jx.SumBoundary(gm.periodicity());
    jy.SumBoundary(gm.periodicity());
    jz.SumBoundary(gm.periodicity());
}

