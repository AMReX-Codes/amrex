
#include <limits>

#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpX_f.H>
#include <ParticleContainer.H>
#include <ParticleIterator.H>

using namespace amrex;

LaserParticleContainer::LaserParticleContainer (AmrCore* amr_core, int ispecies)
    : WarpXParticleContainer(amr_core, ispecies)
{
    charge = PhysConst::q_e; // note that q_e is defined to be positive.
    mass = std::numeric_limits<Real>::max();
}

void
LaserParticleContainer::AllocData ()
{
    // have to resize here, not in the constructor because GDB was not
    // ready in constructor.
    m_particles.resize(GDB().finestLevel()+1);
}

void
LaserParticleContainer::InitData ()
{
    
}

void
LaserParticleContainer::Evolve (int lev,
				const MultiFab&, const MultiFab&, const MultiFab&,
				const MultiFab&, const MultiFab&, const MultiFab&,
				MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt)
{
    BL_PROFILE("Laser::Evolve()");
    BL_PROFILE_VAR_NS("Laser::Evolve::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PICSAR::LaserParticlePush", blp_pxr_pp);
    BL_PROFILE_VAR_NS("PICSAR::LaserCurrentDepo", blp_pxr_cd);

    const Geometry& gm  = GDB().Geom(lev);
    const BoxArray& ba  = jx.boxArray();

#if (BL_SPACEDIM == 3)
    const Real* dx = gm.CellSize();
#elif (BL_SPACEDIM == 2)
    Real dx[3] = { gm.CellSize(0), std::numeric_limits<Real>::quiet_NaN(), gm.CellSize(1) };
#endif

#if (BL_SPACEDIM == 3)
    long ngx_j  = jx.nGrow();
    long ngy_j  = ngx_j;
    long ngz_j  = ngx_j;
#elif (BL_SPACEDIM == 2)
    long ngx_j  = jx.nGrow();;
    long ngy_j  = 0;
    long ngz_j  = ngx_j;
#endif
    
    BL_ASSERT(OnSameGrids(lev,jx));

    {
	Array<Real> xp, yp, zp, wp, uxp, uyp, uzp, giv;

	PartIterInfo info {lev, do_tiling, tile_size};
	for (PartIter pti(*this, info); pti.isValid(); ++pti)
	{
	    const int  gid = pti.index();
	    const Box& tbx = pti.tilebox();
	    const Box& vbx = pti.validbox();
	    const long np  = pti.numParticles();

	    // Data on the grid
	    FArrayBox& jxfab = jx[gid];
	    FArrayBox& jyfab = jy[gid];
	    FArrayBox& jzfab = jz[gid];

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


	    const Box& box = amrex::enclosedCells(ba[gid]);
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
	    // Particle Push
	    //
	    BL_PROFILE_VAR_START(blp_pxr_pp);
	    warpx_laser_pusher(&np, xp.data(), yp.data(), zp.data(),
			       uxp.data(), uyp.data(), uzp.data(), giv.data(),
			       &this->charge, &dt, 
			       &WarpX::laser_pusher_algo);
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
