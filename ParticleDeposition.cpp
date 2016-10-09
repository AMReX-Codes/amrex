
#include <memory>

#include <ParticleContainer.H>
#include <WarpX.H>
#include <PICSAR_f.H>

//
// This funciton is assumed to be called after ParticlePush, and therefore
// gaminv is update to date.
//

void
MyParticleContainer::CurrentDeposition (MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt) const
{
    BL_PROFILE("MyPC::CurrentDeposition()");
    BL_PROFILE_VAR_NS("MyPC::CurrentDeposition::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PICSAR::CurrentDeposition", blp_pxr);

    const int lev = 0;
    const Geometry& gm   = m_gdb->Geom(lev);
    const BoxArray& ba   = jx.boxArray();
    const Real*     dx   = gm.CellSize();
    const PMap&     pmap = m_particles[lev];

    BL_ASSERT(OnSameGrids(lev,jx));

    jx.setVal(0.0);
    jy.setVal(0.0);
    jz.setVal(0.0);

    // Charge
    Real q = 1.;

    Array<Real> xp, yp, zp, wp, uxp, uyp, uzp;

    // Loop over the grids containing particles
    for (const auto& kv : pmap)
    {
        const  int  gid = kv.first;
        const PBox& pbx = kv.second;

        const long np = pbx.size();

	BL_PROFILE_VAR_START(blp_copy);

	// 1D Arrays of particle attributes
	 xp.resize(np);
	 yp.resize(np);
	 zp.resize(np);
	 wp.resize(np);
	uxp.resize(np);
	uyp.resize(np);
	uzp.resize(np);

	// If we save xp, yp, zp, ... like Exp in partdata, we could save the copy here.
	
	// Loop over particles in that box 
        for (auto i = 0; i < np; ++i)
        {
	    const auto& p = pbx[i];
	    BL_ASSERT(p.m_id > 0);
	    xp[i]  = p.m_pos[0];
	    yp[i]  = p.m_pos[1];
	    zp[i]  = p.m_pos[2];
 	    wp[i]  = p.m_data[PIdx::w]; 
 	    uxp[i] = p.m_data[PIdx::ux]; 
 	    uyp[i] = p.m_data[PIdx::uy]; 
 	    uzp[i] = p.m_data[PIdx::uz]; 
        }

	BL_PROFILE_VAR_STOP(blp_copy);

	const auto& pdata = partdata.find(gid)->second;
	const auto& gaminv = pdata[PIdx::gaminv];

	const Box& box = BoxLib::enclosedCells(ba[gid]);
	long nx = box.length(0);
	long ny = box.length(1);
	long nz = box.length(2); 

	long nox = 1;
	long noy = 1;
	long noz = 1;

	long lvect = 8;
	long curr_depo_algo = 3;

	RealBox grid_box = RealBox( box, dx, gm.ProbLo() );
	const Real* xyzmin = grid_box.lo();

	long ng = jx.nGrow();
	
	BL_ASSERT(ng == jx.nGrow());
	BL_ASSERT(ng == jy.nGrow());
	BL_ASSERT(ng == jz.nGrow());

	BL_PROFILE_VAR_START(blp_pxr);

#if 0
	depose_jxjyjz_vecHVv2_1_1_1 (jx[gid].dataPtr(), jy[gid].dataPtr(), jz[gid].dataPtr(),
				     &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
				     uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), 
				     gaminv->dataPtr(), wp.dataPtr(), &q, 
				     &xyzmin[0], &xyzmin[1], &xyzmin[2], 
				     &dt, &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				     &ng, &ng, &ng);
#endif

	warpx_current_deposition( jx[gid].dataPtr(), jy[gid].dataPtr(), jz[gid].dataPtr(),
                                    &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
                                    uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), 
                                    gaminv->dataPtr(), wp.dataPtr(), &q, 
                                    &xyzmin[0], &xyzmin[1], &xyzmin[2], 
                                    &dt, &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				    &ng, &ng, &ng, &nox,&noy,&noz,&lvect,&curr_depo_algo);

	// depose_jxjyjz_scalar_1_1_1( jx[gid].dataPtr(), jy[gid].dataPtr(), jz[gid].dataPtr(),
 //                                    &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
 //                                    uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), 
 //                                    gaminv->dataPtr(), wp.dataPtr(), &q, 
 //                                    &xyzmin[0], &xyzmin[1], &xyzmin[2], 
 //                                    &dt, &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
	// 			    &ng, &ng, &ng);

	BL_PROFILE_VAR_STOP(blp_pxr);
    }

    jx.SumBoundary(gm.periodicity());
    jy.SumBoundary(gm.periodicity());
    jz.SumBoundary(gm.periodicity());
}
