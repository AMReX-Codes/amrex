
#include <memory>

#include <ParticleContainer.H>
#include <PICSAR_f.H>

void
MyParticleContainer::ChargeDeposition (MultiFab& mf_to_be_filled, int lev, int order) const
{
    BL_PROFILE("MyPC::ChargeDeposition")
    BL_PROFILE_VAR_NS("MyPC::ChargeDeposition::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PXR::ChargeDeposition", blp_pxr_cd);

    MultiFab* mf_pointer;

    // We are only deposing one quantity
    int ncomp = 1;

    std::unique_ptr<MultiFab> raii;
    if (OnSameGrids(lev, mf_to_be_filled))
    {
        // If we are already working with the internal mf defined on the 
        // particle_box_array, then we just work with this.
        mf_pointer = &mf_to_be_filled;
    }
    else
    {
        // If mf_to_be_filled is not defined on the particle_box_array, then we need 
        // to make a temporary here and copy into mf_to_be_filled at the end.
	raii = std::unique_ptr<MultiFab>
	    (new MultiFab(m_gdb->ParticleBoxArray(lev), ncomp, mf_to_be_filled.nGrow(),
			  m_gdb->ParticleDistributionMap(lev), Fab_allocate));
	mf_pointer = raii.get();
    }

    // Putting the density to 0 before depositing the charge
    mf_pointer->setVal(0.0);

    const Geometry& gm          = m_gdb->Geom(lev);
    const BoxArray& ba          = mf_pointer->boxArray();
    const Real*     dx          = gm.CellSize();

    const PMap&     pmap        = m_particles[lev];

    // Charge
    Real q = 1.;
    Array<Real> xp, yp, zp, wp;

    // Loop over the grids containing particles
    for (auto& kv : pmap)
    {
        const  int  pgr = kv.first;
        const PBox& pbx = kv.second;

        FArrayBox&  fab = (*mf_pointer)[pgr];

	xp.resize(0);
	yp.resize(0);
	zp.resize(0);
	wp.resize(0);

	xp.reserve( pbx.size() );
	yp.reserve( pbx.size() );
	zp.reserve( pbx.size() );
	wp.reserve( pbx.size() );

	const Box & bx = BoxLib::enclosedCells(ba[pgr]);
	long nx = bx.length(0);
	long ny = bx.length(1);
	long nz = bx.length(2); 

	RealBox grid_box = RealBox( bx, dx, gm.ProbLo() );
	const Real* xyzmin = grid_box.lo();

	long ng = mf_pointer->nGrow();
	long lvect = 8;

	BL_PROFILE_VAR_START(blp_copy);
	
	// Loop over particles in that box (to change array layout)
        for (const auto& p : pbx)
        {
            if (p.m_id <= 0) {
	      continue;
	    }
	    xp.push_back( p.m_pos[0] );
	    yp.push_back( p.m_pos[1] );
	    zp.push_back( p.m_pos[2] );
 	    wp.push_back( 1. ); 
        }
	long np = wp.size();

	BL_PROFILE_VAR_STOP(blp_copy);

	BL_PROFILE_VAR_START(blp_pxr_cd);

#if 0
	depose_rho_vecHVv2_1_1_1( fab.dataPtr(), &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
				  wp.dataPtr(), &q, &xyzmin[0], &xyzmin[1],
				  &xyzmin[2], &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				  &ng, &ng, &ng, &lvect );
#endif

	if (order == 1) 
	{
	    depose_rho_scalar_1_1_1( fab.dataPtr(), &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
	    	                     wp.dataPtr(), &q, &xyzmin[0], &xyzmin[1],
	    	                     &xyzmin[2], &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				     &ng, &ng, &ng, &lvect);
	} 
        else if (order == 2) 
	{
	    depose_rho_scalar_2_2_2( fab.dataPtr(), &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
	    	                     wp.dataPtr(), &q, &xyzmin[0], &xyzmin[1],
	    	                     &xyzmin[2], &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				     &ng, &ng, &ng, &lvect);
	} 
        else if (order == 3) 
	{
	    depose_rho_scalar_3_3_3( fab.dataPtr(), &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
	    	                     wp.dataPtr(), &q, &xyzmin[0], &xyzmin[1],
	    	                     &xyzmin[2], &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				     &ng, &ng, &ng, &lvect);
	} 
  
	BL_PROFILE_VAR_STOP(blp_pxr_cd);
    }

    mf_pointer->SumBoundary(gm.periodicity());

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.  Ghost cells are not copied.
    if (mf_pointer != &mf_to_be_filled)
    {
        mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
    }
}

//
// This is the single-level version for current deposition on faces
//
void
MyParticleContainer::CurrentDeposition (Array< std::unique_ptr<MultiFab> >& mf_to_be_filled,
					int lev, Real dt) const
{
    BL_PROFILE("MyPC::CurrentDeposition");
    BL_PROFILE_VAR_NS("MyPC::CurrentDeposition::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PXR::CurrentDeposition", blp_pxr_cd);

    Array<MultiFab*> mf_pointer(BL_SPACEDIM);

    // We are only deposing one quantity on each face
    int ncomp  = 1;
    int nghost = mf_to_be_filled[0]->nGrow();

    Array< std::unique_ptr<MultiFab> > raii;
    if (OnSameGrids(lev, *mf_to_be_filled[0]))
    {
        // If we are already working with the internal mf defined on the 
        // particle_box_array, then we just work with this.
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    mf_pointer[i] = mf_to_be_filled[i].get();
	}
    }
    else
    {
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    raii[i] = std::unique_ptr<MultiFab>
		(new MultiFab(BoxLib::convert(m_gdb->ParticleBoxArray(lev),
					      mf_to_be_filled[0]->boxArray().ixType()),
			      ncomp, nghost, m_gdb->ParticleDistributionMap(lev), Fab_allocate));
	    mf_pointer[i] = raii[i].get();
	}
    }

    const Geometry& gm          = m_gdb->Geom(lev);
    const BoxArray& ba          = mf_pointer[0]->boxArray();
    const Real*     dx          = gm.CellSize();
    const PMap&     pmap        = m_particles[lev];

    // Setting the current to 0 before depositing the charge
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	mf_pointer[i]->setVal(0.0);
    }
   
    // Charge
    Real q = 1.;

    Array<Real> xp, yp, zp, wp, uxp, uyp, uzp, gip;
    FArrayBox fab;

    // Loop over the grids containing particles
    for (auto& kv : pmap)
    {
        const  int  pgr = kv.first;
        const PBox& pbx = kv.second;

        FArrayBox&  fabx_orig = (*mf_pointer[0])[pgr];
        FArrayBox&  faby_orig = (*mf_pointer[1])[pgr];
        FArrayBox&  fabz_orig = (*mf_pointer[2])[pgr];

	const Box& xbx_orig = fabx_orig.box();
	const Box& ybx_orig = faby_orig.box();
	const Box& zbx_orig = fabz_orig.box();

	Box ndbx = BoxLib::surroundingNodes(xbx_orig);

	BL_PROFILE_VAR_START(blp_copy);

	fab.resize(ndbx, 3);
	fab.setVal(0.0);

	 xp.resize(0);
	 yp.resize(0);
	 zp.resize(0);
	 wp.resize(0);
	uxp.resize(0);
	uyp.resize(0);
	uzp.resize(0);
	gip.resize(0);

	 xp.reserve( pbx.size() );
	 yp.reserve( pbx.size() );
	 zp.reserve( pbx.size() );
	 wp.reserve( pbx.size() );
	uxp.reserve( pbx.size() );
	uyp.reserve( pbx.size() );
	uzp.reserve( pbx.size() );
	gip.reserve( pbx.size() );

	const Box & bx = BoxLib::enclosedCells(ba[pgr]);
	RealBox grid_box = RealBox( bx, dx, gm.ProbLo() );
	const Real* xyzmin = grid_box.lo();
	long nx = bx.length(0);
	long ny = bx.length(1);
	long nz = bx.length(2); 
	long ng = mf_pointer[0]->nGrow();
	
	// Loop over particles in that box (to change array layout)
        for (const auto& p : pbx)
        {
            if (p.m_id <= 0) {
	      continue;
	    }

            // (x,y,z) position
	    xp.push_back( p.m_pos[0] );
	    yp.push_back( p.m_pos[1] );
	    zp.push_back( p.m_pos[2] );

            // weights
 	    wp.push_back( p.m_data[PIdx::w] ); 
	    
 	    uxp.push_back( p.m_data[PIdx::ux] ); 
 	    uyp.push_back( p.m_data[PIdx::uy] ); 
 	    uzp.push_back( p.m_data[PIdx::uz] ); 

            // gaminv 
 	    gip.push_back( p.m_data[PIdx::gaminv] ); 
        }
	long np = gip.size();

	BL_PROFILE_VAR_STOP(blp_copy);

	BL_PROFILE_VAR_START(blp_pxr_cd);

#if 0
	depose_jxjyjz_vecHVv2_1_1_1(fab.dataPtr(0), fab.dataPtr(1), fab.dataPtr(2),
                                    &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
                                    uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), 
                                    gip.dataPtr(), wp.dataPtr(), &q, 
                                    &xyzmin[0], &xyzmin[1], &xyzmin[2], 
                                    &dt, &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				    &ng, &ng, &ng);
#endif

	depose_jxjyjz_scalar_1_1_1( fab.dataPtr(0), fab.dataPtr(1), fab.dataPtr(2),
                                    &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
                                    uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), 
                                    gip.dataPtr(), wp.dataPtr(), &q, 
                                    &xyzmin[0], &xyzmin[1], &xyzmin[2], 
                                    &dt, &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				    &ng, &ng, &ng);

	BL_PROFILE_VAR_STOP(blp_pxr_cd);

	BL_PROFILE_VAR_START(blp_copy);

	fab.SetBoxType(xbx_orig.ixType());
	fabx_orig.copy(fab, xbx_orig, 0, xbx_orig, 0, 1);

	fab.SetBoxType(ybx_orig.ixType());
	faby_orig.copy(fab, ybx_orig, 0, ybx_orig, 0, 1);

	fab.SetBoxType(zbx_orig.ixType());
	fabz_orig.copy(fab, zbx_orig, 0, zbx_orig, 0, 1);

	BL_PROFILE_VAR_STOP(blp_copy);
    }

    for (int i = 0; i < BL_SPACEDIM; ++i) {
	mf_pointer[i]->SumBoundary(gm.periodicity());
    }

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.
    if (mf_pointer[0] != mf_to_be_filled[0].get())
    {
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    mf_to_be_filled[i]->copy(*mf_pointer[i],0,0,ncomp);
	}
    }
}

