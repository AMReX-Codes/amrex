#include "ParticleContainer.H"

extern "C" {
 void depose_rho_vecHVv2_1_1_1( Real* rho, const long* np,
	   const Real* xp, const Real* yp, const Real* zp,
	   const Real* w, const Real* q,
	   const Real* xmin, const Real* ymin, const Real* zmin,
	   const Real* dx, const Real* dy, const Real* dz,
	   const long* nx, const long* ny, const long* nz,
	   const long* nxguard, const long* nyguard, const long* nzguard, const long* lvect );

 void depose_rho_scalar_1_1_1( Real* rho, const long* np,
	   const Real* xp, const Real* yp, const Real* zp,
	   const Real* w, const Real* q,
	   const Real* xmin, const Real* ymin, const Real* zmin,
	   const Real* dx, const Real* dy, const Real* dz,
	   const long* nx, const long* ny, const long* nz,
	   const long* nxguard, const long* nyguard, const long* nzguard);

 void depose_rho_scalar_2_2_2( Real* rho, const long* np,
	   const Real* xp, const Real* yp, const Real* zp,
	   const Real* w, const Real* q,
	   const Real* xmin, const Real* ymin, const Real* zmin,
	   const Real* dx, const Real* dy, const Real* dz,
	   const long* nx, const long* ny, const long* nz,
	   const long* nxguard, const long* nyguard, const long* nzguard);

 void depose_rho_scalar_3_3_3( Real* rho, const long* np,
	   const Real* xp, const Real* yp, const Real* zp,
	   const Real* w, const Real* q,
	   const Real* xmin, const Real* ymin, const Real* zmin,
	   const Real* dx, const Real* dy, const Real* dz,
	   const long* nx, const long* ny, const long* nz,
	   const long* nxguard, const long* nyguard, const long* nzguard);

 void depose_jxjyjz_vecHVv2_1_1_1(Real* jx, Real* jy, Real* jz, const long* np,
           const Real* xp, const Real* yp, const Real* zp,
	   const Real* uxp, const Real* uyp,const Real* uzp,
	   const Real* gip, const Real* w, const Real* q,
	   const Real* xmin, const Real* ymin, const Real* zmin,
	   const Real* dt, 
	   const Real* dx, const Real* dy, const Real* dz,
	   const long* nx, const long* ny, const long* nz,
	   const long* nxguard, const long* nyguard, const long* nzguard);

 void depose_jxjyjz_scalar_1_1_1(Real* jx, Real* jy, Real* jz, const long* np,
           const Real* xp, const Real* yp, const Real* zp,
	   const Real* uxp, const Real* uyp,const Real* uzp,
	   const Real* gip, const Real* w, const Real* q,
	   const Real* xmin, const Real* ymin, const Real* zmin,
	   const Real* dt, 
	   const Real* dx, const Real* dy, const Real* dz,
	   const long* nx, const long* ny, const long* nz,
	   const long* nxguard, const long* nyguard, const long* nzguard);
}

void
MyParticleContainer::ChargeDeposition(MultiFab& mf_to_be_filled, int lev, int order) const
{
    MultiFab* mf_pointer;

    // We are only deposing one quantity
    int ncomp = 1;

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
        mf_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev), ncomp, mf_to_be_filled.nGrow(),
				  m_gdb->ParticleDistributionMap(lev), Fab_allocate);
    }

    // Putting the density to 0 before depositing the charge
    for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi)  
        (*mf_pointer)[mfi].setVal(0);

    const Real      strttime    = ParallelDescriptor::second();

    const Geometry& gm          = m_gdb->Geom(lev);
    const BoxArray& ba          = mf_pointer->boxArray();
    const Real*     dx          = gm.CellSize();

    const PMap&     pmap        = m_particles[lev];

    // Charge
    Real q = 1.;

    // Loop over the grids containing particles
    for (auto& kv : pmap)
    {
        const  int  pgr = kv.first;
        const PBox& pbx = kv.second;

        FArrayBox&  fab = (*mf_pointer)[pgr];

	Array<Real> xp, yp, zp, wp;
	xp.reserve( pbx.size() );
	yp.reserve( pbx.size() );
	zp.reserve( pbx.size() );
	wp.reserve( pbx.size() );

	const Box & bx = ba[pgr];
	RealBox grid_box = RealBox( bx, dx, gm.ProbLo() );
	const Real* xyzmin = grid_box.lo();
	long nx = bx.length(0)-1, ny = bx.length(1)-1, nz = bx.length(2)-1; 
	long ng = mf_pointer->nGrow();
	long lvect = 8;

        Real strt_copy = ParallelDescriptor::second();
	
	// Loop over particles in that box (to change array layout)
	long np = 0;
        for (const auto& p : pbx)
        {
            if (p.m_id <= 0) {
	      continue;
	    }
	    ++np;
	    xp.push_back( p.m_pos[0] );
	    yp.push_back( p.m_pos[1] );
	    zp.push_back( p.m_pos[2] );
 	    wp.push_back( 1. ); 
        }

        Real end_copy = ParallelDescriptor::second() - strt_copy;

        Real strt_depose = ParallelDescriptor::second();

        if (ParallelDescriptor::IOProcessor()) 
            std::cout << "Time in PicsarChargeDeposition : Copy " << end_copy << '\n';

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
				     &ng, &ng, &ng);
	} 
        else if (order == 2) 
	{
	    depose_rho_scalar_2_2_2( fab.dataPtr(), &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
	    	                     wp.dataPtr(), &q, &xyzmin[0], &xyzmin[1],
	    	                     &xyzmin[2], &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				     &ng, &ng, &ng);
	} 
        else if (order == 3) 
	{
	    depose_rho_scalar_3_3_3( fab.dataPtr(), &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
	    	                     wp.dataPtr(), &q, &xyzmin[0], &xyzmin[1],
	    	                     &xyzmin[2], &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				     &ng, &ng, &ng);
	} 
  
        Real end_depose = ParallelDescriptor::second() - strt_depose;
        if (ParallelDescriptor::IOProcessor()) 
            std::cout << "Time in PicsarChargeDeposition : Depose " << end_depose << '\n';
    }

    Real strt_sumb = ParallelDescriptor::second();

    mf_pointer->SumBoundary(gm.periodicity());

    Real end_sumb = ParallelDescriptor::second() - strt_sumb;

    if (ParallelDescriptor::IOProcessor()) 
        std::cout << "Time in PicsarChargeDeposition : SumBoundary " << end_sumb << '\n';

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled)
    {
        mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
	delete mf_pointer;
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::ChargeDeposition time: " << stoptime << '\n';
        }
    }
}

//
// This is the single-level version for current deposition on faces
//
void
MyParticleContainer::CurrentDeposition(PArray<MultiFab>& mf_to_be_filled, int lev, Real dt) const
{
    MultiFab *mf_pointer_x, *mf_pointer_y, *mf_pointer_z;

    // We are only deposing one quantity on each face
    int ncomp  = 1;
    int nghost = mf_to_be_filled[0].nGrow();

    if (OnSameGrids(lev, mf_to_be_filled[0]))
    {
        // If we are already working with the internal mf defined on the 
        // particle_box_array, then we just work with this.
        mf_pointer_x = &mf_to_be_filled[0];
        mf_pointer_y = &mf_to_be_filled[1];
        mf_pointer_z = &mf_to_be_filled[2];
    }
    else
    {
        // If mf_to_be_filled is not defined on the particle_box_array, then we need 
        // to make a temporary here and copy into mf_to_be_filled at the end.

        IntVect nodal(1,1,1);

        mf_pointer_x = new MultiFab(m_gdb->ParticleBoxArray(lev), ncomp, nghost, 
				    m_gdb->ParticleDistributionMap(lev), Fab_allocate, nodal);
        mf_pointer_y = new MultiFab(m_gdb->ParticleBoxArray(lev), ncomp, nghost, 
				    m_gdb->ParticleDistributionMap(lev), Fab_allocate, nodal);
        mf_pointer_z = new MultiFab(m_gdb->ParticleBoxArray(lev), ncomp, nghost, 
				    m_gdb->ParticleDistributionMap(lev), Fab_allocate, nodal);
    }

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = m_gdb->Geom(lev);
    const BoxArray& ba          = mf_pointer_x->boxArray();
    const Real*     dx          = gm.CellSize();
    const PMap&     pmap        = m_particles[lev];

    // Setting the current to 0 before depositing the charge
    for (MFIter mfi(*mf_pointer_x); mfi.isValid(); ++mfi) (*mf_pointer_x)[mfi].setVal(0);
    for (MFIter mfi(*mf_pointer_y); mfi.isValid(); ++mfi) (*mf_pointer_y)[mfi].setVal(0);
    for (MFIter mfi(*mf_pointer_z); mfi.isValid(); ++mfi) (*mf_pointer_z)[mfi].setVal(0);
   
    // Charge
    Real q = 1.;

    // Loop over the grids containing particles
    for (auto& kv : pmap)
    {
        const  int  pgr = kv.first;
        const PBox& pbx = kv.second;

	Array<Real> xp, yp, zp, wp, uxp, uyp, uzp, gip;
        FArrayBox&  fabx = (*mf_pointer_x)[pgr];
        FArrayBox&  faby = (*mf_pointer_y)[pgr];
        FArrayBox&  fabz = (*mf_pointer_z)[pgr];

	 xp.reserve( pbx.size() );
	 yp.reserve( pbx.size() );
	 zp.reserve( pbx.size() );
	 wp.reserve( pbx.size() );
	uxp.reserve( pbx.size() );
	uyp.reserve( pbx.size() );
	uzp.reserve( pbx.size() );
	gip.reserve( pbx.size() );

	const Box & bx = ba[pgr];
	RealBox grid_box = RealBox( bx, dx, gm.ProbLo() );
	const Real* xyzmin = grid_box.lo();
	long nx = bx.length(0)-1, ny = bx.length(1)-1, nz = bx.length(2)-1; 
	long ng = mf_pointer_x->nGrow();
	
	// Loop over particles in that box (to change array layout)
	long np = 0;
        for (const auto& p : pbx)
        {
            if (p.m_id <= 0) {
	      continue;
	    }
	    ++np;

            // (x,y,z) position
	    xp.push_back( p.m_pos[0] );
	    yp.push_back( p.m_pos[1] );
	    zp.push_back( p.m_pos[2] );

            // weights
 	    wp.push_back( 1. ); 

 	    uxp.push_back( 1. ); 
 	    uyp.push_back( 1. ); 
 	    uzp.push_back( 1. ); 

            // gaminv 
 	    gip.push_back( 1. ); 
        }

#if 0
	depose_jxjyjz_vecHVv2_1_1_1(fabx.dataPtr(), faby.dataPtr(), fabz.dataPtr(),
                                    &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
                                    uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), 
                                    gip.dataPtr(), wp.dataPtr(), &q, 
                                    &xyzmin[0], &xyzmin[1], &xyzmin[2], 
                                    &dt, &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				    &ng, &ng, &ng);
#endif

	depose_jxjyjz_scalar_1_1_1( fabx.dataPtr(), faby.dataPtr(), fabz.dataPtr(),
                                    &np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(), 
                                    uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), 
                                    gip.dataPtr(), wp.dataPtr(), &q, 
                                    &xyzmin[0], &xyzmin[1], &xyzmin[2], 
                                    &dt, &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				    &ng, &ng, &ng);

    }

    mf_pointer_x->SumBoundary(gm.periodicity());
    mf_pointer_y->SumBoundary(gm.periodicity());
    mf_pointer_z->SumBoundary(gm.periodicity());

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer_x != &mf_to_be_filled[0])
    {
        mf_to_be_filled[0].copy(*mf_pointer_x,0,0,ncomp);
	delete mf_pointer_x;

        mf_to_be_filled[0].copy(*mf_pointer_y,0,0,ncomp);
	delete mf_pointer_y;

        mf_to_be_filled[0].copy(*mf_pointer_z,0,0,ncomp);
	delete mf_pointer_z;
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::ChargeDeposition time: " << stoptime << '\n';
        }
    }
}

