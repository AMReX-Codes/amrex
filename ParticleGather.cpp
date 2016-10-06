
#include <ParticleContainer.H>
#include <PICSAR_f.H>

void
MyParticleContainer::FieldGather(MultiFab& Ex, MultiFab& Ey, MultiFab& Ez,
                                 MultiFab& Bx, MultiFab& By, MultiFab& Bz,
                                 long order, long field_gathe_algo)
{
    int             lev         = 0; 
    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = m_gdb->Geom(lev);
    const BoxArray& ba          = Ex.boxArray();
    const Real*     dx          = gm.CellSize();

    const PMap&     pmap        = m_particles[lev];

    // Loop over the grids containing particles
    for (auto& kv : pmap)
    {
        const  int  pgr = kv.first;
        const PBox& pbx = kv.second;

	long np = 0;
	Array<Real>  xp,  yp,  zp, wp;
	Array<Real> exp, eyp, ezp;
	Array<Real> bxp, byp, bzp;

	// 1D Arrays of particle attributes
	 xp.reserve( pbx.size() );
	 yp.reserve( pbx.size() );
	 zp.reserve( pbx.size() );
	 wp.reserve( pbx.size() );
	exp.reserve( pbx.size() );
	eyp.reserve( pbx.size() );
	ezp.reserve( pbx.size() );
	bxp.reserve( pbx.size() );
	byp.reserve( pbx.size() );
	bzp.reserve( pbx.size() );

	// Data on the grid
        FArrayBox& exfab = Ex[pgr];
        FArrayBox& eyfab = Ey[pgr];
        FArrayBox& ezfab = Ez[pgr];
        FArrayBox& bxfab = Bx[pgr];
        FArrayBox& byfab = By[pgr];
        FArrayBox& bzfab = Bz[pgr];

	const Box & bx = ba[pgr];
	RealBox grid_box = RealBox( bx, dx, gm.ProbLo() );
	const Real* xyzmin = grid_box.lo();
	long nx = bx.length(0)-1, ny = bx.length(1)-1, nz = bx.length(2)-1; 
	long ng = Ex.nGrow();

        Real strt_copy = ParallelDescriptor::second();
	
	// Loop over particles in that box 
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

 	    // We put dummy values in here just to make sure these values are getting filled below.
 	    exp.push_back( 1.e20 ); 
 	    eyp.push_back( 1.e20 ); 
 	    ezp.push_back( 1.e20 ); 
 	    bxp.push_back( 1.e20 ); 
 	    byp.push_back( 1.e20 ); 
 	    bzp.push_back( 1.e20 ); 
        }

        Real end_copy = ParallelDescriptor::second() - strt_copy;

        if (ParallelDescriptor::IOProcessor()) 
            std::cout << "Time in FieldGather : Copy " << end_copy << '\n';

        bool ll4symtry          = false;
        bool l_lower_order_in_v = true;

        Real strt_gather = ParallelDescriptor::second();

        geteb3d_energy_conserving(&np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(),
                                      exp.dataPtr(),eyp.dataPtr(),ezp.dataPtr(),
                                      bxp.dataPtr(),byp.dataPtr(),bzp.dataPtr(),
				      &xyzmin[0], &xyzmin[1], &xyzmin[2],
				      &dx[0], &dx[1], &dx[2],
				      &nx, &ny, &nz, &ng, &ng, &ng, 
				      &order, &order, &order, 
				      exfab.dataPtr(), eyfab.dataPtr(), ezfab.dataPtr(),
				      bxfab.dataPtr(), byfab.dataPtr(), bzfab.dataPtr(),
				      &ll4symtry, &l_lower_order_in_v, &field_gathe_algo);

       Real end_gather = ParallelDescriptor::second() - strt_gather;
       if (ParallelDescriptor::IOProcessor()) 
           std::cout << "Time in PicsarFieldGather : Gather " << end_gather << '\n';
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::FieldGather time: " << stoptime << '\n';
        }
    }
}
