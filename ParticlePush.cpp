#include "ParticleContainer.H"

extern "C" {
void pxr_epush_v(const long* np, 
           const Real* xp, const Real* yp, const Real* zp,
	   const Real* up, const Real* vp, const Real* wp,
	   const Real* exp, const Real* eyp,const Real* ezp,
	   const Real* charge, const Real* mass, Real* dt);

void pxr_bpush_v(const long* np, 
           const Real* xp, const Real* yp, const Real* zp,
	   const Real* up, const Real* vp, const Real* wp,
	   const Real* gaminv, 
	   const Real* bxp, const Real* byp,const Real* bzp,
	   const Real* charge, const Real* mass, Real* dt);

void pxr_pushxyz(const long* np, 
           const Real* xp, const Real* yp, const Real* zp,
	   const Real* up, const Real* vp, const Real* wp,
	   const Real* gaminv, Real* dt);

void pxr_set_gamma(const long* np, 
	   const Real* up, const Real* vp, const Real* wp,
	   const Real* gaminv);
}

void
MyParticleContainer::ParticlePush(Real dt)
{
    int             lev         = 0; 
    const Real      strttime    = ParallelDescriptor::second();
    const PMap&     pmap        = m_particles[lev];
    const int       ngrids      = pmap.size();

    //
    // This is a little funky.  What in effect this'll do is force
    // each thread to work on a single (separate) grid at a time.  That
    // way no thread will step on any other.  If there's only one grid per CPU,
    // then oh well ....
    //
    // TODO: implement tiling with OpenMP in this grid loop.
    Array<int>         pgrd(ngrids);
    Array<const PBox*> pbxs(ngrids);

    int j = 0;
    for (typename PMap::const_iterator pmap_it = pmap.begin(), pmapEnd = pmap.end();
         pmap_it != pmapEnd;
         ++pmap_it, ++j)
    {
        pgrd[j] =   pmap_it->first;
        pbxs[j] = &(pmap_it->second);
    }

    Real half_dt = 0.5 * dt;

    Real mass = 1.0;

    // Loop over boxes
    for (int j = 0; j < ngrids; j++)
    {
        const PBox& pbx = *pbxs[j];
	long np = 0;
	Array<Real>  xp,  yp,  zp, wp;
	Array<Real> uxp, uyp, uzp;
	Array<Real> exp, eyp, ezp;
	Array<Real> bxp, byp, bzp;
	Array<Real> cp;
	Array<Real> gaminv;

	// 1D Arrays of particle attributes
	 xp.reserve( pbx.size() );
	 yp.reserve( pbx.size() );
	 zp.reserve( pbx.size() );
	 wp.reserve( pbx.size() );
	uxp.reserve( pbx.size() );
	uyp.reserve( pbx.size() );
	uzp.reserve( pbx.size() );
	exp.reserve( pbx.size() );
	eyp.reserve( pbx.size() );
	ezp.reserve( pbx.size() );
	bxp.reserve( pbx.size() );
	byp.reserve( pbx.size() );
	bzp.reserve( pbx.size() );

	gaminv.reserve( pbx.size() );

        Real strt_copy = ParallelDescriptor::second();
	
	// Loop over particles in that box (to change array layout)
        for (typename PBox::const_iterator it = pbx.begin(); it < pbx.end(); ++it)
        {
            const ParticleType& p = *it;
	    
            if (p.m_id <= 0) {
	      continue;
	    }
	    ++np;

            // Position
	     xp.push_back( p.m_pos[0] );
	     yp.push_back( p.m_pos[1] );
	     zp.push_back( p.m_pos[2] );

            // Charge
 	     wp.push_back( p.m_data[0]);

            // Velocity
 	    uxp.push_back( p.m_data[1] ); 
 	    uyp.push_back( p.m_data[2] ); 
 	    uzp.push_back( p.m_data[3] ); 

            // E-field
 	    exp.push_back( p.m_data[5] ); 
 	    eyp.push_back( p.m_data[6] ); 
 	    ezp.push_back( p.m_data[7] ); 

            // B-field
 	    bxp.push_back( p.m_data[8] ); 
 	    byp.push_back( p.m_data[9] ); 
 	    bzp.push_back( p.m_data[10] ); 

            // (1 / Gamma)
 	    gaminv.push_back( 1.e20 );
        }

        Real end_copy = ParallelDescriptor::second() - strt_copy;

        if (ParallelDescriptor::IOProcessor()) 
            std::cout << "Time in ParticlePush : Copy " << end_copy << '\n';

        Real strt_push = ParallelDescriptor::second();

        pxr_epush_v(&np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(),
                        uxp.dataPtr(),uyp.dataPtr(),uzp.dataPtr(),
                        exp.dataPtr(),eyp.dataPtr(),ezp.dataPtr(),
                         wp.dataPtr(),&mass,&half_dt);

        pxr_set_gamma(&np, uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), gaminv.dataPtr());

        pxr_bpush_v(&np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(),
                        uxp.dataPtr(),uyp.dataPtr(),uzp.dataPtr(),
                     gaminv.dataPtr(),
                        bxp.dataPtr(),byp.dataPtr(),bzp.dataPtr(),
                         wp.dataPtr(),&mass,&dt);

        pxr_epush_v(&np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(),
                        uxp.dataPtr(),uyp.dataPtr(),uzp.dataPtr(),
                        exp.dataPtr(),eyp.dataPtr(),ezp.dataPtr(),
                         wp.dataPtr(),&mass,&half_dt);

        pxr_set_gamma(&np, uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), gaminv.dataPtr());

        pxr_pushxyz(&np,  xp.dataPtr(), yp.dataPtr(), zp.dataPtr(),
                        uxp.dataPtr(),uyp.dataPtr(),uzp.dataPtr(),
                      gaminv.dataPtr(), &dt);

        Real end_push = ParallelDescriptor::second() - strt_push;

        if (ParallelDescriptor::IOProcessor()) 
            std::cout << "Time in PicsarPush : Push " << end_push << '\n';
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::ParticlePush time: " << stoptime << '\n';
        }
    }
}
