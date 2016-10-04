#include "ParticleContainer.H"

extern "C" {
void pxr_set_gamma(const long* np,
           const Real* up, const Real* vp, const Real* wp,
           const Real* gaminv);
}

void
MyParticleContainer::Init(MultiFab& dummy_mf)
{
    Real x_off = 0.5;
    Real y_off = 0.5;
    Real z_off = 0.5;

    Real charge = 1.;
    Real  mass  = 1.;

    // Initialize one particle per cell with local position (x_off,y_off,z_off) relative to nondimensional cell size [0:1]
    InitOnePerCell(x_off,y_off,z_off,charge,dummy_mf);

    // Randomly initialize "num_particles" number of particles, each with charge "charge"
    // bool serialize = false;
    // int iseed   = 10;
    // MyPC->InitRandom(num_particles,iseed,charge,serialize);

    int             lev         = 0; 
    const Real      strttime    = ParallelDescriptor::second();
    const PMap&     pmap        = m_particles[lev];
    const int       ngrids      = pmap.size();

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

    // Loop over boxes
    for (int j = 0; j < ngrids; j++)
    {
        const PBox& pbx = *pbxs[j];
	long np = 0;
	Array<Real> uxp, uyp, uzp, gaminv;

	// 1D Arrays of particle attributes
	uxp.reserve( pbx.size() );
	uyp.reserve( pbx.size() );
	uzp.reserve( pbx.size() );

	gaminv.reserve( pbx.size() );

	// Loop over particles in that box (to change array layout)
        for (typename PBox::const_iterator it = pbx.begin(); it < pbx.end(); ++it)
        {
            const ParticleType& p = *it;
	    
            if (p.m_id <= 0) {
	      continue;
	    }
	    ++np;
	    uxp.push_back( p.m_data[0] );
	    uyp.push_back( p.m_data[1] );
	    uzp.push_back( p.m_data[2] );

            // Initialize this to 1e20 to make sure it gets filled below
 	    gaminv.push_back( 1.e20 ); 
        }

        pxr_set_gamma(&np, uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), gaminv.dataPtr());

#if 0
        // Loop over particles in that box (to change array layout)
        int n = 0;
        for (typename PBox::const_iterator it = pbx.begin(); it < pbx.end(); ++it)
        {
            const ParticleType& p = *it;

            if (p.m_id <= 0) {
              continue;
            }
            p.m_data[ 0] = charge;
            p.m_data[10] = gaminv.dataPtr()[n];
            n++;
        }
#endif
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::Init time: " << stoptime << '\n';
        }
    }
}
