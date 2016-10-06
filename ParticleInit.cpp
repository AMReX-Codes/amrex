
#include <BLProfiler.H>

#include <ParticleContainer.H>
#include <PICSAR_f.H>

void
MyParticleContainer::Init(MultiFab& dummy_mf)
{
    BL_PROFILE("MyPC::Init");

    Real x_off = 0.5;
    Real y_off = 0.5;
    Real z_off = 0.5;

    Real charge = 1.;

    // Initialize one particle per cell with local position (x_off,y_off,z_off) relative to nondimensional cell size [0:1]
    InitOnePerCell(x_off,y_off,z_off,charge,dummy_mf);

    // Randomly initialize "num_particles" number of particles, each with charge "charge"
    // bool serialize = false;
    // int iseed   = 10;
    // MyPC->InitRandom(num_particles,iseed,charge,serialize);

    int             lev         = 0; 
    PMap&           pmap        = m_particles[lev];

    // 1D Arrays of particle attributes
    Array<Real> uxp, uyp, uzp, gaminv;

    for (auto& kv : pmap) 
    {
	PBox& pbx = kv.second;

	uxp.resize( pbx.size() );
	uyp.resize( pbx.size() );
	uzp.resize( pbx.size() );

	// Loop over particles in that box (to change array layout)
	for (const auto& p : pbx)
        {
            BL_ASSERT(p.m_id > 0);
	    uxp.push_back( p.m_data[0] );
	    uyp.push_back( p.m_data[1] );
	    uzp.push_back( p.m_data[2] );
        }
    }
}
