#include <AMReX_AmrParticleMeshInSituBridge.H>

#ifdef AMREX_PARTICLES
namespace amrex {
// explicit template specialization for use with AmrParticleContainer class
template<> class AmrParticleMeshInSituBridge<3,0,0,0>;
}
#endif