#include <AMReX_AmrParticleInSituBridge.H>

#ifdef AMREX_PARTICLES
namespace amrex {
// explicit template specialization for use with AmrParticleContainer class
template<> class AmrParticleInSituBridge<3,0,0,0>;
}
#endif
