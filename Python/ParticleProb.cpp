
//
// Each problem must have its own version of PhysicalParticleContainer::InitData()
// to initialize the particle data.  It must also initialize charge and mass.
//

#include <ParticleContainer.H>
#include <WarpXConst.H>

using namespace amrex;

void
PhysicalParticleContainer::InitData()
{
    charge = -PhysConst::q_e;
    mass = PhysConst::m_e;
}
