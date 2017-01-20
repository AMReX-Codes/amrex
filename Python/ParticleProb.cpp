
//
// Each problem must have its own version of SingleSpeciesContainer::InitData()
// to initialize the particle data on this level
//

#include <cmath>

#include <BLProfiler.H>

#include <ParticleContainer.H>
#include <WarpXConst.H>

void
SingleSpeciesContainer::InitData()
{
    charge = -PhysConst::q_e;
    mass = PhysConst::m_e;
}
