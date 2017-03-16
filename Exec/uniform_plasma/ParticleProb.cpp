
//
// Each problem must have its own version of PhysicalParticleContainer::InitData()
// to initialize the particle data.  It must also initialize charge and mass.
//

#include <cmath>

#include <AMReX_BLProfiler.H>

#include <ParticleContainer.H>
#include <WarpXConst.H>
#include <random>

using namespace amrex;

void
PhysicalParticleContainer::InitData()
{
    InitNPerCell();
}
