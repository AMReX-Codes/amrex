
//
// Each problem must have its own version of MyParticleContainer::InitData()
// to initialize the particle data on this level
//

#include <cmath>

#include <BLProfiler.H>

#include <ParticleContainer.H>
#include <PICSAR_f.H>
#include <WarpXConst.H>

void
MyParticleContainer::InitData()
{
    static_assert(false,
		  "Each problem must have its own version of MyParticleContainer::InitData()");
}
