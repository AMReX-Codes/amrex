#ifndef AMREX_PARTICLEASYNCIO_H_
#define AMREX_PARTICLEASYNCIO_H_

#include <AMReX_TypeTraits.H>

namespace amrex
{

template <class PC, EnableIf_t<IsParticleContainer<PC>::value, int> foo = 0>
void
writeParticlePlotFileAsync (PC const& pc, int lev_min, int lev_max, int nGrow)
{
    BL_PROFILE("writeParticlePlotFileAsync");
}

}

#endif
