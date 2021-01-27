#ifndef AMREX_PARTICLECONTAINERBASE_H_
#define AMREX_PARTICLECONTAINERBASE_H_
#include <AMReX_Config.H>

#include <string>

#include <AMReX_INT.H>

namespace amrex {

class ParticleContainerBase {

public:

    const std::string& Version ();
    const std::string& DataPrefix ();
    int MaxReaders ();
    Long MaxParticlesPerRead ();
};

} // namespace amrex

#endif // AMREX_PARTICLECONTAINERBASE_H_
