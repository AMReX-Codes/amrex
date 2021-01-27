#ifndef AMREX_PARTICLECONTAINERBASE_H_
#define AMREX_PARTICLECONTAINERBASE_H_
#include <AMReX_Config.H>

#include <string>

#include <AMReX_INT.H>

namespace amrex {

class ParticleContainerBase {

public:

    static const std::string& Version ();
    static const std::string& DataPrefix ();
    static int MaxReaders ();
    static Long MaxParticlesPerRead ();
};

} // namespace amrex

#endif // AMREX_PARTICLECONTAINERBASE_H_
