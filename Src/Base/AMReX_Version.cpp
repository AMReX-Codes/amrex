//include <AMReX.H>  // contains declaration
#include <AMReX_Version.H>

#include <string>


namespace amrex
{
    std::string Version ()
    {
#ifdef AMREX_GIT_VERSION
        return std::string(AMREX_GIT_VERSION);
#else
        return std::string("Unknown");
#endif
    }
} // namespace amrex
