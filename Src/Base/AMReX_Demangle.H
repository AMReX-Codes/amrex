#ifndef AMREX_DEMANGLE_H_
#define AMREX_DEMANGLE_H_
#include <AMReX_Config.H>

#include <string>

#if __has_include(<cxxabi.h>)
# include <cxxabi.h>
# include <cstdlib>
# define AMREX_USE_CXXABI
#endif

namespace amrex {

/**
 * \brief Demangle C++ name
 *
 * Demange C++ name if possible. For example
 \verbatim
     amrex::Box box;
     std::cout << amrex::demangle(typeid(box).name());
 \endverbatim
 * Demangling turns "N5amrex3BoxE" into "amrex::Box".
 */
inline std::string demangle (const char* name)
{
#ifdef AMREX_USE_CXXABI
    int status;
    char* p = abi::__cxa_demangle(name, nullptr, nullptr, &status);
    if (p) {
        std::string s(p);
        std::free(p);
        return s;
    } else
#endif
    {
        return name;
    }
}

}

#endif
