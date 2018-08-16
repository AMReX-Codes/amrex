//#include <string>
//#include <limits>
//#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif
//#include <AMReX_AmrCore.H>
//#include <AMReX_FluxRegister.H>
//#include <AMReX_AmrMesh.H>
//#include <AMReX_BCRec.H>

extern "C"
{
    // utility to skip to next line in Header
    static void GotoNextLine (std::istream& is);
    std::string restart_chkfile = "";
    std::string chk_file {"chk"};
};
