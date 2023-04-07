#include "verbosity.h"
#include "AMReX.H"

extern "C" int verbosity () {

    return amrex::Verbose();

}
