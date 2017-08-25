#ifndef _EBRedist_H
#define _EBRedist_H

#include <AMReX_REAL.H>

#ifndef BL_LANG_FORT

typedef std::array<int,3>  intDIM;

struct NbrSten
{
    amrex::Real data[3][3][3];
    intDIM iv;
};

#endif

#endif
