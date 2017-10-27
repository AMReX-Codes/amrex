#ifndef MY_PROB_PAR_H_
#define MY_PROB_PAR_H_

#include <AMReX_MLLinOp.H>

namespace prob {

    extern amrex::Real a;
    extern amrex::Real b;
    extern amrex::Real sigma;
    extern amrex::Real w;

    extern amrex::MLLinOp::BCType bc_type;
}

#endif
