#ifndef AMREX_ML_ABECLAPLACIAN_H_
#define AMREX_ML_ABECLAPLACIAN_H_

#include <AMReX_MLLinOp.H>

namespace amrex {

class MLABecLaplacian
    : public MLLinOp
{
public:

    MLABecLaplacian ();
    ~MLABecLaplacian ();

};

}

#endif
