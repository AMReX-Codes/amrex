#ifndef AMREX_EB_UTILS_H_
#define AMREX_EB_UTILS_H_

#include <AMReX.H>


namespace amrex {

    void FillEBNormals(MultiFab & normals, const EBFArrayBoxFactory & eb_factory,
                       const Geometry & geom);


}


#endif
