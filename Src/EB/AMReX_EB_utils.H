#ifndef AMREX_EB_UTILS_H_
#define AMREX_EB_UTILS_H_

#include <AMReX.H>

namespace amrex {

#if (AMREX_SPACEDIM > 1)
    void single_level_redistribute ( int lev, MultiFab& div_tmp_in, MultiFab& div_out,
                                     int div_comp, int ncomp, const Vector<Geometry>& geom);

    void single_level_weighted_redistribute ( int lev, MultiFab& div_tmp_in, MultiFab& div_out, const MultiFab& weights, 
                                              int div_comp, int ncomp, const Vector<Geometry>& geom);
#endif

    void FillEBNormals(MultiFab & normals, const EBFArrayBoxFactory & eb_factory,
                       const Geometry & geom);
}

#endif
