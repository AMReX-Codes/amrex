#ifndef AMREX_EB_UTILS_H_
#define AMREX_EB_UTILS_H_

#include <AMReX.H>

namespace amrex {

#if (AMREX_SPACEDIM > 1)

    void single_level_redistribute (MultiFab& div_tmp_in, MultiFab& div_out,
                                    int div_comp, int ncomp, const Geometry& geom);

    void single_level_weighted_redistribute (MultiFab& div_tmp_in, MultiFab& div_out, const MultiFab& weights,
                                             int div_comp, int ncomp, const Geometry& geom);

#if (__cplusplus >= 201402L)
    [[deprecated]]
#endif
    inline void single_level_redistribute (int lev, MultiFab& div_tmp_in, MultiFab& div_out,
                                           int div_comp, int ncomp, const Vector<Geometry>& geom)
    { single_level_redistribute(div_tmp_in, div_out, div_comp, ncomp, geom[lev]); }

#if (__cplusplus >= 201402L)
    [[deprecated]]
#endif
    inline void single_level_weighted_redistribute (int lev, MultiFab& div_tmp_in, MultiFab& div_out, const MultiFab& weights,
                                                    int div_comp, int ncomp, const Vector<Geometry>& geom)
    { single_level_weighted_redistribute(div_tmp_in, div_out, weights, div_comp, ncomp, geom[lev]); }
#endif

#ifndef BL_NO_FORT
    void FillEBNormals(MultiFab & normals, const EBFArrayBoxFactory & eb_factory,
                       const Geometry & geom);
#endif
}

#endif
