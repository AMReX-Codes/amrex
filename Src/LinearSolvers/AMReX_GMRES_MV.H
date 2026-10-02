#ifndef AMREX_GMRES_MV_H_
#define AMREX_GMRES_MV_H_

#include <AMReX_GMRES.H>
#include <AMReX_KrylovMV.H>

namespace amrex {

/**
 * \file AMReX_GMRES_MV.H
 *
 * GMRES front-end specialized for AlgVector/SpMatrix combinations.
 */

//! GMRES for SpMatrix<T> systems. See KrylovMV.
template <typename T>
using GMRES_MV = KrylovMV<GMRES, T>;

}
#endif
