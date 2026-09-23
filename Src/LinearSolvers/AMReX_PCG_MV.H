#ifndef AMREX_PCG_MV_H_
#define AMREX_PCG_MV_H_

#include <AMReX_PCG.H>
#include <AMReX_KrylovMV.H>

namespace amrex {

/**
 * \file AMReX_PCG_MV.H
 *
 * PCG front-end specialized for AlgVector/SpMatrix combinations.
 */

//! PCG for SpMatrix<T> systems. See KrylovMV.
template <typename T>
using PCG_MV = KrylovMV<PCG, T>;

}
#endif
