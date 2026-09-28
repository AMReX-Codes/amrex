#ifndef AMREX_BICGSTAB_MV_H_
#define AMREX_BICGSTAB_MV_H_

#include <AMReX_BiCGStab.H>
#include <AMReX_KrylovMV.H>

namespace amrex {

/**
 * \file AMReX_BiCGStab_MV.H
 *
 * BiCGStab front-end specialized for AlgVector/SpMatrix combinations.
 */

//! BiCGStab for SpMatrix<T> systems. See KrylovMV.
template <typename T>
using BiCGStab_MV = KrylovMV<BiCGStab, T>;

}
#endif
