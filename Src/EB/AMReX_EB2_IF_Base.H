#ifndef AMREX_EB2_IF_BASE_H_
#define AMREX_EB2_IF_BASE_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <type_traits>

/**
 * \file AMReX_EB2_IF_Base.H
 *
 * Base helpers shared by all EB2 implicit functions, including GPU-aware markers
 * and compile-time utilities to identify GPU-callable functors.
 */

namespace amrex {

//! Marker base that denotes an implicit function callable on device kernels.
struct GPUable {};

namespace EB2 {

//! Type trait that reports whether a functor derives from GPUable.
template <class D> struct IsGPUable : std::false_type {};

template <class D>
requires(std::is_base_of_v<GPUable,D>)
struct IsGPUable<D> : std::true_type {};

}
}

#endif
