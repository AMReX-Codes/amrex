#ifndef AMREX_EB2_IF_BASE_H_
#define AMREX_EB2_IF_BASE_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <type_traits>

namespace amrex {

struct GPUable {};

namespace EB2 {

template <class D, class Enable = void> struct IsGPUable : std::false_type {};

template <class D>
struct IsGPUable<D, typename std::enable_if<std::is_base_of<GPUable,D>::value>::type>
    : std::true_type {};

}
}

#endif
