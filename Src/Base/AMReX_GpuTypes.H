#ifndef AMREX_GPU_TYPES_H_
#define AMREX_GPU_TYPES_H_

#ifdef AMREX_USE_GPU

#ifdef AMREX_USE_DPCPP
#include <CL/sycl.hpp>
namespace sycl = cl::sycl;
#endif

namespace amrex {

#ifdef AMREX_USE_DPCPP

struct dim3 {
    unsigned int x = 1;
    unsigned int y = 1;
    unsigned int z = 1;
    dim3 () = default;
    constexpr dim3 (unsigned int x_, unsigned int y_=1, unsigned int z_=1) : x(x_),y(y_),z(z_) {}
};

struct gpuStream_t {
    sycl::queue* queue = nullptr;
    bool operator== (gpuStream_t const& rhs) noexcept { return queue == rhs.queue; }
};

namespace Gpu {
struct Handler {
    Handler (sycl::nd_item<1> const& a_item, void* a_local = nullptr)
        : item(a_item), local(a_local) {}
    sycl::nd_item<1> const& item;
    void* local; // DPC++ local memory
};
}

#endif

}

#endif
#endif
