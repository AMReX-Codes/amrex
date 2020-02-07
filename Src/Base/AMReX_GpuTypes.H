#ifndef AMREX_GPU_TYPES_H_
#define AMREX_GPU_TYPES_H_

#ifdef AMREX_USE_GPU

namespace amrex {
#ifdef AMREX_USE_DPCPP
    struct dim3 {
        unsigned int x = 1;
        unsigned int y = 1;
        unsigned int z = 1;
        dim3 () = default;
        constexpr dim3 (int x_, int y_=1, int z_=1) : x(x_),y(y_),z(z_) {}
    };
#else
    using dim3 = ::dim3;
#endif
}

#endif

#endif
