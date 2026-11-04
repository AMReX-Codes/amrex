#include <AMReX_FPMP.H>

#if defined(AMREX_USE_FP32MP2) || defined(AMREX_FP32MP2_PARTICLES)

#include <istream>

namespace amrex {

std::istream& operator>> (std::istream& is, cuda::fp32mp2& x)
{
    double d;
    if (is >> d) {
        x = cuda::fp32mp2(d);
    }
    return is;
}

}

#endif
