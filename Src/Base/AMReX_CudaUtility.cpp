
#include <AMReX_CudaUtility.H>

namespace amrex {

#if AMREX_USE_CUDA

std::ostream&
operator<< (std::ostream& os, const dim3& d)
{
    os << '(' << d.x << ' ' << d.y << ' ' << d.z << ')';
    return os;
}

#endif

}
