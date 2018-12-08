
#include <MyKernel.H>

using namespace amrex;

AMREX_GPU_DEVICE void plusone_cudacpp (Box const& bx, FArrayBox& fab)
{
    const auto len = length(bx);  // length of box
    const auto lo  = lbound(bx);  // lower bound of box
    const auto data = fab.view(lo);  // a view starting from lo

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            // We know this is safe for simd on cpu.  So let's give compiler some help.
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                data(i,j,k) += 1.0;
            }
        }
    }
}
