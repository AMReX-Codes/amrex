#include <AMReX_CudaLaunch.H>

namespace amrex {

// Return intersection of the cell for this thread and the entire domain.
// If more threads are assigned than mesh cells in the domain, intersection will return an empty box.
// If box is empty, skip the work in the MFIter loop for that thread.
// If no CUDA, return the entire box. 
AMREX_CUDA_HOST_DEVICE
Box getThreadBox (const Box& bx)
{
#if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
    IntVect iv{AMREX_D_DECL(static_cast<int>(threadIdx.x + blockDim.x*(blockIdx.x)),
                            static_cast<int>(threadIdx.y + blockDim.y*(blockIdx.y)),
                            static_cast<int>(threadIdx.z + blockDim.z*(blockIdx.z)))};
    iv += bx.smallEnd();
    return (bx & Box(iv,iv,bx.type()));
#else
    return bx;
#endif
}


// Extension of getThreadBox that accounts for change of box type.
// Assumes current number of threads & blocks is enough for one cell per nodal box.
// Simply recalculate total  
// Any empty or broken boxes should be ignored by the GPU code threadBox.ok() check.
AMREX_CUDA_HOST_DEVICE
Box getThreadBox (const Box& bx, const IntVect& typ)
{
#if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
    IntVect iv{AMREX_D_DECL(static_cast<int>(threadIdx.x + blockDim.x*(blockIdx.x)),
                            static_cast<int>(threadIdx.y + blockDim.y*(blockIdx.y)),
                            static_cast<int>(threadIdx.z + blockDim.z*(blockIdx.z)))};
    iv += bx.smallEnd();
    return bx.convert(typ) & Box(iv,iv,typ);
#else
    Box threadBox(bx);
    return threadBox.convert(typ);
#endif
}


// Version of threadBox that also breaks threads across by components.
// Each call returns the appropriate box and component for the given thread.
// Components are laid out in the 0-direction.
// e.g. [            0 <= x <  box.length(0)  ]---> scomp
//      [    bx.length <= x < 2*bx.length(0)  ]---> scomp + 1
//      [ 2* bx.length <= x < 3*bx.length(0)  ]---> scomp + 2
AMREX_CUDA_HOST_DEVICE
void getThreadComponentBox (const Box& bx, Box& threadBox, int& scomp, int& ncomp)
{
#if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
     // Don't use getThreadBox because these will initially be outside bx, which represents the
     // component of each thread. So, don't want to do (bx & threadBox) yet until that information
     // is extracted.
    IntVect iv{AMREX_D_DECL(static_cast<int>(threadIdx.x + blockDim.x*(blockIdx.x)),
                            static_cast<int>(threadIdx.y + blockDim.y*(blockIdx.y)),
                            static_cast<int>(threadIdx.z + blockDim.z*(blockIdx.z)))};
    threadBox = Box(iv,iv,bx.type());

    int dcomp =  iv[0] / bx.length(0);
    if (dcomp >= ncomp) {
        threadBox = Box();
    } else {
        scomp += dcomp;
        threadBox.shift(0, -(bx.length(0)*dcomp));  // Shift x-index back into the box.
        threadBox += bx.smallEnd();
        threadBox &= bx;
    }
    ncomp = 1;  // On GPU, each thread works on one component. (Update later for very big boxes.)
#else
    threadBox = bx;
    // On CPUs, scomp and ncomp don't change.
#endif
}


// If on CUDA, return index of particle being calculated locally.
// If index is over N, return size of 0 to skip loop.
// If not on CUDA, return values to run entire loop. (1 to N) 
AMREX_CUDA_HOST_DEVICE
void getThreadIndex (int &index, int &size, const int num_particles)
{
#if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
     index = blockDim.x*blockIdx.x + threadIdx.x;
     size  = (index > num_particles) ? 0 : 1;
#else
     index = 1;
     size = num_particles;
#endif
}

}  // namespace amrex
