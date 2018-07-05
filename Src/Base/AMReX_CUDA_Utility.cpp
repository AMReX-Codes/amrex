#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>

namespace amrex {

// Return intersection of the cell for this thread and the entire domain.
// If more threads are assigned than mesh cells in the domain, intersection will return an empty box.
// If box is empty, skip the work in the MFIter loop for that thread.
// If no CUDA, return the entire box. 
AMREX_CUDA_DEVICE
Box getThreadBox(const Box & bx)
{
#ifdef AMREX_USE_CUDA
     return (bx & Box(IntVect(AMREX_D_DECL(bx.smallEnd()[0] + (threadIdx.x) + blockDim.x*(blockIdx.x),
                                           bx.smallEnd()[1] + (threadIdx.y) + blockDim.y*(blockIdx.y),
                                           bx.smallEnd()[2] + (threadIdx.z) + blockDim.z*(blockIdx.z))),
                      IntVect(AMREX_D_DECL(bx.smallEnd()[0] + (threadIdx.x) + blockDim.x*(blockIdx.x),
                                           bx.smallEnd()[1] + (threadIdx.y) + blockDim.y*(blockIdx.y),
                                           bx.smallEnd()[2] + (threadIdx.z) + blockDim.z*(blockIdx.z)))));
#else
     return bx;
#endif
}

// If on CUDA, return index of particle being calculated locally.
// If index is over N, return size of 0 to skip loop.
// If not on CUDA, return values to run entire loop. (1 to N) 
AMREX_CUDA_DEVICE
void getThreadIndex(int &index, int &size, const int num_particles)
{
#ifdef AMREX_USE_CUDA
     index = blockDim.x * (blockIdx.x - 1) + threadIdx.x;
     size = (index > num_particles) ? 0 : 1;
#else
     index = 1;
     size = num_particles;
#endif
}

}
