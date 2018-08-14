#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_IntVect.H>

namespace amrex {

// Return intersection of the cell for this thread and the entire domain.
// If more threads are assigned than mesh cells in the domain, intersection will return an empty box.
// If box is empty, skip the work in the MFIter loop for that thread.
// If no CUDA, return the entire box. 
AMREX_CUDA_HOST AMREX_CUDA_DEVICE
Box getThreadBox(const Box& bx)
{
#if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
     return (bx & Box(IntVect(AMREX_D_DECL(bx.smallEnd()[0] + (threadIdx.x) + blockDim.x*(blockIdx.x),
                                           bx.smallEnd()[1] + (threadIdx.y) + blockDim.y*(blockIdx.y),
                                           bx.smallEnd()[2] + (threadIdx.z) + blockDim.z*(blockIdx.z))),
                      IntVect(AMREX_D_DECL(bx.smallEnd()[0] + (threadIdx.x) + blockDim.x*(blockIdx.x),
                                           bx.smallEnd()[1] + (threadIdx.y) + blockDim.y*(blockIdx.y),
                                           bx.smallEnd()[2] + (threadIdx.z) + blockDim.z*(blockIdx.z))),
                      bx.type()));
#else
     return bx;
#endif
}

AMREX_CUDA_HOST AMREX_CUDA_DEVICE
Box getThreadBox(const Box& bx, const IntVect& typ)
{
#if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
     Box threadBox; 
     threadBox = getThreadBox(bx);

     if (threadBox.ok())
     {
         const IntVect& Big = bx.bigEnd();
         for (int d=0; d<AMREX_SPACEDIM; ++d) {
             if (bx.type()[d] < typ[d]) {        // Cell to Nodal (0 < 1)
                 if (threadBox.bigEnd(d) == Big[d]) {
                     threadBox.growHi(d, 1);     // Box on the big edge gets the extra work.
                 }
             }
             else if (bx.type()[d] > typ[d]) {   // Nodal to Cell (1 > 0)
                 if (threadBox.bigEnd(d) == Big[d]) {
                     threadBox.growHi(d, -1);    // Thread on the edge looses the nodal work. 
                     break;
                 }
             }
         }
     }
/*
     // GPU output for debugging correctness by hand.
     if (!(threadBox.ok()))
     { 
        IntVect small = threadBox.smallEnd();
        IntVect big   = threadBox.bigEnd();
        IntVect type  = threadBox.type();
        printf(" TBC -> (%i, %i, %i):(%i, %i %i) = (%i, %i, %i), (%i, %i, %i), (%i, %i, %i)\n", threadIdx.x, threadIdx.y, threadIdx.z, blockIdx.x, blockIdx.y, blockIdx.z, small[0], small[1], small[2], big[0], big[1], big[2], type[0], type[1], type[2] );
     }
*/
     return threadBox;
#else
     Box threadBox(bx);
     return threadBox.convert(typ);
#endif
}

// If on CUDA, return index of particle being calculated locally.
// If index is over N, return size of 0 to skip loop.
// If not on CUDA, return values to run entire loop. (1 to N) 
AMREX_CUDA_HOST AMREX_CUDA_DEVICE
void getThreadIndex(int &index, int &size, const int num_particles)
{
#if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
     index = blockDim.x*blockIdx.x + threadIdx.x;
     size  = (index > num_particles) ? 0 : 1;
#else
     index = 1;
     size = num_particles;
#endif
}

}
