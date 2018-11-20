#include <AMReX_CudaLaunch.H>
#include <AMReX_FabArrayBase.H>
#include <AMReX_LayoutData.H>

namespace amrex {

namespace Cuda {

#ifdef AMREX_USE_CUDA
void getGridSize (FabArrayBase const& fa, int ngrow, LayoutData<GridSize>& gs, int& ntotblocks)
{
    gs = LayoutData<GridSize>(fa.boxArray(),fa.DistributionMap());
    ntotblocks = 0;
    for (MFIter mfi(gs); mfi.isValid(); ++mfi) {
        const auto& bx = amrex::grow(mfi.validbox(),ngrow);
        auto ec = ExecutionConfig(bx);
        gs[mfi].numBlocks = ec.numBlocks.x;
        gs[mfi].numThreads = ec.numThreads.x;
        gs[mfi].globalBlockId = ntotblocks;
        ntotblocks += ec.numBlocks.x;
    }
}
#endif

// Return intersection of the cell for this thread and the entire domain.
// If more threads are assigned than mesh cells in the domain, intersection will return an empty box.
// If box is empty, skip the work in the MFIter loop for that thread.
// If no CUDA, return the entire box. 
AMREX_GPU_HOST_DEVICE
Box getThreadBox (const Box& bx)
{
#if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
    long begin, junk;
    getThreadIndex(begin, junk, bx.numPts());
    auto len = bx.length3d();
    long k = begin / (len[0]*len[1]);
    long j = (begin - k*(len[0]*len[1])) / len[0];
    long i = (begin - k*(len[0]*len[1])) - j*len[0];
    IntVect iv{AMREX_D_DECL(static_cast<int>(i),
                            static_cast<int>(j),
                            static_cast<int>(k))};
    iv += bx.smallEnd();
    return (bx & Box(iv,iv,bx.type()));
#else
    return bx;
#endif
    

#if 0
#if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
    IntVect iv{AMREX_D_DECL(static_cast<int>(threadIdx.x + blockDim.x*(blockIdx.x)),
        AMREX_CUDA_Y_STRIDE*static_cast<int>(threadIdx.y + blockDim.y*(blockIdx.y)),
        AMREX_CUDA_Z_STRIDE*static_cast<int>(threadIdx.z + blockDim.z*(blockIdx.z)))};
    iv += bx.smallEnd();
    IntVect iv2 = iv;
#if (AMREX_SPACEDIM >= 2)
    iv2[1] += (AMREX_CUDA_Y_STRIDE-1);
#endif
#if (AMREX_SPACEDIM == 3)
    iv2[2] += (AMREX_CUDA_Z_STRIDE-1);
#endif
    return (bx & Box(iv,iv2,bx.type()));
#else
    return bx;
#endif
#endif
}

// // Extension of getThreadBox that accounts for change of box type.
// // If growing, add extra index on the big edge.
// // If shrinking, remove extra index from the big edge.
// // Any empty or broken boxes should be ignored by the GPU code threadBox.ok() check.
// AMREX_GPU_HOST_DEVICE
// Box getThreadBox (const Box& bx, const IntVect& typ)
// {
// #if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
//      Box threadBox; 
//      threadBox = getThreadBox(bx);

//      if (threadBox.ok())
//      {
//          const IntVect& Big = bx.bigEnd();
//          for (int d=0; d<AMREX_SPACEDIM; ++d) {
//              if (bx.type()[d] < typ[d]) {        // Cell to Nodal (0 < 1)
//                  if (threadBox.bigEnd(d) == Big[d]) {
//                      threadBox.growHi(d, 1);     // Box on the big edge gets the extra work.
//                  }
//              }
//              else if (bx.type()[d] > typ[d]) {   // Nodal to Cell (1 > 0)
//                  if (threadBox.bigEnd(d) == Big[d]) {
//                      threadBox.growHi(d, -1);    // Thread on the edge looses the nodal work. 
//                  }
//              }
//          }
//      }
// /*
//      // GPU output for debugging correctness by hand.
//      if (!(threadBox.ok()))
//      { 
//         IntVect small = threadBox.smallEnd();
//         IntVect big   = threadBox.bigEnd();
//         IntVect type  = threadBox.type();
//         printf(" TBC -> (%i, %i, %i):(%i, %i %i) = (%i, %i, %i), (%i, %i, %i), (%i, %i, %i)\n", threadIdx.x, threadIdx.y, threadIdx.z, blockIdx.x, blockIdx.y, blockIdx.z, small[0], small[1], small[2], big[0], big[1], big[2], type[0], type[1], type[2] );
//      }
// */
//      return threadBox;
// #else
//      Box threadBox(bx);
//      return threadBox.convert(typ);
// #endif
// }


// Version of threadBox that also breaks threads across by components.
// Each call returns the appropriate box and component for the given thread.
// Components are laid out in the 0-direction.
// e.g. [            0 <= x <  box.length(0)  ]---> scomp
//      [    bx.length <= x < 2*bx.length(0)  ]---> scomp + 1
//      [ 2* bx.length <= x < 3*bx.length(0)  ]---> scomp + 2
AMREX_GPU_HOST_DEVICE
ComponentBox
getThreadComponentBox (const Box& bx, int ncomp)
{
#if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
     // Don't use getThreadBox because these will initially be outside bx, which represents the
     // component of each thread. So, don't want to do (bx & threadBox) yet until that information
     // is extracted.
    AMREX_D_TERM(long i0 = threadIdx.x + blockDim.x*(blockIdx.x);,
                 int  i1 = threadIdx.y + blockDim.y*(blockIdx.y);,
                 int  i2 = threadIdx.z + blockDim.z*(blockIdx.z));

    long compDim = static_cast<long>(blockDim.x*gridDim.x)/static_cast<long>(ncomp);
    long quot = i0 / compDim;
    long rem = i0 - quot*compDim;
    int icomp = quot;
    IntVect iv{AMREX_D_DECL(static_cast<int>(rem), i1, i2)};

    ComponentBox r{Box(iv,iv,bx.type()), icomp, 1};

    if (icomp >= ncomp) {
        r.box = Box();
    } else {
        r.box += bx.smallEnd();
        r.box &= bx;
    }

    return r;
#else
    // On CPUs, ncomp don't change.
    return {bx, 0, ncomp};
#endif
}


// If on CUDA, return index of particle being calculated locally.
// If index is over N, return size of 0 to skip loop.
// If not on CUDA, return values to run entire loop. (1 to N) 
AMREX_GPU_HOST_DEVICE
void getThreadIndex (long &index, long &size, const long num_particles)
{
#if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
     index = blockDim.x*blockIdx.x + threadIdx.x;
     size  = (index > num_particles) ? 0 : 1;
#else
     index = 0;
     size = num_particles;
#endif
}

}  // namespace Cuda
}  // namespace amrex
