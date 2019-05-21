#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

AMREX_GPU_GLOBAL
void copy (int size, amrex::Box* bx, amrex::Dim3* offset,
           amrex::Array4<Real>* src, amrex::Array4<Real>* dst, 
           int scomp, int dcomp, int ncomp)
{
    for (int l=0; l<size; ++l)
    {
        int ncells = bx[l].numPts();
        const auto lo = amrex::lbound(bx[l]);
        const auto len = amrex::length(bx[l]);
 
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 icell < ncells; icell += stride) {
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            for (int n = 0; n < ncomp; ++n) {
                (dst[l])(i,j,k,dcomp+n) = (src[l])(i+offset[l].x,j+offset[l].y,k+offset[l].z,scomp+n); 
            }
        }
    }
}

AMREX_GPU_GLOBAL
void copy (amrex::Dim3 lo, amrex::Dim3 len, int ncells,
           amrex::Dim3 offset, amrex::Array4<Real> src, amrex::Array4<Real> dst, 
           int scomp, int dcomp, int ncomp)
{
    for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
        int k =  icell /   (len.x*len.y);
        int j = (icell - k*(len.x*len.y)) /   len.x;
        int i = (icell - k*(len.x*len.y)) - j*len.x;
        i += lo.x;
        j += lo.y;
        k += lo.z;
        for (int n = 0; n < ncomp; ++n) {
            dst(i,j,k,dcomp+n) = src(i+offset.x,j+offset.y,k+offset.z,scomp+n); 
        }
    }
}

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {

        // AMREX_SPACEDIM: number of dimensions
        int n_cell, max_grid_size;
        Vector<int> is_periodic(AMREX_SPACEDIM,1);  // periodic in all direction by default

        // inputs parameters
        {
            // ParmParse is way of reading inputs from the inputs file
            ParmParse pp;

            // We need to get n_cell from the inputs file - this is the number of cells on each side of 
            //   a square (or cubic) domain.
            pp.get("n_cell",n_cell);

            // The domain is broken into boxes of size max_grid_size
            pp.get("max_grid_size",max_grid_size);
        }

        // make BoxArray and Geometry
        BoxArray ba;
        {
            IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
            IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
            Box domain(dom_lo, dom_hi);

            // Initialize the boxarray "ba" from the single box "bx"
            ba.define(domain);
            // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
            ba.maxSize(max_grid_size);
        }

        // Nghost = number of ghost cells for each array 
        int Nghost = 1;
    
        // Ncomp = number of components for each array
        int Ncomp  = 1;
  
        // How Boxes are distrubuted among MPI processes
        DistributionMapping dm(ba);

        // Malloc value for setval testing.
        Real* val;
        cudaMallocManaged(&val, sizeof(Real));

        // Create the MultiFab and touch the data.
        // Ensures the data in on the GPU for all further testing.
        MultiFab x(ba, dm, Ncomp, Nghost);
        MultiFab y(ba, dm, Ncomp, Nghost);
        x.setVal(0.0);
        y.setVal(0.0);

        amrex::Print() << "Testing on " << n_cell << "^3 boxes with max grid size " << max_grid_size << std::endl 
                       << "Number of boxes per MultiFab: " << x.size() << std::endl << std::endl;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Launch without graphs
/*
        {
            x.setVal(4.5);
            y.setVal(5.3);

            BL_PROFILE("HtoD Copy Clean-up");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                Array4<Real> src = x.array(mfi);
                Array4<Real> dst = y.array(mfi);
                Dim3 offset {0,0,0};
                int dcomp = 0;
                int scomp = 0;

                AMREX_HOST_DEVICE_FOR_4D (bx, Ncomp, i, j, k, n,
                {
                    dst(i,j,k,dcomp+n) = src(i+offset.x,j+offset.y,k+offset.z,scomp+n); 
                });
            }

            amrex::Print() << "Clean up sum = " << y.sum() << "; Expected value = " << x.sum() << std::endl;
        }

// ---------------------------------------

        {
            x.setVal(4.0);
            y.setVal(5.0);

            BL_PROFILE("Lambda Copy - Array4");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                Array4<Real> src = x.array(mfi);
                Array4<Real> dst = y.array(mfi);
                Dim3 offset {0,0,0};
                int dcomp = 0;
                int scomp = 0;

                AMREX_HOST_DEVICE_FOR_4D ( bx, Ncomp, i, j, k, n,
                {
                    dst(i,j,k,dcomp+n) = src(i+offset.x,j+offset.y,k+offset.z,scomp+n); 
                });
            }

            amrex::Print() << "No Graph Lambda sum = " << y.sum() << "; Expected value = " << x.sum() << std::endl;
        }
*/
// ---------------------------------------

        {
            x.setVal(0.004);
            y.setVal(0.0027);

            BL_PROFILE_VAR("First", first);

            const int numStreams = 16;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                Gpu::Device::setStreamIndex(mfi.LocalIndex() % numStreams);
                const Box bx = mfi.validbox();

                const auto src = x.array(mfi);
                const auto dst = y.array(mfi);

                const int ncells = bx.numPts();
                const auto lo  = amrex::lbound(bx);
                const auto len = amrex::length(bx);
                const Dim3 offset = {0,0,0};

                const auto ec = Cuda::ExecutionConfig(bx.numPts());
                AMREX_CUDA_LAUNCH_GLOBAL(ec, copy,
                                         lo, len, ncells,
                                         offset, src, dst,
                                         0, 0, Ncomp);
            }
            BL_PROFILE_VAR_STOP(first);

            amrex::Print() << "First = " << y.sum() << "; Expected value = " << x.sum() << std::endl;
        }


        {
            x.setVal(0.75);
            y.setVal(0.25);

            BL_PROFILE_VAR("STANDARD: Function Copy - Array4 (16)", sixteen);

            const int numStreams = 16;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                Gpu::Device::setStreamIndex(mfi.LocalIndex() % numStreams);
                const Box bx = mfi.validbox();

                const auto src = x.array(mfi);
                const auto dst = y.array(mfi);

                const int ncells = bx.numPts();
                const auto lo  = amrex::lbound(bx);
                const auto len = amrex::length(bx);
                const Dim3 offset = {0,0,0};

                const auto ec = Cuda::ExecutionConfig(bx.numPts());
                AMREX_CUDA_LAUNCH_GLOBAL(ec, copy,
                                         lo, len, ncells,
                                         offset, src, dst,
                                         0, 0, Ncomp);
            }
            BL_PROFILE_VAR_STOP(sixteen);

            amrex::Print() << "Standard Function Array sum = " << y.sum() << "; Expected value = " << x.sum() << std::endl;
        }

        {
            x.setVal(0.74);
            y.setVal(0.47);

            BL_PROFILE_VAR("STANDARD: Function Copy - Array4 (1)", one);

            const int numStreams = 1;

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                Gpu::Device::setStreamIndex(mfi.LocalIndex() % numStreams);
                const Box bx = mfi.validbox();

                const auto src = x.array(mfi);
                const auto dst = y.array(mfi);

                const int ncells = bx.numPts();
                const auto lo  = amrex::lbound(bx);
                const auto len = amrex::length(bx);
                const Dim3 offset = {0,0,0};

                const auto ec = Cuda::ExecutionConfig(bx.numPts());
                AMREX_CUDA_LAUNCH_GLOBAL(ec, copy,
                                         lo, len, ncells,
                                         offset, src, dst,
                                         0, 0, Ncomp);
            }
            BL_PROFILE_VAR_STOP(one);


            amrex::Print() << "Standard Function Array sum = " << y.sum() << "; Expected value = " << x.sum() << std::endl;
        }


        {
            x.setVal(0.0001);
            y.setVal(0.0);

            BL_PROFILE_VAR("FUSED: Function Copy - Array4", fusedfull);

            BL_PROFILE_VAR("FUSED: ALLOC", fusedalloc);

            int lsize = x.local_size();
            Box* bx_h =           static_cast<Box*>(std::malloc(lsize*sizeof(Box)));
            Dim3* offset_h =      static_cast<Dim3*>(std::malloc(lsize*sizeof(Dim3)));
            Array4<Real>* src_h = static_cast<Array4<Real>*>(std::malloc(lsize*sizeof(Array4<Real>)));
            Array4<Real>* dst_h = static_cast<Array4<Real>*>(std::malloc(lsize*sizeof(Array4<Real>)));

            Box* bx_d;
            Dim3* offset_d;
            Array4<Real>* src_d, *dst_d;

            AMREX_GPU_SAFE_CALL(cudaMalloc(&bx_d,     lsize*sizeof(Box)));
            AMREX_GPU_SAFE_CALL(cudaMalloc(&offset_d, lsize*sizeof(Dim3)));
            AMREX_GPU_SAFE_CALL(cudaMalloc(&src_d,    lsize*sizeof(Array4<Real>)));
            AMREX_GPU_SAFE_CALL(cudaMalloc(&dst_d,    lsize*sizeof(Array4<Real>)));

            BL_PROFILE_VAR_STOP(fusedalloc);

            BL_PROFILE_VAR("FUSED: SETUP", fusedsetup);

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                int idx = mfi.LocalIndex(); 

                bx_h[idx] = mfi.validbox();
                offset_h[idx] = {0,0,0};
                src_h[idx] = x.array(mfi);
                dst_h[idx] = y.array(mfi);
//                std::memcpy(&src_h[idx], &src, sizeof(Array4<Real>));
//                std::memcpy(&dst_h[idx], &dst, sizeof(Array4<Real>));
            }
            BL_PROFILE_VAR_STOP(fusedsetup);
            BL_PROFILE_VAR("FUSED: MEMCPY", fusedmem);

            AMREX_GPU_SAFE_CALL(cudaMemcpy(bx_d,     bx_h,     lsize*sizeof(Box),          cudaMemcpyHostToDevice));
            AMREX_GPU_SAFE_CALL(cudaMemcpy(offset_d, offset_h, lsize*sizeof(Dim3),         cudaMemcpyHostToDevice));
            AMREX_GPU_SAFE_CALL(cudaMemcpy(src_d,    src_h,    lsize*sizeof(Array4<Real>), cudaMemcpyHostToDevice));
            AMREX_GPU_SAFE_CALL(cudaMemcpy(dst_d,    dst_h,    lsize*sizeof(Array4<Real>), cudaMemcpyHostToDevice));

            amrex::Gpu::Device::synchronize();
            BL_PROFILE_VAR_STOP(fusedmem);

            BL_PROFILE_VAR("FUSED: LAUNCH", fusedl);

            const auto ec = Cuda::ExecutionConfig(bx_h[0].numPts());

            AMREX_CUDA_LAUNCH_GLOBAL(ec, copy, lsize, 
                                     bx_d, offset_d,
                                     src_d, dst_d,
                                     0, 0, Ncomp);

            amrex::Gpu::Device::synchronize();
            BL_PROFILE_VAR_STOP(fusedl);

            BL_PROFILE_VAR_STOP(fusedfull);

            amrex::Print() << "Fused Function Array sum = " << y.sum() << "; Expected value = " << x.sum() << std::endl;

            std::free(bx_h);
            std::free(offset_h);
            std::free(src_h);
            std::free(dst_h);

            cudaFree(bx_d);
            cudaFree(offset_d);
            cudaFree(src_d);
            cudaFree(dst_d);

        }
    }

    amrex::Finalize();
}
