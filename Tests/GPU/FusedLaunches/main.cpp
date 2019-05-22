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

/*
AMREX_GPU_GLOBAL
void copy (int size, amrex::Box* bx, amrex::Dim3* offset,
           amrex::Array4<Real>* src, amrex::Array4<Real>* dst, 
           int scomp, int dcomp, int ncomp, int simul)
{
    for (int l = 0; l < size; l+=simul)
    {
        for (int m = 0; m < simul; m++)
        {

        }

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
*/

AMREX_GPU_GLOBAL
void copy (amrex::Box bx,
           amrex::Dim3 offset, amrex::Array4<Real> src, amrex::Array4<Real> dst, 
           int scomp, int dcomp, int ncomp)
{
    int ncells = bx.numPts();
    const auto lo = amrex::lbound(bx);
    const auto len = amrex::length(bx);

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

void standardLaunch(MultiFab& src_fab, MultiFab& dst_fab, Real src_val, Real dst_val,
                                      int Ncomp, int num_streams = Gpu::Device::numCudaStreams())
{
    src_fab.setVal(src_val);
    dst_fab.setVal(dst_val);

    if (num_streams > Gpu::Device::numCudaStreams())
    {
        num_streams = Gpu::Device::numCudaStreams();
        std::cout << "Too many streams requested. Using maximum value of " << num_streams << std::endl;
    }

    BL_PROFILE_VAR(std::string("STANDARD " + std::to_string(num_streams)), standard);
    for (MFIter mfi(src_fab); mfi.isValid(); ++mfi)
    {
        Gpu::Device::setStreamIndex(mfi.LocalIndex() % num_streams);
        const Box bx = mfi.validbox();

        const auto src = src_fab.array(mfi);
        const auto dst = dst_fab.array(mfi);

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
    BL_PROFILE_VAR_STOP(standard);

    amrex::Print() << "Standard sum = " << dst_fab.sum() 
                   << "; Expected value = " << src_fab.sum() << std::endl;
}

void standardBoxLaunch(MultiFab& src_fab, MultiFab& dst_fab, Real src_val, Real dst_val,
                                         int Ncomp, int num_streams = Gpu::Device::numCudaStreams())
{
    src_fab.setVal(src_val);
    dst_fab.setVal(dst_val);

    if (num_streams > Gpu::Device::numCudaStreams())
    {
        num_streams = Gpu::Device::numCudaStreams();
        std::cout << "Too many streams requested. Using maximum value of " << num_streams << std::endl;
    }

    BL_PROFILE_VAR(std::string("STANDARD BOX " + std::to_string(num_streams)), standard);
    for (MFIter mfi(src_fab); mfi.isValid(); ++mfi)
    {
        Gpu::Device::setStreamIndex(mfi.LocalIndex() % num_streams);
        const Box bx = mfi.validbox();

        const auto src = src_fab.array(mfi);
        const auto dst = dst_fab.array(mfi);

        const Dim3 offset = {0,0,0};

        const auto ec = Cuda::ExecutionConfig(bx.numPts());
        AMREX_CUDA_LAUNCH_GLOBAL(ec, copy, bx,
                                 offset, src, dst,
                                 0, 0, Ncomp);
    }
    BL_PROFILE_VAR_STOP(standard);

    amrex::Print() << "Standard sum = " << dst_fab.sum() 
                   << "; Expected value = " << src_fab.sum() << std::endl;
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

        standardLaunch(x, y, 0.004, 0.0027, 1, 16);
        standardLaunch(x, y, 0.867, 0.5309, 1, 2);
        standardLaunch(x, y, 0.123, 0.4560, 1, 1);

        standardBoxLaunch(x, y, 0.004, 0.0027, 1, 16);
        standardBoxLaunch(x, y, 0.867, 0.5309, 1, 2);
        standardBoxLaunch(x, y, 0.123, 0.4560, 1, 1);

        {
            x.setVal(0.123);
            y.setVal(0.456);

            BL_PROFILE_VAR("FUSED 1", fused);

            BL_PROFILE_VAR("FUSED: ALLOC", fusedalloc);

            int lsize = x.local_size();
            Box* bx_h =           static_cast<Box*>         (The_Pinned_Arena()->alloc(lsize*sizeof(Box)));
            Dim3* offset_h =      static_cast<Dim3*>        (The_Pinned_Arena()->alloc(lsize*sizeof(Dim3)));
            Array4<Real>* src_h = static_cast<Array4<Real>*>(The_Pinned_Arena()->alloc(lsize*sizeof(Array4<Real>)));
            Array4<Real>* dst_h = static_cast<Array4<Real>*>(The_Pinned_Arena()->alloc(lsize*sizeof(Array4<Real>)));

            Box* bx_d =           static_cast<Box*>         (The_Device_Arena()->alloc(lsize*sizeof(Box)));
            Dim3* offset_d =      static_cast<Dim3*>        (The_Device_Arena()->alloc(lsize*sizeof(Dim3)));
            Array4<Real>* src_d = static_cast<Array4<Real>*>(The_Device_Arena()->alloc(lsize*sizeof(Array4<Real>)));
            Array4<Real>* dst_d = static_cast<Array4<Real>*>(The_Device_Arena()->alloc(lsize*sizeof(Array4<Real>))); 

            BL_PROFILE_VAR_STOP(fusedalloc);

            BL_PROFILE_VAR("FUSED: SETUP", fusedsetup);

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                int idx = mfi.LocalIndex(); 

                bx_h[idx] = mfi.validbox();
                offset_h[idx] = {0,0,0};
                src_h[idx] = x.array(mfi);
                dst_h[idx] = y.array(mfi);
            }
            BL_PROFILE_VAR_STOP(fusedsetup);
            BL_PROFILE_VAR("FUSED: MEMCPY", fusedmem);

            AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(bx_d,     bx_h,     lsize*sizeof(Box),          cudaMemcpyHostToDevice));
            AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(offset_d, offset_h, lsize*sizeof(Dim3),         cudaMemcpyHostToDevice));
            AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(src_d,    src_h,    lsize*sizeof(Array4<Real>), cudaMemcpyHostToDevice));
            AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(dst_d,    dst_h,    lsize*sizeof(Array4<Real>), cudaMemcpyHostToDevice));

            amrex::Gpu::Device::synchronize();
            BL_PROFILE_VAR_STOP(fusedmem);

            BL_PROFILE_VAR("FUSED: LAUNCH", fusedl);

            int cpt = 1;
            const auto ec = Cuda::ExecutionConfig((bx_h[0].numPts()+cpt-1)/cpt);

            AMREX_CUDA_LAUNCH_GLOBAL(ec, copy, lsize, 
                                     bx_d, offset_d,
                                     src_d, dst_d,
                                     0, 0, Ncomp);

            amrex::Gpu::Device::synchronize();
            BL_PROFILE_VAR_STOP(fusedl);

            BL_PROFILE_VAR_STOP(fused);

            amrex::Print() << "Fused Function Array sum = " << y.sum() << "; Expected value = " << x.sum() << std::endl;

            The_Pinned_Arena()->free(bx_h);
            The_Pinned_Arena()->free(offset_h);
            The_Pinned_Arena()->free(src_h);
            The_Pinned_Arena()->free(dst_h);

            The_Device_Arena()->free(bx_d);
            The_Device_Arena()->free(offset_d);
            The_Device_Arena()->free(src_d);
            The_Device_Arena()->free(dst_d);
        }

    }

    amrex::Finalize();
}
