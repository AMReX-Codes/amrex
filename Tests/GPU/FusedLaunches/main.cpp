#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

// Single launch version
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

// Single launch, simutaneous box calc version
AMREX_GPU_GLOBAL
void copy (int size, amrex::Box* bx, amrex::Dim3* offset,
           amrex::Array4<Real>* src, amrex::Array4<Real>* dst, 
           int scomp, int dcomp, int ncomp, int simul)
{
    // Works because boxes are currently of equal size.
    //    When box size varies: calc or pass in maximum box size.
    int ncells = bx[0].numPts();
    int tid = blockDim.x*blockIdx.x+threadIdx.x;
    int bidx = tid/ncells;

    if (bidx < simul)
    {
        for (int l = 0; l < size; l+=simul)
        {
            int bid = l+bidx;
            const auto lo  = amrex::lbound(bx[bid]);
            const auto len = amrex::length(bx[bid]);
 
            for (int icell = (tid%ncells), stride = (blockDim.x*gridDim.x)/simul; icell < ncells; icell += stride) {
                int k =  icell /   (len.x*len.y);
                int j = (icell - k*(len.x*len.y)) /   len.x;
                int i = (icell - k*(len.x*len.y)) - j*len.x;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                for (int n = 0; n < ncomp; ++n) {
                    (dst[bid])(i,j,k,dcomp+n) = (src[bid])(i+offset[bid].x,j+offset[bid].y,k+offset[bid].z,scomp+n); 
                }
            }
        }
    }
}

// Launch within MFIter loop version
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
                                      int Ncomp, int cpt = 1, int num_streams = Gpu::Device::numCudaStreams())
{
    src_fab.setVal(src_val);
    dst_fab.setVal(dst_val);

    if (num_streams > Gpu::Device::numCudaStreams())
    {
        num_streams = Gpu::Device::numCudaStreams();
        std::cout << "Too many streams requested. Using maximum value of " << num_streams << std::endl;
    }

    std::string timer_name = "STANDARD: " +
                               std::to_string(num_streams) + " streams, " + 
                               std::to_string(cpt) + " CPT";

    BL_PROFILE_VAR(timer_name, standard);
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

        const auto ec = Cuda::ExecutionConfig((bx.numPts()+cpt-1)/cpt);
        AMREX_CUDA_LAUNCH_GLOBAL(ec, copy,
                                 lo, len, ncells,
                                 offset, src, dst,
                                 0, 0, Ncomp);
    }
    BL_PROFILE_VAR_STOP(standard);

    Real src_sum = src_fab.sum();
    Real dst_sum = dst_fab.sum();

    if (src_sum == dst_sum)
    {
        amrex::Print() << timer_name << ": no error." << std::endl;
    }
    else
    {
        amrex::Print() << timer_name << " error found." << std::endl;
        amrex::Print() << " ---- dst = " << dst_sum 
                       << "; src = " << src_sum << std::endl;
    }
}

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

void fusedLaunch(MultiFab& src_fab, MultiFab& dst_fab, Real src_val, Real dst_val,
                              int Ncomp, int cpt = 1, int simul = 1, int num_launches = 1)
// NUM LAUNCHES.
{
    src_fab.setVal(src_val);
    dst_fab.setVal(dst_val);

    int lsize = src_fab.local_size();

    // Note: Not a parallel copy, so src_ba = dst_ba.
/*
    if (num_launches > lsize)
    {
        num_launches = lsize; 
        std::cout << "Too many launches requested. Using one box per launch: " << dst_fab.local_size() << std::endl;
    }
*/
    if (num_launches != 1)
    {
        num_launches = 1;
        std::cout << "Multiple fused launches feature not yet implemented. num_launces set to 1" << std::endl; 
    }

    std::string timer_name = "FUSED: " +
                               std::to_string(num_launches) + " launches, " + 
                               std::to_string(cpt) + " CPT, " +
                               std::to_string(simul) + " simul";

    BL_PROFILE_VAR(timer_name, fused);
    BL_PROFILE_VAR("FUSED: ALLOC", fusedalloc);

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

    for (MFIter mfi(src_fab); mfi.isValid(); ++mfi)
    {
	int idx = mfi.LocalIndex(); 

	bx_h[idx] = mfi.validbox();
	offset_h[idx] = {0,0,0};
	src_h[idx] = src_fab.array(mfi);
	dst_h[idx] = dst_fab.array(mfi);
    }
 
    BL_PROFILE_VAR_STOP(fusedsetup);
    BL_PROFILE_VAR("FUSED: COPY+LAUNCH", fusedl);
                     
    AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(bx_d,     bx_h,     lsize*sizeof(Box),          cudaMemcpyHostToDevice));
    AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(offset_d, offset_h, lsize*sizeof(Dim3),         cudaMemcpyHostToDevice));
    AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(src_d,    src_h,    lsize*sizeof(Array4<Real>), cudaMemcpyHostToDevice));
    AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(dst_d,    dst_h,    lsize*sizeof(Array4<Real>), cudaMemcpyHostToDevice));

    const auto ec = Cuda::ExecutionConfig(((bx_h[0].numPts()*simul)+cpt-1)/cpt);

    AMREX_CUDA_LAUNCH_GLOBAL(ec, copy, lsize, 
			     bx_d, offset_d,
			     src_d, dst_d,
			     0, 0, Ncomp, simul);

    amrex::Gpu::Device::synchronize();

    BL_PROFILE_VAR_STOP(fusedl);
    BL_PROFILE_VAR_STOP(fused);

    Real src_sum = src_fab.sum();
    Real dst_sum = dst_fab.sum();
    if (src_sum == dst_sum)
    {
        amrex::Print() << timer_name << ": no error." << std::endl;
    }
    else
    {
        amrex::Print() << timer_name << " error found." << std::endl;
        amrex::Print() << " ---- dst = " << dst_sum 
                       << "; src = " << src_sum << std::endl;
    }

    The_Pinned_Arena()->free(bx_h);
    The_Pinned_Arena()->free(offset_h);
    The_Pinned_Arena()->free(src_h);
    The_Pinned_Arena()->free(dst_h);

    The_Device_Arena()->free(bx_d);
    The_Device_Arena()->free(offset_d);
    The_Device_Arena()->free(src_d);
    The_Device_Arena()->free(dst_d);
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
        // NOTE: THIS IS NOT A PARALLEL COPY, SO BA & DM MUST BE IDENTICAL HERE.
        //   (IT'S EITHER THIS, OR SETVAL / LOCAL SUM TYPE STUFF).
        MultiFab x(ba, dm, Ncomp, Nghost);
        MultiFab y(ba, dm, Ncomp, Nghost);
        x.setVal(0.0);
        y.setVal(0.0);

        amrex::Print() << "Testing on " << n_cell << "^3 boxes with max grid size " << max_grid_size << std::endl 
                       << "Number of boxes per MultiFab: " << x.size() << std::endl << std::endl;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        // standardLaunch (src, dst, srcVal, dstVal, Ncomp, cellsPerThread, numStreams);
        // fusedLaunch (src, dst, srcVal, dstVal, Ncomp, cellsPerThread, simultaneousBoxes, numLaunches );

        standardLaunch(x, y, 0.004, 0.0027, 1, 16);
        standardLaunch(x, y, 0.867, 0.5309, 1, 2);
        standardLaunch(x, y, 0.123, 0.4560, 1, 1);
        fusedLaunch(x, y, 0.123, 0.456, 1, 1, 1, 1);      
        fusedLaunch(x, y, 1.0, 0.0, 1, 8, 1, 1);

    }

    amrex::Finalize();
}
