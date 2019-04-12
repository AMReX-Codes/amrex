#include <cuda_runtime.h>

#include <iostream>
#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>

using namespace amrex;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

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

AMREX_GPU_GLOBAL
void copy (amrex::Dim3 lo, amrex::Dim3 len, int ncells,
           amrex::Dim3 offset, amrex::FArrayBox* src_fab, amrex::FArrayBox* dst_fab,
           int scomp, int dcomp, int ncomp)
{
    Array4<Real> src = src_fab->array();
    Array4<Real> dst = dst_fab->array();

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
           amrex::Dim3 offset, amrex::FArrayBox** src_fab, amrex::FArrayBox** dst_fab,
           int scomp, int dcomp, int ncomp)
{
    Array4<Real> src = (*src_fab)->array();
    Array4<Real> dst = (*dst_fab)->array();

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
        amrex::Print() << "amrex::Initialize complete." << "\n";

        // ===================================
        // Simple cuda action to make sure all tests have cuda.
        // Allows nvprof to return data.
        int devices = 0;
#ifdef AMREX_USE_CUDA
        cudaGetDeviceCount(&devices);
#endif
        amrex::Print() << "Hello world from AMReX version " << amrex::Version() << ". GPU devices: " << devices << "\n";
        amrex::Print() << "**********************************\n\n"; 
        // ===================================

        // What time is it now?  We'll use this to compute total run time.
        Real strt_time = amrex::second();

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
        Geometry geom;
        {
            IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
            IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
            Box domain(dom_lo, dom_hi);

            // Initialize the boxarray "ba" from the single box "bx"
            ba.define(domain);
            // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
            ba.maxSize(max_grid_size);

            // This defines the physical box, [-1,1] in each direction.
            RealBox real_box({AMREX_D_DECL(-1.0,-1.0,-1.0)},
                             {AMREX_D_DECL( 1.0, 1.0, 1.0)});

            // This defines a Geometry object
            geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
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
        MultiFab v(ba, dm, Ncomp, Nghost);
        MultiFab w(ba, dm, Ncomp, Nghost);
        x.setVal(0.0);
        y.setVal(0.0);
        v.setVal(0.0);
        w.setVal(0.0);

        int size = x.local_size();
/*
        // Array pointers on host
        Array4<Real>* src_arrs_host;
        Array4<Real>* dst_arrs_host;
        src_arrs_host = static_cast<Array4<Real>*>(malloc(sizeof(Array4<Real>)*size));
        dst_arrs_host = static_cast<Array4<Real>*>(malloc(sizeof(Array4<Real>)*size));
*/
        Arena* the_arena = The_Pinned_Arena();


        FArrayBox** src_fab;
        FArrayBox** dst_fab;

        cudaMallocHost(&src_fab, sizeof(FArrayBox*)*size);
        cudaMallocHost(&dst_fab, sizeof(FArrayBox*)*size);

/*
        FArrayBox* src_fab[size];
        FArrayBox* dst_fab[size];
        cudaHostAlloc(&src_fab[0], sizeof(FArrayBox*)*size, 0);
        cudaHostAlloc(&dst_fab[0], sizeof(FArrayBox*)*size, 0);
*/
/*
        src_fab = static_cast<FArrayBox*>(the_arena->alloc(sizeof(FArrayBox*)*size));
        dst_fab = static_cast<FArrayBox*>(the_arena->alloc(sizeof(FArrayBox*)*size));
*/
        Real points = ba.numPts();

        amrex::Print() << "Testing on " << n_cell << "^3 boxes with max grid size " << max_grid_size
                       << std::endl << std::endl;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Launch without graphs

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

                AMREX_HOST_DEVICE_FOR_4D ( bx, Ncomp, i, j, k, n,
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

// ---------------------------------------

        {
            x.setVal(0.75);
            y.setVal(0.25);

            BL_PROFILE("Function Copy - FAB");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.validbox();

                int idx = mfi.LocalIndex();
                const auto src = x.fabPtr(mfi);
                const auto dst = y.fabPtr(mfi);

                int ncells = bx.numPts();
                const auto lo  = amrex::lbound(bx);
                const auto len = amrex::length(bx);
                const auto ec = Cuda::ExecutionConfig(ncells);
                const Dim3 offset = {0,0,0};

                AMREX_CUDA_LAUNCH_GLOBAL(ec, copy,
                                         lo, len, ncells,
                                         offset, src, dst,
                                         0, 0, 1);
            }

            amrex::Print() << "No Graph Function FAB sum = " << y.sum() << "; Expected value = " << x.sum() << std::endl;
        }

// ---------------------------------------

        {
            x.setVal(0.867);
            y.setVal(0.5309);

            BL_PROFILE("Lambda Copy - FAB");

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.validbox();

                int idx = mfi.LocalIndex();
                const auto src_fab = x.fabPtr(mfi);
                const auto dst_fab = y.fabPtr(mfi);
                const Dim3 offset = {0,0,0};
                int dcomp = 0;
                int scomp = 0;

                AMREX_HOST_DEVICE_FOR_4D ( bx, Ncomp, i, j, k, n,
                {
                    Array4<Real> src = src_fab->array();
                    Array4<Real> dst = dst_fab->array();
                    dst(i,j,k,dcomp+n) = src(i+offset.x,j+offset.y,k+offset.z,scomp+n); 
                });
            }

            amrex::Print() << "No Graph Lambda FAB sum = " << y.sum() << "; Expected value = " << x.sum() << std::endl;
        }

        amrex::Print() << "=============" << std::endl;

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create a single graph for the MFIter loop using a function:
//        an empty node at the start linked to each individually captured stream graph.

        cudaGraphExec_t graphExec;

        {
            x.setVal(2e-5);
            y.setVal(0.0);

            BL_PROFILE("cudaGraphFunction");

// --------- Capture each stream in the MFIter loop ----------

            BL_PROFILE_VAR("cudaGraphFunction-create", cgfc);

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                if (mfi.LocalIndex() == 0)
                {
                    amrex::Gpu::Device::startGraphRecording();
                } 

                const Box bx = mfi.validbox();

                int idx = mfi.LocalIndex();
                src_fab[idx] = x.fabPtr(mfi);
                dst_fab[idx] = y.fabPtr(mfi);

                int ncells = bx.numPts();
                const auto lo  = amrex::lbound(bx);
                const auto len = amrex::length(bx);
                const auto ec = Cuda::ExecutionConfig(ncells);
                const Dim3 offset = {0,0,0};

                AMREX_CUDA_LAUNCH_GLOBAL(ec, copy,
                                         lo, len, ncells,
                                         offset, &(src_fab[idx]), &(dst_fab[idx]), 0, 0, 1); 

                if (mfi.LocalIndex() == (x.local_size() - 1) )
                {
                    graphExec = amrex::Gpu::Device::stopGraphRecording(); 
                }
            }

            BL_PROFILE_VAR_STOP(cgfc);

// --------- Launch the graph  ----------

            BL_PROFILE_VAR("cudaGraphFunction-launch", cgfl);

            amrex::Gpu::Device::executeGraph(graphExec);

            BL_PROFILE_VAR_STOP(cgfl);

            amrex::Print() << "Graphed = " << y.sum() << "; Expected value = " << x.sum() << std::endl;

// --------- Relaunch the graph with a different result  ----------

            x.setVal(4337654e-9);
            y.setVal(0.0);

            BL_PROFILE_VAR("cudaGraphFunction-relaunch", cgfrl);

            amrex::Gpu::Device::executeGraph(graphExec);

            BL_PROFILE_VAR_STOP(cgfrl);

            amrex::Print() << "Regraphed = " << y.sum() << "; Expected value = " << x.sum() << std::endl;

// --------- Relaunch the graph on different MFIters  ----------
// --------- Doesn't work changing the Array4 in CPU memory, even with function. ----------
// --------- Trying with Array4 in device memory defined by Arena. ----------

            x.setVal(0.238761);
            v.setVal(0.5e-5);
            w.setVal(0.0);

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                int idx = mfi.LocalIndex();
                src_fab[idx] = v.fabPtr(mfi);
                dst_fab[idx] = w.fabPtr(mfi);
            }

            BL_PROFILE_VAR("cudaGraphFunction-diff", cgfdiff);

            amrex::Gpu::Device::executeGraph(graphExec);

            BL_PROFILE_VAR_STOP(cgfdiff);

            amrex::Print() << "Diff Graph Function = " << v.sum() << "; Expected value = " << w.sum() << std::endl;
            amrex::Print() << " x = " << x.sum() << "; y = " << y.sum() << std::endl;
        }

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//      Create a single graph for the MFIter loop using a lambda:
//        an empty node at the start linked to each individually captured stream graph.
        amrex::Print() << "=============" << std::endl;

        {
            x.setVal(4e-5);
            y.setVal(0.0);

            BL_PROFILE("cudaGraphLambda");

// --------- Capture each stream in the MFIter loop ----------

            BL_PROFILE_VAR("cudaGraphLambda-create", cgfc);

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                if (mfi.LocalIndex() == 0)
                {
                    amrex::Gpu::Device::startGraphRecording();
                } 

                const Box bx = mfi.validbox();

                int idx = mfi.LocalIndex();
                src_fab[idx] = x.fabPtr(mfi);
                dst_fab[idx] = y.fabPtr(mfi);
                const Dim3 offset = {0,0,0};
                int dcomp = 0;
                int scomp = 0;

                AMREX_HOST_DEVICE_FOR_4D ( bx, Ncomp, i, j, k, n,
                {
                    Array4<Real> src = src_fab[idx]->array();
                    Array4<Real> dst = dst_fab[idx]->array();
                    dst(i,j,k,dcomp+n) = src(i+offset.x,j+offset.y,k+offset.z,scomp+n); 
                });

                if (mfi.LocalIndex() == (x.local_size() - 1) )
                {
                    graphExec = amrex::Gpu::Device::stopGraphRecording(); 
                }
            }

            BL_PROFILE_VAR_STOP(cgfc);

// --------- Launch the graph  ----------

            BL_PROFILE_VAR("cudaGraphLambda-launch", cgfl);

            amrex::Gpu::Device::executeGraph(graphExec);

            BL_PROFILE_VAR_STOP(cgfl);

            amrex::Print() << "Graphed = " << y.sum() << "; Expected value = " << x.sum() << std::endl;

// --------- Relaunch the graph with a different result  ----------

            x.setVal(6.5784e-9);
            y.setVal(0.0);

            BL_PROFILE_VAR("cudaGraphLambda-relaunch", cgfrl);

            amrex::Gpu::Device::executeGraph(graphExec);

            BL_PROFILE_VAR_STOP(cgfrl);

            amrex::Print() << "Regraphed = " << y.sum() << "; Expected value = " << x.sum() << std::endl;

// --------- Relaunch the graph on different MFIters  ----------
// --------- Doesn't work changing the Array4 in CPU memory, even with function. ----------
// --------- Trying with Array4 in device memory defined by Arena. ----------

            x.setVal(0.167852);
            v.setVal(0.15e-5);
            w.setVal(0.0);

            for (MFIter mfi(x); mfi.isValid(); ++mfi)
            {
                int idx = mfi.LocalIndex();
                src_fab[idx] = v.fabPtr(mfi);
                dst_fab[idx] = w.fabPtr(mfi);
            }

            BL_PROFILE_VAR("cudaGraphLambda-diff", cgfdiff);

            amrex::Gpu::Device::executeGraph(graphExec);

            BL_PROFILE_VAR_STOP(cgfdiff);

            amrex::Print() << "Diff Graph Function = " << v.sum() << "; Expected value = " << w.sum() << std::endl;
            amrex::Print() << " x = " << x.sum() << "; y = " << y.sum() << std::endl;
        }

        cudaFreeHost(src_fab);
        cudaFreeHost(dst_fab);

        amrex::Print() << "Test Completed." << std::endl;
    }

    amrex::Finalize();
}
