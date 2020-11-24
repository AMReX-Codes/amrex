
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Arena.H>
#include <AMReX_KernelTimer.H>


using namespace amrex;

void main_main();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();
}


void main_main ()
{
    // Initialize box array
    BoxArray ba;
    {
        // This parameter can control GPU compute work
        int n_cell = 256;
        int max_grid_size = 64;
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.query("max_grid_size", max_grid_size);
        Box domain_box(IntVect(0), IntVect(n_cell-1));
        ba.define(domain_box);
        ba.maxSize(max_grid_size);
    }

    MultiFab mf(ba,DistributionMapping{ba},1,0);

    // Now compared different 'timer' instrumentation

    //////////////////////
    // -1) Dummy launch //
    //////////////////////
    mf.setVal(0.0);
    {
        amrex::Gpu::Device::synchronize();
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            amrex::ParallelFor(bx,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k)
                               {
                                   fab(i,j,k) += 1.;
                               });
        }
        
        amrex::Gpu::Device::synchronize();
    }
    
    /////////////////
    // 0) Baseline //
    /////////////////
    mf.setVal(0.0);
    {
        BL_PROFILE("baseline");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            amrex::ParallelFor(bx,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k)
                               {
                                   fab(i,j,k) += 1.;
                               });
        }
        
        amrex::Gpu::Device::synchronize();
    }

    //////////////////
    // 1) GPU clock //
    //////////////////
    // Sample KernelTimer instrumentation of amrex::ParallelFor function;
    // we pass this to KernelTimer to store accumulated thread cycles,
    // as a proxy for GPU compute work; for good performance, we allocate
    // pinned host memory
    mf.setVal(0.0);
    {
        BL_PROFILE("gpu_clock");
        amrex::Real* cost = (amrex::Real*) The_Managed_Arena()->alloc(sizeof(amrex::Real));
        *cost = 0.;

        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            amrex::ParallelFor(bx,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k)
                               {
                                   // GPU timer instrumentation; constructor takes two
                                   // arguments, first controls whether timer is active
                                   // and second is a pointer to the Real that stores 
                                   // accumulated GPU thread cycles
                                   amrex::KernelTimer KnlTimer(true, cost);
                                   fab(i,j,k) += 1.;
                               });
        }
        // Now cost is filled with thread-wise summed cycles
        The_Managed_Arena()->free(cost);

        amrex::Gpu::Device::synchronize();
    }

    ///////////////////
    // 2) Reduce sum //
    ///////////////////
    mf.setVal(0.0);
    {
        BL_PROFILE("reduce_sum");

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            reduce_op.eval(bx, reduce_data,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                           {
                               auto t0 = clock64();
                               fab(i,j,k) += 1.;   
                               auto t1 = clock64();
                               return {amrex::Real(t1-t0)};
                           });
        }
        ReduceTuple hv = reduce_data.value();
        auto result = amrex::get<0>(hv);

        amrex::Gpu::Device::synchronize();
    }

    ///////////////////////////////
    // 3) PTX timer + reduce sum //
    ///////////////////////////////
    mf.setVal(0.0);
    {
        BL_PROFILE("reduce_sum_with_PTX");

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            reduce_op.eval(bx, reduce_data,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                           {
                               unsigned cost = 0;
                               asm("mov.u32 %0, %%clock;" : "=r"(cost));
                               fab(i,j,k) += 1.;
                               asm(".reg .u32 t1;\n\t"         // temp reg t1
                                   " mov.u32 t1, %%clock;\n\t" // t1 = clock
                                   " sub.u32 %0, t1, %1;"  // cost = t1 - cost
                                   : "=r"(cost) : "r" (cost));
                               return {amrex::Real(cost)};
                           });
        }
        ReduceTuple hv = reduce_data.value();
        auto result = amrex::get<0>(hv);
        
        amrex::Gpu::Device::synchronize();
    }

    ////////////////////////
    // 4) Synchronization //
    ////////////////////////
    mf.setVal(0.0);        
    {
        BL_PROFILE("synchronization");
        auto t0 = clock();
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            amrex::ParallelFor(bx,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k)
                               {
                                   fab(i,j,k) += 1.;
                               });
            amrex::Gpu::Device::synchronize();
        }
        auto t1 = clock();

        amrex::Gpu::Device::synchronize();
    }

    ////////////////////
    // 5) CUDA events //
    ////////////////////
    mf.setVal(0.0);        
    {
        BL_PROFILE("CUDA_events");
        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);

        cudaEventRecord(start);
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            amrex::ParallelFor(bx,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k)
                               {
                                   fab(i,j,k) += 1.;
                               });
        }
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);

        amrex::Gpu::Device::synchronize();
    }
    
}
