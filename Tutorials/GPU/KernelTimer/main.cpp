
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Arena.H>
#include <AMReX_KernelTimer.H>


using namespace amrex;

void main_main();

void CUDART_CB getTimer (cudaStream_t stream, cudaError_t status, void* data);

ReduceOps<ReduceOpSum> a_reduce_op;
using ReduceTuple = typename ReduceData<Real>::Type;


struct myData {
    
    ReduceTuple hv;
    Real cost;
    
    myData (Real& a_cost, ReduceTuple& a_hv) : cost{a_cost}, hv{a_hv} {}
    
};


int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();
}

void CUDART_CB getTimer (cudaStream_t stream, cudaError_t status, void *data)
{
    //((myData*)data)->cost = amrex::get<0>( ((myData*)data)->hv);
}

void main_main ()
{
    // Initialize box array
    BoxArray ba;
    {
        // This parameter can control GPU compute work
        int n_cell = 256;
        int max_grid_size = 16;
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.query("max_grid_size", max_grid_size);
        Box domain_box(IntVect(0), IntVect(n_cell-1));
        ba.define(domain_box);
        ba.maxSize(max_grid_size);
    }

    MultiFab mf(ba,DistributionMapping{ba},1,0);

    // Costs for each box of boxarray
    std::unique_ptr<amrex::Real> costs( (amrex::Real*) The_Managed_Arena()->alloc(ba.size()*sizeof(amrex::Real)) );

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

    ////////////////////
    // 1.0) GPU clock //
    ////////////////////
    // Sample KernelTimer instrumentation of amrex::ParallelFor function;
    // we pass this to KernelTimer to store accumulated thread cycles,
    // as a proxy for GPU compute work; for good performance, we allocate
    // pinned host memory
    mf.setVal(0.0);
    {
        BL_PROFILE("gpu_clock");

        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            amrex::Real* cost = &costs.get()[mfi.index()];
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
        // Now costs is filled with thread-wise summed cycles
        
        amrex::Gpu::Device::synchronize();
    }

    //////////////////////////////////
    // 1.1) GPU clock shared memory //
    //////////////////////////////////
    // Sample KernelTimer instrumentation of amrex::ParallelFor function;
    // we pass this to KernelTimer to store accumulated thread cycles,
    // as a proxy for GPU compute work; for good performance, we allocate
    // pinned host memory
    mf.setVal(0.0);
    {
        BL_PROFILE("gpu_clock_shared");

        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            amrex::Real* cost = &costs.get()[mfi.index()];
            amrex::ParallelFor(bx,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k)
                               {
                                   // GPU timer instrumentation; constructor takes two
                                   // arguments, first controls whether timer is active
                                   // and second is a pointer to the Real that stores 
                                   // accumulated GPU thread cycles
                                   amrex::KernelTimerShared KnlTimer(true, cost);
                                   fab(i,j,k) += 1.;
                               });
        }
        // Now costs is filled with thread-wise summed cycles
                
        amrex::Gpu::Device::synchronize();
    }

    //////////////////////////////////////
    // 1.2) GPU clock shared memory PTX //
    //////////////////////////////////////
    // Sample KernelTimer instrumentation of amrex::ParallelFor function;
    // we pass this to KernelTimer to store accumulated thread cycles,
    // as a proxy for GPU compute work; for good performance, we allocate
    // pinned host memory
    mf.setVal(0.0);
    {
        BL_PROFILE("gpu_clock_shared_PTX");

        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            amrex::Real* cost = &costs.get()[mfi.index()];
            amrex::ParallelFor(bx,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k)
                               {
                                   // GPU timer instrumentation; constructor takes two
                                   // arguments, first controls whether timer is active
                                   // and second is a pointer to the Real that stores 
                                   // accumulated GPU thread cycles
                                   amrex::KernelTimerSharedPTX KnlTimer(true, cost);
                                   fab(i,j,k) += 1.;
                               });
        }
        // Now costs is filled with thread-wise summed cycles
        
        amrex::Gpu::Device::synchronize();
    }

    /////////////////////
    // 2.0) Reduce sum //
    /////////////////////
    mf.setVal(0.0);
    {
        BL_PROFILE("reduce_sum_callback");

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


            ReduceTuple hv = reduce_data.value();
            myData data(costs.get()[mfi.index()], hv);
            cudaStreamAddCallback(amrex::Gpu::gpuStream(), getTimer, (void*)&data, 0);
        }

        // amrex::Gpu::Device::synchronize();

        // for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        // {
        //     auto x = costs.get()[mfi.index()];
        //     amrex::Print() << mfi.index() << ": cost is " << x << "\n";
        // }
        
        amrex::Gpu::Device::synchronize();
    }

    /////////////////////////////////
    // 2.1) Reduce sum no callback //
    /////////////////////////////////
    mf.setVal(0.0);
    {
        BL_PROFILE("reduce_sum_no_callback");

        ReduceOps<ReduceOpSum> reduce_op;
        using ReduceTuple = typename ReduceData<Real>::Type;
        Vector<ReduceTuple*> reduce_data;

        //auto x = ReduceTuple(reduce_op);
        
        // for (int i=0; i<mf.size(); ++i)
        // {
        //     reduce_data.push_back(new ReduceTuple(reduce_op));
        // }
 
        // for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        // {
        //     const Box& bx = mfi.tilebox();
        //     Array4<Real> const& fab = mf.array(mfi);
        //     reduce_op.eval(bx, *reduce_data[mfi.index()],
        //                    [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        //                    {
        //                        auto t0 = clock64();
        //                        fab(i,j,k) += 1.;   
        //                        auto t1 = clock64();
        //                        return {amrex::Real(t1-t0)};
        //                    });


        //     //ReduceTuple hv = reduce_data.value();
        //     //myData data(costs.get()[mfi.index()], hv);
        //     //cudaStreamAddCallback(amrex::Gpu::gpuStream(), getTimer, (void*)&data, 0);
        // }

        for (int i=0; i<mf.size(); ++i)
        {
            //ReduceTuple hv = reduce_data.value();
            //costs.get()[i] = amrex::get<0>(hv);
        }
        
        amrex::Gpu::Device::synchronize();
    }

    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        auto x = costs.get()[mfi.index()];
        amrex::Print() << mfi.index() << ": cost is " << x << "\n";
    }
    
    /////////////////////////////////
    // 2.2) PTX timer + reduce sum //
    /////////////////////////////////
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
                                   " sub.u32 %0, t1, %1;"      // cost = t1 - cost
                                   : "=r"(cost) : "r" (cost));
                               return {amrex::Real(cost)};
                           });
        }
        ReduceTuple hv = reduce_data.value();
        auto result = amrex::get<0>(hv);
        
        amrex::Gpu::Device::synchronize();
    }

    ////////////////////////
    // 3) Synchronization //
    ////////////////////////
    mf.setVal(0.0);
    {
        BL_PROFILE("synchronization");

        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            auto t0 = clock();
            
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            amrex::ParallelFor(bx,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k)
                               {
                                   fab(i,j,k) += 1.;
                               });
            amrex::Gpu::Device::synchronize();
            
            auto t1 = clock();

            costs.get()[mfi.index()] = amrex::Real(t1 - t0);
        }

        amrex::Gpu::Device::synchronize();
    }

    ////////////////////
    // 4) CUDA events //
    ////////////////////
    mf.setVal(0.0);
    {
        BL_PROFILE("CUDA_events");

        cudaEvent_t starts[ba.size()], stops[ba.size()];
        for (int i=0; i<ba.size(); i++)
        {
            cudaEventCreate(&starts[i]);
            cudaEventCreate(&stops[i]);
        }
        
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            
            cudaEventRecord(starts[mfi.index()]);
        
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            amrex::ParallelFor(bx,
                               [=] AMREX_GPU_DEVICE (int i, int j, int k)
                               {
                                   fab(i,j,k) += 1.;
                               });

            cudaEventRecord(stops[mfi.index()]);
        }

        for (int i=0; i<ba.size(); i++)
        {
            cudaEventSynchronize(stops[i]);
            float milliseconds;
            cudaEventElapsedTime(&milliseconds, starts[i], stops[i]);
            costs.get()[i] = amrex::Real(milliseconds);
        }

        amrex::Gpu::Device::synchronize();
    }

    // Managed memory cleanup
    costs.release();
}
