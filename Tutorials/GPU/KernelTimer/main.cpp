#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Arena.H>
#include <AMReX_KernelTimer.H>


using namespace amrex;

using ReduceTuple = typename ReduceData<Real>::Type;

void main_main();

// struct myData {
    
//     ReduceTuple hv;
//     Real cost;
    
//     myData (Real& a_cost, ReduceTuple& a_hv) : cost{a_cost}, hv{a_hv} {}
    
// };


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

    // Now compare different 'timer' instrumentation

    //////////////////////
    // -1) Dummy launch //
    //////////////////////
    mf.setVal(0.0);
    {
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
    }

    //////////////////////////////////
    // 1.1) GPU clock shared memory //
    //////////////////////////////////
    // Same as GPU clock test but uses a shared memory optimization;
    // this is experimental and not guaranteed to be thread safe
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
                                   amrex::KernelTimerShared KnlTimer(true, cost);
                                   fab(i,j,k) += 1.;
                               });
        }
    }

    //////////////////////////////////////
    // 1.2) GPU clock shared memory PTX //
    //////////////////////////////////////
    // Same as GPU clock shared memory test but with PTX clock
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
                                   amrex::KernelTimerSharedPTX KnlTimer(true, cost);
                                   fab(i,j,k) += 1.;
                               });
        }        
    }

    ///////////////////
    // 2) Reduce sum //
    ///////////////////
    // Cost measured with reduction operation
    mf.setVal(0.0);
    {
        BL_PROFILE("reduce_sum");

        ReduceOps<ReduceOpSum> reduce_op;
        Vector<std::unique_ptr<ReduceData<Real>>> reduce_data(mf.size());

        for (int i=0; i<mf.size(); ++i)
        {
            reduce_data[i] = std::unique_ptr<ReduceData<Real>>(new ReduceData<Real>(reduce_op));
        }
 
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            reduce_op.eval(bx, *reduce_data[mfi.index()],
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                           {
                               auto t0 = clock64();
                               fab(i,j,k) += 1.;   
                               auto t1 = clock64();
                               return {amrex::Real(t1-t0)};
                           });
        }

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            ReduceTuple hv = reduce_data[mfi.index()]->value();
            costs.get()[mfi.index()] = amrex::get<0>(hv);
        }        
    }
    
    ////////////////////////
    // 3) Synchronization //
    ////////////////////////
    // MFIter loop is synchronized
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
    }

    ////////////////////
    // 4) CUDA Events //
    ////////////////////
    // Use CUDA Events API to get box time
    mf.setVal(0.0);
    {
        BL_PROFILE("CUDA_Events");

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
    }

    // Managed memory cleanup
    costs.release();
}
