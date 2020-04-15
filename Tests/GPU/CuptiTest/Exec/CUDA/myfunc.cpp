#include "myfunc.H"
#include "mykernel.H"

#ifdef AMREX_USE_CUPTI
#ifndef AMREX_TINY_PROFILING
#include <AMReX_CuptiTrace.H> // No need to import if using only the CUPTI profiling macros
#endif
#endif

void doDeviceSleep (MultiFab& mf, int& n) {
    // Call device sleep function at each iteration of MFIter loop

    // Code is instrumented here to capture CUPTI
    // kernel activity from outside MFIter loop
#ifdef AMREX_USE_CUPTI
#ifdef AMREX_TINY_PROFILING
    BL_PROFILE_VAR_NS("CPU::deviceSleep", blpCpuSleep);
    BL_PROFILE_VAR_NS_CUPTI("GPU::deviceSleep", blpGpuSleep);
    BL_PROFILE_VAR_START(blpCpuSleep);      
    BL_PROFILE_VAR_START_CUPTI(blpGpuSleep);
#else
    CuptiTrace cuptiTrace = CuptiTrace();
    cuptiTrace.start();
#endif
#else
#ifdef AMREX_TINY_PROFILING
    BL_PROFILE_VAR_NS("CPU::deviceSleep", blpCpuSleep);
    BL_PROFILE_VAR_START(blpCpuSleep);
#endif
#endif

    for ( MFIter mfi(mf); mfi.isValid(); ++mfi )
    {
        // Test 1: launch kernel
        // deviceSleep<<<1, 1, 0, amrex::Gpu::Device::gpuStream()>>>();

        // Test 2: launch kernel as inlined lambda      
        const Box& bx = mfi.tilebox();
        const Dim3 lo = amrex::lbound(bx);
        const Dim3 hi = amrex::ubound(bx);
        amrex::ParallelFor(bx, 1,
                           [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                           {
                               if (i==lo.x & j==lo.y & k==lo.z) {
                                   // Test 2.1: single sleep
                                   // deviceSleep(1e8);
                                   
                                   // Test 2.2: nested sleep
                                   // Expected to be similar to above
                                   deviceNestedSleep(1e8/3, 3);
                               }
                           });
    }
#ifdef AMREX_USE_CUPTI
#ifdef AMREX_TINY_PROFILING
    BL_PROFILE_VAR_STOP_CUPTI_ID(blpGpuSleep, (int) 1729*amrex::Random());
    // BL_PROFILE_VAR_STOP(blpGpuSleep); // Alternatively, may call without unsigned flag
    BL_PROFILE_VAR_STOP(blpCpuSleep);
#else
    cuptiTrace.stop(1729*amrex::Random());
    // cuptiTrace.stop();  // Alternatively, may call without unsigned flag
#endif

    if (true)
    {
        amrex::Print() << "Average time per box (s): "
                       << computeElapsedTimeUserdata(activityRecordUserdata) << "\n";

        unsigned long long t_start = 0;
        unsigned long long t_stop = 0;
        for (auto& record : activityRecordUserdata) {
            t_start = (unsigned long long) record->getStartTime();
            t_stop = (unsigned long long) record->getEndTime();

            unsigned long long dt = 0;
            dt = (((unsigned long long)t_stop) - ((unsigned long long)t_start));

            // Kernel data captured in a vector `activityRecordUserdata`
            // Print here some information from the captured records
            amrex::AllPrint() << "  t_elapsed (ns):  " << dt
                              << "; n_step: " << n
                              << "; Proc.: " << ParallelContext::MyProcSub()
                              << "; Act. Size: " << activityRecordUserdata.size()
                // In this example, random unsigned flag is assigned
                // for all records captured outside same MFIter loop,
                // but if instrumenting kernel launch within an
                // MFIter loop, records may be flagged with, e.g.,
                // an ID corresponding to the loop iteration
                              << "; Unsigned flag: " << record->getUintID()
                              << "\n";
        }
    }
#endif
}
