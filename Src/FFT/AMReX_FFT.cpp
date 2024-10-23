#include <AMReX_FFT_Helper.H>

namespace amrex::FFT::detail
{

DistributionMapping make_iota_distromap (Long n)
{
    AMREX_ASSERT(n <= ParallelContext::NProcsSub());
    Vector<int> pm(n);
    for (int i = 0; i < n; ++i) {
        pm[i] = ParallelContext::local_to_global_rank(i);
    }
    return DistributionMapping(std::move(pm));
}

#ifdef AMREX_USE_HIP
void hip_execute (rocfft_plan plan, void **in, void **out)
{
    rocfft_execution_info execinfo = nullptr;
    AMREX_ROCFFT_SAFE_CALL(rocfft_execution_info_create(&execinfo));

    std::size_t buffersize = 0;
    AMREX_ROCFFT_SAFE_CALL(rocfft_plan_get_work_buffer_size(plan, &buffersize));

    auto* buffer = (void*)amrex::The_Arena()->alloc(buffersize);
    AMREX_ROCFFT_SAFE_CALL(rocfft_execution_info_set_work_buffer(execinfo, buffer, buffersize));

    AMREX_ROCFFT_SAFE_CALL(rocfft_execution_info_set_stream(execinfo, amrex::Gpu::gpuStream()));

    AMREX_ROCFFT_SAFE_CALL(rocfft_execute(plan, in, out, execinfo));

    amrex::Gpu::streamSynchronize();
    amrex::The_Arena()->free(buffer);

    AMREX_ROCFFT_SAFE_CALL(rocfft_execution_info_destroy(execinfo));
}
#endif

}
