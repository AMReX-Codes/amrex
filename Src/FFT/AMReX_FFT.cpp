#include <AMReX_FFT.H>

namespace amrex::FFT
{

#ifdef AMREX_USE_HIP
namespace detail
{
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
}
#endif

}
