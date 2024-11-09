#include <AMReX_FFT.H>
#include <AMReX_FFT_Helper.H>

#include <map>

namespace amrex::FFT
{

namespace
{
    bool s_initialized = false;
    std::map<Key, PlanD> s_plans_d;
    std::map<Key, PlanF> s_plans_f;
}

void Initialize ()
{
    if (!s_initialized)
    {
        s_initialized = true;

#if defined(AMREX_USE_HIP) && defined(AMREX_USE_FFT)
        AMREX_ROCFFT_SAFE_CALL(rocfft_setup());
#endif
    }

    amrex::ExecOnFinalize(amrex::FFT::Finalize);
}

void Finalize ()
{
    if (s_initialized)
    {
        s_initialized = false;

        Clear();

#if defined(AMREX_USE_HIP) && defined(AMREX_USE_FFT)
        AMREX_ROCFFT_SAFE_CALL(rocfft_cleanup());
#endif
    }
}

void Clear ()
{
    for (auto& [k, p] : s_plans_d) {
        Plan<double>::destroy_vendor_plan(p);
    }

    for (auto& [k, p] : s_plans_f) {
        Plan<float>::destroy_vendor_plan(p);
    }
}

PlanD* get_vendor_plan_d (Key const& key)
{
    if (auto found = s_plans_d.find(key); found != s_plans_d.end()) {
        return &(found->second);
    } else {
        return nullptr;
    }
}

PlanF* get_vendor_plan_f (Key const& key)
{
    if (auto found = s_plans_f.find(key); found != s_plans_f.end()) {
        return &(found->second);
    } else {
        return nullptr;
    }
}

void add_vendor_plan_d (Key const& key, PlanD plan)
{
    s_plans_d[key] = plan;
}

void add_vendor_plan_f (Key const& key, PlanF plan)
{
    s_plans_f[key] = plan;
}

}

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
