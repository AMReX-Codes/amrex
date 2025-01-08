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

SubHelper::SubHelper (Box const& domain)
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(domain);
#elif (AMREX_SPACEDIM == 2)
    if (domain.length(0) == 1) {
        m_case = case_1n;
    }
#else
    if (domain.length(0) == 1 && domain.length(1) == 1) {
        m_case = case_11n;
    } else if (domain.length(0) == 1 && domain.length(2) == 1) {
        m_case = case_1n1;
    } else if (domain.length(0) == 1) {
        m_case = case_1nn;
    } else if (domain.length(1) == 1) {
        m_case = case_n1n;
    }
#endif
}

Box SubHelper::make_box (Box const& box) const
{
    return Box(make_iv(box.smallEnd()), make_iv(box.bigEnd()), box.ixType());
}

Periodicity SubHelper::make_periodicity (Periodicity const& period) const
{
    return Periodicity(make_iv(period.intVect()));
}

bool SubHelper::ghost_safe (IntVect const& ng) const
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(ng,this);
    return true;
#elif (AMREX_SPACEDIM == 2)
    if (m_case == case_1n) {
        return (ng[0] == 0);
    } else {
        return true;
    }
#else
    if (m_case == case_11n) {
        return (ng[0] == 0) && (ng[1] == 0);
    } else if (m_case == case_1n1) {
        return (ng[0] == 0);
    } else if (m_case == case_1nn) {
        return (ng[0] == 0);
    } else if (m_case == case_n1n) {
        return (ng[1] == 0);
    } else {
        return true;
    }
#endif
}

IntVect SubHelper::make_iv (IntVect const& iv) const
{
    return this->make_array(iv);
}

IntVect SubHelper::make_safe_ghost (IntVect const& ng) const
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(this);
    return ng;
#elif (AMREX_SPACEDIM == 2)
    if (m_case == case_1n) {
        return IntVect{0,ng[1]};
    } else {
        return ng;
    }
#else
    if (m_case == case_11n) {
        return IntVect{0,0,ng[2]};
    } else if (m_case == case_1n1) {
        return IntVect{0,ng[1],ng[2]};
    } else if (m_case == case_1nn) {
        return IntVect{0,ng[1],ng[2]};
    } else if (m_case == case_n1n) {
        return IntVect{ng[0],0,ng[2]};
    } else {
        return ng;
    }
#endif
}

BoxArray SubHelper::inverse_boxarray (BoxArray const& ba) const
{ // sub domain order -> original domain order
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(this);
    return ba;
#elif (AMREX_SPACEDIM == 2)
    AMREX_ALWAYS_ASSERT(m_case == case_1n);
    BoxList bl = ba.boxList();
    // sub domain order: y, x
    for (auto& b : bl) {
        auto const& lo = b.smallEnd();
        auto const& hi = b.bigEnd();
        b.setSmall(IntVect(lo[1],lo[0]));
        b.setBig  (IntVect(hi[1],hi[0]));
    }
    return BoxArray(std::move(bl));
#else
    BoxList bl = ba.boxList();
    if (m_case == case_11n) {
        // sub domain order: z, x, y
        for (auto& b : bl) {
            auto const& lo = b.smallEnd();
            auto const& hi = b.bigEnd();
            b.setSmall(IntVect(lo[1],lo[2],lo[0]));
            b.setBig  (IntVect(hi[1],hi[2],hi[0]));
        }
    } else if (m_case == case_1n1) {
        // sub domain order: y, x, z
        for (auto& b : bl) {
            auto const& lo = b.smallEnd();
            auto const& hi = b.bigEnd();
            b.setSmall(IntVect(lo[1],lo[0],lo[2]));
            b.setBig  (IntVect(hi[1],hi[0],hi[2]));
        }
    } else if (m_case == case_1nn) {
        // sub domain order: y, z, x
        for (auto& b : bl) {
            auto const& lo = b.smallEnd();
            auto const& hi = b.bigEnd();
            b.setSmall(IntVect(lo[2],lo[0],lo[1]));
            b.setBig  (IntVect(hi[2],hi[0],hi[1]));
        }
    } else if (m_case == case_n1n) {
        // sub domain order: x, z, y
        for (auto& b : bl) {
            auto const& lo = b.smallEnd();
            auto const& hi = b.bigEnd();
            b.setSmall(IntVect(lo[0],lo[2],lo[1]));
            b.setBig  (IntVect(hi[0],hi[2],hi[1]));
        }
    } else {
        amrex::Abort("SubHelper::inverse_boxarray: how did this happen?");
    }
    return BoxArray(std::move(bl));
#endif
}

IntVect SubHelper::inverse_order (IntVect const& order) const
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(this);
    return order;
#elif (AMREX_SPACEDIM == 2)
    amrex::ignore_unused(this);
    return IntVect(order[1],order[0]);
#else
    auto translate = [&] (int index) -> int
    {
        int r = index;
        if (m_case == case_11n) {
            // sub domain order: z, x, y
            if (index == 0) {
                r = 2;
            } else if (index == 1) {
                r = 0;
            } else {
                r = 1;
            }
        } else if (m_case == case_1n1) {
            // sub domain order: y, x, z
            if (index == 0) {
                r = 1;
            } else if (index == 1) {
                r = 0;
            } else {
                r = 2;
            }
        } else if (m_case == case_1nn) {
            // sub domain order: y, z, x
            if (index == 0) {
                r = 1;
            } else if (index == 1) {
                r = 2;
            } else {
                r = 0;
            }
        } else if (m_case == case_n1n) {
            // sub domain order: x, z, y
            if (index == 0) {
                r = 0;
            } else if (index == 1) {
                r = 2;
            } else {
                r = 1;
            }
        }
        return r;
    };

    IntVect iv;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        iv[idim] = translate(order[idim]);
    }
    return iv;
#endif
}

GpuArray<int,3> SubHelper::xyz_order () const
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(this);
    return GpuArray<int,3>{0,1,2};
#elif (AMREX_SPACEDIM == 2)
    if (m_case == case_1n) {
        return GpuArray<int,3>{1,0,2};
    } else {
        return GpuArray<int,3>{0,1,2};
    }
#else
    if (m_case == case_11n) {
        return GpuArray<int,3>{1,2,0};
    } else if (m_case == case_1n1) {
        return GpuArray<int,3>{1,0,2};
    } else if (m_case == case_1nn) {
        return GpuArray<int,3>{2,0,1};
    } else if (m_case == case_n1n) {
        return GpuArray<int,3>{0,2,1};
    } else {
        return GpuArray<int,3>{0,1,2};
    }
#endif
}

}
