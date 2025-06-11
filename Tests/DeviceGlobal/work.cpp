
#include "global_vars.H"

void work ()
{
    amrex::Gpu::PinnedVector<amrex::Long> pv;
    pv.resize(5,0);
    auto* p = pv.data();
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int)
    {
        p[0] = dg_x;
        for (int n = 0; n < 4; ++n) {
            p[1+n] = dg_y[n];
        }
    });
    amrex::Gpu::streamSynchronize();
    AMREX_ALWAYS_ASSERT(pv[0] == 1 &&
                        pv[1] == 100 &&
                        pv[2] == 101 &&
                        pv[3] == 102 &&
                        pv[4] == 103);
}

void work2 ()
{
    amrex::Gpu::PinnedVector<amrex::Long> pv;
    pv.resize(5,0);
    amrex::Gpu::memcpy_from_device_global_to_host_async
        (pv.data(), dg_x, sizeof(amrex::Long));
    amrex::Gpu::memcpy_from_device_global_to_host_async
        (pv.data()+1, dg_y, sizeof(amrex::Long));
    amrex::Gpu::memcpy_from_device_global_to_host_async
        (pv.data()+2, dg_y, sizeof(amrex::Long)*3, sizeof(amrex::Long));
    amrex::Gpu::streamSynchronize();
    AMREX_ALWAYS_ASSERT(pv[0] == 2 &&
                        pv[1] == 200 &&
                        pv[2] == 201 &&
                        pv[3] == 202 &&
                        pv[4] == 203);
}
