
#include "global_vars.H"

void init ()
{
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int)
    {
        dg_x = 1;
        for (int n = 0; n < 4; ++n) {
            dg_y[n] = 100 + n;
        }
    });

    amrex::Gpu::streamSynchronize();
}

void init2 ()
{
    amrex::Gpu::PinnedVector<amrex::Long> pv{2,200,201,202,203};
    amrex::Gpu::memcpy_from_host_to_device_global_async
        (dg_x, pv.data(), sizeof(amrex::Long));
    amrex::Gpu::memcpy_from_host_to_device_global_async
        (dg_y, pv.data()+1, sizeof(amrex::Long));
    amrex::Gpu::memcpy_from_host_to_device_global_async
        (dg_y, pv.data()+2, sizeof(amrex::Long)*3, sizeof(amrex::Long));
    amrex::Gpu::streamSynchronize();
}
