#include <AMReX_GpuElixir.H>
#include <AMReX_GpuDevice.H>

namespace amrex::Gpu {

void
Elixir::clear ()
{
#if defined(AMREX_USE_GPU)
    if (Gpu::inLaunchRegion())
    {
        for (auto const& pa : m_pa) {
            if (pa.first != nullptr) { Gpu::freeAsync(pa.second, pa.first); }
        }
    }
    else
#endif
    {
        for (auto const& pa : m_pa) {
            pa.second->free(pa.first);
        }
    }
    m_pa.clear();
}

}
